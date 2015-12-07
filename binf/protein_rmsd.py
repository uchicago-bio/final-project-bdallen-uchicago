#!/usr/bin/env python

import os
import argparse
import itertools
import random
import urllib
import multiprocessing
import math

from pdbremix import pdbatoms, rmsd, v3
from pdbremix.data import res_name_to_char

import _pathfix

from binf.conserved_residues import iter_residues_for_conserved_alignment
from binf.pymol_batch import run_pymol_script

PDB_SITE = "http://www.rcsb.org/pdb/files"


class StructuralAlign(object):
    """
    Class to encapsulate different strategies for structural alignment of two
    pdbremix.pdbatoms.Soup objects using RMSD. Keeps track of the best
    solution found (the one with the lowest RMSD score).

    The main challenge is selecting which residues to align, since RMSD
    requires aligning residue chains of equal length. Trying all possible
    subsequences of the requested length is generally computationally
    infeasible.

    @param local_segment: align subsets of this length
    @param residues1: use this subset of residues from the soup instead of
                      using the standard C atom for each residue
    """
    def __init__(self, soup1, soup2, local_segment=None,
                 atoms1=None, atoms2=None):
        self.soup1 = soup1
        self.soup2 = soup2
        self.local_segment = local_segment

        self.atoms1 = atoms1 or get_ca_atom_positions(soup1)
        self.atoms2 = atoms2 or get_ca_atom_positions(soup2)

        self.len1 = len(self.atoms1)
        self.len2 = len(self.atoms2)

        if self.local_segment is not None:
            self.align_len = self.local_segment
            if self.len1 < self.align_len:
                raise ValueError("local_segment larger than size of soup1 %d"
                                 % self.len1)
            if self.len2 < self.align_len:
                raise ValueError("local_segment larger than size of soup2 %d"
                                 % self.len2)
        else:
            self.align_len = min(self.len1, self.len2)

        self.best = None

    def _rmsd_rot(self, idx_set1, idx_set2):
        """
        Order preserving RMSD by residue numbers.
        """
        res1 = [r for r in self.atoms1 if r.i in idx_set1]
        res2 = [r for r in self.atoms2 if r.i in idx_set2]
        return self._rmsd_rot_res(res1, res2)

    def _rmsd_rot_res(self, res1, res2):
        """
        RMSD of arbitrary residue lists, may be out of order or a subset
        of total.
        """
        assert_eq(len(res1), len(res2))
        coords1 = [r.position for r in res1]
        coords2 = [r.position for r in res2]
        center1, coords1 = center_vlist(coords1)
        center2, coords2 = center_vlist(coords2)
        score, rot_matrix = rmsd.calc_rmsd_rot(coords1, coords2)
        return StructuralAlignSolution(self,
                                       score=score,
                                       rot_matrix=rot_matrix,
                                       center1=center1,
                                       center2=center2,
                                       res1=res1,
                                       res2=res2)

    def align_contiguous(self):
        """
        Try all contiguous alignments, and keep track of the solution with
        the lowest rmsd score. Returns the best solutions found by this
        method, and saves it to instance best if it's the best found overall.

        The number of calc_rmsd_rot calls is
        (len1 - align_len) * (len2 - align_len), which is O(n^2).
        """
        local_best = None
        for offset1 in xrange(self.len1 - self.align_len + 1):
            idx_set1 = set(range(offset1, offset1 + self.align_len))
            for offset2 in xrange(self.len2 - self.align_len + 1):
                idx_set2 = set(range(offset2, offset2 + self.align_len))
                s = self._rmsd_rot(idx_set1, idx_set2)
                if local_best is None or s.score < local_best.score:
                    s.strategy = ("contigous", offset1, offset2)
                    local_best = s
        if self.best is None or local_best.score < self.best.score:
            self.best = local_best
        return local_best

    def align_random(self, attempts=1000):
        """
        Try aligning a random subsequence of each soup.
        """
        local_best = None
        for i in xrange(attempts):
            idx_set1 = random_subseq_idx(self.len1, self.align_len)
            idx_set2 = random_subseq_idx(self.len2, self.align_len)
            s = self._rmsd_rot(idx_set1, idx_set2)
            if local_best is None or s.score < local_best.score:
                s.strategy = ("random", i)
                local_best = s
        if self.best is None or local_best.score < self.best.score:
            self.best = local_best
        return local_best

    def align_conserved(self, start=0, stop=None):
        """
        Try all possible conserved alignments, or optionally a slice of the
        total possibilities for use with multiprocessing. Does not preserve
        residue order, and not all atoms can be matched so it will typical
        exclude many atoms.

        Note: ignores local_segment.
        """
        # It's inefficient to calculate this separately for each
        # process, but reduces the data that needs to be sent via IPC
        # and this is not a significant portion of the total work for
        # soups of typical size.
        total, it = iter_residues_for_conserved_alignment(self.atoms1,
                                                          self.atoms2)
        if start > 0 or stop is not None:
            it = itertools.islice(it, start, stop)

        local_best = None
        for res1bytype, res2bytype in it:
            # flatten tuples, they are list of list by type
            res1 = [r for rstype in res1bytype for r in rstype]
            res2 = [r for rstype in res2bytype for r in rstype]
            s = self._rmsd_rot_res(res1, res2)
            if local_best is None or s.score < local_best.score:
                s.strategy = ("conserved", start, stop)
                local_best = s
        if self.best is None or local_best.score < self.best.score:
            self.best = local_best
        return local_best

    def get_transformed_soup1(self, solution=None):
        """
        Get a transformed copy of soup1 based on the best alignment,
        so it will align with unchanged soup2.
        """
        if solution is None:
            if self.best is None:
                raise ValueError("must be called after finding solution")
            solution = self.best
        soup1x = self.soup1.copy()
        soup1x.transform(v3.translation(-solution.center1))
        soup1x.transform(solution.rot_matrix)
        soup1x.transform(v3.translation(solution.center2))
        return soup1x

    def get_transformed_soup2(self, solution=None):
        """
        Get a transformed copy of soup2 based on the best alignment,
        so it will align with unchanged soup1.
        """
        if solution is None:
            if self.best is None:
                raise ValueError("must be called after finding solution")
            solution = self.best
        soup2x = self.soup2.copy()
        soup2x.transform(v3.translation(-solution.center2))
        soup2x.transform(v3.left_inverse(solution.rot_matrix))
        soup2x.transform(v3.translation(solution.center2))
        return soup2x

    def get_combined_soup(self, solution=None):
        """
        Get a transformed copy of soup2 based on the best alignment,
        so it will align with unchanged soup1, combined with soup1.

        Hacks the chain ids of soup2 from A, B, C... to Z, Y, X..., assuming
        that both soups will be using the same chain id and it will be
        useful to distinguish them in visualizations.
        """
        soup2x = self.get_transformed_soup2(solution)
        soup1x = self.soup1.copy()
        for res in soup1x.residues():
            # A->Z, B->Y, etc
            res.set_chain_id(chr(ord("Z") - (ord(res.chain_id) - ord("A"))))
        soup2x.insert_soup(soup2x.n_residue(), soup1x)
        return soup2x


class AtomPosition(object):
    """
    @ivar i: index in sequence of all atoms
    @ivar residue_name: three letter code for the residue
    @ivar atom_name: three letter code for the residue
    @ivar position: position vector of center of mass (CA atom) of the
                    residue, as a v3 vector
    """
    def __init__(self, i, residue_type, atom_name, position):
        self.i = i
        self.residue_type = residue_type
        self.atom_name = atom_name
        self.position = position
        self.letter = res_name_to_char.get(residue_type)
        if self.letter in (None, "<", ">"):
            raise ValueError("Unsupported residue type: %s" % residue_type)

    def get_position(self, matrix):
        """
        Get position of residue transformed by matrix, as a new v3 vector.
        Does not modify instance variable.
        """
        return v3.transform(matrix, self.position)

    def __repr__(self):
        return ("<AtomPosition %d %s %s>" % (self.i, self.residue_type,
                                             self.atom_name))


class StructuralAlignSolution(object):
    def __init__(self, structural_align, score,
                 rot_matrix, center1, center2,
                 res1, res2, strategy=None):
        self.sa = structural_align
        self.score = score
        self.rot_matrix = rot_matrix
        self.center1 = center1
        self.center2 = center2
        self.res1 = res1
        self.res2 = res2
        self.strategy = strategy

    def print_multiline(self):
        print "Solution", repr(self.strategy)
        print "   score:", self.score
        print "     rot:", self.rot_matrix
        print " center1:", self.center1
        print " center2:", self.center2
        ar1, ar2 = self.get_aligned_residues()
        print "  %ident:", percent_identity(ar1, ar2)
        print "", ar1
        print "", ar2
        ar1, ar2 = self.get_aligned_residues_numbered()
        print "", ar1
        print "", ar2

    def get_percent_identity(self):
        ar1, ar2 = self.get_aligned_residues()
        return percent_identity(ar1, ar2)

    def get_aligned_residues(self):
        aligned1 = "".join(r.letter for r in self.res1)
        aligned2 = "".join(r.letter for r in self.res2)
        return (aligned1, aligned2)

    def get_aligned_residues_numbered(self):
        aligned1 = ",".join("%s%03d" % (r.letter, r.i) for r in self.res1)
        aligned2 = ",".join("%s%03d" % (r.letter, r.i) for r in self.res2)
        return (aligned1, aligned2)


def percent_identity(s1, s2):
    assert len(s1) == len(s2)
    ident_count = 0
    for a, b in itertools.izip(s1, s2):
        if a == b:
            ident_count += 1
    return 100.0 * ident_count / len(s1)


def _job_align_conserved(args):
    sa, start, stop = args
    return sa.align_conserved(start, stop)


def load_soups(*pdbs):
    paths = [_pdbpath(pdb) for pdb in pdbs]
    download_pdbs(*paths)
    return [pdbatoms.Soup(path) for path in paths]


def _main():
    args = parse_args()

    pdb1_path = _pdbpath(args.pdb1)
    pdb2_path = _pdbpath(args.pdb2)
    download_pdbs(pdb1_path, pdb2_path)

    soup1 = pdbatoms.Soup(pdb1_path)
    soup2 = pdbatoms.Soup(pdb2_path)
    sa = StructuralAlign(soup1, soup2, args.local_segment)
    best_contiguous = sa.align_contiguous()
    best_contiguous.print_multiline()
    if args.random_attempts:
        best_random = sa.align_random(args.random_attempts)
        best_random.print_multiline()
    if args.conserved:
        ncpus = multiprocessing.cpu_count()
        total, it = iter_residues_for_conserved_alignment(sa.atoms1,
                                                          sa.atoms2)
        per_job = int(math.ceil(float(total) / ncpus))
        jobs = [(sa, start, start + per_job)
                for start in xrange(0, total, per_job)]
        assert len(jobs) == ncpus
        pool = multiprocessing.Pool(ncpus)
        results = pool.map(_job_align_conserved, jobs)
        best_conserved = None
        for result in results:
            result.print_multiline()
            if best_conserved is None or result.score < best_conserved.score:
                best_conserved = result
        if best_conserved.score < sa.best.score:
            sa.best = best_conserved

    print "======== best =========="
    sa.best.print_multiline()

    if args.pdb1_transform_out:
        soup1x = sa.get_transformed_soup1()
        pdb1_outpath = _pdbpath(args.pdb1_transform_out)
        soup1x.write_pdb(pdb1_outpath)
        if args.png_out:
            pymol_cartoon(pdb1_outpath, pdb2_path, args.png_out)
    elif args.pdb2_transform_out:
        soup2x = sa.get_transformed_soup2()
        pdb2_outpath = _pdbpath(args.pdb2_transform_out)
        soup2x.write_pdb(pdb2_outpath)
        if args.png_out:
            pymol_cartoon(pdb1_path, pdb2_outpath, args.png_out)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=
        "Align two pdb's using RMSD"
    )
    parser.add_argument("pdb1", help="code of first pdb")
    parser.add_argument("pdb2", help="code of second pdb")
    def _pos_int_or_none(value):
        value = int(value)
        if value <= 0:
            return None
        return value
    parser.add_argument("-l", "--local-segment", type=_pos_int_or_none,
                        help="segment length to align. if not specified, "
                             "do global alignment, i.e. align all of "
                             "shorter sequence")
    parser.add_argument("-c", "--conserved", action="store_true",
                        help="do parallel conserved alignment")
    parser.add_argument("-r", "--random-attempts", type=int,
                        help="attempt random alignments")
    parser.add_argument("--pdb1-transform-out",
                        help="output a transformed copy of pdb1 to this path")
    parser.add_argument("--pdb2-transform-out",
                        help="output a transformed copy of pdb2 to this path")
    parser.add_argument("--png-out",
                        help="output a visualization of the two proteins")

    args = parser.parse_args(argv)

    if args.conserved and args.local_segment is not None:
        parser.error("options --conserved and --local-segment are exclusive")
    if args.pdb1_transform_out and args.pdb2_transform_out:
        parser.error("specify only one of --pdbN-transform-out")
    if (args.png_out and not args.pdb1_transform_out
    and not args.pdb2_transform_out):
        parser.error("--png-out requires one of --pdbN-transform-out")
    return args


def _pdbpath(path_or_code):
    """
    Get absolute path to pdb argument. Relative paths are assumed to be
    relative to the project data dir. Bare codes are also accepted, and
    saved in (data_path)/(code).pdb.
    """
    if not path_or_code.endswith(".pdb"):
        path_or_code += ".pdb"
    if not path_or_code.startswith("/"):
        path_or_code = os.path.join(_pathfix.data_path, path_or_code)
    return path_or_code


def _pdbname(path):
    return os.path.splitext(os.path.basename(path))[0]


def get_ca_atom_positions(soup):
    """
    Get CA atom positions for all residues in a soup. Ignores residues that
    don't have a CA atom.
    """
    rps = []
    for r in soup.residues():
        try:
            ca_atom = r.atom("CA")
            if (r.type not in res_name_to_char
            or r.type in ("NME", "ACE")):
                raise ValueError("ERR: unexpected residue with CA (%s %s)",
                                 (r.type, r.tag))
            rp = AtomPosition(len(rps), r.type, "CA", ca_atom.pos)
            rps.append(rp)
        except KeyError:
            pass
    return rps


def _get_alignment_idx(tb, fasta1, fasta2):
    """
    Given a traceback from an S-W alignment of the fasta sequence, build
    a list of indexes corresponding to the alignment.
    """
    n = m = 0
    idx1 = set()
    idx2 = set()
    for i in xrange(len(tb)):
        c1 = tb.aligna[i]
        c2 = tb.alignb[i]
        if c1 == "_":
            # gap in sequence 1, skip one letter in fasta2
            m += 1
            continue
        elif c2 == "_":
            # gap in sequence 2, skip one letter in fasta1
            n += 1
            continue
        idx1.add(n)
        idx2.add(m)
        n += 1
        m += 1
    return idx1, idx2


def center_vlist(vlist):
    """
    Center a list of v3 vectors and return (center, centered_vlist_copy)
    """
    center = v3.get_center(vlist)
    center_matrix = v3.translation(-center)
    return (center, [v3.transform(center_matrix, v) for v in vlist])


def pymol_cartoon(pdb_path1, pdb_path2, outpath, surface=False):
    data = dict(
        pdb_name1=_pdbname(pdb_path1),
        pdb_name2=_pdbname(pdb_path2),
        pdb_path1=pdb_path1,
        pdb_path2=pdb_path2,
        outpath=outpath,
        bgcolor="white"
    )

    if surface:
        data["show_type"] = "surface"
    else:
        data["show_type"] = "cartoon"

    script = """
load %(pdb_path1)s
load %(pdb_path2)s
zoom
hide everything
show %(show_type)s
bg_color %(bgcolor)s
color yellow, %(pdb_name1)s
color cyan, %(pdb_name2)s
set antialias,1
ray
png %(outpath)s
""" % data

    return run_pymol_script(script)


def random_subseq_idx(seq_len, sub_len):
    """
    Construct a random set of indexes of length @sub_len for indexes in
    range(0, seq_len)
    """
    option_list = range(seq_len)
    chosen_set = set()

    while len(chosen_set) < sub_len:
        idx = random.choice(option_list)
        option_list.remove(idx)
        chosen_set.add(idx)

    return chosen_set


def download_pdbs(*pdb_paths):
    """
    Fetches PDB files using HTTP from a list of pdb-codes/text-files.
    """
    for pdb_path in pdb_paths:
        if os.path.isfile(pdb_path):
            continue
        pdb_name = _pdbname(pdb_path)
        url = os.path.join(PDB_SITE, pdb_name + ".pdb")
        _, info = urllib.urlretrieve(url, pdb_path)
        if info.gettype() != "text/plain":
            #os.unlink(pdb_path)
            raise Exception("error downloading pdb '%s' to '%s', got type %s"
                            % (pdb_name, pdb_path, info.gettype()))


def assert_eq(a, b):
    assert a == b, "%r != %r" % (a, b)


if __name__ == "__main__":
    _main()
