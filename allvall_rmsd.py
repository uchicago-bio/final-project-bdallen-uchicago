#!/usr/bin/env python
"""
Compute pairwise RMSD of all pdbs listed in the input file, restricting
the calculation to the atoms provided in the second argument, which must
be present in all pdbs. Designed to be used in conjunction with
the "ALL" atom list from calc_presenting_surface.py.

Usage:

 find /path/to/project-data -name '*_complex_0001.pdb' > pdb_list.txt
 ./allvall_rmsd.py pdb_list.txt atom_list.txt output/

"""

import sys
import multiprocessing
import os.path
import traceback

from pdbremix import pdbatoms, rmsd, v3


def rmsd_all(pdb_path, pdb_list, atom_list, outdir):
    """
    Compute RMSD between pdb_path and every pdb in pdb_list, writing output
    to outpath/rmsd_{basename(pdb_path)}.txt, and restricting the comparison
    to atoms in atom_list, which are assumed to be present in all pdbs.
    """
    pdb_name1 = get_pdb_name(pdb_path)
    soup1 = pdbatoms.Soup(pdb_path)
    coords1 = get_atom_coords(soup1, atom_list)
    center1, coords1 = center_vlist(coords1)
    outpath = os.path.join(outdir, "rmsd_%s.txt" % pdb_name1)
    with open(outpath, "w") as out:
        for pdb_path2 in pdb_list:
            # leave self v self in as a sanity check, it's only
            # 1/3600 of total work anyway
            #if pdb_path2 == pdb_path:
            #    continue
            pdb_name2 = get_pdb_name(pdb_path2)
            soup2 = pdbatoms.Soup(pdb_path2)
            coords2 = get_atom_coords(soup2, atom_list)
            center2, coords2 = center_vlist(coords2)
            score, rot_matrix = rmsd.calc_rmsd_rot(coords1, coords2)
            out.write("%0.6f %s %s\n" % (score, pdb_name1, pdb_name2))


def get_atom_coords(soup, atom_list):
    atom_positions = []
    for atom in atom_list:
        chain, resi, resn, name = atom
        residue = soup.residue_by_tag("%s:%d" % (chain, resi))
        assert residue.num == resi, "%s != %s" % (residue.num, resi)
        if resn != "X":
            assert residue.type == resn, "%s != %s" % (residue.type, resn)
        atom = residue.atom(name)
        atom_positions.append(atom.pos)
    return atom_positions


def center_vlist(vlist):
    """
    Center a list of v3 vectors and return (center, centered_vlist_copy)
    """
    center = v3.get_center(vlist)
    center_matrix = v3.translation(-center)
    return (center, [v3.transform(center_matrix, v) for v in vlist])


def _job_rmsd_all(args):
    try:
        return rmsd_all(*args)
    except Exception:
        traceback.print_exc()
        raise


def get_pdb_name(pdb_path):
    return os.path.splitext(os.path.basename(pdb_path))[0]


def _parse_atom_list(infile):
    atom_list = []
    for line in infile:
        line = line.strip()
        chain, resi, resn, atom_name = line.split()
        atom_list.append((chain, int(resi), resn, atom_name))
    return atom_list


def _main():
    if len(sys.argv) != 4:
        print "Usage: %s pdb_list.txt atom_list.txt outdir" % sys.argv[0]
        sys.exit(1)

    with open(sys.argv[1]) as f:
        pdb_list = [line.strip() for line in f.readlines()]

    with open(sys.argv[2]) as f:
        atom_list = _parse_atom_list(f)

    outdir = sys.argv[3]

    jobs = [(pdb_path, pdb_list, atom_list, outdir) for pdb_path in pdb_list]

    ncpus = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(ncpus)
    pool.map(_job_rmsd_all, jobs)


if __name__ == '__main__':
    _main()
