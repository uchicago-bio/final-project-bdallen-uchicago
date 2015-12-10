#!/usr/bin/env python
"""
Find a common presenting surface for all the hlac complexes, allowing for
different residue names in mutated positions and in the docked peptide chain.
Uses X in the output to denote when any residue name is allowed.

Also calculates the union of positions and statistics about the presenting
surfaces extracted from each pdb.

Usage:

 find /path/to/project-data -name '*_complex_0001.pdb' > pdb_list.txt
 ./calc_presenting_surface.py pdb_list.txt | tee psurface_stats.txt | less

"""
import sys

from binf.pymol_batch import pymol_get_presenting_surface_atoms


MUTATED_POSITIONS = set([91, 97, 156])

def get_restricted_psurface(pdb_path, mutated_positions, within=15):
    """
    Get psurface list of (chain, resi, resn, atom_name), but with mutated
    positions replaced resn replaced with X to allow the residue to vary
    and still find the largest set of common atoms.

    Note that when using within=10, residue 91 wasn't included in all, which
    suggests that a wider area is needed.

    Also replaces the peptide amino acids with X, since they also vary among
    samples.
    """
    psurface = pymol_get_presenting_surface_atoms(pdb_path, "c", within=within)
    restricted = []
    mutation_atoms_added = set()
    for atom in psurface:
        chain, resi, resn, atom_name = atom
        if resi in mutated_positions:
            resn = "X"
        if chain == "C":
            resn = "X"
        restricted.append((chain, resi, resn, atom_name))
    return restricted


def compare_psurfaces(pdb_list):
    first_path = pdb_list[0]
    psurface = get_restricted_psurface(first_path, MUTATED_POSITIONS)
    all_psurface = set(psurface)
    any_psurface = set(psurface)
    lengths = [len(psurface)]
    all_len = len(all_psurface)
    for i, pdb_path in enumerate(pdb_list[1:]):
        psurface = get_restricted_psurface(pdb_path, MUTATED_POSITIONS)
        lengths.append(len(psurface))
        psurface = set(psurface)
        all_psurface &= psurface
        any_psurface |= psurface
        new_all_len = len(all_psurface)
        if all_len - new_all_len > 10:
            print ("WARN: blacksheep pdb[%d] '%s', %d -> %d (%d)"
                   % (i, pdb_path, all_len, new_all_len, len(psurface)))
        all_len = new_all_len
    all_psurface_sorted = list(all_psurface)
    all_psurface_sorted.sort()
    any_psurface_sorted = list(any_psurface)
    any_psurface_sorted.sort()
    print "------ STATS -------"
    print "pdb count:", len(pdb_list)
    print "avg len:", float(sum(lengths)) / len(lengths)
    print "med len:", median(lengths)
    print "max len:", max(lengths)
    print "min len:", min(lengths)
    print "all len:", len(all_psurface)
    print "any len:", len(any_psurface)
    print
    print "------ ALL -------"
    for atom in all_psurface_sorted:
        print " ".join(str(x) for x in atom)
    print
    print "------ ANY -------"
    for atom in any_psurface_sorted:
        print " ".join(str(x) for x in atom)


def median(alist):
    alen = len(alist)
    if alen % 2 == 1:
        return float(alist[alen/2])
    return (alist[alen/2 - 1] + alist[alen/2]) / 2.0


def _main():
    if len(sys.argv) != 2:
        print "Usage: %s pdb_list.txt" % sys.argv[0]
        sys.exit(1)

    with open(sys.argv[1]) as f:
        pdb_list = [line.strip() for line in f.readlines()]

    compare_psurfaces(pdb_list)


if __name__ == '__main__':
    _main()
