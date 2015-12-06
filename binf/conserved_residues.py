#!/usr/bin/env python
"""
Split a pair of PDBs into multiple pairs of PDBs based on aligning conserved
residues of certain types, as input into parallel alignment.
"""

import os
import argparse
import itertools
import random
import subprocess
import urllib
from collections import defaultdict
import itertools
import math

from pdbremix import pdbatoms, rmsd, v3
from pdbremix.data import res_name_to_char

import _pathfix


def iter_residues_for_conserved_alignment(res1, res2):
    """
    Get an iterator for all possible conserved alignments between two
    lists of residue positions, res1 and res2, where each residue is
    paired with another residue of the same type. Does not preserve
    the order of residues.
    """
    residues_by_type = [defaultdict(list), defaultdict(list)]
    for i, res in enumerate((res1, res2)):
        for r in res:
            residues_by_type[i][r.residue_type].append(r)
    min_counts = defaultdict(int)
    for rtype in residues_by_type[0].iterkeys():
        min_counts[rtype] = min(len(residues_by_type[i][rtype])
                                for i in (0, 1))

    #for rtype in sorted(min_counts.keys(), key=lambda k: min_counts[k],
    #                    reverse=True):
    #    print rtype, min_counts[rtype],
    #    print len(residues_by_type[0][rtype]),
    #    print len(residues_by_type[1][rtype])

    all_types = [rtype for rtype in min_counts.iterkeys()
                 if min_counts[rtype] > 0]

    #print all_types

    count1, iter1 = _iter_residues(all_types, residues_by_type[0], min_counts)
    count2, iter2 = _iter_residues(all_types, residues_by_type[1], min_counts,
                                   permute=True)
    return (count1 * count2, itertools.product(iter1, iter2))


def _iter_residues(types, residues_by_type,  max_counts, permute=False):
    """
    Get an iterator for all possible combinations of residues of the specified
    types, up to the the specified maximums for each. Optionally also permute
    each residue set. Also returns the total number of items in the iterator,
    useful for estimating the cost of checking all possibilities.

    @param types: list of residue types, three letter codes, all caps
    @param residues_by_type: dict mapping residue type to list of
                             ResiduePosition objects
    @param max_counts: dict mapping type to int, max count for residues of that
                       type
    @param permute: if specified, use permutations instead of combinations.
                    If combining with another soup, permute only one of them.

    Returns (total_count, iterator)
    """
    type_iters = []
    total_count = 1
    for t in types:
        rs = residues_by_type[t]
        max_count = max_counts[t]
        if permute:
            type_iter = itertools.permutations(rs, max_count)
            count = permutation_count(len(rs), max_count)
        else:
            type_iter = itertools.combinations(rs, max_count)
            count = combination_count(len(rs), max_count)
        type_iters.append(type_iter)
        total_count *= count
    return total_count, itertools.product(*type_iters)


def permutation_count(n, k):
    """
    Calculate the number of permutations of size k from a set of size n,
    i.e. n! / (n-k)!.
    """
    assert k <= n
    total = 1
    for i in xrange(n, n-k, -1):
        total *= i
    return total


def combination_count(n, k):
    """
    Calculate nCk, i.e. n! / [k!(n-k)!]
    """
    assert k <= n
    return permutation_count(n, k) / math.factorial(k)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=
        "Generate new pdbs for parallel paired residue alignment"
    )
    parser.add_argument("pdb1", help="code of first pdb")
    parser.add_argument("pdb2", help="code of second pdb")

    args = parser.parse_args(argv)
    return args


def _test_conserved_align_iter_main():
    args = parse_args()

    # avoid circular import
    from binf.protein_rmsd import get_residue_positions, _pdbpath

    pdb1_path = _pdbpath(args.pdb1)
    pdb2_path = _pdbpath(args.pdb2)

    soup1 = pdbatoms.Soup(pdb1_path)
    soup2 = pdbatoms.Soup(pdb2_path)

    res1 = get_residue_positions(soup1)
    res2 = get_residue_positions(soup2)

    count, it = iter_residues_for_conserved_alignment(res1, res2)
    count2 = len(list(it))
    assert count == count2
    print "count", count


if __name__ == '__main__':
    _test_conserved_align_iter_main()
