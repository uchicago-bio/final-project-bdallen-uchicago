#!/usr/bin/env python
"""
Script to do structure RMSD alignment of the presenting surfaces of two
HLA molecules.
"""
import os.path

from pdbremix import pdbatoms

import _pathfix

from binf.protein_rmsd import StructuralAlign, AtomPosition
from pymol_batch import pymol_get_presenting_surface_atoms, run_pymol_script

BACKBONE_ATOMS = ["N", "CA", "C", "O"]

def psurface_rmsd(pdb_path1, pdb_path2):
    psurface1 = pymol_get_presenting_surface_atoms(pdb_path1, "c")
    #psurface2 = pymol_get_presenting_surface_atoms(pdb_path2, "c")
    #print len(psurface1), len(psurface2)
    #min_len = min(len(psurface1), len(psurface2))
    #for i in xrange(min_len):
    #    print psurface1[i], psurface2[i]

    mutated_pos = set([_parse_mutated_position(pdb_path1),
                       _parse_mutated_position(pdb_path2)])

    soup1 = pdbatoms.Soup(pdb_path1)
    soup2 = pdbatoms.Soup(pdb_path2)
    aps1 = _get_atom_positions(soup1, psurface1, mutated_pos)
    aps2 = _get_atom_positions(soup2, psurface1, mutated_pos)

    sa = StructuralAlign(soup1, soup2, atoms1=aps1, atoms2=aps2)
    sa.align_contiguous()
    return sa


def _parse_mutated_position(pdb_path):
    # Example hlac-156-A-AAADAAAAL_complex_0001.pdb
    fname = os.path.basename(pdb_path)
    parts = fname.split("-", 3)
    return int(parts[1])


def _get_atom_positions(soup, psurface, mutated_positions=None):
    if mutated_positions is None:
        mutated_positions = set()
    else:
        mutated_positions = set(mutated_positions)
    atom_positions = []
    seen_mutated_positions = set()
    for atom in psurface:
        chain, resi, resn, name = atom
        residue = soup.residue_by_tag("%s:%d" % (chain, resi))
        assert residue.num == resi, "%s != %s" % (residue.num, resi)
        if resi in mutated_positions:
            if resi not in seen_mutated_positions:
                # add entire backbone
                for name in BACKBONE_ATOMS:
                    atom = residue.atom(name)
                    ap = AtomPosition(len(atom_positions), residue.type, name,
                                      atom.pos)
                    atom_positions.append(ap)
                # ignore it if we see it again
                seen_mutated_positions.add(resi)
            continue
        assert residue.type == resn, "%s != %s" % (residue.type, resn)
        atom = residue.atom(name)
        ap = AtomPosition(len(atom_positions), residue.type, name, atom.pos)
        atom_positions.append(ap)
    return atom_positions


def pymol_psurface_png(pdb_paths, outpath):
    """
    Create an autozoomed image of a list of pdb files.
    """
    lines = ["load %s" % path for path in pdb_paths]

    lines.append("""
zoom
hide everything
show surface
bg_color %(bgcolor)s
select psurface, all within %(within)s of chain %(chain)s
color purple, psurface
set antialias,1
ray
png %(outpath)s
""" % dict(
        outpath=outpath,
        bgcolor="white"))
    script = "\n".join(lines)
    return run_pymol_script(script)



if __name__ == '__main__':
    import sys
    sa = psurface_rmsd(sys.argv[1], sys.argv[2])
    print "RMSD:", sa.best.score
    if len(sys.argv) > 3:
        pdb_outpath = os.path.abspath(sys.argv[3])
        soupx = sa.get_combined_soup()
        print "writing combined pdb", pdb_outpath
        soupx.write_pdb(pdb_outpath)
