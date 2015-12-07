"""
Module for running pymol scripts in batch mode use subprocess.
"""
import subprocess


def run_pymol_script(script):
    """
    The script should be a new line separated list of pymol commands.

    Note that for the png command to work, the ray command must be called
    first.

    Returns (stdout_text, stderr_text)
    Raises PymolError if pymol has a nonzero return code.
    """
    p = subprocess.Popen(["pymol", "-cpq"],
                         stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout, stderr = p.communicate(script)
    if p.returncode != 0:
        raise PymolError(p.returncode, stdout, stderr)
    return stdout, stderr


class PymolError(Exception):
    def __init__(self, returncode, stdout, stderr):
        Exception.__init__(self, "error running pymol: %d %s"
                           % (returncode, stderr))
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr


def pymol_get_presenting_surface_atoms(pdb_path, chain, within=10):
    """
    Get a list of all atoms within a certain distance of the specified
    chain. Useful for finding the T-cell presenting surface of HLA molecules.

    Returns a list of tuples (chain, residue_number, residue_name, atom_name)
    """
    script = """
load %(pdb_path)s
select psurface, all within %(within)s of chain %(chain)s
iterate psurface, print 'OUTPUT',chain,resi,resn,name
""" % dict(
        pdb_path=pdb_path,
        chain=chain,
        within=within
    )
    stdout, stderr = run_pymol_script(script)
    atom_tuples = []
    in_iterate_output = False
    for line in stdout.split("\n"):
        line = line.strip()
        if not line.startswith("OUTPUT"):
            continue
        _, chain, resi, resn, name = line.split(" ")
        atom_tuples.append((chain, int(resi), resn, name))
    return atom_tuples


def pymol_png(pdb_paths, outpath, show_type="cartoon"):
    """
    Create an autozoomed image of a list of pdb files.
    """
    lines = ["load %s" % path for path in pdb_paths]

    lines.append("""
zoom
hide everything
show %(show_type)s
bg_color %(bgcolor)s
set antialias,1
ray
png %(outpath)s
""" % dict(
        pdb_path1=pdb_path1,
        pdb_path2=pdb_path2,
        outpath=outpath,
        bgcolor="white",
        show_type=show_type))
    script = "\n".join(lines)
    return run_pymol_script(script)


def _test_psurface(pdb_path):
    atoms = pymol_get_presenting_surface_atoms(pdb_path, "c", 10)
    for atom in atoms:
        print " ".join([str(x) for x in atom])


if __name__ == '__main__':
    import sys
    _test_psurface(sys.argv[1])
