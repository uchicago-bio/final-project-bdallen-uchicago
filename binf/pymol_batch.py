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
