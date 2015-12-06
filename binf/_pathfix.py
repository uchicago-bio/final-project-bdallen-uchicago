"""
Hack to get scripts to run from source checkout without having to set
PYTHONPATH.
"""

import sys
from os.path import dirname, join, abspath

bin_path = dirname(__file__)
project_path = abspath(join(bin_path, ".."))
sys.path.insert(0, project_path)

# Data is stored in a separate directory along side the project
# directory; too large to put into source control.
data_path = abspath(join(project_path, "..", "project-data"))
