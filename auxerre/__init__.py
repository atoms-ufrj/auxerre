"""
Auxerre
Transport Properties from Molecular Dynamics Using Reciprocal-Space Correlation Analysis
"""

# Add imports here
from .auxerre import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
