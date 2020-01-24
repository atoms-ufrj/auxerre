"""
Unit and regression test for the auxerre package.
"""

# Import package, test suite, and other packages as needed
import auxerre
import pytest
import sys

def test_auxerre_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "auxerre" in sys.modules
