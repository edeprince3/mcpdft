import psi4
import re
import os
import inputparser
import math
import warnings
from driver import *
from wrappers import *
from molutil import *
import p4util
from psiexceptions import *


def run_bpsdp(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    bpsdp can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('bpsdp')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    # Your plugin's psi4 run sequence goes here
    scf_helper(name, **kwargs)
    returnvalue = psi4.plugin('bpsdp.so')

    return psi4.get_variable('CURRENT ENERGY')


# Integration with driver routines
procedures['energy']['bpsdp'] = run_bpsdp


def exampleFN():
    # Your Python code goes here
    pass
