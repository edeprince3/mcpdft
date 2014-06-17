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


def run_tdhf_cqed(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    tdhf_cqed can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('tdhf_cqed')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    # Your plugin's psi4 run sequence goes here
    scf_helper(name, **kwargs)
    returnvalue = psi4.plugin('tdhf_cqed.so')
    psi4.set_variable('CURRENT ENERGY', returnvalue)


# Integration with driver routines
procedures['energy']['tdhf_cqed'] = run_tdhf_cqed
procedures['energy']['cqed'] = run_tdhf_cqed


def exampleFN():
    # Your Python code goes here
    pass
