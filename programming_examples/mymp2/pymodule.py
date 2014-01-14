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


def run_mymp2(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    mymp2 can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('mymp2')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    # Your plugin's psi4 run sequence goes here
    scf_helper(name, **kwargs)
    returnvalue = psi4.plugin('mymp2.so')
    psi4.set_variable('CURRENT ENERGY', returnvalue)


# Integration with driver routines
procedures['energy']['mymp2'] = run_mymp2


def exampleFN():
    # Your Python code goes here
    pass
