import psi4
import re
import os
import inputparser
import math
import warnings
import driver
#from wrappers import *
from molutil import *
import p4util
#from psiexceptions import *

def run_mean_field_cqed(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    mean_field_cqed can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('mean_field_cqed')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    # Your plugin's psi4 run sequence goes here
    scf_wfn = driver.scf_helper(name, **kwargs)
    mean_field_cqed_wfn = psi4.plugin('mean_field_cqed.so',scf_wfn)
    #psi4.set_variable('CURRENT ENERGY', returnvalue)
    return mean_field_cqed_wfn


# Integration with driver routines
driver.procedures['energy']['mean_field_cqed'] = run_mean_field_cqed
driver.procedures['energy']['cqed'] = run_mean_field_cqed

def exampleFN():
    # Your Python code goes here
    pass
