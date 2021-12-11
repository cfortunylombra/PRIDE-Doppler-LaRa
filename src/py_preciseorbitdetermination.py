"""
Description: Environment Setup for the Precise Orbit Determination

Author: C. Fortuny-Lombraña
"""

########################################################################################################################
################################################## IMPORT PACKAGES ##################################################
########################################################################################################################

import numpy as np
from tudatpy.kernel import constants
from tudatpy.kernel import numerical_simulation
from tudatpy.kernel.astro import element_conversion
from tudatpy.kernel.interface import spice_interface
from tudatpy.kernel.numerical_simulation import environment_setup
from tudatpy.kernel.numerical_simulation import propagation_setup
from tudatpy.kernel.numerical_simulation import propagation

spice_interface.load_standard_kernels()