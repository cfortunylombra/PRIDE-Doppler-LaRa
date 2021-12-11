"""
Description: Environment Setup for the Precise Orbit Determination

Author: C. Fortuny-Lombraña
"""

########################################################################################################################
################################################## IMPORT PACKAGES #####################################################
########################################################################################################################

import numpy as np
from tudatpy.kernel import constants
from tudatpy.kernel.interface import spice_interface
from tudatpy.kernel.numerical_simulation import environment_setup


#from tudatpy.kernel import numerical_simulation
#from tudatpy.kernel.astro import element_conversion
#from tudatpy.kernel.numerical_simulation import propagation_setup
#from tudatpy.kernel.numerical_simulation import propagation

########################################################################################################################
################################################## CONSTANTS AND VARIABLES #############################################
########################################################################################################################

# J2000 epoch
J2000_in_Julian_days = 2451545.0

# Initial date of the simulation
start_date = 2459215.5 #in Julian days (J2000) = 01/01/2021 00:00:00

# Duration of the simulation
simulation_duration = 700*constants.JULIAN_DAY #seconds

########################################################################################################################
################################################## CREATE ENVIRONMENT ##################################################
########################################################################################################################

# Load spice kernels
spice_interface.load_standard_kernels()

# Initial and end time of the simulation

simulation_start_epoch = (start_date-J2000_in_Julian_days)*constants.JULIAN_DAY #seconds
simulation_end_epoch = simulation_start_epoch+simulation_duration #seconds

# Define bodies in the simulation
bodies_to_create = ["Saturn","Jupyter","Mars","Moon","Earth","Venus","Mercury","Sun"]


global_frame_origin = "SSB" #Barycenter of Solar System
global_frame_orientation = "ECLIPJ2000"
body_settings = environment_setup.get_default_body_settings(
    bodies_to_create, global_frame_origin, global_frame_orientation )

bodies = environment_setup.create_system_of_bodies(body_settings)