"""
Description: Environment Setup for the Precise Orbit Determination

Author: C. Fortuny-Lombraña
"""

if __name__=="__main__":
    ########################################################################################################################
    ################################################## IMPORT PACKAGES #####################################################
    ########################################################################################################################

    import os
    import numpy as np
    from tudatpy.kernel import constants
    from tudatpy.kernel.interface import spice_interface
    from tudatpy.kernel.numerical_simulation import environment_setup

    #import tudatpy as tudat
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
    bodies_to_create = ["Saturn","Jupiter","Mars","Moon","Earth","Venus","Mercury","Sun"]


    global_frame_origin = "SSB" #Barycenter of Solar System
    global_frame_orientation = "ECLIPJ2000"
    body_settings = environment_setup.get_default_body_settings(bodies_to_create, global_frame_origin, global_frame_orientation )

    #body_settings=environment_setup.get_default_body_settings_time_limited(bodies_to_create, initial_time, final_time, global_frame_origin, global_frame_orientation, time_step)

    #In Documentation and does not work
    #body_settings.get( "Sun" ).gravity_field_settings = environment_setup.gravity_field.point_mass( 1.32712440042E20 )
    #https://tudat-space.readthedocs.io/en/latest/_src_user_guide/state_propagation/environment_setup/create_bodies/create_body_settings.html#create-celestial-body-settings

    #Functions not defined
    #v1
    #https://github.com/tudat-team/tudatpy/blob/master/tudatpy/kernel/expose_simulation/expose_environment_setup.cpp
    #body_settings.get("Moon").ephemeris_settings = environment_setup.ephemeris.reset_frame_origin("Sun")
    #body_settings.get("Mars").rotation_model_settings = environment_setup.rotation_model.getHighAccuracyMarsRotationModel(simulation_start_epoch,simulation_end_epoch)
    #v2
    #https://tudat-space.readthedocs.io/en/latest/_src_user_guide/state_propagation/environment_setup/create_models/available.html?highlight=simple%20rotational%20model#rotational
    #bodySettings[ "Mars" ]->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >("ECLIPJ2000", "IAU_Mars", spice_interface::computeRotationQuaternionBetweenFrames("ECLIPJ2000", "IAU_Mars", initialEphemerisTime ),initialEphemerisTime, 2.0 * mathematical_constants::PI /( physical_constants::JULIAN_DAY + 40.0 * 60.0 ) );
    #body_settings.get("Mars").rotation_model_settings = environment_setup.rotation_model.simple("ECLIPJ2000", "IAU_Mars",spice_interface.compute_rotation_quaternion_between_frames("ECLIPJ2000", "IAU_Mars", simulation_start_epoch),simulation_end_epoch, 2*math.pi/(constants.JULIAN_DAY+40*60))

    bodies = environment_setup.create_system_of_bodies(body_settings)

    ########################################################################################################################
    ################################################## CREATE GROUND STATIONS AND LANDER ###################################
    ########################################################################################################################

    # Creation Dictionary for Ground Stations
    ground_station_dict = {}

    # Read the text file containing the name and cartesian coordinates of the ground stations
    with open(os.path.dirname(os.path.realpath(__file__))+'\gs_locations.dat') as file:
        lines = file.read().splitlines()
        
        # Variables
        skiplines = 29 #lines to be removed from the description at the beginning of the text file
        eachgroundstationlines = 6 #lines of specs that contains each ground stations

        lines = lines[skiplines:]
        number_ground_stations_file = int(len(lines)/eachgroundstationlines) #total number of ground stations 

        for pointer_ground_station in range(0,number_ground_stations_file):
            name_line_ground_station = lines[pointer_ground_station*eachgroundstationlines+1]
            coordinates_line_ground_station = lines[pointer_ground_station*eachgroundstationlines+2]
            
            if len(name_line_ground_station.split("DBNAME=",1)) == 2:
                name_ground_station = name_line_ground_station.split("DBNAME=",1)[1].split()[0]
            elif len(name_line_ground_station.split("DBCODE=",1)) == 2:
                name_ground_station = name_line_ground_station.split("DBCODE=",1)[1].split()[0]
            x_coordinate_ground_station = float(coordinates_line_ground_station.split("X=",1)[1].split()[0])
            y_coordinate_ground_station = float(coordinates_line_ground_station.split("Y=",1)[1].split()[0])
            z_coordinate_ground_station = float(coordinates_line_ground_station.split("Z=",1)[1].split()[0])

            ground_station_dict[name_ground_station] = np.array([x_coordinate_ground_station,y_coordinate_ground_station,z_coordinate_ground_station])
            
    # Ground Station Creation
    for pointer_ground_station in range(0,len(ground_station_dict.keys())):
        environment_setup.add_ground_station(bodies.get_body("Earth"),list(ground_station_dict.keys())[pointer_ground_station],ground_station_dict[list(ground_station_dict.keys())[pointer_ground_station]])
    print(environment_setup.get_ground_station_list(bodies.get_body("Earth")))