"""
Description: Ground Station Analysis

Author: C. Fortuny-Lombraña
"""
import time
run_time = time.time()
if __name__=="__main__":
    ########################################################################################################################
    ################################################## IMPORT PACKAGES #####################################################
    ########################################################################################################################

    import sys
    sys.path.insert(0, "/home/cfortunylombra/tudat-bundle/cmake-build-release-wsl/tudatpy/")
    import os
    import copy
    import numpy as np
    from tudatpy.kernel import constants, numerical_simulation
    from tudatpy.kernel.astro import element_conversion, frame_conversion
    from tudatpy.kernel.interface import spice_interface
    from tudatpy.kernel.numerical_simulation import environment_setup,environment,propagation_setup,propagation,estimation_setup,estimation
    from tudatpy.kernel.numerical_simulation.estimation_setup import parameter,observation

    ########################################################################################################################
    ################################################## CONSTANTS AND VARIABLES #############################################
    ########################################################################################################################

    # J2000 epoch
    J2000_in_Julian_days = 2451545.0

    # days in a week
    days_in_a_week = 7 #days

    # Days of observations per week
    observation_days_per_week = 2 

    # Initial date of the simulation
    start_date = 2459215.5 #in Julian days (J2000) = 01/01/2021 00:00:00

    # Duration of the simulation
    simulation_duration_days = 700 #days
    simulation_duration_weeks = simulation_duration_days/days_in_a_week #weeks
    simulation_duration = simulation_duration_days*constants.JULIAN_DAY #seconds

    # LaRa landing site
    reflector_name = "LaRa"
    reflector_latitude_deg = 18.4 #North degrees
    reflector_longitude_deg = 335.37 #East degrees

    # Earth-based transmitter
    transmitter_name = "DSS63"
    transmitter_position_cartesian = np.array([4849092.6814,-360180.5350,4115109.1298]) #Taken from https://www.aoc.nrao.edu/software/sched/catalogs/locations.dat

    ########################################################################################################################
    ################################################## CREATE ENVIRONMENT ##################################################
    ########################################################################################################################

    # Load spice kernels
    spice_interface.load_standard_kernels()

    # Initial and end time of the simulation
    simulation_start_epoch = (start_date-J2000_in_Julian_days)*constants.JULIAN_DAY #seconds
    simulation_end_epoch = simulation_start_epoch+simulation_duration #seconds

    # Define bodies in the simulation
    bodies_to_create = ["Earth","Mars"]


    global_frame_origin = "SSB" #Barycenter of Solar System
    global_frame_orientation = "ECLIPJ2000"
    body_settings = environment_setup.get_default_body_settings_time_limited(
        bodies_to_create,
        simulation_start_epoch-constants.JULIAN_DAY,
        simulation_end_epoch+constants.JULIAN_DAY,
        global_frame_origin,
        global_frame_orientation)

    # Simple rotation model before moving to the more realistic Mars rotation model
    body_settings.get("Mars").rotation_model_settings = environment_setup.rotation_model.simple_from_spice("ECLIPJ2000","IAU_Mars","IAU_Mars",simulation_start_epoch)

    # Complex rotation model 
    #body_settings.get("Mars").rotation_model_settings = environment_setup.rotation_model.getHighAccuracyMarsRotationModel(simulation_start_epoch,simulation_end_epoch)

    bodies = environment_setup.create_system_of_bodies(body_settings)

    ########################################################################################################################
    ################################################## CREATE GROUND STATIONS AND LANDER ###################################
    ########################################################################################################################

    # Creation dictionary for ground stations
    ground_station_dict = {}

    # Adding transmitter 
    ground_station_dict [transmitter_name] = transmitter_position_cartesian

    # Read the text file containing the name and cartesian coordinates of the ground stations
    with open(os.path.dirname(os.path.realpath(__file__))+'/gs_locations.dat') as file:
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
            
    # Earth-based ground station creation
    for pointer_ground_station in range(0,len(ground_station_dict.keys())):
        environment_setup.add_ground_station(
            bodies.get_body("Earth"),
            list(ground_station_dict.keys())[pointer_ground_station],
            ground_station_dict[list(ground_station_dict.keys())[pointer_ground_station]])
    
    Earth_ground_station_list = environment_setup.get_ground_station_list(bodies.get_body("Earth"))

    # Mars-based ground station creation
    environment_setup.add_ground_station(
        bodies.get_body("Mars"),
        reflector_name,
        np.array([spice_interface.get_average_radius("Mars"),np.deg2rad(reflector_latitude_deg),np.deg2rad(reflector_longitude_deg)]),
         element_conversion.spherical_position_type)

    Mars_ground_station_list = environment_setup.get_ground_station_list(bodies.get_body("Mars"))

    ########################################################################################################################
    ################################################## GROUND STATIONS ELEVATION HISTORY ###################################
    ########################################################################################################################

    # Define time of first observation
    observation_start_epoch = simulation_start_epoch + constants.JULIAN_DAY
    
    # Define time between two observations
    observation_interval = 60 #seconds

    # Define observation simulation times for each link
    observation_times_list = list()
    for pointer_weeks in range(0,int(simulation_duration_weeks)):
        for pointer_days_per_week in range(0,int(observation_days_per_week)):
            for pointer_interval in range(0,int(constants.JULIAN_DAY/observation_interval)):
                observation_times_list.append(observation_start_epoch+pointer_weeks*days_in_a_week*constants.JULIAN_DAY \
                    +pointer_days_per_week*(days_in_a_week/2)*constants.JULIAN_DAY \
                        +pointer_interval*observation_interval)

    # Specifications of the reflector
    reflector_station = bodies.get_body("Mars").get_ground_station(reflector_name)
    reflector_nominal_state_object = reflector_station.station_state
    reflector_pointing_angle_calculator_object = reflector_station.pointing_angles_calculator
    rotation_from_Mars_body_frame_to_inertial_frame = bodies.get_body("Mars").body_fixed_to_inertial_frame

    # Specifications of the transmitter
    transmitter_station = bodies.get_body("Earth").get_ground_station(transmitter_name)
    transmitter_nominal_state_object = transmitter_station.station_state
    transmitter_pointing_angle_calculator_object = transmitter_station.pointing_angles_calculator
    rotation_from_Earth_body_frame_to_inertial_frame = bodies.get_body("Earth").body_fixed_to_inertial_frame

    ########################################################################################################################
    ################################################## PROVIDE OUTPUT TO CONSOLE AND FILES #################################
    ########################################################################################################################

    output_folder_path = os.path.dirname(os.path.realpath(__file__)).replace('/src','/output')
    os.makedirs(output_folder_path,exist_ok=True)

    np.savetxt(output_folder_path+"/observation_time.dat",observation_time,fmt='%.16e')
    np.savetxt(output_folder_path+"/DSS63_observation_time.dat",DSS63_observation_time,fmt='%.16e')
    np.savetxt(output_folder_path+"/DSS63_elevation.dat",DSS63_elevation,fmt='%.16e')
    np.savetxt(output_folder_path+"/earth_elevation.dat",earth_elevation,fmt='%.16e')
    np.savetxt(output_folder_path+"/earth_azimuth.dat",earth_azimuth,fmt='%.16e')   
    np.savetxt(output_folder_path+"/ground_station_observation_time.dat",ground_station_observation_time,fmt='%.16e')
    np.savetxt(output_folder_path+"/ground_station_elevation.dat",ground_station_elevation,fmt='%.16e')
    np.savetxt(output_folder_path+"/ground_station_ids.dat",ground_station_ids,fmt='%.16e')         

print("--- %s seconds ---" % (time.time() - run_time))