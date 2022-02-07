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

    # days in a week
    days_in_a_week = 7 #days

    # Days of observations per week
    observation_days_per_week = 2 

    # Initial date of the simulation
    start_date = 2460004.5 #in Julian days (J2000) = 01/03/2023 00:00:00 # Two years later than March 2021 (taken from "LaRa after RISE: Expected improvement in the Mars rotation and interior models")

    # Duration of the simulation
    simulation_duration_days = 700 #days
    simulation_duration_weeks = simulation_duration_days/days_in_a_week #weeks
    simulation_duration = simulation_duration_days*constants.JULIAN_DAY #seconds

    # LaRa landing site, taken from "LaRa after RISE: Expected improvement in the Mars rotation and interior models"
    reflector_name = "LaRa"
    reflector_latitude_deg = 18.3 #North degrees
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
    simulation_start_epoch = (start_date-constants.JULIAN_DAY_ON_J2000)*constants.JULIAN_DAY #seconds
    simulation_end_epoch = simulation_start_epoch+simulation_duration #seconds
    print(simulation_start_epoch)

    # Define bodies in the simulation
    bodies_to_create = ["Earth","Mars"]


    global_frame_origin = "SSB" #Barycenter of Solar System
    global_frame_orientation = "ECLIPJ2000"
    body_settings = environment_setup.get_default_body_settings_time_limited(
        bodies_to_create,
        simulation_start_epoch-constants.JULIAN_DAY,
        simulation_end_epoch+constants.JULIAN_DAY,
        global_frame_origin,
        global_frame_orientation,
        60)

    # Mars rotation model 
    body_settings.get("Mars").rotation_model_settings = environment_setup.rotation_model.mars_high_accuracy()

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
    print(observation_start_epoch)
    # Define time between two observations
    observation_interval = 60 #seconds

    # Define observation simulation times for each link
    observation_times_list = list()
    for pointer_weeks in range(0,int(np.ceil(simulation_duration_weeks))):
        for pointer_days_per_week in range(0,int(observation_days_per_week)):
            for pointer_interval in range(0,int(np.ceil(constants.JULIAN_DAY/observation_interval))):
                observation_times_list.append(observation_start_epoch+pointer_weeks*days_in_a_week*constants.JULIAN_DAY \
                    +pointer_days_per_week*np.floor(days_in_a_week/observation_days_per_week)*constants.JULIAN_DAY \
                        +pointer_interval*observation_interval)

    # Specifications of the reflector
    reflector_station = bodies.get_body("Mars").get_ground_station(reflector_name)
    reflector_nominal_state_object = reflector_station.station_state
    reflector_pointing_angle_calculator_object = reflector_station.pointing_angles_calculator

    # Specifications of the transmitter
    transmitter_station = bodies.get_body("Earth").get_ground_station(transmitter_name)
    transmitter_nominal_state_object = transmitter_station.station_state
    transmitter_pointing_angle_calculator_object = transmitter_station.pointing_angles_calculator

    earth_elevation = list()
    earth_azimuth = list()
    DSS63_observation_time = list()
    DSS63_elevation = list()

    for pointer_time in observation_times_list:
        rotation_from_Mars_body_frame_to_inertial_frame = bodies.get_body("Mars").rotation_model.body_fixed_to_inertial_rotation(pointer_time)
        rotation_from_Earth_body_frame_to_inertial_frame = bodies.get_body("Earth").rotation_model.body_fixed_to_inertial_rotation(pointer_time)

        # Elevation and azimuth as seen by LaRa 
        earth_elevation.append(reflector_pointing_angle_calculator_object.calculate_elevation_angle(
            bodies.get_body("Earth").state_in_base_frame_from_ephemeris(pointer_time)[:3] \
                -np.matmul(rotation_from_Mars_body_frame_to_inertial_frame,reflector_nominal_state_object.get_cartesian_position(pointer_time)),
                pointer_time))
        earth_azimuth.append(reflector_pointing_angle_calculator_object.calculate_azimuth_angle(
            bodies.get_body("Earth").state_in_base_frame_from_ephemeris(pointer_time)[:3] \
                -np.matmul(rotation_from_Mars_body_frame_to_inertial_frame,reflector_nominal_state_object.get_cartesian_position(pointer_time)),
                pointer_time))

        # Angle viability
        if earth_elevation[-1] >= np.deg2rad(35) and earth_elevation[-1] <= np.deg2rad(45): 
            DSS63_observation_time.append(pointer_time)
            
            # As seen by DSS63
            DSS63_elevation.append(transmitter_pointing_angle_calculator_object.calculate_elevation_angle(
                bodies.get_body("Mars").state_in_base_frame_from_ephemeris(pointer_time)[:3] \
                    -np.matmul(rotation_from_Earth_body_frame_to_inertial_frame,transmitter_nominal_state_object.get_cartesian_position(pointer_time)),
                    pointer_time))

    ground_station_ids = list()
    ground_station_observation_time = list()
    ground_station_elevation = list()

    id = 0

    for Earth_ground_station_pointer in Earth_ground_station_list:
        if Earth_ground_station_pointer[1] == transmitter_name:
            continue

        current_station = bodies.get("Earth").get_ground_station(Earth_ground_station_pointer[1])
        current_nominal_state_object = current_station.station_state
        current_ground_station_pointing_angle_calculator_object = current_station.pointing_angles_calculator
        
        ind = 0

        # Angle viability
        for pointer_transmitter_time in DSS63_observation_time:
            if DSS63_elevation[ind] >= np.deg2rad(20):
                rotation_from_Earth_body_frame_to_inertial_frame = bodies.get_body("Earth").rotation_model.body_fixed_to_inertial_rotation(pointer_transmitter_time)
                
                # As seen by the ground station
                ground_station_observation_time.append(pointer_transmitter_time)
                ground_station_elevation.append(current_ground_station_pointing_angle_calculator_object.calculate_elevation_angle(
                    bodies.get_body("Mars").state_in_base_frame_from_ephemeris(pointer_transmitter_time)[:3] \
                        -np.matmul(rotation_from_Earth_body_frame_to_inertial_frame,current_nominal_state_object.get_cartesian_position(pointer_transmitter_time)),
                        pointer_transmitter_time))
                ground_station_ids.append(id)
            ind+=1
        id+=1

    ########################################################################################################################
    ################################################## PROVIDE OUTPUT TO CONSOLE AND FILES #################################
    ########################################################################################################################

    output_folder_path = os.path.dirname(os.path.realpath(__file__)).replace('/src','/output/GS')
    os.makedirs(output_folder_path,exist_ok=True)

    np.savetxt(output_folder_path+"/observation_time.dat",observation_times_list,fmt='%.15e')
    np.savetxt(output_folder_path+"/DSS63_observation_time.dat",DSS63_observation_time,fmt='%.15e')
    np.savetxt(output_folder_path+"/DSS63_elevation.dat",DSS63_elevation,fmt='%.15e')
    np.savetxt(output_folder_path+"/earth_elevation.dat",earth_elevation,fmt='%.15e')
    np.savetxt(output_folder_path+"/earth_azimuth.dat",earth_azimuth,fmt='%.15e')   
    np.savetxt(output_folder_path+"/ground_station_observation_time.dat",ground_station_observation_time,fmt='%.15e')
    np.savetxt(output_folder_path+"/ground_station_elevation.dat",ground_station_elevation,fmt='%.15e')
    np.savetxt(output_folder_path+"/ground_station_ids.dat",ground_station_ids,fmt='%.15e')                 

print("--- %s seconds ---" % (time.time() - run_time))