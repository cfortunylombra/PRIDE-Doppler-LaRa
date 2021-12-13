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
    from tudatpy.kernel.astro import element_conversion
    from tudatpy.kernel.interface import spice_interface
    from tudatpy.kernel.numerical_simulation import environment_setup,propagation_setup,propagation
    from tudatpy.kernel.numerical_simulation.estimation_setup import observations

    ########################################################################################################################
    ################################################## CONSTANTS AND VARIABLES #############################################
    ########################################################################################################################

    # J2000 epoch
    J2000_in_Julian_days = 2451545.0

    # Initial date of the simulation
    start_date = 2459215.5 #in Julian days (J2000) = 01/01/2021 00:00:00

    # Duration of the simulation
    simulation_duration = 700*constants.JULIAN_DAY #seconds

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
    bodies_to_create = ["Saturn","Jupiter","Mars","Moon","Earth","Venus","Mercury","Sun"]


    global_frame_origin = "SSB" #Barycenter of Solar System
    global_frame_orientation = "ECLIPJ2000"
    body_settings = environment_setup.get_default_body_settings_time_limited(
        bodies_to_create,
        simulation_start_epoch-constants.JULIAN_DAY,
        simulation_end_epoch+constants.JULIAN_DAY,
        global_frame_origin,
        global_frame_orientation)

    # Reset frame origin
    environment_setup.ephemeris.frame_origin = "Sun"

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
    ################################################## CREATE ACCELERATION MODELS ##########################################
    ########################################################################################################################

    # Define accelelerations
    accelerations_settings_Mars = dict(
        Saturn = [propagation_setup.acceleration.point_mass_gravity()],
        Jupiter = [propagation_setup.acceleration.point_mass_gravity()],
        Earth = [propagation_setup.acceleration.point_mass_gravity()],
        Venus = [propagation_setup.acceleration.point_mass_gravity()],
        Mercury = [propagation_setup.acceleration.point_mass_gravity()],
        Sun = [propagation_setup.acceleration.point_mass_gravity()],
    )

    acceleration_settings = {"Mars": accelerations_settings_Mars}

    # Define list of bodies to propagate
    bodies_to_propagate = ["Mars"]

    # Define central bodies to use in the propagation
    central_bodies = ["SSB"]

    # Create acceleration models and propagation settings
    acceleration_models = propagation_setup.create_acceleration_models(bodies,acceleration_settings,bodies_to_propagate,central_bodies)

    ########################################################################################################################
    ################################################## CREATE PROPAGATION SETTINGS #########################################
    ########################################################################################################################   

    # Define intitial state system state
    initial_state = propagation.get_state_of_bodies(bodies_to_propagate,central_bodies,bodies,simulation_start_epoch)

    # Define termination condition
    termination_condition = propagation_setup.propagator.time_termination(simulation_end_epoch)

    # Define propagator settings
    propagator_settings = propagation_setup.propagator.translational(central_bodies,acceleration_models,bodies_to_propagate,initial_state,termination_condition)

    # Define integrator settings
    initial_time_step = 1 #second
    minimum_step_size = initial_time_step #seconds
    maximum_step_size = 60 #seconds
    relative_error_tolerance = 1E-14
    absolute_error_tolerance = 1E-14

    integrator_settings = propagation_setup.integrator.runge_kutta_variable_step_size(
        simulation_start_epoch,
        initial_time_step,
        propagation_setup.integrator.rkf_78,
        minimum_step_size,
        maximum_step_size,
        relative_error_tolerance,
        absolute_error_tolerance)

    ########################################################################################################################
    ################################################## DEFINE LINK ENDS FOR OBSERVATIONS ###################################
    ######################################################################################################################## 

    # Create list of observation settings
    observation_settings_list = list()

    # Define link ends
    for pointer_ground_station_receiver in range(0,len(ground_station_dict.keys())):
        receiver_name = list(ground_station_dict.keys())[pointer_ground_station_receiver]
        if receiver_name != transmitter_name:
            two_way_link_ends = dict()
            two_way_link_ends[observations.transmitter] = ("Earth", transmitter_name )
            two_way_link_ends[observations.reflector1] = ( "Mars", reflector_name );
            two_way_link_ends[observations.receiver] = ( "Earth", receiver_name );

            observation_settings_list.append(two_way_link_ends)