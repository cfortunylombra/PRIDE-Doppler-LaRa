"""
Description: Environment Setup for the Precise Orbit Determination

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
    from tudatpy.kernel.astro import element_conversion
    from tudatpy.kernel.interface import spice_interface
    from tudatpy.kernel.numerical_simulation import environment_setup,propagation_setup,propagation,estimation_setup,estimation
    from tudatpy.kernel.numerical_simulation.estimation_setup import observation

    np.set_printoptions(suppress=True,precision=15)
    ########################################################################################################################
    ################################################## CONSTANTS AND VARIABLES #############################################
    ########################################################################################################################

    # days in a week
    days_in_a_week = np.double(7) #days

    # Days of observations per week
    observation_days_per_week = np.double(2)

    # Initial date of the simulation
    start_date = np.double(2459215.5) #in Julian days (J2000) = 01/01/2021 00:00:00

    # Duration of the simulation
    simulation_duration_days = 700 #days #NOTE
    simulation_duration_weeks = simulation_duration_days/days_in_a_week #weeks
    simulation_duration = simulation_duration_days*constants.JULIAN_DAY #seconds

    # LaRa landing site
    reflector_name = "LaRa"
    reflector_latitude_deg = 18.20 #18.4 #North degrees
    reflector_longitude_deg = 335.45 #335.37 #East degrees

    # Earth-based transmitter
    transmitter_name = "DSS63"
    transmitter_position_cartesian = np.array([4849092.6814,-360180.5350,4115109.1298], dtype='d') #Taken from https://www.aoc.nrao.edu/software/sched/catalogs/locations.dat

    ########################################################################################################################
    ################################################## CREATE ENVIRONMENT ##################################################
    ########################################################################################################################

    # Load spice kernels
    spice_interface.load_standard_kernels()

    # Initial and end time of the simulation
    simulation_start_epoch = (start_date-constants.JULIAN_DAY_ON_J2000)*constants.JULIAN_DAY #seconds
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
        global_frame_orientation,
        time_step = 60)

    # Reset frame origin
    #environment_setup.ephemeris.frame_origin = "Sun"

    #Mars rotation model
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

            ground_station_dict[name_ground_station] = np.array([x_coordinate_ground_station,y_coordinate_ground_station,z_coordinate_ground_station],dtype='d')

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
        np.array([spice_interface.get_average_radius("Mars"),np.deg2rad(reflector_latitude_deg,dtype='d'),np.deg2rad(reflector_longitude_deg,dtype='d')]),
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
    maximum_step_size = 60.0 #seconds
    relative_error_tolerance = 1.0E-14
    absolute_error_tolerance = 1.0E-14

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
            two_way_link_ends[observation.transmitter] = ("Earth", transmitter_name )
            two_way_link_ends[observation.reflector1] = ( "Mars", reflector_name )
            two_way_link_ends[observation.receiver] = ( "Earth", receiver_name )

            observation_settings_list.append(two_way_link_ends)
    
    # Create the uplink list
    observation_settings_uplink_list = copy.deepcopy(observation_settings_list)
    observation_settings_uplink_list = list([observation_settings_uplink_list[0]])

    # Copy the entire list of dictionaries for downlink
    observation_settings_downlink_list = copy.deepcopy(observation_settings_list)

    # Remove receiver for uplink, and rename the reflector to receiver
    del observation_settings_uplink_list[0][observation.receiver]
    observation_settings_uplink_list[0][observation.receiver] =  observation_settings_uplink_list[
        0].pop(observation.reflector1)

    for pointer_link_ends in range(0,len(observation_settings_list)):
        # Remove transmitter for downlink, and rename the reflector to transmitter
        del observation_settings_downlink_list[pointer_link_ends][observation.transmitter]
        observation_settings_downlink_list[pointer_link_ends][observation.transmitter] =  observation_settings_downlink_list[
            pointer_link_ends].pop(observation.reflector1)

    ########################################################################################################################
    ################################################## DEFINE PARAMETERS TO ESTIMATE #######################################
    ######################################################################################################################## 

    # Create list of parameters that are to be estimated
    parameter_settings = estimation_setup.parameter.initial_states(propagator_settings,bodies)
    parameter_settings.append(estimation_setup.parameter.core_factor("Mars"))
    parameter_settings.append(estimation_setup.parameter.free_core_nutation_rate("Mars"))
    parameter_settings.append(estimation_setup.parameter.ground_station_position("Mars", reflector_name))
    parameter_settings.append(estimation_setup.parameter.periodic_spin_variations("Mars"))
    parameter_settings.append(estimation_setup.parameter.polar_motion_amplitudes("Mars"))

    parameters_set = estimation_setup.create_parameter_set(parameter_settings,bodies,propagator_settings)

    truth_parameter = parameters_set.parameter_vector
    # Print identifiers and indices of parameters to terminal
    estimation_setup.print_parameter_names(parameters_set)
    print('Initial parameter estimate is: ')
    print(truth_parameter) 

    ########################################################################################################################
    ################################################## CREATE OBSERVATION SETTINGS #########################################
    ######################################################################################################################## 
    
    # Define twoway Doppler observation settings
    two_way_doppler_observation_settings = list()
    for pointer_link_ends in range(0,len(observation_settings_list)):
        two_way_doppler_observation_settings.append(observation.two_way_open_loop_doppler(
            observation_settings_list[pointer_link_ends]))

    ########################################################################################################################
    ################################################## INITIALIZE OD  ######################################################
    ######################################################################################################################## 

    # Create observation simulators
    observation_simulators = estimation_setup.create_observation_simulators(two_way_doppler_observation_settings,bodies) 

    # Create physical environment (as set of physical bodies)
    estimator = numerical_simulation.Estimator(bodies,parameters_set,two_way_doppler_observation_settings,
        integrator_settings,propagator_settings)

    # Variational equations and dynamics
    #variational_equations_simulator = estimator.variational_solver
    #dynamics_simulator = variational_equations_simulator.dynamics_simulator

    # Extract observation simulators
    observation_simulators = estimator.observation_simulators

    ########################################################################################################################
    ################################################## SIMULATE OBSERVATIONS ###############################################
    ########################################################################################################################

    # Define time of first observation
    observation_start_epoch = simulation_start_epoch + constants.JULIAN_DAY
    
    # Define time between two observations
    observation_interval = 60*60 #seconds #NOTE

    # Define observation simulation times for each link
    observation_times_list = list()
    for pointer_weeks in range(0,int(np.ceil(simulation_duration_weeks))):
        for pointer_interval in range(0,int(np.ceil(constants.JULIAN_DAY/observation_interval))):
            for pointer_days_per_week in range(0,int(observation_days_per_week)):
                observation_times_list.append(observation_start_epoch+pointer_weeks*days_in_a_week*constants.JULIAN_DAY \
                    +pointer_days_per_week*3.25*constants.JULIAN_DAY \
                        +pointer_interval*observation_interval)

    # Create observation viability settings and calculators
    #viability_settings_list = list()
    #viability_settings_list.append(observation.elevation_angle_viability(["Earth",""],np.deg2rad(20)))
    #viability_settings_list.append(observation.elevation_angle_viability(["Mars",""],np.deg2rad(35)))
    #viability_settings_list.append(observations.elevation_angle_viability(["Mars",""],np.deg2rad(45))) #NOTE maximum elevation angle viability
    #viability_settings_list.append(observation.body_avoidance_viability(["Earth",""],"Sun",np.deg2rad(20)))
    #viability_settings_list.append(observation.body_occultation_viability(("Earth",""),"Moon"))

    # Create measurement simulation input
    observation_simulation_settings = list()
    observation_simulation_settings.extend(observation.tabulated_simulation_settings_list(
        dict({observation.two_way_doppler_type:observation_settings_list}),observation_times_list,reference_link_end_type = observation.transmitter))
    #for link_end_pointer in range(0,len(observation_settings_list)):
    #    observation_simulation_settings.extend(observation.tabulated_simulation_settings_list(
    #        dict({observation.two_way_doppler_type:[observation_settings_list[link_end_pointer]]}),observation_times_list,
            #viability_settings = viability_settings_list,
    #         reference_link_end_type = observation.transmitter))
    
    #observation.add_viability_check_to_settings(observation_simulation_settings,viability_settings_list) 

    # Simulate required observation
    simulated_observations = estimation.simulate_observations(observation_simulation_settings, observation_simulators, bodies)

    # Perturbation
    parameter_perturbation = np.zeros(parameters_set.parameter_set_size) #parameters_set.parameter_vector*0.01
    parameter_perturbation[0:3]=1000*np.ones(3)
    parameter_perturbation[3:6]=10*np.ones(3)
    parameter_perturbation[6]=1*10**(-2)
    parameter_perturbation[7]=1*10**(-7)
    parameter_perturbation[8:11]=100*np.ones(3)
    parameter_perturbation[11:]=1*10**(-9)*np.ones(28)

    print("Perturbation vector is:")
    print(parameter_perturbation)

    # Perturb estimate
    initial_parameter_deviation = parameter_perturbation

    # Define a priori covariance
    #inverse_a_priori_covariance = np.zeros((parameters_set.parameter_set_size,parameters_set.parameter_set_size))

    # Estimate parameters
    pod_input = estimation.PodInput(simulated_observations,parameters_set.parameter_set_size, apriori_parameter_correction = initial_parameter_deviation)
    #pod_input.define_estimation_settings(reintegrate_variational_equations = False)

    # Define noise levels
    doppler_noise = 0.05e-3/constants.SPEED_OF_LIGHT 
    weights_per_observable = dict({observation.two_way_doppler_type:doppler_noise**(-2)}) #observations.two_way_doppler_type NOTE
    pod_input.set_constant_weight_per_observable(weights_per_observable)

    # Perform estimation
    pod_output = estimator.perform_estimation(pod_input)

    # Create noise functions
    #observations.add_gaussian_noise_to_settings(observation_simulation_settings,doppler_noise,observations.one_way_doppler_type) #observations.two_way_doppler_type NOTE

    ########################################################################################################################
    ################################################## PROVIDE OUTPUT TO CONSOLE AND FILES #################################
    ########################################################################################################################

    output_folder_path = os.path.dirname(os.path.realpath(__file__)).replace('/src','/output/POD')
    os.makedirs(output_folder_path,exist_ok=True)

    estimation_error = np.subtract(pod_output.parameter_estimate,truth_parameter)
    print('True estimation error is:')
    print(estimation_error.astype('d'))
    formal_error = pod_output.formal_errors
    print('Formal estimation error is:')
    print(formal_error.astype('d'))
    true_to_form_estimation_error_ratio = formal_error/estimation_error 
    print('True to form estimation error is:')
    print(true_to_form_estimation_error_ratio)
    estimation_information_matrix = pod_output.normalized_design_matrix
    estimation_information_matrix_normalization = pod_output.normalization_terms
    concatenated_times = simulated_observations.concatenated_times 
    concatenated_observations = simulated_observations.concatenated_observations

    np.savetxt(output_folder_path+"/estimation_information_matrix.dat",estimation_information_matrix,fmt='%.15e')
    np.savetxt(output_folder_path+"/estimation_information_matrix_normalization.dat",
        estimation_information_matrix_normalization,fmt='%.15e')
    np.savetxt(output_folder_path+"/concatenated_times.dat",concatenated_times,fmt='%.15e')
    np.savetxt(output_folder_path+"/concatenated_observations.dat",concatenated_observations,fmt='%.15e')

print("--- %s seconds ---" % (time.time() - run_time))