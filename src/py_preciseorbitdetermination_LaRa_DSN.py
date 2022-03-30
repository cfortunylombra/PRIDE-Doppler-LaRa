"""
Description: Environment Setup for the Precise Orbit Determination (LaRa)

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
    import datetime
    import numpy as np
    import matplotlib.pyplot as plt
    from tudatpy.kernel import constants, numerical_simulation
    from tudatpy.kernel.astro import element_conversion,time_conversion
    from tudatpy.kernel.interface import spice_interface
    from tudatpy.kernel.numerical_simulation import environment_setup,propagation_setup,propagation,estimation_setup,estimation
    from tudatpy.kernel.numerical_simulation.estimation_setup import observation

    np.set_printoptions(suppress=True,precision=15)
    ########################################################################################################################
    ################################################## CONSTANTS AND VARIABLES #############################################
    ########################################################################################################################

    # days in a week
    days_in_a_week = 7 #days

    # Initial date of the simulation
    #start_date = 2460096.5 #in Julian days = 01/06/2023 00:00:00 # Two years later than March 2021 (taken from "LaRa after RISE: Expected improvement in the Mars rotation and interior models")
    start_date = 2459580.5 #in Julian days = 01/01/2022
    
    # Duration of the simulation
    simulation_duration_days = 1000 #days
    simulation_duration_weeks = simulation_duration_days/days_in_a_week #weeks
    simulation_duration = simulation_duration_days*constants.JULIAN_DAY #seconds

    # LaRa landing site, taken from "LaRa after RISE: Expected improvement in the Mars rotation and interior models"
    reflector_name = "LaRa"
    reflector_latitude_deg = 18.3 #North degrees
    reflector_longitude_deg = 335.37 #East degrees
    reflector_altitude = -2000 #m

    # Earth-based transmitter
    transmitter_names = ['DSS 43','DSS 63','DSS 14']

    transmitter_positions_cartesian = list()  #Taken from JPL web site
    transmitter_positions_cartesian.append(np.array([-4460894.9170,2682361.5070,-3674748.1517])) # DSS 43
    transmitter_positions_cartesian.append(np.array([4849092.5175,-360180.3480,4115109.2506])) # DSS 63
    transmitter_positions_cartesian.append(np.array([-2353621.4197,-4641341.4717,3677052.3178])) # DSS 14

    # Viability settings
    earth_min = 35 #deg
    earth_max = 45 #deg
    antenna_min_elevation = 10 #deg
    body_avoidance_angle = 10 #deg

    ########################################################################################################################
    ################################################## CREATE ENVIRONMENT ##################################################
    ########################################################################################################################

    # Load spice kernels
    spice_interface.load_standard_kernels()

    # Initial and end time of the simulation
    simulation_start_epoch = (start_date-constants.JULIAN_DAY_ON_J2000)*constants.JULIAN_DAY #seconds
    simulation_end_epoch = simulation_start_epoch+simulation_duration #seconds

    # Define bodies in the simulation
    bodies_to_create = ["Saturn","Jupiter","Mars","Moon","Earth","Venus","Mercury","Sun"] #Phobos Deimos

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
    environment_setup.ephemeris.frame_origin = "Sun"

    #Mars rotation model
    body_settings.get("Mars").rotation_model_settings = environment_setup.rotation_model.mars_high_accuracy()

    bodies = environment_setup.create_system_of_bodies(body_settings)

    ########################################################################################################################
    ################################################## CREATE GROUND STATIONS AND LANDER ###################################
    ########################################################################################################################

    # Creation dictionary for ground stations
    ground_station_dict = {}

    # Adding transmitter 
    for transmitter_index in range(0,len(transmitter_names)):
        ground_station_dict[transmitter_names[transmitter_index]] = transmitter_positions_cartesian[transmitter_index]

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
        np.array([reflector_altitude,np.deg2rad(reflector_latitude_deg,dtype='d'),np.deg2rad(reflector_longitude_deg,dtype='d')]),
         element_conversion.geodetic_position_type)

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
    maximum_step_size = 60*60*24 #seconds
    relative_error_tolerance = 1.0E-12
    absolute_error_tolerance = 1.0E-12

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
        two_way_link_ends = dict()
        two_way_link_ends[observation.transmitter] = ("Earth", receiver_name )
        two_way_link_ends[observation.reflector1] = ( "Mars", reflector_name )
        two_way_link_ends[observation.receiver] = ( "Earth", receiver_name )

        observation_settings_list.append(two_way_link_ends)

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
    variational_equations_simulator = estimator.variational_solver
    dynamics_simulator = variational_equations_simulator.dynamics_simulator

    # Extract observation simulators
    observation_simulators = estimator.observation_simulators

    ########################################################################################################################
    ################################################## SIMULATE OBSERVATIONS ###############################################
    ########################################################################################################################
    
    # Define time between two observations
    observation_interval = 60 #seconds 

    # Define observation simulation times for each link
    observation_times_list = list()
    
    # Read the epoch times
    with open(os.path.dirname(os.path.realpath(__file__))+'/lsor_export_forTUDelft.txt') as file:
        lines = file.read().splitlines()
        observation_start_epoch = np.inf
        for pointer_pass in range(1,len(lines)):
            line = lines[pointer_pass]
            line_info = line.split()

            startpass_year = int(line_info[0].split('-')[0])
            startpass_day_of_year = int(line_info[0].split('-')[1].split('T')[0])
            startpass_hour = int(line_info[0].split('-')[1].split('T')[1].split(':')[0])
            startpass_min = int(line_info[0].split('-')[1].split('T')[1].split(':')[1])
            startpass_sec = int(line_info[0].split('-')[1].split('T')[1].split(':')[2])
            startpass_date = datetime.datetime(startpass_year,1,1)+datetime.timedelta(days=startpass_day_of_year-1,hours=startpass_hour,minutes=startpass_min,seconds=startpass_sec)
            startpass_epoch = (startpass_date - datetime.datetime(2000,1,1,12,0,0,0)).total_seconds()
            
            endpass_year = int(line_info[1].split('-')[0])
            endpass_day_of_year = int(line_info[1].split('-')[1].split('T')[0])
            endpass_hour = int(line_info[1].split('-')[1].split('T')[1].split(':')[0])
            endpass_min = int(line_info[1].split('-')[1].split('T')[1].split(':')[1])
            endpass_sec = int(line_info[1].split('-')[1].split('T')[1].split(':')[2])
            endpass_date = datetime.datetime(startpass_year,1,1)+datetime.timedelta(days=endpass_day_of_year-1,hours=endpass_hour,minutes=endpass_min,seconds=endpass_sec)
            endpass_epoch = (endpass_date - datetime.datetime(2000,1,1,12,0,0,0)).total_seconds()

            observation_times_list.extend(np.arange(startpass_epoch,endpass_epoch+observation_interval,observation_interval))

    # Create observation viability settings and calculators
    viability_settings_list = list()
    viability_settings_list.append(observation.elevation_angle_viability(["Earth",""],np.deg2rad(antenna_min_elevation)))
    viability_settings_list.append(observation.elevation_angle_viability(["Mars",""],np.deg2rad(earth_min)))
    #viability_settings_list.append(observations.elevation_angle_viability(["Mars",""],np.deg2rad(earth_max)))
    viability_settings_list.append(observation.body_avoidance_viability(["Earth",""],"Sun",np.deg2rad(body_avoidance_angle)))
    viability_settings_list.append(observation.body_occultation_viability(("Earth",""),"Moon"))

    # Create measurement simulation input
    observation_simulation_settings = observation.tabulated_simulation_settings_list(
        dict({observation.two_way_doppler_type:observation_settings_list}),observation_times_list,
        viability_settings = viability_settings_list,reference_link_end_type = observation.receiver)

    # Define noise levels
    doppler_noise = 0.05e-3/constants.SPEED_OF_LIGHT_LONG # Taken from the Radioscience LaRa instrument onboard ExoMars to investigate the rotation and interior of Mars
    weights_per_observable = dict({observation.two_way_doppler_type:doppler_noise**(-2)})

    # Create noise functions
    observation.add_gaussian_noise_to_settings(observation_simulation_settings,doppler_noise,observation.two_way_doppler_type)
    
    # Simulate required observation
    simulated_observations = estimation.simulate_observations(observation_simulation_settings, observation_simulators, bodies)

    # Perturbation
    parameter_perturbation = np.zeros(parameters_set.parameter_set_size) 
    mas =np.pi/(180.0*1000.0*3600.0) # Conversion from milli arc seconds to seconds 
    # Position of Mars
    parameter_perturbation[0:3]=1000*np.ones(3) # meters; Taken from Improving the Accuracy of the Martian Ephemeris Short-Term Prediction
    # Velocity of Mars
    parameter_perturbation[3:6]=0.0002*np.ones(3) #meters; Taken from Improving the Accuracy of the Martian Ephemeris Short-Term Prediction
    # Core factor of the celestial body of Mars
    parameter_perturbation[6]=0.014 # Unitless; Taken from A global solution for the Mars static and seasonal gravity, Mars orientation, Phobos and Deimos masses, and Mars ephemeris
    # Free core nutation rate of the celestial body of Mars
    parameter_perturbation[7]=np.deg2rad(0.075)/constants.JULIAN_DAY #rad/s; Taken from A global solution for the Mars static and seasonal gravity, Mars orientation, Phobos and Deimos masses, and Mars ephemeris
    # Ground station position of Mars    
    parameter_perturbation[8:11]=30*np.ones(3) # meters; Taken from Position Determination of a Lander and Rover at Mars With Warth-Based Differential Tracking
    # Periodic spin variation for full planetary rotational model of Mars
    # First order - cosine term
    parameter_perturbation[11]=23*mas # seconds; Taken from the PhD from Sebastien LeMaistre
    # First order - sine term
    parameter_perturbation[12]=26*mas # seconds; Taken from the PhD from Sebastien LeMaistre
    # Second order - cosine term
    parameter_perturbation[13]=22*mas # seconds; Taken from the PhD from Sebastien LeMaistre
    # Second order - sine term
    parameter_perturbation[14]=22*mas # seconds; Taken from the PhD from Sebastien LeMaistre
    # Third order - cosine term
    parameter_perturbation[15]=18*mas # seconds; Taken from the PhD from Sebastien LeMaistre
    # Third order - sine term
    parameter_perturbation[16]=19*mas # seconds; Taken from the PhD from Sebastien LeMaistre
    # Fourth order - cosine term
    parameter_perturbation[17]=16*mas # seconds; Taken from the PhD from Sebastien LeMaistre
    # Fourth order - sine term
    parameter_perturbation[18]=16*mas # seconds; Taken from the PhD from Sebastien LeMaistre
    # Polar motion amplitude for full planetary rotational model of Mars
    parameter_perturbation[19:]=2*mas*np.ones(20) # seconds; Taken from UNCERTAINTIES ON MARS INTERIOR PARAMETERS DEDUCED FROM ORIENTATION PARAMETERS USING DIFFERENT RADIOLINKS: ANALYTICAL SIMULATIONS.

    print("Perturbation vector is:")
    print(parameter_perturbation)

    # Define a priori covariance
    inverse_a_priori_covariance = np.diag(1/parameter_perturbation**2)

    # Estimate parameters
    pod_input = estimation.PodInput(simulated_observations,parameters_set.parameter_set_size, inverse_apriori_covariance = inverse_a_priori_covariance, apriori_parameter_correction = parameter_perturbation)
    pod_input.set_constant_weight_per_observable(weights_per_observable)
    #pod_input.define_estimation_settings(reintegrate_variational_equations = False)
    
    # Perform estimation
    pod_output = estimator.perform_estimation(pod_input)
    
    ########################################################################################################################
    ################################################## PROVIDE OUTPUT TO CONSOLE AND FILES #################################
    ########################################################################################################################

    output_folder_path = os.path.dirname(os.path.realpath(__file__)).replace('/src','/output/POD_LaRa')
    os.makedirs(output_folder_path,exist_ok=True)

    estimation_error = np.subtract(pod_output.parameter_estimate,truth_parameter)
    print('True estimation error is:')
    print(estimation_error.astype('d'))
    formal_error = pod_output.formal_errors
    print('Formal estimation error is:')
    print(formal_error.astype('d'))
    true_to_form_estimation_error_ratio = estimation_error/formal_error 
    print('True to form estimation error is:')
    print(true_to_form_estimation_error_ratio)
    estimation_information_matrix = pod_output.normalized_design_matrix
    estimation_information_matrix_normalization = pod_output.normalization_terms
    concatenated_times = simulated_observations.concatenated_times 
    concatenated_link_ends = simulated_observations.concatenated_link_ends
    doppler_residuals = pod_output.residual_history

    np.savetxt(output_folder_path+"/estimation_information_matrix.dat",estimation_information_matrix,fmt='%.15e')
    np.savetxt(output_folder_path+"/estimation_information_matrix_normalization.dat",
        estimation_information_matrix_normalization,fmt='%.15e')
    np.savetxt(output_folder_path+"/concatenated_times.dat",concatenated_times,fmt='%.15e')
    np.savetxt(output_folder_path+"/concatenated_link_ends.dat",concatenated_link_ends,fmt='%.15e')
    np.savetxt(output_folder_path+"/doppler_residuals.dat",doppler_residuals,fmt='%.15e')
    
    ########################################################################################################################
    ################################################## PLOTTING TRUE TO FORM RATIO #########################################
    ########################################################################################################################

    plt.hist(np.abs(true_to_form_estimation_error_ratio), bins = 8)
    plt.ylabel('Frequency [-]')
    plt.xlabel('True to form ratio [-]')
    plt.grid()
    plt.savefig(output_folder_path+"/true_to_form_ratio.pdf",bbox_inches="tight")
    plt.show()

print("--- %s seconds ---" % (time.time() - run_time))