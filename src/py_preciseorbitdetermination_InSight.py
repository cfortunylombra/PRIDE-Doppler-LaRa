"""
Description: Environment Setup for the Precise Orbit Determination (RISE)

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
    import numpy as np
    import scipy.interpolate
    import matplotlib.pyplot as plt
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
    days_in_a_week = 7 #days

    # Initial date of the simulation
    start_date = 2458449.5 #in Julian days = 27/11/2018 00:00:00; Taken from the txt file sent by Sebastien

    # Duration of the simulation
    simulation_duration_days = 998 #days 
    simulation_duration = simulation_duration_days*constants.JULIAN_DAY #seconds

    # RISE landing site, taken from "LaRa after RISE: Expected improvement in the Mars rotation and interior models"
    reflector_name = "RISE"
    reflector_latitude_deg = 4.5 #North degrees
    reflector_longitude_deg = 135.62 #East degrees
    reflector_altitude = -2600 #m

    # Earth-based transmitters
    transmitter_names = ['DSS 43','DSS 34','DSS 35','DSS 36','DSS 65','DSS 63','DSS 55','DSS 54','DSS 56','DSS 14','DSS 26', 'DSS 24', 'DSS 25']

    transmitter_positions_cartesian = list()  #Taken from https://www.aoc.nrao.edu/software/sched/catalogs/locations.dat
    transmitter_positions_cartesian.append(np.array([-4460894.7273,2682361.5296,-3674748.4238])) # DSS 43
    transmitter_positions_cartesian.append(np.array([-4461147.4205,2682439.2423,-3674392.5623])) # DSS 34
    transmitter_positions_cartesian.append(np.array([-4461273.4175,2682568.9283,-3674151.5223])) # DSS 35
    transmitter_positions_cartesian.append(np.array([-4461168.7425,2682814.6603,-3674083.3303])) # DSS 36
    transmitter_positions_cartesian.append(np.array([4849339.5378,-360427.4851,4114750.8520])) # DSS 65A
    transmitter_positions_cartesian.append(np.array([4849092.6814,-360180.5350,4115109.1298])) # DSS 63
    transmitter_positions_cartesian.append(np.array([4849525.256,-360606.09,4114495.08])) # DSS 55 #http://astrogeo.org/aplo/vlbi.inp
    transmitter_positions_cartesian.append(np.array([4849434.4880,-360723.8999,4114618.8350])) # DSS 54
    transmitter_positions_cartesian.append(np.array([4849421.500903,-360549.2280048,4114647.264832])) # DSS 56 #https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/stations/earth_topo_201023.tf
    transmitter_positions_cartesian.append(np.array([-2353621.2459,-4641341.5369,3677052.2305])) # DSS 14
    transmitter_positions_cartesian.append(np.array([-2354890.967,-4647166.93,3668872.21])) # DSS 26
    transmitter_positions_cartesian.append(np.array([-2354906.495,-4646840.13,3669242.317])) # DSS 24
    transmitter_positions_cartesian.append(np.array([-2355022.066,-4646953.64,3669040.90])) # DSS 25

    # Viability settings
    earth_min = 10 #deg
    earth_max = 30 #deg
    antenna_min_elevation = 20 #deg

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
        # Only if shadow tracking
        if receiver_name in transmitter_names:
            continue
        for transmitter_name in transmitter_names:
            two_way_link_ends = dict()
            two_way_link_ends[observation.transmitter] = ("Earth", transmitter_name )
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
    with open(os.path.dirname(os.path.realpath(__file__))+'/InSight_mes_upto_31122021.forCarlos') as file:
        lines = file.read().splitlines()
        observation_start_epoch = np.inf
        # Iterate along each transmitter
        for transmitter_pointer in transmitter_names:
            transmitter_ground_station_number =  [int(s) for s in transmitter_pointer.split() if s.isdigit()][0]
            # Iterate along each observation time
            for pointer_time in range(0,len(lines)):
                line = lines[pointer_time]
                line_info = line.split()
                # Condition to save the observation time
                if float(line_info[0]) == transmitter_ground_station_number and float(line_info[1]) == transmitter_ground_station_number \
                    and float(line_info[2])<=simulation_end_epoch:
                    observation_times_list.append(float(line_info[2]))
                    # Save the minimum epoch
                    if float(line_info[2])<observation_start_epoch:
                        observation_start_epoch = float(line_info[2])

    observation_times_list.sort()

    # Create observation viability settings and calculators
    viability_settings_list = list()
    viability_settings_list.append(observation.elevation_angle_viability(["Earth",""],np.deg2rad(antenna_min_elevation)))
    viability_settings_list.append(observation.elevation_angle_viability(["Mars",""],np.deg2rad(earth_min)))
    #viability_settings_list.append(observations.elevation_angle_viability(["Mars",""],np.deg2rad(earth_max))) #NOTE maximum elevation angle viability
    viability_settings_list.append(observation.body_avoidance_viability(["Earth",""],"Sun",np.deg2rad(antenna_min_elevation)))
    viability_settings_list.append(observation.body_occultation_viability(("Earth",""),"Moon"))

    # Change directory in order to read ResStatPerPass_ForCarlos.txt
    output_folder_path = os.path.dirname(os.path.realpath(__file__))
    os.makedirs(output_folder_path,exist_ok=True)

    time_days_mHz_pass = list()
    std_mHz = list()

    # Append the standard deviations and times to the empty lists
    with open(output_folder_path+'/ResStatPerPass_ForCarlos.txt') as f:
        lines = f.readlines()
        for line in lines[1:]:
            line_split = line.split()
            if not (np.isnan(float(line_split[1])) and np.isnan(float(line_split[2]))):
                time_days_mHz_pass.append(float(line_split[0]))
                std_mHz.append(float(line_split[2]))
    
    # Nearest interpolation
    std_mHz_function = scipy.interpolate.interp1d(time_days_mHz_pass, std_mHz, fill_value='extrapolate', kind='nearest')

    # Insert seed
    np.random.seed(42)

    # Function to compute the standard deviation
    def std_mHz_callable(t):
        return np.array([np.random.normal(0,std_mHz_function((t-observation_times_list[0])/constants.JULIAN_DAY)*10**(-3)/constants.SPEED_OF_LIGHT_LONG)])

    observation_simulation_settings = list()
    for pointer_link_ends in range(0,len(observation_settings_list)):
        observation_simulation_settings.append(observation.tabulated_simulation_settings(observation.two_way_doppler_type,
            observation_settings_list[pointer_link_ends],observation_times_list,
            viability_settings = viability_settings_list,reference_link_end_type = observation.receiver,
            noise_function = std_mHz_callable))

    # Simulate required observation
    simulated_observations = estimation.simulate_observations(observation_simulation_settings, observation_simulators, bodies)

    # Perturbation
    parameter_perturbation = np.zeros(parameters_set.parameter_set_size) 
    mas = np.pi/(180.0*1000.0*3600.0) # Conversion from milli arc seconds to seconds 
    # Position of Mars
    parameter_perturbation[0:3]=1000*np.ones(3) # meters; Taken from Improving the Accuracy of the Martian Ephemeris Short-Term Prediction
    # Velocity of Mars
    parameter_perturbation[3:6]=0.0002*np.ones(3) # meters; Taken from Improving the Accuracy of the Martian Ephemeris Short-Term Prediction
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
    #pod_input.define_estimation_settings(reintegrate_variational_equations = False)
    
    # Define noise levels for weights
    vector_weights = std_mHz_function((simulated_observations.concatenated_times-observation_times_list[0]*np.ones(len(simulated_observations.concatenated_times)))/constants.JULIAN_DAY)*10**(-3)/constants.SPEED_OF_LIGHT_LONG
    pod_input.set_weight(1/vector_weights**2)

    # Just to check the viability of the sun, remove the last per pass
    plt.figure()
    plt.scatter((simulated_observations.concatenated_times-observation_times_list[0]*np.ones(len(simulated_observations.concatenated_times)))/constants.JULIAN_DAY,std_mHz_function((simulated_observations.concatenated_times-observation_times_list[0]*np.ones(len(simulated_observations.concatenated_times)))/constants.JULIAN_DAY))
    plt.ylabel('Weight [-]')
    plt.xlabel('Time [-]')
    plt.grid()
    plt.ylim([-3,3])
    plt.show()
    plt.close('all')    

    # Perform estimation
    pod_output = estimator.perform_estimation(pod_input)
    
    ########################################################################################################################
    ################################################## PROVIDE OUTPUT TO CONSOLE AND FILES #################################
    ########################################################################################################################

    output_folder_path = os.path.dirname(os.path.realpath(__file__)).replace('/src','/output/POD_RISE')
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

    np.savetxt(output_folder_path+"/weights_diagonal.dat",pod_output.weights_matrix_diagonal,fmt='%.15e')
    np.savetxt(output_folder_path+"/estimation_information_matrix.dat",estimation_information_matrix,fmt='%.15e')
    np.savetxt(output_folder_path+"/estimation_information_matrix_normalization.dat",
        estimation_information_matrix_normalization,fmt='%.15e')
    np.savetxt(output_folder_path+"/concatenated_times.dat",concatenated_times,fmt='%.15e')
    np.savetxt(output_folder_path+"/concatenated_link_ends.dat",concatenated_link_ends,fmt='%.15e')
    np.savetxt(output_folder_path+"/doppler_residuals.dat",doppler_residuals,fmt='%.15e')
    
    ########################################################################################################################
    ################################################## PLOTTING TRUE TO FORM RATIO #########################################
    ########################################################################################################################

    plt.figure()
    plt.hist(np.abs(true_to_form_estimation_error_ratio), bins = 8)
    plt.ylabel('Frequency [-]')
    plt.xlabel('True to form ratio [-]')
    plt.grid()
    plt.savefig(output_folder_path+"/true_to_form_ratio.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

print("--- %s seconds ---" % (time.time() - run_time))