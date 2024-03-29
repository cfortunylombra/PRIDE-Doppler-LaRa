"""
Description: Environment Setup for the Precise Orbit Determination (RISE and LaRa with DSN-PRIDE) without intial state, but including the cross-correlation metric

Author: C. Fortuny-Lombraña
"""
import time
run_time = time.time()
if __name__=="__main__":
    ########################################################################################################################
    ################################################## IMPORT PACKAGES #####################################################
    ########################################################################################################################

    import sys
    sys.path.insert(0, "/home/carlos/tudat-bundle/build_release/tudatpy/")#sys.path.insert(0, "/home/cfortunylombra/tudat-bundle/cmake-build-release-wsl/tudatpy/")
    import os
    import numpy as np
    import datetime
    import scipy.interpolate
    import scipy.sparse
    import scipy.sparse.linalg
    import matplotlib.pyplot as plt
    from multiprocessing import Pool
    from tudatpy.kernel import constants, numerical_simulation
    from tudatpy.kernel.astro import element_conversion
    from tudatpy.kernel.interface import spice_interface
    from tudatpy.kernel.numerical_simulation import environment_setup,propagation_setup,propagation,estimation_setup,estimation
    from tudatpy.kernel.numerical_simulation.estimation_setup import observation

    np.set_printoptions(suppress=True,precision=15)
    ########################################################################################################################
    ################################################## CONSTANTS AND VARIABLES #############################################
    ########################################################################################################################

    # Days in a week
    days_in_a_week = 7 #days

    # CPU number for parallel computing
    CPU_par = 14

    # Booleans to understand whether we want to simulate together RISE and LaRa missions, or separetely - These booleans can be swapped to FALSE
    RISE_boolean = True
    LaRa_boolean = True
    
    if LaRa_boolean:
        # PRIDE stations boolean - This boolean can be swapped to FALSE
        PRIDE_boolean = True
        
        # Boolean for removing PRIDE stations only after the POD - Just for verification of the W matrix
        remove_PRIDE_weight_boolean = False

    else:
        # PRIDE stations boolean - This boolean can be swapped to TRUE (Since RISE does not have any PRIDE stations=NO DIFFERENCE)
        PRIDE_boolean = False
        
        # Boolean for removing PRIDE stations only after the POD - Jus for verification of the W matrix
        remove_PRIDE_weight_boolean = False

    # Evaluation step size for the plots of the formal errors as function of time - This int can be increased for decreasing the computation time
    step_eval = 1

    # Output folder - Creates new folders for the new plots
    if LaRa_boolean:
        output_folder_path = os.path.dirname(os.path.realpath(__file__)).replace('/src','/output/POD_RISE'+str(RISE_boolean)+'_LaRa'+str(LaRa_boolean)+'_LaRanointcomplex'+str(PRIDE_boolean)+str(remove_PRIDE_weight_boolean))
    else:
        output_folder_path = os.path.dirname(os.path.realpath(__file__)).replace('/src','/output/POD_RISE'+str(RISE_boolean)+'_LaRanointcomplex'+str(LaRa_boolean))
    os.makedirs(output_folder_path,exist_ok=True)

    if RISE_boolean:
        # Initial date of the simulation
        start_date = 2458449.5 #in Julian days = 27/11/2018 00:00:00; Taken from the txt file sent by Sebastien

        # Duration of the simulation
        RISE_simulation_duration_days = 998 #days 
        RISE_simulation_duration = RISE_simulation_duration_days*constants.JULIAN_DAY #seconds
        if LaRa_boolean:
            simulation_duration_days = 2956 #days; until 31/12/2026
        else:
            simulation_duration_days = 998 #days
        simulation_duration = simulation_duration_days*constants.JULIAN_DAY #seconds
    
    elif LaRa_boolean:
        # Initial date of the simulation
        start_date = 2459580.5 #in Julian days = 01/01/2022; Taken from the txt file sent by Sebastien
        
        # Duration of the simulation
        RISE_simulation_duration_days = 0 #days
        RISE_simulation_duration = RISE_simulation_duration_days*constants.JULIAN_DAY #seconds
        simulation_duration_days = 2956-998 #days; until 31/12/2026
        simulation_duration = simulation_duration_days*constants.JULIAN_DAY #seconds

    # RISE landing site, taken from "LaRa after RISE: Expected improvement in the Mars rotation and interior models"
    RISE_reflector_name = "RISE"
    RISE_reflector_latitude_deg = 4.5 #North degrees
    RISE_reflector_longitude_deg = 135.62 #East degrees
    RISE_reflector_altitude = -2600 #m

    # LaRa landing site, taken from "LaRa after RISE: Expected improvement in the Mars rotation and interior models"
    LaRa_reflector_name = "LaRa"
    LaRa_reflector_latitude_deg = 18.3 #North degrees
    LaRa_reflector_longitude_deg = 335.37 #East degrees
    LaRa_reflector_altitude = -2000 #m

    base_frequency = 8400.5e6 #Hz; taken from Giuseppe's files

    # Earth-based transmitters
    transmitters_dict = dict() #Taken from JPL web site
    transmitters_dict['DSS 43']=np.array([-4460894.9170,2682361.5070,-3674748.1517]) # DSS 43
    transmitters_dict['DSS 34']=np.array([-4461147.0925,2682439.2385,-3674393.1332]) # DSS 34
    transmitters_dict['DSS 35']=np.array([-4461273.4175,2682568.9283,-3674151.5223]) # DSS 35 (https://www.aoc.nrao.edu/software/sched/catalogs/locations.dat)
    transmitters_dict['DSS 36']=np.array([-4461168.7425,2682814.6603,-3674083.3303]) # DSS 36 (https://www.aoc.nrao.edu/software/sched/catalogs/locations.dat)
    transmitters_dict['DSS 65']=np.array([4849339.6448,-360427.6560,4114750.7428]) # DSS 65
    transmitters_dict['DSS 63']=np.array([4849092.5175,-360180.3480,4115109.2506]) # DSS 63
    transmitters_dict['DSS 55']=np.array([4849525.2561,-360606.0932,4114495.0843]) # DSS 55
    transmitters_dict['DSS 54']=np.array([4849434.4877,-360723.8999,4114618.8354]) # DSS 54
    transmitters_dict['DSS 56']=np.array([4849421.500903,-360549.2280048,4114647.264832]) # DSS 56 (https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/stations/earth_topo_201023.tf)
    transmitters_dict['DSS 14']=np.array([-2353621.4197,-4641341.4717,3677052.3178]) # DSS 14
    transmitters_dict['DSS 26']=np.array([-2354890.7996,-4647166.3182,3668871.7546]) # DSS 26
    transmitters_dict['DSS 24']=np.array([-2354906.7087,-4646840.0834,3669242.3207]) # DSS 24
    transmitters_dict['DSS 25']=np.array([-2355022.0140,-4646953.2040,3669040.5666]) # DSS 25

    # Earth-based transmitter for RISE
    RISE_transmitter_names = ['DSS 43','DSS 34','DSS 35','DSS 36','DSS 65','DSS 63','DSS 55','DSS 54','DSS 56','DSS 14','DSS 26', 'DSS 24', 'DSS 25']

    # Earth-based transmitter for LaRa
    LaRa_transmitter_names = ['DSS 43','DSS 63','DSS 14']

    # Viability settings for RISE
    RISE_earth_min = 10 #deg
    RISE_earth_max = 30 #deg
    RISE_antenna_min_elevation = 10 #deg
    RISE_body_avoidance_angle = 10 #deg

    # Viability settings for LaRa
    LaRa_earth_min = 35 #deg
    LaRa_earth_max = 45 #deg
    LaRa_antenna_min_elevation = 10 #deg
    LaRa_body_avoidance_angle = 10 #deg

    # Noise sources (X-band at tau=60s)
    frequency_standard=10**(np.log10(8e-16)+1/2*(np.log10(1000)-np.log10(60)))
    antenna_mechanical=1.6e-14
    ground_electronics=10**(np.log10(2e-16)+1/2*(np.log10(1000)-np.log10(60)))
    stochastic_spacecraft_motion=10**(np.log10(2e-16)+1/2*(np.log10(1000)-np.log10(60)))
    receiver_thermal_noise=10**(np.log10(1e-16)+1/2*(np.log10(1000)-np.log10(60)))
    spacecraft_transponder=1.8e-14
    tropospheric_scintillation=6.5e-14

    # Plasma noise function (tau=60s)
    def plasma_noise_function(SEP_angle_value):
        if SEP_angle_value>=np.deg2rad(0) and SEP_angle_value<=np.deg2rad(90):
            return (1.76*10**(-14)*(np.sin(SEP_angle_value))**(-1.98)+6.25*10**(-14)*(np.sin(SEP_angle_value))**(0.06))
        elif SEP_angle_value>np.deg2rad(90) and SEP_angle_value<=np.deg2rad(170):
            return (1.76*10**(-14)+6.25*10**(-14))*(np.sin(SEP_angle_value))**(1.05)
        elif SEP_angle_value>np.deg2rad(170) and SEP_angle_value<=np.deg2rad(180):
            return 1.27*10**(-14)

    # Common noise percentage
    def common_Doppler_noise_percentage_function(SEP_angle_value):
        # Noise sources (X-band at tau=60s)
        frequency_standard=10**(np.log10(8e-16)+1/2*(np.log10(1000)-np.log10(60)))
        antenna_mechanical=1.6e-14
        ground_electronics=10**(np.log10(2e-16)+1/2*(np.log10(1000)-np.log10(60)))
        stochastic_spacecraft_motion=10**(np.log10(2e-16)+1/2*(np.log10(1000)-np.log10(60)))
        receiver_thermal_noise=10**(np.log10(1e-16)+1/2*(np.log10(1000)-np.log10(60)))
        spacecraft_transponder=1.8e-14
        tropospheric_scintillation=6.5e-14
        plasma_noise = plasma_noise_function(SEP_angle_value)
        
        common_noise = np.sqrt(frequency_standard**2+plasma_noise**2+stochastic_spacecraft_motion**2+spacecraft_transponder**2+(tropospheric_scintillation/np.sqrt(2))**2)
        total_noise = np.sqrt(frequency_standard**2+antenna_mechanical**2+ground_electronics**2+stochastic_spacecraft_motion**2+receiver_thermal_noise**2+spacecraft_transponder**2+tropospheric_scintillation**2+plasma_noise**2)
        return common_noise**2/total_noise**2

    def weather_percentage_function(SEP_angle_value):
        frequency_standard=10**(np.log10(8e-16)+1/2*(np.log10(1000)-np.log10(60)))
        antenna_mechanical=1.6e-14
        ground_electronics=10**(np.log10(2e-16)+1/2*(np.log10(1000)-np.log10(60)))
        stochastic_spacecraft_motion=10**(np.log10(2e-16)+1/2*(np.log10(1000)-np.log10(60)))
        receiver_thermal_noise=10**(np.log10(1e-16)+1/2*(np.log10(1000)-np.log10(60)))
        spacecraft_transponder=1.8e-14
        tropospheric_scintillation=6.5e-14
        plasma_noise = plasma_noise_function(SEP_angle_value)

        weather_noise = np.sqrt(tropospheric_scintillation**2)
        total_noise = np.sqrt(frequency_standard**2+antenna_mechanical**2+ground_electronics**2+stochastic_spacecraft_motion**2+receiver_thermal_noise**2+spacecraft_transponder**2+tropospheric_scintillation**2+plasma_noise**2)
        return weather_noise**2/total_noise**2

    # Min Doppler Noise from ED045 observations
    min_Doppler_noise_real_data = dict()
    min_Doppler_noise_real_data['MEDICINA'] = 1.1893669211845945e-13
    min_Doppler_noise_real_data['WETTZELL'] = 5.523841298583862e-13
    min_Doppler_noise_real_data['ONSALA60'] = 1.460047359922946e-13
    min_Doppler_noise_real_data['EFLSBERG'] = 6.601389242000477e-14
    min_Doppler_noise_real_data['WRT0'] = 2.8182590417184753e-13
    min_Doppler_noise_real_data['YEBES40M'] = 1.0347841930949488e-13
    min_Doppler_noise_real_data['TIANMA65'] = 3.799449475342986e-14
    min_Doppler_noise_real_data['CEDUNA'] = 2.4417833773786286e-12
    min_Doppler_noise_real_data['BADARY'] = 1.9464314322469196e-13
    min_Doppler_noise_real_data['HARTRAO'] = 7.121852569380006e-13
    min_Doppler_noise_real_data['IRBENE'] = 3.467891060872128e-10
    for LaRa_transmitter in LaRa_transmitter_names:
        min_Doppler_noise_real_data[LaRa_transmitter] = 2.564768287601929e-14

    # Parameters chosen from an analysis performed
    distance_ionospheric = 413e3

    def indicator_Doppler_noise_function(SEP_angle_value):
        common_Doppler_noise = common_Doppler_noise_percentage_function(SEP_angle_value)*min_Doppler_noise_real_data[LaRa_transmitter_names[0]]
        
        indicator_Doppler_noise = dict()
        for station_1 in list(min_Doppler_noise_real_data.keys()):
            indicator_Doppler_noise[station_1] = dict()
            for station_2 in list(min_Doppler_noise_real_data.keys()):
                if station_1 != station_2:
                    linear_value = ((min_Doppler_noise_real_data[station_1]-common_Doppler_noise)+(min_Doppler_noise_real_data[station_2]-common_Doppler_noise))/(2*common_Doppler_noise+2*(plasma_noise_function(SEP_angle_value)-plasma_noise_function(np.pi)))
                    indicator_Doppler_noise[station_1][station_2] = 2-2/(1+np.exp(-10**(-1/2)*linear_value))
        return indicator_Doppler_noise

    ########################################################################################################################
    ################################################## CREATE ENVIRONMENT ##################################################
    ########################################################################################################################

    # Load spice kernels
    spice_interface.load_standard_kernels()

    # Initial and end time of the simulation
    simulation_start_epoch = (start_date-constants.JULIAN_DAY_ON_J2000)*constants.JULIAN_DAY #seconds
    simulation_end_epoch = simulation_start_epoch+simulation_duration #seconds
    RISE_simulation_end_epoch = simulation_start_epoch+RISE_simulation_duration #seconds

    # Define bodies in the simulation
    bodies_to_create = ["Saturn","Jupiter","Mars","Moon","Earth","Venus","Mercury","Sun"] 

    # Define frame
    global_frame_origin = "SSB" #Barycenter of Solar System
    global_frame_orientation = "ECLIPJ2000"
    
    # Define body settings
    body_settings = environment_setup.get_default_body_settings_time_limited(
        bodies_to_create,
        simulation_start_epoch-constants.JULIAN_DAY,
        simulation_end_epoch+constants.JULIAN_DAY,
        global_frame_origin,
        global_frame_orientation,
        time_step = 60)

    # Reset frame origin
    environment_setup.ephemeris.frame_origin = "Sun"

    #Mars rotation model (High-accuracy)
    body_settings.get("Mars").rotation_model_settings = environment_setup.rotation_model.mars_high_accuracy()

    bodies = environment_setup.create_system_of_bodies(body_settings)

    ########################################################################################################################
    ################################################## CREATE GROUND STATIONS AND LANDER ###################################
    ########################################################################################################################

    # Empty dictionary for the radio telescopes coordinates
    radio_telescopes_dict = dict()

    # Read the text file containing the name and cartesian coordinates of the radio telescopes
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
            
            # Since the Sebastien files do not have Hart15M and Hobart12 radio telescopes, they are not included in the simulation
            if name_ground_station=="HART15M" or name_ground_station=="HOBART12" or name_ground_station=="WARK30M" or name_ground_station=="HOBART26" or name_ground_station=="IRBENE":
                continue
            else:
                x_coordinate_ground_station = float(coordinates_line_ground_station.split("X=",1)[1].split()[0])
                y_coordinate_ground_station = float(coordinates_line_ground_station.split("Y=",1)[1].split()[0])
                z_coordinate_ground_station = float(coordinates_line_ground_station.split("Z=",1)[1].split()[0])

                radio_telescopes_dict[name_ground_station] = np.array([x_coordinate_ground_station,y_coordinate_ground_station,z_coordinate_ground_station])
    
    # Earth-based ground station creation
    for pointer_ground_station in range(0,len(transmitters_dict.keys())):
        environment_setup.add_ground_station(
            bodies.get_body("Earth"),
            list(transmitters_dict.keys())[pointer_ground_station],
            transmitters_dict[list(transmitters_dict.keys())[pointer_ground_station]])

    # Earth-based radio telescope creation
    for pointer_radio_telescope in range(0,len(radio_telescopes_dict.keys())):
        environment_setup.add_ground_station(
            bodies.get_body("Earth"),
            list(radio_telescopes_dict.keys())[pointer_radio_telescope],
            radio_telescopes_dict[list(radio_telescopes_dict.keys())[pointer_radio_telescope]])
    
    Earth_ground_station_list = environment_setup.get_ground_station_list(bodies.get_body("Earth"))
    
    # Mars-based ground stations creation
    environment_setup.add_ground_station(
        bodies.get_body("Mars"),
        RISE_reflector_name,
        np.array([RISE_reflector_altitude,np.deg2rad(RISE_reflector_latitude_deg,dtype='d'),np.deg2rad(RISE_reflector_longitude_deg,dtype='d')]),
         element_conversion.geodetic_position_type)

    environment_setup.add_ground_station(
        bodies.get_body("Mars"),
        LaRa_reflector_name,
        np.array([LaRa_reflector_altitude,np.deg2rad(LaRa_reflector_latitude_deg,dtype='d'),np.deg2rad(LaRa_reflector_longitude_deg,dtype='d')]),
         element_conversion.geodetic_position_type)

    Mars_ground_station_list = environment_setup.get_ground_station_list(bodies.get_body("Mars"))

    # Distance between stations
    distance_stations = dict()
    Earth_ground_station_dict = {**transmitters_dict,**radio_telescopes_dict}
    for station_1 in Earth_ground_station_dict:
        distance_stations[station_1] = dict()
        for station_2 in Earth_ground_station_dict:
            if station_1 != station_2:
                distance_stations[station_1][station_2] = np.sqrt((Earth_ground_station_dict[station_1][0]-Earth_ground_station_dict[station_2][0])**2+(Earth_ground_station_dict[station_1][1]-Earth_ground_station_dict[station_2][1])**2+(Earth_ground_station_dict[station_1][2]-Earth_ground_station_dict[station_2][2])**2)

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
    ################################################## DEFINE OBSERVATIONS TIMES ###########################################
    ######################################################################################################################## 

    # Define time between two observations
    observation_interval = 60 #seconds

    # Define observation simulation times for each link (RISE)
    RISE_observation_times_list = list()
    RISE_observation_times_dict = dict()
    RISE_observation_start_epoch = None
    
    # Read the epoch times for RISE
    with open(os.path.dirname(os.path.realpath(__file__))+'/InSight_mes_upto_31122021.forCarlos') as file:
        lines = file.read().splitlines()
        # Iterate along each transmitter
        for RISE_transmitter_pointer in RISE_transmitter_names:
            RISE_observation_times_dict[RISE_transmitter_pointer] = list()
            transmitter_ground_station_number =  [int(s) for s in RISE_transmitter_pointer.split() if s.isdigit()][0]
            # Iterate along each observation time
            for pointer_time in range(0,len(lines)):
                line = lines[pointer_time]
                line_info = line.split()
                if pointer_time==0:
                    RISE_observation_start_epoch_reference_noise = float(line_info[2])

                # Condition to save the observation time
                if float(line_info[0]) == transmitter_ground_station_number and float(line_info[1]) == transmitter_ground_station_number \
                    and float(line_info[2])<=RISE_simulation_end_epoch and float(line_info[2])>=simulation_start_epoch and float(line_info[2])<=simulation_end_epoch:
                    RISE_observation_times_dict[RISE_transmitter_pointer].append(float(line_info[2]))
                    RISE_observation_times_list.append(float(line_info[2]))

        if len(RISE_observation_times_list)!=0:
            RISE_observation_start_epoch = min(RISE_observation_times_list)

    # Define observation simulation times for each link (LaRa)
    LaRa_observation_times_list = list()
    LaRa_observation_start_epoch = None
    
    # Read the epoch times for LaRa
    with open(os.path.dirname(os.path.realpath(__file__))+'/lsor_export_forTUDelft.txt') as file:
        lines = file.read().splitlines()
        for pointer_pass in range(1,len(lines)):
            line = lines[pointer_pass]
            line_info = line.split()

            # Take the start time of the observation pass
            startpass_year = int(line_info[0].split('-')[0])
            startpass_day_of_year = int(line_info[0].split('-')[1].split('T')[0])
            startpass_hour = int(line_info[0].split('-')[1].split('T')[1].split(':')[0])
            startpass_min = int(line_info[0].split('-')[1].split('T')[1].split(':')[1])
            startpass_sec = int(line_info[0].split('-')[1].split('T')[1].split(':')[2])
            startpass_date = datetime.datetime(startpass_year,1,1)+datetime.timedelta(days=startpass_day_of_year-1,hours=startpass_hour,minutes=startpass_min,seconds=startpass_sec)
            startpass_epoch = (startpass_date - datetime.datetime(2000,1,1,12,0,0,0)).total_seconds()
            
            # Take the end time of the observation pass
            endpass_year = int(line_info[1].split('-')[0])
            endpass_day_of_year = int(line_info[1].split('-')[1].split('T')[0])
            endpass_hour = int(line_info[1].split('-')[1].split('T')[1].split(':')[0])
            endpass_min = int(line_info[1].split('-')[1].split('T')[1].split(':')[1])
            endpass_sec = int(line_info[1].split('-')[1].split('T')[1].split(':')[2])
            endpass_date = datetime.datetime(startpass_year,1,1)+datetime.timedelta(days=endpass_day_of_year-1,hours=endpass_hour,minutes=endpass_min,seconds=endpass_sec)
            endpass_epoch = (endpass_date - datetime.datetime(2000,1,1,12,0,0,0)).total_seconds()
            
            # Create the list with observation from start to the end of the pass
            if simulation_end_epoch>=endpass_epoch:
                if simulation_start_epoch<=startpass_epoch:
                    LaRa_observation_times_list.extend(np.arange(startpass_epoch,endpass_epoch+observation_interval,observation_interval))
                elif simulation_start_epoch>startpass_epoch and simulation_start_epoch<=endpass_epoch:
                    LaRa_observation_times_list.extend(np.arange(simulation_start_epoch,endpass_epoch+observation_interval,observation_interval))
            elif simulation_end_epoch>startpass_epoch and simulation_end_epoch<=endpass_epoch:
                if simulation_start_epoch<=startpass_epoch:
                    LaRa_observation_times_list.extend(np.arange(startpass_epoch,simulation_end_epoch+observation_interval,observation_interval))
                elif simulation_start_epoch>startpass_epoch and simulation_start_epoch<=endpass_epoch:
                    LaRa_observation_times_list.extend(np.arange(simulation_start_epoch,simulation_end_epoch+observation_interval,observation_interval))
        
        if len(LaRa_observation_times_list)!=0:
            LaRa_observation_start_epoch = min(LaRa_observation_times_list)

    # Sortht the observation times, and find the start epoch
    RISE_observation_times_list.sort()
    LaRa_observation_times_list.sort()
    observation_start_epoch = min(RISE_observation_times_list+LaRa_observation_times_list)

    # Computing solar plasma noise for ExoMars mission
    LaRa_SEP_angle = list()
    LaRa_solar_noise = list()
    LaRa_time_solar_SEP = list()
    LaRa_time_solar_noise = list()

    # Iterate along all the LaRa observation times
    for LaRa_time_pointer in LaRa_observation_times_list:

        # Compute the the distances from the Sun to Mars and Earth
        r_earthsun = np.linalg.norm(spice_interface.get_body_cartesian_position_at_epoch("Earth","Sun","ECLIPJ2000","NONE",LaRa_time_pointer))
        r_marssun = np.linalg.norm(spice_interface.get_body_cartesian_position_at_epoch("Mars","Sun","ECLIPJ2000","NONE",LaRa_time_pointer))
        r_earthmars = np.linalg.norm(spice_interface.get_body_cartesian_position_at_epoch("Mars","Earth","ECLIPJ2000","NONE",LaRa_time_pointer))

        # Cosine rule in order to compute the SEP angle
        #LaRa_SEP_angle.append((np.pi-np.arccos(np.dot(r1,r2)/(np.linalg.norm(r1)*np.linalg.norm(r2))))/2)
        LaRa_SEP_angle.append(np.arccos((r_earthsun**2+r_earthmars**2-r_marssun**2)/(2*r_earthsun*r_earthmars)))

        # Compute the solar plasma noise using the equations shown in IMPROVED DOPPLER TRACKING SYSTEMS FOR DEEP SPACE NAVIGATION
        if LaRa_SEP_angle[-1]>=np.deg2rad(0) and LaRa_SEP_angle[-1]<=np.deg2rad(90):
            LaRa_solar_noise.append(0.21545336*10**(-3)/base_frequency+1.76*10**(-14)*(np.sin(LaRa_SEP_angle[-1]))**(-1.98)+6.25*10**(-14)*(np.sin(LaRa_SEP_angle[-1]))**(0.06))
            LaRa_time_solar_noise.append(LaRa_time_pointer)
        elif LaRa_SEP_angle[-1]>np.deg2rad(90) and LaRa_SEP_angle[-1]<=np.deg2rad(170):
            LaRa_solar_noise.append(0.21545336*10**(-3)/base_frequency+(1.76*10**(-14)+6.25*10**(-14))*(np.sin(LaRa_SEP_angle[-1]))**(1.05))
            LaRa_time_solar_noise.append(LaRa_time_pointer)
        elif LaRa_SEP_angle[-1]>np.deg2rad(170) and LaRa_SEP_angle[-1]<=np.deg2rad(180):
            LaRa_solar_noise.append(0.21545336*10**(-3)/base_frequency+1.27*10**(-14))
            LaRa_time_solar_noise.append(LaRa_time_pointer)

        LaRa_time_solar_SEP.append(LaRa_time_pointer)

    ########################################################################################################################
    ################################################## DEFINE LINK ENDS FOR OBSERVATIONS ###################################
    ######################################################################################################################## 

    # Create list of observation settings
    observation_settings_list = list()

    # Define link ends for RISE
    if len(RISE_observation_times_list)!=0:
        for pointer_RISE_transmitter in RISE_transmitter_names:
            two_way_link_ends = dict()
            two_way_link_ends[observation.transmitter] = ("Earth",pointer_RISE_transmitter)
            two_way_link_ends[observation.reflector1] = ("Mars",RISE_reflector_name)
            two_way_link_ends[observation.receiver] = ("Earth",pointer_RISE_transmitter)

            observation_settings_list.append(two_way_link_ends)
        
    # Number of link ends for InSight mission
    RISE_link_ends_length = len(observation_settings_list)

    # Define link ends for LaRa
    if len(LaRa_observation_times_list)!=0:
        for pointer_LaRa_transmitter in LaRa_transmitter_names:
            two_way_link_ends = dict()
            two_way_link_ends[observation.transmitter] = ("Earth",pointer_LaRa_transmitter)
            two_way_link_ends[observation.reflector1] = ("Mars",LaRa_reflector_name)
            two_way_link_ends[observation.receiver] = ("Earth",pointer_LaRa_transmitter)

            observation_settings_list.append(two_way_link_ends)

        # Define link ends for LaRa-PRIDE (adding PRIDE stations)
        if PRIDE_boolean:
            for pointer_LaRa_transmitter in LaRa_transmitter_names:
                for pointer_LaRa_pride_station in list(radio_telescopes_dict.keys()):
                    two_way_link_ends = dict()
                    two_way_link_ends[observation.transmitter] = ("Earth",pointer_LaRa_transmitter)
                    two_way_link_ends[observation.reflector1] = ("Mars",LaRa_reflector_name)
                    two_way_link_ends[observation.receiver] = ("Earth",pointer_LaRa_pride_station)

                    observation_settings_list.append(two_way_link_ends)
    
    # Number of link ends for ExoMars mission
    LaRa_link_ends_length = len(observation_settings_list)-RISE_link_ends_length

    # Since simulated_observations orders the observations by alphabetical order, the sorted link ends list is determined
    link_ends_sort = list()
    DSN_link_ends_number = 0

    # Iterate along all the Earth DSN stations
    for transmitter_pointer in Earth_ground_station_list:
        if transmitter_pointer[1].startswith("DSS") and len(LaRa_observation_times_list)!=0:
            # Alphabetical order for the LaRa mission
            if transmitter_pointer[1] in LaRa_transmitter_names:
                boolean_DSN_receiver = True
                for receiver_pointer in Earth_ground_station_list:
                    # DSN receiving stations
                    if len(transmitter_pointer[1]) == len(receiver_pointer[1]) and transmitter_pointer[1] < receiver_pointer[1] and boolean_DSN_receiver:
                        two_way_link_ends = dict()
                        two_way_link_ends[observation.transmitter] = ("Earth",transmitter_pointer[1])
                        two_way_link_ends[observation.reflector1] = ("Mars",LaRa_reflector_name)
                        two_way_link_ends[observation.receiver] = ("Earth",transmitter_pointer[1])
                        link_ends_sort.append(two_way_link_ends)
                        DSN_link_ends_number+=1
                        boolean_DSN_receiver = False

                    # PRIDE stations
                    if not(receiver_pointer[1].startswith("DSS")) and PRIDE_boolean:
                        two_way_link_ends = dict()
                        two_way_link_ends[observation.transmitter] = ("Earth",transmitter_pointer[1])
                        two_way_link_ends[observation.reflector1] = ("Mars",LaRa_reflector_name)
                        two_way_link_ends[observation.receiver] = ("Earth",receiver_pointer[1])
                        link_ends_sort.append(two_way_link_ends)

        # Alphabetical order for the RISE mission
        if transmitter_pointer[1].startswith("DSS") and len(RISE_observation_times_list)!=0:
            if transmitter_pointer[1] in RISE_transmitter_names:
                two_way_link_ends = dict()
                two_way_link_ends[observation.transmitter] = ("Earth",transmitter_pointer[1])
                two_way_link_ends[observation.reflector1] = ("Mars",RISE_reflector_name)
                two_way_link_ends[observation.receiver] = ("Earth",transmitter_pointer[1])
                link_ends_sort.append(two_way_link_ends)
                DSN_link_ends_number+=1

    # Find the index position, transmitter and receiver for each link-end
    link_ends_numbers = list()
    link_ends_transmitter = list()
    link_ends_receiver = list()
    for link_end_pointer in observation_settings_list:
        link_ends_numbers.append(link_ends_sort.index(link_end_pointer))
        link_ends_transmitter.append(link_end_pointer[observation.transmitter][1])
        link_ends_receiver.append(link_end_pointer[observation.receiver][1])

    ########################################################################################################################
    ################################################## DEFINE PARAMETERS TO ESTIMATE #######################################
    ######################################################################################################################## 

    # Create list of parameters that are to be estimated
    parameter_settings = estimation_setup.parameter.initial_states(propagator_settings,bodies)
    parameter_settings.append(estimation_setup.parameter.core_factor("Mars"))
    parameter_settings.append(estimation_setup.parameter.free_core_nutation_rate("Mars"))
    if len(RISE_observation_times_list)!=0:
        parameter_settings.append(estimation_setup.parameter.ground_station_position("Mars", RISE_reflector_name))
    if len(LaRa_observation_times_list)!=0:
        parameter_settings.append(estimation_setup.parameter.ground_station_position("Mars", LaRa_reflector_name))
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
    
    # Create observation viability settings and calculators for RISE
    if len(RISE_observation_times_list)!=0:
        RISE_viability_settings_list = list()
        RISE_viability_settings_list.append(observation.minimum_elevation_angle_viability(["Earth",""],np.deg2rad(RISE_antenna_min_elevation)))
        RISE_viability_settings_list.append(observation.body_avoidance_viability(["Earth",""],"Sun",np.deg2rad(RISE_antenna_min_elevation)))
        RISE_viability_settings_list.append(observation.body_occultation_viability(("Earth",""),"Moon"))

    # Create observation viability settings and calculators for LaRa
    if len(LaRa_observation_times_list)!=0:
        LaRa_viability_settings_list = list()
        LaRa_viability_settings_list.append(observation.minimum_elevation_angle_viability(["Earth",""],np.deg2rad(LaRa_antenna_min_elevation)))
        LaRa_viability_settings_list.append(observation.minimum_elevation_angle_viability(["Mars",""],np.deg2rad(LaRa_earth_min)))
        LaRa_viability_settings_list.append(observation.maximum_elevation_angle_viability(["Mars",""],np.deg2rad(LaRa_earth_max)))
        LaRa_viability_settings_list.append(observation.body_avoidance_viability(["Earth",""],"Sun",np.deg2rad(LaRa_antenna_min_elevation)))
        LaRa_viability_settings_list.append(observation.body_occultation_viability(("Earth",""),"Moon"))
    
    #Change directory in order to read ResStatPerPass_ForCarlos.txt
    noise_folder_path = os.path.dirname(os.path.realpath(__file__))
    os.makedirs(noise_folder_path,exist_ok=True)

    RISE_time_days_mHz_pass = list()
    RISE_std_mHz = list()

    # Append the standard deviations and times to the empty lists
    with open(noise_folder_path+'/ResStatPerPass_ForCarlos.txt') as f:
        lines = f.readlines()
        for line in lines[1:]:
            line_split = line.split()
            if not (np.isnan(float(line_split[1])) and np.isnan(float(line_split[2]))):
                RISE_time_days_mHz_pass.append(float(line_split[0]))
                RISE_std_mHz.append(float(line_split[2]))

    # Nearest interpolation for RISE
    if len(RISE_observation_times_list)!=0:
        RISE_std_mHz_function = scipy.interpolate.interp1d(RISE_time_days_mHz_pass, RISE_std_mHz, fill_value='extrapolate', kind='nearest')

    # Nearest interpolation for LaRa
    if len(LaRa_observation_times_list)!=0:
        LaRa_std_noise_function = scipy.interpolate.interp1d(LaRa_time_solar_noise,LaRa_solar_noise,fill_value='extrapolate', kind='nearest')

    # Insert seed
    np.random.seed(42)

    # Function to compute the standard deviation for RISE
    def RISE_std_mHz_callable(t):
        return np.array([np.random.normal(0,RISE_std_mHz_function((t-RISE_observation_start_epoch_reference_noise)/constants.JULIAN_DAY)*10**(-3)/base_frequency)])

    # Function to compute the standard deviation for LaRa
    def LaRa_std_mHz_callable(t):
        return np.array([np.random.normal(0,LaRa_std_noise_function(t))])

    # Create global observation simulation settings
    observation_simulation_settings = list()

    # Create observation simulation settings for RISE
    if len(RISE_observation_times_list)!=0:
        for RISE_pointer_link_ends in range(0,RISE_link_ends_length):
            if len(RISE_observation_times_dict[RISE_transmitter_names[RISE_pointer_link_ends]])!=0:
                observation_simulation_settings.append(observation.tabulated_simulation_settings(observation.two_way_doppler_type,
                    observation_settings_list[RISE_pointer_link_ends],RISE_observation_times_dict[RISE_transmitter_names[RISE_pointer_link_ends]],
                    viability_settings = RISE_viability_settings_list,reference_link_end_type = observation.receiver,
                    noise_function = RISE_std_mHz_callable))

    # Create observation simulation settings for LaRa
    if len(LaRa_observation_times_list)!=0:
        for LaRa_pointer_link_ends in range(0,LaRa_link_ends_length):
            observation_simulation_settings.append(observation.tabulated_simulation_settings(observation.two_way_doppler_type,
                observation_settings_list[RISE_link_ends_length+LaRa_pointer_link_ends],LaRa_observation_times_list,
                viability_settings = LaRa_viability_settings_list,reference_link_end_type = observation.transmitter,
                noise_function = LaRa_std_mHz_callable))
    
    # Simulate required observation
    simulated_observations = estimation.simulate_observations(observation_simulation_settings, observation_simulators, bodies)

    # A priori
    apriori_vector = np.zeros(parameters_set.parameter_set_size) 
    mas =np.pi/(180.0*1000.0*3600.0) # Conversion from milli arc seconds to seconds 
    # Position of Mars
    apriori_vector[0:3]=1000*np.ones(3) # meters; Taken from Improving the Accuracy of the Martian Ephemeris Short-Term Prediction
    # Velocity of Mars
    apriori_vector[3:6]=0.0002*np.ones(3) #meters per sec; Taken from Improving the Accuracy of the Martian Ephemeris Short-Term Prediction
    # Core factor of the celestial body of Mars
    apriori_vector[6]=0.07 # Unitless; Taken from A global solution for the Mars static and seasonal gravity, Mars orientation, Phobos and Deimos masses, and Mars ephemeris
    # Free core nutation rate of the celestial body of Mars
    apriori_vector[7]=-np.deg2rad(1.5)/constants.JULIAN_DAY #rad/s; Taken from A global solution for the Mars static and seasonal gravity, Mars orientation, Phobos and Deimos masses, and Mars ephemeris
    # Ground station position of Mars    
    if len(RISE_observation_times_list)!=0 and len(LaRa_observation_times_list)!=0:
        add_par = 3
    else:
        add_par = 0
    apriori_vector[8:11+add_par]=30*np.ones(3+add_par) # meters; Taken from Position Determination of a Lander and Rover at Mars With Warth-Based Differential Tracking
    # Periodic spin variation for full planetary rotational model of Mars # Mars Pathfinder model
    # First order - cosine term
    apriori_vector[11+add_par]=23*mas # seconds; Taken from the PhD from Sebastien LeMaistre
    # First order - sine term
    apriori_vector[12+add_par]=26*mas # seconds; Taken from the PhD from Sebastien LeMaistre
    # Second order - cosine term
    apriori_vector[13+add_par]=22*mas # seconds; Taken from the PhD from Sebastien LeMaistre
    # Second order - sine term
    apriori_vector[14+add_par]=22*mas # seconds; Taken from the PhD from Sebastien LeMaistre
    # Third order - cosine term
    apriori_vector[15+add_par]=18*mas # seconds; Taken from the PhD from Sebastien LeMaistre
    # Third order - sine term
    apriori_vector[16+add_par]=19*mas # seconds; Taken from the PhD from Sebastien LeMaistre
    # Fourth order - cosine term
    apriori_vector[17+add_par]=16*mas # seconds; Taken from the PhD from Sebastien LeMaistre
    # Fourth order - sine term
    apriori_vector[18+add_par]=16*mas # seconds; Taken from the PhD from Sebastien LeMaistre
    # Polar motion amplitude for full planetary rotational model of Mars
    apriori_vector[19+add_par:]=50*mas*np.ones(20) # seconds; Taken from UNCERTAINTIES ON MARS INTERIOR PARAMETERS DEDUCED FROM ORIENTATION PARAMETERS USING DIFFERENT RADIOLINKS: ANALYTICAL SIMULATIONS.
    
    print("Apriori vector is:")
    print(apriori_vector)

    # Define a priori covariance
    inverse_a_priori_covariance = np.diag(1/apriori_vector**2)

    # Estimate parameters
    pod_input = estimation.PodInput(simulated_observations,parameters_set.parameter_set_size, inverse_apriori_covariance = inverse_a_priori_covariance)
    pod_input.define_estimation_settings(reintegrate_equations_on_first_iteration = False,reintegrate_variational_equations = False)
    
    # Observations visibles
    concatenated_times_array = np.array(simulated_observations.concatenated_times)
    if len(RISE_observation_times_list)!=0 and LaRa_observation_start_epoch!=None:
        RISE_concatenated_times = concatenated_times_array[concatenated_times_array<LaRa_observation_start_epoch]
    elif len(RISE_observation_times_list)!=0:
        RISE_concatenated_times = concatenated_times_array
    if len(LaRa_observation_times_list)!=0:
        LaRa_concatenated_times = concatenated_times_array[concatenated_times_array>=LaRa_observation_start_epoch]

    # Define noise levels for weights
    vector_weights = list()
    link_ends = simulated_observations.concatenated_link_end_names
    for index in range(0,len(concatenated_times_array)):
        link_end_name = link_ends[index]
        if link_end_name[observation.reflector1][1] == RISE_reflector_name:
            vector_weights.extend(list([RISE_std_mHz_function((concatenated_times_array[index]-RISE_observation_start_epoch_reference_noise)/constants.JULIAN_DAY)*10**(-3)/base_frequency]))
        if link_end_name[observation.reflector1][1] == LaRa_reflector_name:
            if link_end_name[observation.transmitter] == link_end_name[observation.receiver]:
                vector_weights.extend(list(LaRa_std_noise_function([concatenated_times_array[index]])))
            else:
                vector_weights.extend(list(LaRa_std_noise_function([concatenated_times_array[index]])))

    # Weights list becomes an array
    vector_weights = np.array(vector_weights)

    # .set_weight function has been created
    pod_input.set_weight(1/vector_weights**2) 

    # Perform estimation
    pod_output = estimator.perform_estimation(pod_input,convergence_checker = estimation.estimation_convergence_checker(maximum_iterations = 1))

    # Understand whether there are any duplicated time values
    if len(RISE_observation_times_list)!=0: 
        print('len(RISE concatenated times):',len(RISE_concatenated_times),' len(RISE obs times):',len(RISE_observation_times_list))
        print('Is there any duplicated RISE time value? :',any(list(RISE_concatenated_times).count(x) > 1 for x in list(RISE_observation_times_list)))
    
    if len(LaRa_observation_times_list)!=0:
        print('len(LaRa concatenated times):',len(LaRa_concatenated_times),' len(LaRa obs times):',len(LaRa_observation_times_list))
        print('Is there any duplicated LaRa time value? :',any(list(LaRa_concatenated_times).count(x) > 1 for x in list(LaRa_observation_times_list)))
     
    ########################################################################################################################
    ################################################## PROVIDE OUTPUT TO CONSOLE  ##########################################
    ########################################################################################################################

    # Computation estimation, formal errors
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
    concatenated_link_ends = np.array(simulated_observations.concatenated_link_ends)
    concatenated_link_end_names = simulated_observations.concatenated_link_end_names
    doppler_residuals = pod_output.residual_history
    observations_list = np.array(simulated_observations.concatenated_observations)

    # Compute how many DSN link ends there are
    DSN_concatenated_link_ends = list()
    for i in range(0,DSN_link_ends_number):
       DSN_concatenated_link_ends.extend(list(concatenated_link_ends[concatenated_link_ends==link_ends_numbers[i]]))
    print("DSN link ends:",DSN_link_ends_number," - DSN accumulated link ends:",len(DSN_concatenated_link_ends))
    
    concatenated_link_end_names_list = list([str(concatenated_link_end_names[i]) for i in range(0,len(concatenated_link_end_names))])

    # Save unsorted data
    np.savetxt(output_folder_path+"/weights_diagonal.dat",pod_output.weights_matrix_diagonal,fmt='%.15e')
    np.savetxt(output_folder_path+"/estimation_information_matrix.dat",estimation_information_matrix,fmt='%.15e')
    np.savetxt(output_folder_path+"/estimation_information_matrix_normalization.dat",
        estimation_information_matrix_normalization,fmt='%.15e')
    np.savetxt(output_folder_path+"/concatenated_times.dat",concatenated_times,fmt='%.15e')
    np.savetxt(output_folder_path+"/concatenated_link_ends.dat",concatenated_link_ends,fmt='%.15e')
    np.savetxt(output_folder_path+"/concatenated_link_end_names.dat",concatenated_link_end_names_list, fmt='%s')
    np.savetxt(output_folder_path+"/doppler_residuals.dat",doppler_residuals,fmt='%.15e')
    np.savetxt(output_folder_path+"/vector_weights.dat",vector_weights,fmt='%.15e')
    np.savetxt(output_folder_path+"/observations_list.dat",observations_list,fmt='%.15e')
    np.savetxt(output_folder_path+"/LaRa_time_solar_SEP.dat",np.rad2deg(LaRa_SEP_angle),fmt='%.15e')

    # Sort index
    index_sort = np.argsort(concatenated_times)

    # Function to sort concatenated times
    def sort_concatenated_times_func(index):
        return concatenated_times[index]

    # Sort concatenated times
    sort_concatenated_times_dict = Pool(CPU_par)
    print("Sorting concatenated times")
    concatenated_times_sort = sort_concatenated_times_dict.map(sort_concatenated_times_func,index_sort)
    sort_concatenated_times_dict.close()
    sort_concatenated_times_dict.join()

    # Remove duplicated concatenated times
    concatenated_times_no_duplicated = list(set(list(concatenated_times_sort)))
    concatenated_times_no_duplicated.sort()

    # Function to sort concatenated times index
    def sort_concatenated_times_index_func(index):
        return list(concatenated_times_sort).index(index)

    # Sort concatenated times index
    sort_concatenated_times_index_dict = Pool(CPU_par)
    print("Sorting time indices")
    concatenated_times_index = sort_concatenated_times_index_dict.map(sort_concatenated_times_index_func,concatenated_times_no_duplicated)
    sort_concatenated_times_index_dict.close()
    sort_concatenated_times_index_dict.join()

    # Function to sort concatenated times count
    def sort_concatenated_times_count_func(index):
        return list(concatenated_times_sort).count(index)

    # Sort concatenated times count
    sort_concatenated_times_count_dict = Pool(CPU_par)
    print("Sorting time counts")
    concatenated_times_count = sort_concatenated_times_count_dict.map(sort_concatenated_times_count_func,concatenated_times_no_duplicated)
    sort_concatenated_times_count_dict.close()
    sort_concatenated_times_count_dict.join()

    # Function to sort vector weights
    def sort_vector_weights_func(index):
        return vector_weights[index]

    # Sort concatenated times
    sort_vector_weights_dict = Pool(CPU_par)
    print("Sorting vector weigths")
    vector_weights_sort = sort_vector_weights_dict.map(sort_vector_weights_func,index_sort)
    sort_vector_weights_dict.close()
    sort_vector_weights_dict.join()

    # Function to sort estimation information matrix
    def sort_estimation_information_matrix_func(index):
        return estimation_information_matrix[index][6:]

    # Sort concatenated times
    sort_estimation_information_matrix_dict = Pool(CPU_par)
    print("Sorting estimation information matrices")
    estimation_information_matrix_sort = sort_estimation_information_matrix_dict.map(sort_estimation_information_matrix_func,index_sort)
    sort_estimation_information_matrix_dict.close()
    sort_estimation_information_matrix_dict.join()

    # Function to sort doppler residuals
    def sort_residuals_func(index):
        return doppler_residuals[index][:]

    # Sort concatenated residuals
    sort_residuals_dict = Pool(CPU_par)
    print("Sorting residuals")
    residuals_sort = sort_residuals_dict.map(sort_residuals_func,index_sort)
    sort_residuals_dict.close()
    sort_residuals_dict.join()

    # Function to sort concatenated link ends
    def sort_concatenated_link_ends_func(index):
        return concatenated_link_ends[index]

    # Sort concatenated link ends
    sort_concatenated_link_ends_dict = Pool(CPU_par)
    print("Sorting concatenated link ends")
    concatenated_link_ends_sort = sort_concatenated_link_ends_dict.map(sort_concatenated_link_ends_func,index_sort)
    sort_concatenated_link_ends_dict.close()
    sort_concatenated_link_ends_dict.join()

    # Function to sort concatenated link ends names
    def sort_concatenated_link_ends_names_func(index):
        return concatenated_link_end_names_list[index]

    # Sort concatenated link ends names
    sort_concatenated_link_ends_names_dict = Pool(CPU_par)
    print("Sorting concatenated link ends names")
    concatenated_link_end_names_list_sort = sort_concatenated_link_ends_names_dict.map(sort_concatenated_link_ends_names_func,index_sort)
    sort_concatenated_link_ends_names_dict.close()
    sort_concatenated_link_ends_names_dict.join()

    # Function to sort observations
    def sort_observations_func(index):
        return observations_list[index]

    # Sort cobservations
    sort_observations_dict = Pool(CPU_par)
    print("Sorting observations")
    observations_list_sort = sort_observations_dict.map(sort_observations_func,index_sort)
    sort_observations_dict.close()
    sort_observations_dict.join()

    # Initialize inverted weighting matrix
    inv_weight_complex = scipy.sparse.coo_matrix((0,0))

    # Save the correlation coefficients
    correlation_list = list()
    weather_list = list()
    indicator_Doppler_list = list()
    indicator_weather_list = list()
    Doppler_total_list = list()
    weather_total_list = list()
    correlation_dict = dict()
    LaRa_stations = ['DSS 43','DSS 63','DSS 14','MEDICINA','WETTZELL','ONSALA60','EFLSBERG','WRT0','YEBES40M','TIANMA65','CEDUNA','BADARY','HARTRAO']
    for i_station in LaRa_stations:
        correlation_dict[i_station] = dict()
        for j_station in LaRa_stations:
            if i_station!=j_station:
                correlation_dict[i_station][j_station] = list()
    
    # Delete arrays where two DSN have observation at the same time
    indices_delete = list()

    # Iteratate along each time value
    for time_index in range(0,len(concatenated_times_no_duplicated)):
        # Understand whether there are PRIDE stations (or two or more DSN link ends)
        count_time_obs = concatenated_times_count[time_index]
        start_index = concatenated_times_index[time_index]
        
        # If one observation, means that there are no PRIDE stations and the matrix added is a single value
        if count_time_obs == 1:
            inv_weight_complex = scipy.sparse.block_diag((inv_weight_complex,(vector_weights_sort[start_index])**(-2)))
        
        # More than observations = PRIDE stations present or two (or more) DSN observations at the same time
        else:
            end_index = start_index + (count_time_obs)
            # Sort the receiver link ends, since there can be several transmitters
            receiver_link_ends = np.argsort(concatenated_link_ends_sort[start_index:end_index])

            # Sorting again part of the arrays (with respect to the order of the link ends)
            estimation_information_matrix_sort_short = [None]*count_time_obs
            residuals_sort_short = [None]*count_time_obs
            for row_index in range(0,np.shape(estimation_information_matrix_sort_short)[0]):
                estimation_information_matrix_sort_short[row_index] = estimation_information_matrix_sort[start_index+receiver_link_ends[row_index]]
                residuals_sort_short[row_index] = residuals_sort[start_index+receiver_link_ends[row_index]]
            estimation_information_matrix_sort[start_index:end_index] = estimation_information_matrix_sort_short
            residuals_sort[start_index:end_index] = residuals_sort_short
            concatenated_link_end_names_list_sort[start_index:end_index] = np.array(concatenated_link_end_names_list_sort[start_index:end_index])[receiver_link_ends]
            vector_weights_sort[start_index:end_index] = np.array(vector_weights_sort[start_index:end_index])[receiver_link_ends]
            observations_list_sort[start_index:end_index] = np.array(observations_list_sort[start_index:end_index])[receiver_link_ends]
            concatenated_link_ends_sort[start_index:end_index] = np.array(concatenated_link_ends_sort[start_index:end_index])[receiver_link_ends]
            
            # Understanding which transmitters and receivers are involved
            transmitters_count = list()
            transmitters_total = list()
            receivers_total = list()
            for receiver_index in range(start_index,end_index):
                transmitter_station = link_ends_transmitter[link_ends_numbers.index(concatenated_link_ends_sort[receiver_index])]
                transmitters_total.append(transmitter_station)
                receivers_total.append(link_ends_receiver[link_ends_numbers.index(concatenated_link_ends_sort[receiver_index])])
                if not transmitter_station in transmitters_count:
                    transmitters_count.append(transmitter_station)
            
            # Counting how many receivers there are for each transmitter
            transmitter_count_number = list()    
            for transmitter_index in transmitters_count:
                transmitter_count_number.append(transmitters_total.count(transmitter_index))

            # Verification
            if len(concatenated_times_sort[start_index:end_index]) != concatenated_times_sort.count(concatenated_times_sort[start_index]):
                sys.exit()

            # Building the weighting matrix block
            start_index_block = start_index
            # When only we have the same transmitter
            if len(transmitter_count_number)==1:
                # When only the closed-loop observations are taken into account
                if remove_PRIDE_weight_boolean:
                    # Number of receivers for a transmitter
                    split_count = transmitter_count_number[0]
                    block_weight = (scipy.sparse.coo_matrix((1,1))).toarray()
                    for row_index in range(0,split_count):
                        # When transmitter = receiver -> add in the inv_weight
                        if transmitters_total[row_index]==receivers_total[row_index]:
                            block_weight[0][0] = vector_weights_sort[start_index_block+row_index]**2
                            inv_weight_complex = scipy.sparse.block_diag((inv_weight_complex,np.linalg.inv(block_weight)))
                        # When transmitter != receiver -> delete index
                        else:
                            indices_delete.append(start_index_block+row_index)

                # When all the PRIDE observations are taken into account
                else:
                    # Number of receivers for a transmitter
                    split_count = transmitter_count_number[0]
                    block_weight = (scipy.sparse.coo_matrix((split_count,split_count))).toarray()

                    # Compute the the distances from the Sun to Mars and Earth
                    r_earthsun = np.linalg.norm(spice_interface.get_body_cartesian_position_at_epoch("Earth","Sun","ECLIPJ2000","NONE",concatenated_times_sort[start_index]))
                    r_marssun = np.linalg.norm(spice_interface.get_body_cartesian_position_at_epoch("Mars","Sun","ECLIPJ2000","NONE",concatenated_times_sort[start_index]))
                    r_earthmars = np.linalg.norm(spice_interface.get_body_cartesian_position_at_epoch("Mars","Earth","ECLIPJ2000","NONE",concatenated_times_sort[start_index]))

                    # Cosine rule in order to compute the SEP angle
                    SEP_angle = np.arccos((r_earthsun**2+r_earthmars**2-r_marssun**2)/(2*r_earthsun*r_earthmars))

                    indicator_Doppler_noise = indicator_Doppler_noise_function(SEP_angle)

                    for row_index in range(0,split_count):
                        for column_index in range(0,split_count):
                            # Diagonal element
                            if row_index == column_index:
                                block_weight[row_index][column_index] = vector_weights_sort[start_index_block+row_index]**2
                            # Non-diagonal element
                            else:
                                station_1 = receivers_total[row_index]
                                station_2 = receivers_total[column_index]
                                correlation_coefficient = (1-weather_percentage_function(SEP_angle))*indicator_Doppler_noise[station_1][station_2]+weather_percentage_function(SEP_angle)*np.exp(-distance_stations[station_1][station_2]/(2*distance_ionospheric))
                                block_weight[row_index][column_index] = correlation_coefficient*vector_weights_sort[start_index_block+row_index]*vector_weights_sort[start_index_block+column_index]
                                if row_index>column_index:
                                    correlation_dict[station_1][station_2].append(correlation_coefficient)
                                    correlation_list.append(correlation_coefficient)
                                    weather_list.append(weather_percentage_function(SEP_angle))
                                    indicator_Doppler_list.append(indicator_Doppler_noise[station_1][station_2])
                                    indicator_weather_list.append(np.exp(-distance_stations[station_1][station_2]/(2*distance_ionospheric)))
                                    Doppler_total_list.append((1-weather_percentage_function(SEP_angle))*indicator_Doppler_noise[station_1][station_2])
                                    weather_total_list.append(weather_percentage_function(SEP_angle)*np.exp(-distance_stations[station_1][station_2]/(2*distance_ionospheric)))
                    inv_weight_complex = scipy.sparse.block_diag((inv_weight_complex,np.linalg.inv(block_weight)))

            # When different transmitters are available
            elif len(transmitter_count_number)>1:
                index_not_delete = 1 # Always choosing DSS 43 or 63 -> Reason: DSS 14 (Goldstone) there are no PRIDE stations close-by
                start_index_block = start_index
                total_split = 0
                for index_split_count in range(0,len(transmitter_count_number)):
                    # Number of receivers for a transmitter
                    split_count = transmitter_count_number[index_split_count]

                    # Only the transmitter with index_not_delete is taken into consideration
                    if index_split_count == index_not_delete:
                        # When only the closed-loop observations are taken into account
                        if remove_PRIDE_weight_boolean:
                            block_weight = (scipy.sparse.coo_matrix((1,1))).toarray()
                            for row_index in range(0,split_count):
                                # When transmitter = receiver -> add in the inv_weight
                                if transmitters_total[total_split+row_index]==receivers_total[total_split+row_index]:
                                    block_weight[0][0] = vector_weights_sort[start_index_block+row_index]**2
                                    inv_weight_complex = scipy.sparse.block_diag((inv_weight_complex,np.linalg.inv(block_weight)))
                                # When transmitter !=r receiver -> delete index
                                else:
                                    indices_delete.append(start_index_block+row_index)

                        # When all the PRIDE observations are taken into account
                        else:
                            block_weight = (scipy.sparse.coo_matrix((split_count,split_count))).toarray()

                            # Compute the the distances from the Sun to Mars and Earth
                            r_earthsun = np.linalg.norm(spice_interface.get_body_cartesian_position_at_epoch("Earth","Sun","ECLIPJ2000","NONE",concatenated_times_sort[start_index]))
                            r_marssun = np.linalg.norm(spice_interface.get_body_cartesian_position_at_epoch("Mars","Sun","ECLIPJ2000","NONE",concatenated_times_sort[start_index]))
                            r_earthmars = np.linalg.norm(spice_interface.get_body_cartesian_position_at_epoch("Mars","Earth","ECLIPJ2000","NONE",concatenated_times_sort[start_index]))

                            # Cosine rule in order to compute the SEP angle
                            SEP_angle = np.arccos((r_earthsun**2+r_earthmars**2-r_marssun**2)/(2*r_earthsun*r_earthmars))

                            indicator_Doppler_noise = indicator_Doppler_noise_function(SEP_angle)
                            for row_index in range(0,split_count):
                                for column_index in range(0,split_count):
                                    # Diagonal element
                                    if row_index == column_index:
                                        block_weight[row_index][column_index] = vector_weights_sort[start_index_block+row_index]**2
                                    # Non-diagonal element
                                    else:
                                        station_1 = receivers_total[total_split+row_index]
                                        station_2 = receivers_total[total_split+column_index]
                                        correlation_coefficient = (1-weather_percentage_function(SEP_angle))*indicator_Doppler_noise[station_1][station_2]+weather_percentage_function(SEP_angle)*np.exp(-distance_stations[station_1][station_2]/(2*distance_ionospheric))
                                        block_weight[row_index][column_index] = correlation_coefficient*vector_weights_sort[start_index_block+row_index]*vector_weights_sort[start_index_block+column_index]
                                        if row_index>column_index:
                                            correlation_dict[station_1][station_2].append(correlation_coefficient)
                                            correlation_list.append(correlation_coefficient)
                                            weather_list.append(weather_percentage_function(SEP_angle))
                                            indicator_Doppler_list.append(indicator_Doppler_noise[station_1][station_2])
                                            indicator_weather_list.append(np.exp(-distance_stations[station_1][station_2]/(2*distance_ionospheric)))
                                            Doppler_total_list.append((1-weather_percentage_function(SEP_angle))*indicator_Doppler_noise[station_1][station_2])
                                            weather_total_list.append(weather_percentage_function(SEP_angle)*np.exp(-distance_stations[station_1][station_2]/(2*distance_ionospheric)))
                            inv_weight_complex = scipy.sparse.block_diag((inv_weight_complex,np.linalg.inv(block_weight)))

                    # Remove all the observations from other transmitters
                    else:
                        indices_delete.extend(range(start_index_block,start_index_block+split_count))

                    # Increasing the start_index_block and total_split
                    start_index_block += split_count
                    total_split+=split_count
    
    # Delete terms
    for index_delete in indices_delete[::-1]:
        estimation_information_matrix_sort.pop(index_delete)
        residuals_sort.pop(index_delete)
        concatenated_times_sort.pop(index_delete)
        concatenated_link_ends_sort.pop(index_delete)
        concatenated_link_end_names_list_sort.pop(index_delete)
        vector_weights_sort.pop(index_delete)
        observations_list_sort.pop(index_delete)

    # Unit test for when remove_PRIDE_weight_boolean==TRUE, check whether the concatenated_times_sort is the same as concatenated_times_no_duplicated
    if remove_PRIDE_weight_boolean:
        if concatenated_times_sort != concatenated_times_no_duplicated:
                sys.exit()

    # Sort concatenated times index
    sort_concatenated_times_index_dict = Pool(CPU_par)
    print("Sorting time indices")
    concatenated_times_index = sort_concatenated_times_index_dict.map(sort_concatenated_times_index_func,concatenated_times_no_duplicated)
    sort_concatenated_times_index_dict.close()
    sort_concatenated_times_index_dict.join()

    # Function to sort concatenated times count
    def sort_concatenated_times_count_func(index):
        return list(concatenated_times_sort).count(index)

    # Sort concatenated times count
    sort_concatenated_times_count_dict = Pool(CPU_par)
    print("Sorting time counts")
    concatenated_times_count = sort_concatenated_times_count_dict.map(sort_concatenated_times_count_func,concatenated_times_no_duplicated)
    sort_concatenated_times_count_dict.close()
    sort_concatenated_times_count_dict.join()
    
    # Convert the inverted weighting matrix to Compressed Sparse Row format
    inv_weight_complex_total = (inv_weight_complex).tocsr()

    # Partial covariance
    if PRIDE_boolean==False:
        partial_cov = np.transpose(estimation_information_matrix_sort)@(inv_weight_complex_total.sqrt())

    # Save sorted data
    np.savetxt(output_folder_path+"/estimation_information_matrix_sort.dat",estimation_information_matrix_sort,fmt='%.15e')
    np.savetxt(output_folder_path+"/concatenated_times_sort.dat",concatenated_times_sort,fmt='%.15e')
    np.savetxt(output_folder_path+"/concatenated_times_sort_no_duplicated.dat",concatenated_times_no_duplicated,fmt='%.15e')
    np.savetxt(output_folder_path+"/concatenated_times_sort_index.dat",concatenated_times_index,fmt='%.15e')
    np.savetxt(output_folder_path+"/concatenated_times_sort_count.dat",concatenated_times_count,fmt='%.15e')
    np.savetxt(output_folder_path+"/concatenated_link_ends_sort.dat",concatenated_link_ends_sort,fmt='%.15e')
    np.savetxt(output_folder_path+"/concatenated_link_end_names_sort.dat",concatenated_link_end_names_list_sort, fmt='%s')
    np.savetxt(output_folder_path+"/doppler_residuals_sort.dat",residuals_sort,fmt='%.15e')
    np.savetxt(output_folder_path+"/vector_weights_sort.dat",vector_weights_sort,fmt='%.15e')
    np.savetxt(output_folder_path+"/observations_list_sort.dat",observations_list_sort,fmt='%.15e')
    np.savetxt(output_folder_path+"/correlation_list.dat",correlation_list,fmt='%.15e')
    np.savetxt(output_folder_path+"/weather_list.dat",weather_list,fmt='%.15e')
    np.savetxt(output_folder_path+"/indicator_Doppler_list.dat",indicator_Doppler_list,fmt='%.15e')
    np.savetxt(output_folder_path+"/indicator_weather_list.dat",indicator_weather_list,fmt='%.15e')
    np.savetxt(output_folder_path+"/Doppler_total_list.dat",Doppler_total_list,fmt='%.15e')
    np.savetxt(output_folder_path+"/weather_total_list.dat",weather_total_list,fmt='%.15e')
    if PRIDE_boolean==False:
        np.savetxt(output_folder_path+"/partial_cov.dat",partial_cov,fmt='%.15e')

    import json
    with open(output_folder_path+'/correlations.json', 'w') as fp:
            json.dump(correlation_dict, fp)  

    # Histogram correlation
    
    plt.figure(figsize=(8, 5))
    plt.hist(correlation_list, bins = 20)
    average_correlation = np.mean(correlation_list)
    plt.axvline(x=average_correlation, color='k', linestyle='--',label='Mean = '+str(average_correlation))
    plt.ylabel('Frequency [-]')
    plt.xlabel('Correlation [-]')
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.grid()
    plt.legend()
    plt.savefig(output_folder_path+"/correlation.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(8, 5))
    plt.hist(weather_list, bins = 20)
    plt.ylabel('Frequency [-]')
    plt.xlabel('Weather Coeff [-]')
    plt.grid()
    plt.savefig(output_folder_path+"/weather_list.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(8, 5))
    plt.hist(indicator_weather_list, bins = 20)
    plt.ylabel('Frequency [-]')
    plt.xlabel('Weather Corr [-]')
    plt.grid()
    plt.savefig(output_folder_path+"/indicator_weather_list.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(8, 5))
    plt.hist(indicator_Doppler_list, bins = 20)
    plt.ylabel('Frequency [-]')
    plt.xlabel('Doppler Coeff [-]')
    plt.grid()
    plt.savefig(output_folder_path+"/indicator_Doppler_list.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    ########################################################################################################################
    ################################################## COVARIANCE ANALYSIS #################################################
    ########################################################################################################################

    # Step evaluation
    arange_eval = np.arange(0,len(concatenated_times_no_duplicated),step_eval) 
    time_eval = [concatenated_times_no_duplicated[i] for i in arange_eval]

    np.savetxt(output_folder_path+"/time_plot.dat",time_eval,fmt='%.15e')
    
    # Normalized inverse a priori covariance
    norm_inverse_a_priori_covariance = np.diag(inverse_a_priori_covariance.diagonal()[6:]/(estimation_information_matrix_normalization[6:]**2))

    # Function for the normalized covariance matrix
    def norm_covariance_matrix_func(time_index):
        return np.transpose(estimation_information_matrix_sort[:concatenated_times_index[time_index]+1])@inv_weight_complex_total[:concatenated_times_index[time_index]+1,:concatenated_times_index[time_index]+1]@estimation_information_matrix_sort[:concatenated_times_index[time_index]+1]\
            +norm_inverse_a_priori_covariance

    # Compute the normalized covariance matrix using several CPUs
    norm_covariance_matrix_dict = Pool(CPU_par)
    print("Calculating the inverse of W and the normalized covariance matrix")
    norm_covariance_values = norm_covariance_matrix_dict.map(norm_covariance_matrix_func,arange_eval)  
    norm_covariance_matrix_dict.close()
    norm_covariance_matrix_dict.join()

    # Function for the invert the covariance matrix
    def inv_covariance_matrix_func(time_index):
        return np.linalg.inv(norm_covariance_values[time_index])

    # Compute the inversion of the covariance matrix using several CPUs
    inv_covariance_matrix_dict = Pool(CPU_par)
    print("Calculating the inverse normalized covariance matrix")
    inv_covariance_matrix_values = inv_covariance_matrix_dict.map(inv_covariance_matrix_func,range(0,len(arange_eval)))  
    inv_covariance_matrix_dict.close()
    inv_covariance_matrix_dict.join()

    # Function for the covariance matrix
    def covariance_matrix_func(time_index):
        covariance_matrix = np.zeros(np.shape(inv_covariance_matrix_values[time_index]))
        for i in range(0,np.shape(inv_covariance_matrix_values[time_index])[0]):
            for j in range(0,np.shape(inv_covariance_matrix_values[time_index])[1]):
                covariance_matrix[i][j] = inv_covariance_matrix_values[time_index][i][j]/\
                   (estimation_information_matrix_normalization[i+6]*estimation_information_matrix_normalization[j+6])
        return covariance_matrix 

    # Compute the unnormalized covariance matrix using several CPUs
    covariance_matrix_dict = Pool(CPU_par)
    print("Calculating the unnormalized covariance matrix")
    covariance_values = covariance_matrix_dict.map(covariance_matrix_func,range(0,len(arange_eval)))
    covariance_matrix_dict.close()
    covariance_matrix_dict.join()

    # Function for the standard deviation (formal)
    def sigma_covariance_matrix_func(time_index):
        return np.sqrt(covariance_values[time_index].diagonal())

    # Take the standard deviation (formal) from the diagonal of the unnormalized covariance matrix using several CPUs
    sigma_covariance_matrix_dict = Pool(CPU_par)
    print("Calculating the standard deviation (formal)")
    sigma_values = sigma_covariance_matrix_dict.map(sigma_covariance_matrix_func,range(0,len(arange_eval)))
    sigma_covariance_matrix_dict.close()
    sigma_covariance_matrix_dict.join()

    # Verification
    #print("A:",inverse_a_priori_covariance.diagonal()**(-0.5))
    #print("F:",sigma_values[1])

    # Function for the correlation matrix
    def correlation_matrix_func(time_index):
        correlation_matrix = np.zeros(np.shape(covariance_values[time_index]))
        for i in range(0,np.shape(covariance_values[time_index])[0]):
            for j in range(0,np.shape(covariance_values[time_index])[1]):
                correlation_matrix[i][j] = covariance_values[time_index][i][j]/\
                    (sigma_values[time_index][i]*sigma_values[time_index][j])
        return correlation_matrix

    # Compute correlation matrix using 4 CPUs
    correlation_matrix_dict = Pool(CPU_par)
    print("Calculating the correlation matrix")
    correlation_values = correlation_matrix_dict.map(correlation_matrix_func,range(0,len(arange_eval)))
    correlation_matrix_dict.close()
    correlation_matrix_dict.join()

    # Compute correlation matrix only for LaRa mission
    if len(LaRa_observation_times_list)!=0:
        index_start_LaRa = list(concatenated_times_sort).index(min(LaRa_concatenated_times))
        
        inv_norm_covariance_LaRa_value = np.linalg.inv(np.transpose(estimation_information_matrix_sort[index_start_LaRa:])@inv_weight_complex_total[index_start_LaRa:,index_start_LaRa:]@estimation_information_matrix_sort[index_start_LaRa:]\
                +norm_inverse_a_priori_covariance)

        covariance_matrix_LaRa = np.zeros(np.shape(inv_norm_covariance_LaRa_value))
        for i in range(0,np.shape(inv_norm_covariance_LaRa_value)[0]):
            for j in range(0,np.shape(inv_norm_covariance_LaRa_value)[1]):
                covariance_matrix_LaRa[i][j] = inv_norm_covariance_LaRa_value[i][j]/\
                    (estimation_information_matrix_normalization[i+6]*estimation_information_matrix_normalization[j+6])

        sigma_LaRa_value = np.sqrt(covariance_matrix_LaRa.diagonal())

        correlation_matrix_LaRa = np.zeros(np.shape(covariance_matrix_LaRa))
        for i in range(0,np.shape(covariance_matrix_LaRa)[0]):
            for j in range(0,np.shape(covariance_matrix_LaRa)[1]):
                correlation_matrix_LaRa[i][j] =covariance_matrix_LaRa[i][j]/\
                    (sigma_LaRa_value[i]*sigma_LaRa_value[j])

    ########################################################################################################################
    ################################################## PLOTS ###############################################################
    ########################################################################################################################

    # True to form ratio histogram
    plt.figure(figsize=(15, 6))
    plt.hist(np.abs(true_to_form_estimation_error_ratio), bins = 8)
    plt.ylabel('Frequency [-]')
    plt.xlabel('True to form ratio [-]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.grid()
    plt.savefig(output_folder_path+"/true_to_form_ratio.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # Plot to check the viability of the Sun
    plt.figure(figsize=(8, 4))
    plt.rcParams.update({'font.size': 18})
    if len(RISE_observation_times_list)!=0: 
        plt.scatter((RISE_concatenated_times-observation_start_epoch)/constants.JULIAN_DAY,RISE_std_mHz_function((RISE_concatenated_times-RISE_observation_start_epoch_reference_noise)/constants.JULIAN_DAY))
    if len(LaRa_observation_times_list)!=0: 
        plt.scatter((LaRa_concatenated_times-observation_start_epoch)/constants.JULIAN_DAY,LaRa_std_noise_function(LaRa_concatenated_times)/10**(-3)*base_frequency)
    plt.ylabel(r'1-$\sigma$ Doppler Residual [mHz]')
    plt.xlabel('Time [days]')
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.grid()
    plt.savefig(output_folder_path+"/std_noise_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all') 
    plt.rcParams.update({'font.size': 10})

    # Plot to check the viability of the Sun v2
    plt.figure(figsize=(15, 6))
    plt.scatter((concatenated_times_sort-observation_start_epoch*np.ones(len(concatenated_times_sort)))/constants.JULIAN_DAY,vector_weights_sort)
    plt.ylabel('Std noise [-]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.grid()
    plt.savefig(output_folder_path+"/weights_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # Plot observations as a function of time
    plt.figure(figsize=(15,6))
    plt.scatter((concatenated_times_sort-observation_start_epoch*np.ones(len(concatenated_times_sort)))/constants.JULIAN_DAY,np.array(observations_list_sort)*constants.SPEED_OF_LIGHT)
    plt.ylabel('Range rate [m/s]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.grid()
    plt.savefig(output_folder_path+"/observations_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close("all")

    # Plot residual as a function of time
    plt.figure(figsize=(15,6))
    plt.scatter((concatenated_times_sort-observation_start_epoch*np.ones(len(concatenated_times_sort)))/constants.JULIAN_DAY,residuals_sort)
    plt.ylabel('Residual [-]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.grid()
    plt.savefig(output_folder_path+"/residuals_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close("all")
    
    np.savetxt(output_folder_path+"/final_mop_correlation.dat",correlation_values[-1],fmt='%.15e')
    
    end_RISE_observation_time = max(RISE_observation_times_list)
    idx_nearest = np.abs(end_RISE_observation_time-np.array(time_eval)).argmin()
    np.savetxt(output_folder_path+"/RISE_mop_correlation.dat",correlation_values[idx_nearest],fmt='%.15e')
        
    np.savetxt(output_folder_path+"/LaRa_mop_correlation.dat",correlation_matrix_LaRa,fmt='%.15e')

    # Formal to apriori ratio
    plt.figure(figsize=(15, 6))
    plt.plot(range(0,len(apriori_vector)-6),np.abs(sigma_values[-1][:]/apriori_vector[6:]),'o--')
    plt.ylabel('Formal to Apriori Ratio')
    plt.xlabel('Estimated Parameters')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    #'''
    if len(LaRa_observation_times_list)==0:
        plt.xticks(range(0,len(apriori_vector)-6),labels=['F',r'$\sigma_{FCN}$',r'$x_{{RISE}}$',r'$y_{{RISE}}$',r'$z_{{RISE}}$',
            r'$\phi^c_1$',r'$\phi^s_1$',r'$\phi^c_2$',r'$\phi^s_2$',r'$\phi^c_3$',r'$\phi^s_3$',r'$\phi^c_4$',r'$\phi^s_4$',
            r'$Xp^c_1$',r'$Xp^s_1$',r'$Yp^c_1$',r'$Yp^s_1$',r'$Xp^c_2$',r'$Xp^s_2$',r'$Yp^c_2$',r'$Yp^s_2$',
            r'$Xp^c_3$',r'$Xp^s_3$',r'$Yp^c_3$',r'$Yp^s_3$',r'$Xp^c_4$',r'$Xp^s_4$',r'$Yp^c_4$',r'$Yp^s_4$',
            r'$Xp^c_5$',r'$Xp^s_5$',r'$Yp^c_5$',r'$Yp^s_5$'])
    elif len(RISE_observation_times_list)==0:
        plt.xticks(range(0,len(apriori_vector)-6),labels=['F',r'$\sigma_{FCN}$',r'$x_{{LaRa}}$',r'$y_{{LaRa}}$',r'$z_{{LaRa}}$',
            r'$\phi^c_1$',r'$\phi^s_1$',r'$\phi^c_2$',r'$\phi^s_2$',r'$\phi^c_3$',r'$\phi^s_3$',r'$\phi^c_4$',r'$\phi^s_4$',
            r'$Xp^c_1$',r'$Xp^s_1$',r'$Yp^c_1$',r'$Yp^s_1$',r'$Xp^c_2$',r'$Xp^s_2$',r'$Yp^c_2$',r'$Yp^s_2$',
            r'$Xp^c_3$',r'$Xp^s_3$',r'$Yp^c_3$',r'$Yp^s_3$',r'$Xp^c_4$',r'$Xp^s_4$',r'$Yp^c_4$',r'$Yp^s_4$',
            r'$Xp^c_5$',r'$Xp^s_5$',r'$Yp^c_5$',r'$Yp^s_5$'])
    else:
    #'''
        plt.xticks(range(0,len(apriori_vector)-6),labels=['F',r'$\sigma_{FCN}$',r'$x_{{RISE}}$',r'$y_{{RISE}}$',r'$z_{{RISE}}$',r'$x_{{LaRa}}$',r'$y_{{LaRa}}$',r'$z_{{LaRa}}$',
            r'$\phi^c_1$',r'$\phi^s_1$',r'$\phi^c_2$',r'$\phi^s_2$',r'$\phi^c_3$',r'$\phi^s_3$',r'$\phi^c_4$',r'$\phi^s_4$',
            r'$Xp^c_1$',r'$Xp^s_1$',r'$Yp^c_1$',r'$Yp^s_1$',r'$Xp^c_2$',r'$Xp^s_2$',r'$Yp^c_2$',r'$Yp^s_2$',
            r'$Xp^c_3$',r'$Xp^s_3$',r'$Yp^c_3$',r'$Yp^s_3$',r'$Xp^c_4$',r'$Xp^s_4$',r'$Yp^c_4$',r'$Yp^s_4$',
            r'$Xp^c_5$',r'$Xp^s_5$',r'$Yp^c_5$',r'$Yp^s_5$'])
    plt.grid()
    plt.savefig(output_folder_path+"/formal_to_a_priori.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # Formal to apriori ratio simple
    plt.figure(figsize=(15, 6))
    plt.plot(range(0,len(apriori_vector)-6),np.abs(sigma_values[-1][:]/apriori_vector[6:]),'o--')
    plt.ylabel('Formal to Apriori Ratio')
    plt.xlabel('Estimated Parameters')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    #'''
    if len(LaRa_observation_times_list)==0:
        plt.xticks(range(0,len(apriori_vector)-6),labels=['F',r'$\sigma_{FCN}$',r'$x_{{RISE}}$',r'$y_{{RISE}}$',r'$z_{{RISE}}$',
            r'$\phi^c_1$',r'$\phi^s_1$',r'$\phi^c_2$',r'$\phi^s_2$',r'$\phi^c_3$',r'$\phi^s_3$',r'$\phi^c_4$',r'$\phi^s_4$',
            r'$Xp^c_1$',r'$Xp^s_1$',r'$Yp^c_1$',r'$Yp^s_1$',r'$Xp^c_2$',r'$Xp^s_2$',r'$Yp^c_2$',r'$Yp^s_2$',
            r'$Xp^c_3$',r'$Xp^s_3$',r'$Yp^c_3$',r'$Yp^s_3$',r'$Xp^c_4$',r'$Xp^s_4$',r'$Yp^c_4$',r'$Yp^s_4$',
            r'$Xp^c_5$',r'$Xp^s_5$',r'$Yp^c_5$',r'$Yp^s_5$'])
    elif len(RISE_observation_times_list)==0:
        plt.xticks(range(0,len(apriori_vector)-6),labels=['F',r'$\sigma_{FCN}$',r'$x_{{LaRa}}$',r'$y_{{LaRa}}$',r'$z_{{LaRa}}$',
            r'$\phi^c_1$',r'$\phi^s_1$',r'$\phi^c_2$',r'$\phi^s_2$',r'$\phi^c_3$',r'$\phi^s_3$',r'$\phi^c_4$',r'$\phi^s_4$',
            r'$Xp^c_1$',r'$Xp^s_1$',r'$Yp^c_1$',r'$Yp^s_1$',r'$Xp^c_2$',r'$Xp^s_2$',r'$Yp^c_2$',r'$Yp^s_2$',
            r'$Xp^c_3$',r'$Xp^s_3$',r'$Yp^c_3$',r'$Yp^s_3$',r'$Xp^c_4$',r'$Xp^s_4$',r'$Yp^c_4$',r'$Yp^s_4$',
            r'$Xp^c_5$',r'$Xp^s_5$',r'$Yp^c_5$',r'$Yp^s_5$'])
    else:
    #'''
        plt.xticks(range(0,len(apriori_vector)-6),labels=['F',r'$\sigma_{FCN}$',r'$x_{{RISE}}$',r'$y_{{RISE}}$',r'$z_{{RISE}}$',r'$x_{{LaRa}}$',r'$y_{{LaRa}}$',r'$z_{{LaRa}}$',
            r'$\phi^c_1$',r'$\phi^s_1$',r'$\phi^c_2$',r'$\phi^s_2$',r'$\phi^c_3$',r'$\phi^s_3$',r'$\phi^c_4$',r'$\phi^s_4$',
            r'$Xp^c_1$',r'$Xp^s_1$',r'$Yp^c_1$',r'$Yp^s_1$',r'$Xp^c_2$',r'$Xp^s_2$',r'$Yp^c_2$',r'$Yp^s_2$',
            r'$Xp^c_3$',r'$Xp^s_3$',r'$Yp^c_3$',r'$Yp^s_3$',r'$Xp^c_4$',r'$Xp^s_4$',r'$Yp^c_4$',r'$Yp^s_4$',
            r'$Xp^c_5$',r'$Xp^s_5$',r'$Yp^c_5$',r'$Yp^s_5$'])
    plt.grid()
    plt.savefig(output_folder_path+"/formal_to_a_priori_simple.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # 1-sigma F as a function of time
    plt.figure(figsize=(15, 6))
    F_values = list()
    for time_index in range(0,len(time_eval)):
        F_values.append(sigma_values[time_index][0])
    np.savetxt(output_folder_path+"/corefactor_plot.dat",F_values,fmt='%.15e')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        F_values,'-o',label='F')
    if RISE_boolean and LaRa_boolean:
        plt.axvline(x=(LaRa_observation_times_list[0]-observation_start_epoch)/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ F [-]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.grid()
    plt.legend()
    plt.savefig(output_folder_path+"/Fvalues_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # 1-sigma sigma_FCN as a function of time
    plt.figure(figsize=(15, 6))
    sigma_FCN_values = list()
    for time_index in range(0,len(time_eval)):
        sigma_FCN_values.append(sigma_values[time_index][1])
    np.savetxt(output_folder_path+"/sigmaFCN_plot.dat",sigma_FCN_values,fmt='%.15e')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        sigma_FCN_values,'-o',label=r'$\sigma_{FCN}$')
    if RISE_boolean and LaRa_boolean:
        plt.axvline(x=(LaRa_observation_times_list[0]-observation_start_epoch)/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $\sigma_{FCN}$ [rad/s]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.grid()
    plt.legend()
    plt.savefig(output_folder_path+"/sigmaFCNvalues_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # 1-sigma x_RISE,y_RISE,z_RISE,x_LaRa,y_LaRa,z_LaRa as a function of time
    plt.figure(figsize=(15, 6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(apriori_vector[8:11+add_par])))))
    if len(RISE_observation_times_list)!=0:
        xRISElander_values = list()
        yRISElander_values = list()
        zRISElander_values = list()
    if len(LaRa_observation_times_list)!=0:
        xLaRalander_values = list()
        yLaRalander_values = list()
        zLaRalander_values = list()
    for time_index in range(0,len(time_eval)):
        if len(RISE_observation_times_list)!=0:
            xRISElander_values.append(sigma_values[time_index][2])
            yRISElander_values.append(sigma_values[time_index][3])
            zRISElander_values.append(sigma_values[time_index][4])
        
        if len(LaRa_observation_times_list)!=0:
            xLaRalander_values.append(sigma_values[time_index][2+add_par])
            yLaRalander_values.append(sigma_values[time_index][3+add_par])
            zLaRalander_values.append(sigma_values[time_index][4+add_par])

    if len(RISE_observation_times_list)!=0:
        np.savetxt(output_folder_path+"/xRISE_plot.dat",xRISElander_values,fmt='%.15e')
        np.savetxt(output_folder_path+"/yRISE_plot.dat",yRISElander_values,fmt='%.15e')
        np.savetxt(output_folder_path+"/zRISE_plot.dat",zRISElander_values,fmt='%.15e')        
        plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
            xRISElander_values,'-o',label='$x_{RISE}$')
        plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
            yRISElander_values,'-o',label='$y_{RISE}$')
        plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
            zRISElander_values,'-o',label='$z_{RISE}$')
    if len(LaRa_observation_times_list)!=0:
        np.savetxt(output_folder_path+"/xLaRa_plot.dat",xLaRalander_values,fmt='%.15e')
        np.savetxt(output_folder_path+"/yLaRa_plot.dat",yLaRalander_values,fmt='%.15e')
        np.savetxt(output_folder_path+"/zLaRa_plot.dat",zLaRalander_values,fmt='%.15e')
        plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
            xLaRalander_values,'-o',label='$x_{LaRa}$')
        plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
            yLaRalander_values,'-o',label='$y_{LaRa}$')
        plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
            zLaRalander_values,'-o',label='$z_{LaRa}$')
    if RISE_boolean and LaRa_boolean:
        plt.axvline(x=(LaRa_observation_times_list[0]-observation_start_epoch)/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ x,y,z [m]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.legend()
    plt.grid()
    plt.savefig(output_folder_path+"/xyzlander_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # 1-sigma spin variations as a function of time
    plt.figure(figsize=(15, 6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(apriori_vector[11+add_par:19+add_par])))))
    cos1spin_values = list()
    sin1spin_values = list()
    cos2spin_values = list()
    sin2spin_values = list()
    cos3spin_values = list()
    sin3spin_values = list()
    cos4spin_values = list()
    sin4spin_values = list()
    for time_index in range(0,len(time_eval)):
        cos1spin_values.append(sigma_values[time_index][5+add_par])
        sin1spin_values.append(sigma_values[time_index][6+add_par])
        cos2spin_values.append(sigma_values[time_index][7+add_par])
        sin2spin_values.append(sigma_values[time_index][8+add_par])
        cos3spin_values.append(sigma_values[time_index][9+add_par])
        sin3spin_values.append(sigma_values[time_index][10+add_par])
        cos4spin_values.append(sigma_values[time_index][11+add_par])
        sin4spin_values.append(sigma_values[time_index][12+add_par])
    np.savetxt(output_folder_path+"/cos1spin_plot.dat",cos1spin_values,fmt='%.15e')
    np.savetxt(output_folder_path+"/sin1spin_plot.dat",sin1spin_values,fmt='%.15e')
    np.savetxt(output_folder_path+"/cos2spin_plot.dat",cos2spin_values,fmt='%.15e')
    np.savetxt(output_folder_path+"/sin2spin_plot.dat",sin2spin_values,fmt='%.15e')
    np.savetxt(output_folder_path+"/cos3spin_plot.dat",cos3spin_values,fmt='%.15e')
    np.savetxt(output_folder_path+"/sin3spin_plot.dat",sin3spin_values,fmt='%.15e')
    np.savetxt(output_folder_path+"/cos4spin_plot.dat",cos4spin_values,fmt='%.15e')
    np.savetxt(output_folder_path+"/sin4spin_plot.dat",sin4spin_values,fmt='%.15e')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(cos1spin_values)/mas,'-o',label=r'$\phi^c_1$')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(sin1spin_values)/mas,'-o',label=r'$\phi^s_1$')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(cos2spin_values)/mas,'-o',label=r'$\phi^c_2$')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(sin2spin_values)/mas,'-o',label=r'$\phi^s_2$')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(cos3spin_values)/mas,'-o',label=r'$\phi^c_3$')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(sin3spin_values)/mas,'-o',label=r'$\phi^s_3$')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(cos4spin_values)/mas,'-o',label=r'$\phi^c_4$')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(sin4spin_values)/mas,'-o',label=r'$\phi^s_4$')
    if RISE_boolean and LaRa_boolean:
        plt.axvline(x=(LaRa_observation_times_list[0]-observation_start_epoch)/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $\phi$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.legend()
    plt.grid()
    plt.savefig(output_folder_path+"/psispin_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # 1-sigma polar motion (order 1) as a function of time
    plt.figure(figsize=(15, 6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(apriori_vector[19+add_par:23+add_par])))))
    xpcos1_values = list()
    xpsin1_values = list()
    ypcos1_values = list()
    ypsin1_values = list()
    for time_index in range(0,len(time_eval)):
        xpcos1_values.append(sigma_values[time_index][13+add_par])
        xpsin1_values.append(sigma_values[time_index][14+add_par])
        ypcos1_values.append(sigma_values[time_index][15+add_par])
        ypsin1_values.append(sigma_values[time_index][16+add_par])
    np.savetxt(output_folder_path+"/xpcos1_plot.dat",xpcos1_values,fmt='%.15e')
    np.savetxt(output_folder_path+"/xpsin1_plot.dat",xpsin1_values,fmt='%.15e')
    np.savetxt(output_folder_path+"/ypcos1_plot.dat",ypcos1_values,fmt='%.15e')
    np.savetxt(output_folder_path+"/ypsin1_plot.dat",ypsin1_values,fmt='%.15e')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(xpcos1_values)/mas,'-o',label=r'$Xp^c_1$')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(xpsin1_values)/mas,'-o',label=r'$Xp^s_1$')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(ypcos1_values)/mas,'-o',label=r'$Yp^c_1$')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(ypsin1_values)/mas,'-o',label=r'$Yp^s_1$')
    if RISE_boolean and LaRa_boolean:
        plt.axvline(x=(LaRa_observation_times_list[0]-observation_start_epoch)/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $Xp, Yp$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.legend()
    plt.grid()
    plt.savefig(output_folder_path+"/polarmotionamp1_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # 1-sigma polar motion (order 2) as a function of time
    plt.figure(figsize=(15, 6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(apriori_vector[23+add_par:27+add_par])))))
    xpcos2_values = list()
    xpsin2_values = list()
    ypcos2_values = list()
    ypsin2_values = list()
    for time_index in range(0,len(time_eval)):
        xpcos2_values.append(sigma_values[time_index][17+add_par])
        xpsin2_values.append(sigma_values[time_index][18+add_par])
        ypcos2_values.append(sigma_values[time_index][19+add_par])
        ypsin2_values.append(sigma_values[time_index][20+add_par])
    np.savetxt(output_folder_path+"/xpcos2_plot.dat",xpcos2_values,fmt='%.15e')
    np.savetxt(output_folder_path+"/xpsin2_plot.dat",xpsin2_values,fmt='%.15e')
    np.savetxt(output_folder_path+"/ypcos2_plot.dat",ypcos2_values,fmt='%.15e')
    np.savetxt(output_folder_path+"/ypsin2_plot.dat",ypsin2_values,fmt='%.15e')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(xpcos2_values)/mas,'-o',label=r'$Xp^c_2$')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(xpsin2_values)/mas,'-o',label=r'$Xp^s_2$')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(ypcos2_values)/mas,'-o',label=r'$Yp^c_2$')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(ypsin2_values)/mas,'-o',label=r'$Yp^s_2$')
    if RISE_boolean and LaRa_boolean:
        plt.axvline(x=(LaRa_observation_times_list[0]-observation_start_epoch)/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $Xp, Yp$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.legend()
    plt.grid()
    plt.savefig(output_folder_path+"/polarmotionamp2_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # 1-sigma polar motion (order 3) as a function of time
    plt.figure(figsize=(15, 6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(apriori_vector[27+add_par:31+add_par])))))
    xpcos3_values = list()
    xpsin3_values = list()
    ypcos3_values = list()
    ypsin3_values = list()
    for time_index in range(0,len(time_eval)):
        xpcos3_values.append(sigma_values[time_index][21+add_par])
        xpsin3_values.append(sigma_values[time_index][22+add_par])
        ypcos3_values.append(sigma_values[time_index][23+add_par])
        ypsin3_values.append(sigma_values[time_index][24+add_par])
    np.savetxt(output_folder_path+"/xpcos3_plot.dat",xpcos3_values,fmt='%.15e')
    np.savetxt(output_folder_path+"/xpsin3_plot.dat",xpsin3_values,fmt='%.15e')
    np.savetxt(output_folder_path+"/ypcos3_plot.dat",ypcos3_values,fmt='%.15e')
    np.savetxt(output_folder_path+"/ypsin3_plot.dat",ypsin3_values,fmt='%.15e') 
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(xpcos3_values)/mas,'-o',label=r'$Xp^c_3$')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(xpsin3_values)/mas,'-o',label=r'$Xp^s_3$')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(ypcos3_values)/mas,'-o',label=r'$Yp^c_3$')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(ypsin3_values)/mas,'-o',label=r'$Yp^s_3$')
    if RISE_boolean and LaRa_boolean:
        plt.axvline(x=(LaRa_observation_times_list[0]-observation_start_epoch)/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $Xp, Yp$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.legend()
    plt.grid()
    plt.savefig(output_folder_path+"/polarmotionamp3_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # 1-sigma polar motion (order 4) as a function of time
    plt.figure(figsize=(15, 6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(apriori_vector[31+add_par:35+add_par])))))
    xpcos4_values = list()
    xpsin4_values = list()
    ypcos4_values = list()
    ypsin4_values = list()
    for time_index in range(0,len(time_eval)):
        xpcos4_values.append(sigma_values[time_index][25+add_par])
        xpsin4_values.append(sigma_values[time_index][26+add_par])
        ypcos4_values.append(sigma_values[time_index][27+add_par])
        ypsin4_values.append(sigma_values[time_index][28+add_par]) 
    np.savetxt(output_folder_path+"/xpcos4_plot.dat",xpcos4_values,fmt='%.15e')
    np.savetxt(output_folder_path+"/xpsin4_plot.dat",xpsin4_values,fmt='%.15e')
    np.savetxt(output_folder_path+"/ypcos4_plot.dat",ypcos4_values,fmt='%.15e')
    np.savetxt(output_folder_path+"/ypsin4_plot.dat",ypsin4_values,fmt='%.15e')   
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(xpcos4_values)/mas,'-o',label=r'$Xp^c_4$')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(xpsin4_values)/mas,'-o',label=r'$Xp^s_4$')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(ypcos4_values)/mas,'-o',label=r'$Yp^c_4$')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(ypsin4_values)/mas,'-o',label=r'$Yp^s_4$')
    if RISE_boolean and LaRa_boolean:
        plt.axvline(x=(LaRa_observation_times_list[0]-observation_start_epoch)/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $Xp, Yp$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.legend()
    plt.grid()
    plt.savefig(output_folder_path+"/polarmotionamp4_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # 1-sigma polar motion (order 5) as a function of time
    plt.figure(figsize=(15, 6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(apriori_vector[35+add_par:])))))
    xpcos5_values = list()
    xpsin5_values = list()
    ypcos5_values = list()
    ypsin5_values = list()
    for time_index in range(0,len(time_eval)):
        xpcos5_values.append(sigma_values[time_index][29+add_par])
        xpsin5_values.append(sigma_values[time_index][30+add_par])
        ypcos5_values.append(sigma_values[time_index][31+add_par])
        ypsin5_values.append(sigma_values[time_index][32+add_par])
    np.savetxt(output_folder_path+"/xpcos5_plot.dat",xpcos5_values,fmt='%.15e')
    np.savetxt(output_folder_path+"/xpsin5_plot.dat",xpsin5_values,fmt='%.15e')
    np.savetxt(output_folder_path+"/ypcos5_plot.dat",ypcos5_values,fmt='%.15e')
    np.savetxt(output_folder_path+"/ypsin5_plot.dat",ypsin5_values,fmt='%.15e')  
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(xpcos5_values)/mas,'-o',label=r'$Xp^c_5$')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(xpsin5_values)/mas,'-o',label=r'$Xp^s_5$')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(ypcos5_values)/mas,'-o',label=r'$Yp^c_5$')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(ypsin5_values)/mas,'-o',label=r'$Yp^s_5$')
    if RISE_boolean and LaRa_boolean:
        plt.axvline(x=(LaRa_observation_times_list[0]-observation_start_epoch)/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $Xp, Yp$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.legend()
    plt.grid()
    plt.savefig(output_folder_path+"/polarmotionamp5_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    np.savetxt(output_folder_path+"/final_formalerrors_plot.dat",sigma_values[-1][:],fmt='%.15e')

    print('Is there any duplicated total time value? :',any(list(concatenated_times).count(x) > 1 for x in list(concatenated_times)))

print("--- %s seconds ---" % (time.time() - run_time))