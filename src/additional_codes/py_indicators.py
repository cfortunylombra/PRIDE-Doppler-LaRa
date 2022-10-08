"""
Description: Save the data of the two indicators (cross-correlations)

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
    transmitters_dict['DSS 63']=np.array([4849092.5175,-360180.3480,4115109.2506]) # DSS 63
    transmitters_dict['DSS 14']=np.array([-2353621.4197,-4641341.4717,3677052.3178]) # DSS 14

    # Earth-based transmitter for LaRa
    LaRa_transmitter_names = ['DSS 43','DSS 63','DSS 14']

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
            if name_ground_station=="HART15M" or name_ground_station=="HOBART12" or name_ground_station=="WARK30M" or name_ground_station=="HOBART26":
                continue
            else:
                x_coordinate_ground_station = float(coordinates_line_ground_station.split("X=",1)[1].split()[0])
                y_coordinate_ground_station = float(coordinates_line_ground_station.split("Y=",1)[1].split()[0])
                z_coordinate_ground_station = float(coordinates_line_ground_station.split("Z=",1)[1].split()[0])

                radio_telescopes_dict[name_ground_station] = np.array([x_coordinate_ground_station,y_coordinate_ground_station,z_coordinate_ground_station])
    
    ########################################################################################################################
    ################################################## CREATE ENVIRONMENT ##################################################
    ########################################################################################################################

    # Load spice kernels
    spice_interface.load_standard_kernels()

    # Initial date of the simulation
    start_date = 2458449.5 #in Julian days = 27/11/2018 00:00:00; Taken from the txt file sent by Sebastien

    # Duration of the simulation
    RISE_simulation_duration_days = 998 #days 
    RISE_simulation_duration = RISE_simulation_duration_days*constants.JULIAN_DAY #seconds
    simulation_duration_days = 2956 #days; until 31/12/2026
    simulation_duration = simulation_duration_days*constants.JULIAN_DAY #seconds

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

    # Determination of the correlation coefficient
    indicator_Doppler_noise = indicator_Doppler_noise_function(np.deg2rad(1))
    correlation_coefficient = dict()
    for station_1 in Earth_ground_station_dict:
        correlation_coefficient[station_1] = dict()
        for station_2 in Earth_ground_station_dict:
            if station_1 != station_2:
                correlation_coefficient[station_1][station_2] = (1-weather_percentage_function(np.deg2rad(1)))*indicator_Doppler_noise[station_1][station_2]+weather_percentage_function(np.deg2rad(1))*np.exp(-distance_stations[station_1][station_2]/(2*distance_ionospheric))

    #print(correlation_coefficient)

    print(weather_percentage_function(np.pi))

    # Computation of correlation coefficient
    #station_1 = receivers_total[total_split+row_index]
    #station_2 = receivers_total[total_split+column_index]
    #correlation_coefficient = weather_percentage*np.exp(-distance_stations[station_1][station_2]/(2*distance_ionospheric))+(1-weather_percentage)*indicator_Doppler_noise[station_1][station_2]
