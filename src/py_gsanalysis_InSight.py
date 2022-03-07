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
    import datetime
    import numpy as np
    from tudatpy.kernel import constants, numerical_simulation
    from tudatpy.kernel.astro import element_conversion, frame_conversion
    from tudatpy.kernel.interface import spice_interface
    from tudatpy.kernel.numerical_simulation import environment_setup,environment,propagation_setup,propagation,estimation_setup,estimation
    from tudatpy.kernel.numerical_simulation.estimation_setup import parameter,observation
    import matplotlib.pyplot as plt

    ########################################################################################################################
    ################################################## CONSTANTS AND VARIABLES #############################################
    ########################################################################################################################

    # days in a week
    days_in_a_week = 7 #days

    # Initial date of the simulation
    start_date = 2458423.5 #in Julian days = 01/11/2018 00:00:00; Taken from "LaRa after RISE: Expected improvement in the Mars rotation and interior models"

    # Duration of the simulation
    simulation_duration_days = 1200 #days
    simulation_duration_weeks = simulation_duration_days/days_in_a_week #weeks
    simulation_duration = simulation_duration_days*constants.JULIAN_DAY #seconds

    # LaRa landing site, taken from "LaRa after RISE: Expected improvement in the Mars rotation and interior models"
    reflector_name = 'RISE'
    reflector_latitude_deg = 4.5 #North degrees
    reflector_longitude_deg = 135.62 #East degrees
    reflector_altitude = -2600 #m

    # Earth-based transmitter
    #transmitter_name = 'DSS 63'
    transmitter_names = ['DSS 43','DSS 65','DSS 63','DSS 55','DSS 14','DSS 26','DSS 34','DSS 35','DSS 36','DSS 54','DSS 24', 'DSS 25', 'DSS 56']

    #transmitter_position_cartesian = np.array([4849092.6814,-360180.5350,4115109.1298]) #Taken from https://www.aoc.nrao.edu/software/sched/catalogs/locations.dat
    transmitter_positions_cartesian = list()  #Taken from https://www.aoc.nrao.edu/software/sched/catalogs/locations.dat
    transmitter_positions_cartesian.append(np.array([-4460894.7273,2682361.5296,-3674748.4238])) # DSS 43
    transmitter_positions_cartesian.append(np.array([4849339.5378,-360427.4851,4114750.8520])) # DSS 65A
    transmitter_positions_cartesian.append(np.array([4849092.6814,-360180.5350,4115109.1298])) # DSS 63
    transmitter_positions_cartesian.append(np.array([4849525.256,-360606.09,4114495.08])) # DSS 55 #http://astrogeo.org/aplo/vlbi.inp
    transmitter_positions_cartesian.append(np.array([-2353621.2459,-4641341.5369,3677052.2305])) # DSS 14
    transmitter_positions_cartesian.append(np.array([-2354890.967,-4647166.93,3668872.21])) # DSS 26
    transmitter_positions_cartesian.append(np.array([-4461147.4205,2682439.2423,-3674392.5623])) # DSS 34
    transmitter_positions_cartesian.append(np.array([-4461273.4175,2682568.9283,-3674151.5223])) # DSS 35
    transmitter_positions_cartesian.append(np.array([-4461168.7425,2682814.6603,-3674083.3303])) # DSS 36
    transmitter_positions_cartesian.append(np.array([4849434.4880,-360723.8999,4114618.8350])) # DSS 54
    transmitter_positions_cartesian.append(np.array([-2354906.495,-4646840.13,3669242.317])) # DSS 24
    transmitter_positions_cartesian.append(np.array([-2355022.066,-4646953.64,3669040.90])) # DSS 25
    transmitter_positions_cartesian.append(np.array([4849421.500903,-360549.2280048,4114647.264832])) # DSS 56 #https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/stations/earth_topo_201023.tf

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
    #ground_station_dict [transmitter_name] = transmitter_position_cartesian
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
        np.array([reflector_altitude,np.deg2rad(reflector_latitude_deg),np.deg2rad(reflector_longitude_deg)]),
         element_conversion.geodetic_position_type)

    Mars_ground_station_list = environment_setup.get_ground_station_list(bodies.get_body("Mars"))

    ########################################################################################################################
    ################################################## GROUND STATIONS ELEVATION HISTORY ###################################
    ########################################################################################################################

    # Define observation simulation times for each link
    #observation_times_list = list()
    data_transmitter = dict()
    
    # Read the epoch times
    with open(os.path.dirname(os.path.realpath(__file__))+'/InSight_mes_upto_31122021.forCarlos') as file:
        lines = file.read().splitlines()
        observation_start_epoch = np.inf
        for transmitter_pointer in transmitter_names:
            data_transmitter[transmitter_pointer] = dict()
            data_transmitter[transmitter_pointer]['Time at receiver'] = list()
            transmitter_ground_station_number =  [int(s) for s in transmitter_pointer.split() if s.isdigit()][0]
            for pointer_time in range(0,len(lines)):
                line = lines[pointer_time]
                line_info = line.split()
                if float(line_info[0]) == transmitter_ground_station_number and float(line_info[1]) == transmitter_ground_station_number:
                    data_transmitter[transmitter_pointer]['Time at receiver'].append(float(line_info[2]))
                    if float(line_info[2])<observation_start_epoch:
                        observation_start_epoch = float(line_info[2])
    
    for transmitter_pointer in transmitter_names:
        data_transmitter[transmitter_pointer]['Time at reflector'] = list()
        data_transmitter[transmitter_pointer]['Earth elevation'] = list()
        data_transmitter[transmitter_pointer]['Earth azimuth'] = list()
        data_transmitter[transmitter_pointer]['Time at transmitter'] = list()
        data_transmitter[transmitter_pointer]['Elevation at transmitter'] = list()
        data_transmitter[transmitter_pointer]['Elevation at receiver'] = list()
        data_transmitter[transmitter_pointer]['Azimuth at transmitter'] = list()
        data_transmitter[transmitter_pointer]['Azimuth at receiver'] = list()
        for receiver_time_pointer in range(0,len(data_transmitter[transmitter_pointer]['Time at receiver'])):
            bool_receiver = estimation.compute_target_angles_and_range(bodies,('Earth',transmitter_pointer),'Mars',
                [data_transmitter[transmitter_pointer]['Time at receiver'][receiver_time_pointer]],False)

            data_transmitter[transmitter_pointer]['Time at reflector'].append(data_transmitter[transmitter_pointer]['Time at receiver'][receiver_time_pointer]\
                -bool_receiver[list(bool_receiver.keys())[0]][2]/constants.SPEED_OF_LIGHT_LONG)
            bool_reflector = estimation.compute_target_angles_and_range(bodies,('Mars',reflector_name),'Earth',[data_transmitter[transmitter_pointer]['Time at reflector'][receiver_time_pointer]],False)

            transmitter_time = data_transmitter[transmitter_pointer]['Time at reflector'][receiver_time_pointer]-bool_reflector[list(bool_reflector.keys())[0]][2]/constants.SPEED_OF_LIGHT_LONG
            bool_transmitter = estimation.compute_target_angles_and_range(bodies,('Earth',transmitter_pointer),'Mars',
                [transmitter_time],True)

            if True:#np.deg2rad(earth_min) <= bool_reflector[list(bool_reflector.keys())[0]][0] <= np.deg2rad(earth_max) and\
                 #bool_transmitter[list(bool_transmitter.keys())[0]][0] >= np.deg2rad(antenna_min_elevation):
                data_transmitter[transmitter_pointer]['Elevation at receiver'].append(bool_receiver[list(bool_receiver.keys())[0]][0])
                data_transmitter[transmitter_pointer]['Azimuth at receiver'].append(bool_receiver[list(bool_receiver.keys())[0]][1])

                data_transmitter[transmitter_pointer]['Earth elevation'].append(bool_reflector[list(bool_reflector.keys())[0]][0])
                data_transmitter[transmitter_pointer]['Earth azimuth'].append(bool_reflector[list(bool_reflector.keys())[0]][1])

                data_transmitter[transmitter_pointer]['Time at transmitter'].append(transmitter_time)
                data_transmitter[transmitter_pointer]['Elevation at transmitter'].append(bool_transmitter[list(bool_transmitter.keys())[0]][0])
                data_transmitter[transmitter_pointer]['Azimuth at transmitter'].append(bool_transmitter[list(bool_transmitter.keys())[0]][1])


    data_receiver = dict()

    # Iterate along receivers
    for Earth_ground_station_pointer in Earth_ground_station_list:
        data_receiver[Earth_ground_station_pointer[1]] = dict()
        data_receiver[Earth_ground_station_pointer[1]]['Total'] = dict()
        data_receiver[Earth_ground_station_pointer[1]]['Total']['Observation time at receiver'] = list()
        data_receiver[Earth_ground_station_pointer[1]]['Total']['Elevation at receiver'] = list()
        data_receiver[Earth_ground_station_pointer[1]]['Total']['Azimuth at receiver'] = list()

        # Shadow tracking
        #if Earth_ground_station_pointer[1] in transmitter_names:
        #    continue

        # Iterate along transmitters
        for transmitter_pointer in transmitter_names:
            data_receiver[Earth_ground_station_pointer[1]][transmitter_pointer] = dict()
            data_receiver[Earth_ground_station_pointer[1]][transmitter_pointer]['Observation time at receiver'] = list()
            data_receiver[Earth_ground_station_pointer[1]][transmitter_pointer]['Elevation at receiver'] = list()
            data_receiver[Earth_ground_station_pointer[1]][transmitter_pointer]['Azimuth at receiver'] = list()

            # Iterate along transmitting times
            for transmitter_time_pointer in range(0,len(data_transmitter[transmitter_pointer]['Time at transmitter'])):
                bool_transmitter = estimation.compute_target_angles_and_range(bodies,('Earth',transmitter_pointer),'Mars',
                [data_transmitter[transmitter_pointer]['Time at transmitter'][transmitter_time_pointer]],True)

                reflector_time = data_transmitter[transmitter_pointer]['Time at transmitter'][transmitter_time_pointer]+bool_transmitter[list(bool_transmitter.keys())[0]][2]/constants.SPEED_OF_LIGHT_LONG
                bool_reflector = estimation.compute_target_angles_and_range(bodies,('Mars',reflector_name),'Earth',[reflector_time],False)

                receiver_time = reflector_time+bool_transmitter[list(bool_transmitter.keys())[0]][2]/constants.SPEED_OF_LIGHT_LONG
                bool_receiver = estimation.compute_target_angles_and_range(bodies,('Earth',Earth_ground_station_pointer[1]),'Mars',[receiver_time],False)

                if bool_receiver[list(bool_receiver.keys())[0]][0] >= np.deg2rad(antenna_min_elevation):
                    data_receiver[Earth_ground_station_pointer[1]][transmitter_pointer]['Observation time at receiver'].append(receiver_time)
                    data_receiver[Earth_ground_station_pointer[1]][transmitter_pointer]['Elevation at receiver'].append(bool_receiver[list(bool_receiver.keys())[0]][0])
                    data_receiver[Earth_ground_station_pointer[1]][transmitter_pointer]['Azimuth at receiver'].append(bool_receiver[list(bool_receiver.keys())[0]][1])

            data_receiver[Earth_ground_station_pointer[1]]['Total']['Observation time at receiver'].extend(data_receiver[Earth_ground_station_pointer[1]][transmitter_pointer]['Observation time at receiver'])
            index_sort = np.argsort(data_receiver[Earth_ground_station_pointer[1]]['Total']['Observation time at receiver'])
            data_receiver[Earth_ground_station_pointer[1]]['Total']['Elevation at receiver'].extend(data_receiver[Earth_ground_station_pointer[1]][transmitter_pointer]['Elevation at receiver'])
            data_receiver[Earth_ground_station_pointer[1]]['Total']['Azimuth at receiver'].extend(data_receiver[Earth_ground_station_pointer[1]][transmitter_pointer]['Azimuth at receiver'])
        
            data_receiver[Earth_ground_station_pointer[1]]['Total']['Observation time at receiver'] = [data_receiver[Earth_ground_station_pointer[1]]['Total']['Observation time at receiver'][i] for i in index_sort]
            data_receiver[Earth_ground_station_pointer[1]]['Total']['Elevation at receiver'] = [data_receiver[Earth_ground_station_pointer[1]]['Total']['Elevation at receiver'][i] for i in index_sort]
            data_receiver[Earth_ground_station_pointer[1]]['Total']['Azimuth at receiver'] = [data_receiver[Earth_ground_station_pointer[1]]['Total']['Azimuth at receiver'][i] for i in index_sort]

    ########################################################################################################################
    ################################################## PLOTS ###############################################################
    ########################################################################################################################
    
    output_folder_path = os.path.dirname(os.path.realpath(__file__)).replace('/src','/output/GS_RISE')
    os.makedirs(output_folder_path,exist_ok=True)

    plt.figure(figsize=(15, 6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(transmitter_names)))))
    for transmitter_pointer in transmitter_names:
        plt.scatter((np.array(data_transmitter[transmitter_pointer]['Time at reflector'])-observation_start_epoch)/constants.JULIAN_DAY,
            np.rad2deg(data_transmitter[transmitter_pointer]['Earth elevation']),s=5,label=transmitter_pointer)
    plt.ylabel('Earth elevation as seen by '+reflector_name+' [deg]')
    plt.xlabel('Time after landing [Earth days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.grid()
    plt.savefig(output_folder_path+'/Earth_elevation_seen_'+reflector_name+'_vs_time.pdf',bbox_inches='tight')
    plt.show()
    plt.close('all')

    plt.figure(figsize=(15,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(transmitter_names)))))
    for transmitter_pointer in transmitter_names:
        plt.scatter((np.array(data_transmitter[transmitter_pointer]['Time at reflector'])-observation_start_epoch)/constants.JULIAN_DAY,
            np.mod(-np.rad2deg(data_transmitter[transmitter_pointer]['Earth azimuth'])+90,360),s=5,label=transmitter_pointer)
    plt.ylabel('Earth azimuth as seen by '+reflector_name+' [deg]')
    plt.xlabel('Time after landing [Earth days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.grid()
    plt.savefig(output_folder_path+'/Earth_azimuth_seen_'+reflector_name+'_vs_time.pdf',bbox_inches='tight')
    plt.show()
    plt.close('all')

    plt.figure(figsize=(15,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(transmitter_names)))))
    for transmitter_pointer in transmitter_names:
        plt.scatter(np.mod(-np.rad2deg(data_transmitter[transmitter_pointer]['Earth azimuth'])+90,360),
            np.rad2deg(data_transmitter[transmitter_pointer]['Earth elevation']),s=5,label=transmitter_pointer)
    plt.ylabel('Earth elevation as seen by '+reflector_name+' [deg]')
    plt.xlabel('Earth azimuth as seen by '+reflector_name+' [deg]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.grid()
    plt.savefig(output_folder_path+'/Earth_elevation_seen_'+reflector_name+'_vs_Earth_azimuth_seen_'+reflector_name+'.pdf',bbox_inches='tight')
    plt.show()
    plt.close('all')

    plt.figure(figsize=(15,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(transmitter_names)))))
    for transmitter_pointer in transmitter_names:
        plt.scatter((np.array(data_transmitter[transmitter_pointer]['Time at transmitter'])-observation_start_epoch)/constants.JULIAN_DAY,
            np.rad2deg(data_transmitter[transmitter_pointer]['Elevation at transmitter']),s=5,label=transmitter_pointer)
    plt.ylabel('Elevation at DSN transmitter [deg]')
    plt.xlabel('Time after landing [Earth days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.grid()
    plt.savefig(output_folder_path+'/Transmitter_elevation_vs_time_transmission.pdf',bbox_inches='tight')
    plt.show()
    plt.close('all')

    plt.figure(figsize=(15,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(transmitter_names)))))
    for transmitter_pointer in transmitter_names:
        plt.scatter((np.array(data_transmitter[transmitter_pointer]['Time at receiver'])-observation_start_epoch)/constants.JULIAN_DAY,
            np.rad2deg(data_transmitter[transmitter_pointer]['Elevation at receiver']),s=5,label=transmitter_pointer)
    plt.ylabel('Elevation at DSN receiver [deg]')
    plt.xlabel('Time after landing [Earth days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.grid()
    plt.savefig(output_folder_path+'/Transmitter_elevation_vs_time_reception.pdf',bbox_inches='tight')
    plt.show()
    plt.close('all')

    
    plt.figure(figsize=(15,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(transmitter_names)))))
    for transmitter_pointer in transmitter_names:
        plt.plot(np.arange(0,(data_transmitter[transmitter_pointer]['Time at transmitter'][-1]-observation_start_epoch)/constants.JULIAN_DAY,1),
            np.polyval(np.polyfit((np.array(data_transmitter[transmitter_pointer]['Time at transmitter'])-observation_start_epoch)/constants.JULIAN_DAY,
                np.rad2deg(data_transmitter[transmitter_pointer]['Elevation at transmitter']),4),
            np.arange(0,(data_transmitter[transmitter_pointer]['Time at transmitter'][-1]-observation_start_epoch)/constants.JULIAN_DAY,1)),label='Mean - '+transmitter_pointer)
    plt.ylabel('Elevation at DSN transmitter [deg]')
    plt.xlabel('Time after landing [Earth days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.ylim([0,90])
    plt.grid()
    plt.savefig(output_folder_path+'/Transmitter_elevation_vs_time_mean_transmission.pdf',bbox_inches='tight')
    plt.show()
    plt.close('all')

    plt.figure(figsize=(15,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(transmitter_names)))))
    for transmitter_pointer in transmitter_names:
        plt.plot(np.arange(0,(data_transmitter[transmitter_pointer]['Time at receiver'][-1]-observation_start_epoch)/constants.JULIAN_DAY,1),
            np.polyval(np.polyfit((np.array(data_transmitter[transmitter_pointer]['Time at receiver'])-observation_start_epoch)/constants.JULIAN_DAY,
                np.rad2deg(data_transmitter[transmitter_pointer]['Elevation at receiver']),4),
            np.arange(0,(data_transmitter[transmitter_pointer]['Time at receiver'][-1]-observation_start_epoch)/constants.JULIAN_DAY,1)),label='Mean - '+transmitter_pointer)
    plt.ylabel('Elevation at DSN receiver [deg]')
    plt.xlabel('Time after landing [Earth days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.ylim([0,90])
    plt.grid()
    plt.savefig(output_folder_path+'/Transmitter_elevation_vs_time_mean_receiver.pdf',bbox_inches='tight')
    plt.show()
    plt.close('all')

    plt.figure(figsize=(15,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(Earth_ground_station_list)))))
    for Earth_ground_station_pointer in Earth_ground_station_list:
        plt.scatter((np.array(data_receiver[Earth_ground_station_pointer[1]]['Total']['Observation time at receiver'])-observation_start_epoch)/constants.JULIAN_DAY,
            np.rad2deg(data_receiver[Earth_ground_station_pointer[1]]['Total']['Elevation at receiver']),s=5,label=Earth_ground_station_pointer[1])
    plt.ylabel('Elevation at receivers [deg]')
    plt.xlabel('Time after landing [Earth days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.grid()
    plt.savefig(output_folder_path+'/Receivers_elevation_vs_time.pdf',bbox_inches='tight')
    plt.show()
    plt.close('all')  

    plt.figure(figsize=(15,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(Earth_ground_station_list)))))
    for Earth_ground_station_pointer in Earth_ground_station_list:
        plt.plot(np.arange(0,(data_receiver[Earth_ground_station_pointer[1]]['Total']['Observation time at receiver'][-1]-observation_start_epoch)/constants.JULIAN_DAY,1),
            np.polyval(np.polyfit((np.array(data_receiver[Earth_ground_station_pointer[1]]['Total']['Observation time at receiver'])-observation_start_epoch)/constants.JULIAN_DAY,
                np.rad2deg(data_receiver[Earth_ground_station_pointer[1]]['Total']['Elevation at receiver']),4),
            np.arange(0,(data_receiver[Earth_ground_station_pointer[1]]['Total']['Observation time at receiver'][-1]-observation_start_epoch)/constants.JULIAN_DAY,1)),
            label='Mean - '+Earth_ground_station_pointer[1])
    plt.ylabel('Elevation at receivers [deg]')
    plt.xlabel('Time after landing [Earth days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.grid()
    plt.savefig(output_folder_path+'/Receivers_elevation_vs_time_mean.pdf',bbox_inches='tight')
    plt.show()
    plt.close('all')
    
    plt.figure(figsize=(15,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(Earth_ground_station_list)))))
    station_names = list()
    receiver_observation_number = list()
    for Earth_ground_station_pointer in Earth_ground_station_list:
        station_names.append(Earth_ground_station_pointer[1])
        receiver_observation_number.append(len(data_receiver[Earth_ground_station_pointer[1]]['Total']['Observation time at receiver']))
        plt.scatter((np.array(data_receiver[Earth_ground_station_pointer[1]]['Total']['Observation time at receiver'])-observation_start_epoch)/constants.JULIAN_DAY,
            station_names.index(Earth_ground_station_pointer[1])*np.ones(len(data_receiver[Earth_ground_station_pointer[1]]['Total']['Observation time at receiver'])),s=5)
    plt.grid()
    plt.yticks(range(0,len(station_names)),labels=station_names)
    plt.xlabel('Time after landing [Earth days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.savefig(output_folder_path+'/Receivers_scatter.pdf',bbox_inches='tight')
    plt.show()
    plt.close('all')

    plt.figure(figsize=(15,6))
    plot1 = plt.barh(np.arange(len(Earth_ground_station_list)),receiver_observation_number)
    plt.yticks(np.arange(len(station_names)),labels=station_names)
    plt.xlabel('Number of Observations')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.bar_label(plot1)
    plt.grid()
    plt.savefig(output_folder_path+'/Receivers_hbar.pdf',bbox_inches='tight')
    plt.show()
    plt.close('all')

print("--- %s seconds ---" % (time.time() - run_time))

'''
    # Initialize useful lists            
    earth_elevation = list()
    earth_azimuth = list()
    transmitter_observation_time = list()
    transmitter_elevation = list()
    transmitter_azimuth = list()
    receiver_ids = list()
    receiver_observation_time = list()
    receiver_elevation = list()
    receiver_azimuth = list()

    # Iterate along each observation time
    for pointer_time in range(0,len(observation_times_list)):
        # Compute elevation and azimuths
        bool_reflector = estimation.compute_target_angles_and_range(bodies,('Mars',reflector_name),'Earth',[observation_times_list[pointer_time]],False)
        bool_transmitter = estimation.compute_target_angles_and_range(bodies,('Earth',transmitter_name),'Mars',[observation_times_list[pointer_time]],True)

        # Append the Earth elevation and azimuth as seen by RISE
        earth_elevation.append(bool_reflector[list(bool_reflector.keys())[0]][0])
        earth_azimuth.append(bool_reflector[list(bool_reflector.keys())[0]][1])

        # Check angle viability for transmitting the signal
        if (earth_elevation[-1] >= np.deg2rad(earth_min) and earth_elevation[-1] <= np.deg2rad(earth_max)) and np.rad2deg(bool_transmitter[list(bool_transmitter.keys())[0]][0]) >= antenna_min_elevation:
            # Appent the transmission time, DSN elevation and azimuth
            transmitter_observation_time.append(observation_times_list[pointer_time])
            transmitter_elevation.append(bool_transmitter[list(bool_transmitter.keys())[0]][0]) 
            transmitter_azimuth.append(bool_transmitter[list(bool_transmitter.keys())[0]][1]) 

    # Iterate along the receivers (shadow tracking)
    for Earth_ground_station_pointer in Earth_ground_station_list:
        if Earth_ground_station_pointer[1] == transmitter_name:
            continue
        
        # Iterate along the transmitting time
        for pointer_time in range(0,len(transmitter_observation_time)):
            # Compute elevation and azimuth
            bool_reflector = estimation.compute_target_angles_and_range(bodies,('Mars',reflector_name),'Earth',[transmitter_observation_time[pointer_time]],True)
            bool_receiver = estimation.compute_target_angles_and_range(bodies,('Earth',Earth_ground_station_pointer[1]),'Mars',[transmitter_observation_time[pointer_time]],False)

            # Check angle viability to receive the signal
            if np.rad2deg(bool_receiver[list(bool_receiver.keys())[0]][0]) >= antenna_min_elevation:
                # Append the receiver ids, time, elevation and azimuth
                receiver_ids.append(Earth_ground_station_list.index(Earth_ground_station_pointer))
                receiver_observation_time.append(transmitter_observation_time[pointer_time])
                receiver_elevation.append(bool_receiver[list(bool_receiver.keys())[0]][0])
                receiver_azimuth.append(bool_receiver[list(bool_receiver.keys())[0]][1])

    
    ########################################################################################################################
    ################################################## PLOTS ###############################################################
    ########################################################################################################################

    output_folder_path = os.path.dirname(os.path.realpath(__file__)).replace('/src','/output/GS_RISE')
    os.makedirs(output_folder_path,exist_ok=True)

    plt.figure()
    plt.scatter((np.array(observation_times_list)-min(observation_times_list)*np.ones(len(observation_times_list)))/constants.JULIAN_DAY,
        np.rad2deg(earth_elevation),s=5,c='r',label='No Tracking Windows of '+reflector_name+' Found')
    plt.scatter((np.array(transmitter_observation_time)-min(observation_times_list)*np.ones(len(transmitter_observation_time)))/constants.JULIAN_DAY,
        list(map(np.rad2deg(earth_elevation).__getitem__,np.nonzero(np.in1d(observation_times_list,transmitter_observation_time))[0])),s=5,c='b',
        label='Tracking Windows of '+reflector_name)
    plt.ylabel('Earth elevation as seen by '+reflector_name+' [deg]')
    plt.xlabel('Time after landing [Earth days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_times_list[0])))
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.grid()
    plt.savefig(output_folder_path+'/Earth_elevation_seen_'+reflector_name+'_vs_time.pdf',bbox_inches='tight')
    plt.show()
    plt.close('all')

    plt.figure()
    plt.scatter((np.array(observation_times_list)-min(observation_times_list)*np.ones(len(observation_times_list)))/constants.JULIAN_DAY,
        np.mod(-np.rad2deg(earth_azimuth)+90*np.ones(len(observation_times_list)),360),s=5,c='r',label='No Tracking Windows of '+reflector_name+' Found')
    plt.scatter((np.array(transmitter_observation_time)-min(observation_times_list)*np.ones(len(transmitter_observation_time)))/constants.JULIAN_DAY,
        list(map(np.mod(-np.rad2deg(earth_azimuth)+90*np.ones(len(observation_times_list)),360).__getitem__,
        np.nonzero(np.in1d(observation_times_list,transmitter_observation_time))[0])),s=5,c='b',label='Tracking Windows of '+reflector_name)
    plt.ylabel('Earth azimuth as seen by '+reflector_name+' [deg]')
    plt.xlabel('Time after landing [Earth days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_times_list[0])))
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.grid()
    plt.savefig(output_folder_path+'/Earth_azimuth_seen_'+reflector_name+'_vs_time.pdf',bbox_inches='tight')
    plt.show()
    plt.close('all')

    plt.figure()
    plt.scatter(np.mod(-np.rad2deg(earth_azimuth)+90*np.ones(len(observation_times_list)),360),np.rad2deg(earth_elevation),s=5,c='r',
        label='No Tracking Windows of '+reflector_name+' Found')
    plt.scatter(list(map(np.mod(-np.rad2deg(earth_azimuth)+90*np.ones(len(observation_times_list)),360).__getitem__,
        np.nonzero(np.in1d(observation_times_list,transmitter_observation_time))[0])),
        list(map(np.rad2deg(earth_elevation).__getitem__,np.nonzero(np.in1d(observation_times_list,transmitter_observation_time))[0])),s=5,c='b',
        label='Tracking Windows of '+reflector_name)
    plt.ylabel('Earth elevation as seen by '+reflector_name+' [deg]')
    plt.xlabel('Earth azimuth as seen by '+reflector_name+' [deg]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_times_list[0])))
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.grid()
    plt.savefig(output_folder_path+'/Earth_elevation_seen_'+reflector_name+'_vs_Earth_azimuth_seen_'+reflector_name+'.pdf',bbox_inches='tight')
    plt.show()
    plt.close('all')

    plt.figure()
    plt.scatter((np.array(transmitter_observation_time)-min(observation_times_list)*np.ones(len(transmitter_observation_time)))/constants.JULIAN_DAY,np.rad2deg(transmitter_elevation),s=5,c='b',
        label='Tracking Windows of '+reflector_name)
    plt.plot(np.arange(0,(observation_times_list[-1]-observation_times_list[0])/constants.JULIAN_DAY,1),
        np.polyval(np.polyfit((np.array(transmitter_observation_time)-min(observation_times_list)*np.ones(len(transmitter_observation_time)))/constants.JULIAN_DAY,
        np.rad2deg(transmitter_elevation),4),np.arange(0,(observation_times_list[-1]-observation_times_list[0])/constants.JULIAN_DAY,1)),'g-',label='Mean Elevation for '+reflector_name+' Tracking')
    plt.ylabel('Elevation at '+transmitter_name+' [deg]')
    plt.xlabel('Time after landing [Earth days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_times_list[0])))
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.grid()
    plt.savefig(output_folder_path+'/Transmitter_elevation_vs_time.pdf',bbox_inches='tight')
    plt.show()
    plt.close('all')

    plt.figure()
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(Earth_ground_station_list)-1))))
    
    for Earth_ground_station_pointer in Earth_ground_station_list:
        if Earth_ground_station_pointer[1] == transmitter_name or \
            len(list(map(receiver_observation_time.__getitem__,np.nonzero(np.in1d(receiver_ids,Earth_ground_station_list.index(Earth_ground_station_pointer)))[0])))==0:
            continue
        plt.plot(np.arange(0,(observation_times_list[-1]-observation_times_list[0])/constants.JULIAN_DAY,7),
        np.polyval(np.polyfit((np.array(list(map(receiver_observation_time.__getitem__,np.nonzero(np.in1d(receiver_ids,Earth_ground_station_list.index(Earth_ground_station_pointer)))[0])))\
            -min(observation_times_list)*np.ones(len(list(map(receiver_observation_time.__getitem__,np.nonzero(np.in1d(receiver_ids,Earth_ground_station_list.index(Earth_ground_station_pointer)))[0])))))\
                /constants.JULIAN_DAY,
        list(map(np.rad2deg(receiver_elevation).__getitem__,np.nonzero(np.in1d(receiver_ids,Earth_ground_station_list.index(Earth_ground_station_pointer)))[0])),4),
            np.arange(0,(observation_times_list[-1]-observation_times_list[0])/constants.JULIAN_DAY,7)),label='Mean Elevation for '+Earth_ground_station_pointer[1])
        
    plt.ylabel('Elevation at receivers [deg]')
    plt.xlabel('Time after landing [Earth days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_times_list[0])))
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.ylim([0,90])
    plt.grid()
    plt.savefig(output_folder_path+'/Receivers_elevation_vs_time.pdf',bbox_inches='tight')
    plt.show()
    plt.close('all')

    plt.figure()
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(Earth_ground_station_list)-1))))
    
    for Earth_ground_station_pointer in Earth_ground_station_list:
        if Earth_ground_station_pointer[1] == transmitter_name or \
            len(list(map(receiver_observation_time.__getitem__,np.nonzero(np.in1d(receiver_ids,Earth_ground_station_list.index(Earth_ground_station_pointer)))[0])))==0:
            continue
        plt.scatter((np.array(list(map(receiver_observation_time.__getitem__,np.nonzero(np.in1d(receiver_ids,Earth_ground_station_list.index(Earth_ground_station_pointer)))[0])))\
            -min(observation_times_list)*np.ones(len(list(map(receiver_observation_time.__getitem__,np.nonzero(np.in1d(receiver_ids,Earth_ground_station_list.index(Earth_ground_station_pointer)))[0])))))/constants.JULIAN_DAY,
        list(map(np.rad2deg(receiver_elevation).__getitem__,np.nonzero(np.in1d(receiver_ids,Earth_ground_station_list.index(Earth_ground_station_pointer)))[0])),
        s=5,label=Earth_ground_station_pointer[1])
    plt.ylabel('Elevation at receivers [deg]')
    plt.xlabel('Time after landing [Earth days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_times_list[0])))
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.grid()
    plt.savefig(output_folder_path+'/Receivers_elevation_vs_time_scatter.pdf',bbox_inches='tight')
    plt.show()
    plt.close('all')

    plt.figure()
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(Earth_ground_station_list)-1))))
    station_names = list()
    receiver_observation_number = list()
    for Earth_ground_station_pointer in Earth_ground_station_list:
        station_names.append(Earth_ground_station_pointer[1])
        receiver_observation_number.append(len(list(map(receiver_observation_time.__getitem__,np.nonzero(np.in1d(receiver_ids,Earth_ground_station_list.index(Earth_ground_station_pointer)))[0]))))
        if Earth_ground_station_pointer[1] == transmitter_name or \
            len(list(map(receiver_observation_time.__getitem__,np.nonzero(np.in1d(receiver_ids,Earth_ground_station_list.index(Earth_ground_station_pointer)))[0])))==0:
            continue
        plt.scatter((np.array(list(map(receiver_observation_time.__getitem__,np.nonzero(np.in1d(receiver_ids,Earth_ground_station_list.index(Earth_ground_station_pointer)))[0])))\
            -min(observation_times_list)*np.ones(len(list(map(receiver_observation_time.__getitem__,np.nonzero(np.in1d(receiver_ids,Earth_ground_station_list.index(Earth_ground_station_pointer)))[0])))))/constants.JULIAN_DAY,
            Earth_ground_station_list.index(Earth_ground_station_pointer)*np.ones(len(list(map(receiver_observation_time.__getitem__,np.nonzero(np.in1d(receiver_ids,Earth_ground_station_list.index(Earth_ground_station_pointer)))[0])))),
            s=5)
    plt.grid()
    plt.yticks(range(0,len(station_names)),labels=station_names)
    plt.xlabel('Time after landing [Earth days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_times_list[0])))
    plt.savefig(output_folder_path+'/Receivers_scatter.pdf',bbox_inches='tight')
    plt.show()
    plt.close('all')

    plt.figure()
    plot1 = plt.barh(np.arange(len(station_names)),receiver_observation_number)
    plt.yticks(np.arange(len(station_names)),labels=station_names)
    plt.xlabel('Number of Observations')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_times_list[0])))
    plt.bar_label(plot1)
    plt.grid()
    plt.savefig(output_folder_path+'/Receivers_hbar.pdf',bbox_inches='tight')
    plt.show()
    plt.close('all')
    
    ########################################################################################################################
    ################################################## PROVIDE OUTPUT TO CONSOLE AND FILES #################################
    ########################################################################################################################

    np.savetxt(output_folder_path+"/observation_time.dat",observation_times_list,fmt='%.15e')
    np.savetxt(output_folder_path+"/DSS63_observation_time.dat",transmitter_observation_time,fmt='%.15e')
    np.savetxt(output_folder_path+"/DSS63_elevation.dat",transmitter_elevation,fmt='%.15e')
    np.savetxt(output_folder_path+"/earth_elevation.dat",earth_elevation,fmt='%.15e')
    np.savetxt(output_folder_path+"/earth_azimuth.dat",earth_azimuth,fmt='%.15e')   
    np.savetxt(output_folder_path+"/ground_station_observation_time.dat",receiver_observation_time,fmt='%.15e')
    np.savetxt(output_folder_path+"/ground_station_elevation.dat",receiver_elevation,fmt='%.15e')
    np.savetxt(output_folder_path+"/ground_station_ids.dat",receiver_ids,fmt='%.15e')  
'''