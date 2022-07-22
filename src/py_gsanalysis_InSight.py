"""
Description: Ground Station Analysis for RISE

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
    start_date = 2458449.5 #in Julian days = 27/11/2018 00:00:00; Taken from the txt file sent by Sebastien

    # Duration of the simulation
    simulation_duration_days = 1000 #days
    simulation_duration_weeks = simulation_duration_days/days_in_a_week #weeks
    simulation_duration = simulation_duration_days*constants.JULIAN_DAY #seconds

    # RISE landing site, taken from "LaRa after RISE: Expected improvement in the Mars rotation and interior models"
    reflector_name = 'RISE'
    reflector_latitude_deg = 4.5 #North degrees
    reflector_longitude_deg = 135.62 #East degrees
    reflector_altitude = -2600 #m

    # Earth-based transmitter
    transmitter_names = ['DSS 43','DSS 34','DSS 35','DSS 36','DSS 65','DSS 63','DSS 55','DSS 54','DSS 56','DSS 14','DSS 26', 'DSS 24', 'DSS 25']

    transmitter_positions_cartesian = list()  #Taken from JPL web site
    transmitter_positions_cartesian.append(np.array([-4460894.9170,2682361.5070,-3674748.1517])) # DSS 43
    transmitter_positions_cartesian.append(np.array([-4461147.0925,2682439.2385,-3674393.1332])) # DSS 34
    transmitter_positions_cartesian.append(np.array([-4461273.4175,2682568.9283,-3674151.5223])) # DSS 35 (https://www.aoc.nrao.edu/software/sched/catalogs/locations.dat)
    transmitter_positions_cartesian.append(np.array([-4461168.7425,2682814.6603,-3674083.3303])) # DSS 36 (https://www.aoc.nrao.edu/software/sched/catalogs/locations.dat)
    transmitter_positions_cartesian.append(np.array([4849339.6448,-360427.6560,4114750.7428])) # DSS 65
    transmitter_positions_cartesian.append(np.array([4849092.5175,-360180.3480,4115109.2506])) # DSS 63
    transmitter_positions_cartesian.append(np.array([4849525.2561,-360606.0932,4114495.0843])) # DSS 55
    transmitter_positions_cartesian.append(np.array([4849434.4877,-360723.8999,4114618.8354])) # DSS 54
    transmitter_positions_cartesian.append(np.array([4849421.500903,-360549.2280048,4114647.264832])) # DSS 56 (https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/stations/earth_topo_201023.tf)
    transmitter_positions_cartesian.append(np.array([-2353621.4197,-4641341.4717,3677052.3178])) # DSS 14
    transmitter_positions_cartesian.append(np.array([-2354890.7996,-4647166.3182,3668871.7546])) # DSS 26
    transmitter_positions_cartesian.append(np.array([-2354906.7087,-4646840.0834,3669242.3207])) # DSS 24
    transmitter_positions_cartesian.append(np.array([-2355022.0140,-4646953.2040,3669040.5666])) # DSS 25

    zones = {'Australia':['DSS 43','DSS 34','DSS 35','DSS 36'],'Spain':['DSS 65','DSS 63','DSS 55','DSS 54','DSS 56'],'USA':['DSS 14','DSS 26','DSS 24','DSS 25']}
    
    # Viability settings
    earth_min = 10 #deg
    earth_max = 30 #deg
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
    bodies_to_create = ["Earth","Mars","Sun"]

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

    data_transmitter = dict()
    
    # Read the epoch times
    with open(os.path.dirname(os.path.realpath(__file__))+'/InSight_mes_upto_31122021.forCarlos') as file:
        lines = file.read().splitlines()
        observation_start_epoch = np.inf
        # Iterate along each transmitter
        for transmitter_pointer in transmitter_names:
            data_transmitter[transmitter_pointer] = dict()
            data_transmitter[transmitter_pointer]['Epoch'] = list()
            transmitter_ground_station_number =  [int(s) for s in transmitter_pointer.split() if s.isdigit()][0]
            # Iterate along each observation time
            for pointer_time in range(0,len(lines)):
                line = lines[pointer_time]
                line_info = line.split()
                # Condition to save the observation time
                if float(line_info[0]) == transmitter_ground_station_number and float(line_info[1]) == transmitter_ground_station_number \
                    and float(line_info[2])<=simulation_end_epoch:
                    data_transmitter[transmitter_pointer]['Epoch'].append(float(line_info[2]))
                    # Save the minimum epoch
                    if float(line_info[2])<observation_start_epoch:
                        observation_start_epoch = float(line_info[2])

    # Iterate along transmitters
    for transmitter_pointer in transmitter_names:
        data_transmitter[transmitter_pointer]['Time at reflector'] = list()
        data_transmitter[transmitter_pointer]['Earth elevation'] = list()
        data_transmitter[transmitter_pointer]['Earth azimuth'] = list()
        data_transmitter[transmitter_pointer]['Time at receiver'] = list()
        data_transmitter[transmitter_pointer]['Time at transmitter'] = list()
        data_transmitter[transmitter_pointer]['Elevation at transmitter'] = list()
        data_transmitter[transmitter_pointer]['Elevation at receiver'] = list()
        data_transmitter[transmitter_pointer]['Azimuth at transmitter'] = list()
        data_transmitter[transmitter_pointer]['Azimuth at receiver'] = list()
        # Iterate along each receiver time
        for receiver_time_pointer in range(0,len(data_transmitter[transmitter_pointer]['Epoch'])):
            # Compute azimuth, elevation angles and range for receiver
            bool_receiver = estimation.compute_target_angles_and_range(bodies,('Earth',transmitter_pointer),'Mars',
                [data_transmitter[transmitter_pointer]['Epoch'][receiver_time_pointer]],False)
            
            # Compute observation time, azimuth, elevation angles and range for reflector
            time_reflector = data_transmitter[transmitter_pointer]['Epoch'][receiver_time_pointer]-bool_receiver[list(bool_receiver.keys())[0]][2]/constants.SPEED_OF_LIGHT_LONG
            bool_reflector = estimation.compute_target_angles_and_range(bodies,('Mars',reflector_name),'Earth',[time_reflector],False)

            # Compute observation time, azimuth, elevation angles and range for transmitter
            time_transmitter = time_reflector-bool_reflector[list(bool_reflector.keys())[0]][2]/constants.SPEED_OF_LIGHT_LONG
            bool_transmitter = estimation.compute_target_angles_and_range(bodies,('Earth',transmitter_pointer),'Mars',
                [time_transmitter],True)

            # Compute azimuth, elevation angles and range for the body avoidance angle
            def polar2cart(r, theta, phi):
                return np.array([r*np.sin(theta)*np.cos(phi),r*np.sin(theta)*np.sin(phi),r*np.cos(theta)])

            bool_receiver_Sun = estimation.compute_target_angles_and_range(bodies,('Earth',transmitter_pointer),'Sun',
                [data_transmitter[transmitter_pointer]['Epoch'][receiver_time_pointer]],False)

            values_receiver_Sun = polar2cart(1,np.pi-bool_receiver_Sun[list(bool_receiver_Sun.keys())[0]][0],bool_receiver_Sun[list(bool_receiver_Sun.keys())[0]][1])
            values_receiver_Mars = polar2cart(1,np.pi-bool_receiver[list(bool_receiver.keys())[0]][0],bool_receiver[list(bool_receiver.keys())[0]][1])
            angle_body_between_receiver = np.arccos(np.dot(values_receiver_Sun,values_receiver_Mars))

            bool_transmitter_Sun = estimation.compute_target_angles_and_range(bodies,('Earth',transmitter_pointer),'Sun',
                [time_transmitter],True)

            values_transmitter_Sun = polar2cart(1,np.pi-bool_transmitter_Sun[list(bool_transmitter_Sun.keys())[0]][0],bool_transmitter_Sun[list(bool_transmitter_Sun.keys())[0]][1])
            values_transmitter_Mars = polar2cart(1,np.pi-bool_transmitter[list(bool_transmitter.keys())[0]][0],bool_transmitter[list(bool_transmitter.keys())[0]][1])
            angle_body_between_transmitter = np.arccos(np.dot(values_transmitter_Sun,values_transmitter_Mars))

            # Viability setting
            if not (angle_body_between_transmitter < np.deg2rad(body_avoidance_angle) and bool_transmitter[list(bool_transmitter.keys())[0]][2]>bool_transmitter_Sun[list(bool_transmitter_Sun.keys())[0]][2]) \
                and not (angle_body_between_receiver < np.deg2rad(body_avoidance_angle) and bool_receiver[list(bool_receiver.keys())[0]][2]>bool_receiver_Sun[list(bool_receiver_Sun.keys())[0]][2]):
                #np.deg2rad(earth_min) <= bool_reflector[list(bool_reflector.keys())[0]][0] <= np.deg2rad(earth_max) and\
                 #bool_transmitter[list(bool_transmitter.keys())[0]][0] >= np.deg2rad(antenna_min_elevation):
                
                data_transmitter[transmitter_pointer]['Time at receiver'].append(data_transmitter[transmitter_pointer]['Epoch'][receiver_time_pointer])
                data_transmitter[transmitter_pointer]['Elevation at receiver'].append(bool_receiver[list(bool_receiver.keys())[0]][0])
                data_transmitter[transmitter_pointer]['Azimuth at receiver'].append(bool_receiver[list(bool_receiver.keys())[0]][1])

                data_transmitter[transmitter_pointer]['Time at reflector'].append(time_reflector)
                data_transmitter[transmitter_pointer]['Earth elevation'].append(bool_reflector[list(bool_reflector.keys())[0]][0])
                data_transmitter[transmitter_pointer]['Earth azimuth'].append(bool_reflector[list(bool_reflector.keys())[0]][1])

                data_transmitter[transmitter_pointer]['Time at transmitter'].append(time_transmitter)
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
        if Earth_ground_station_pointer[1] in transmitter_names:
            continue

        # Iterate along transmitters
        for transmitter_pointer in transmitter_names:
            data_receiver[Earth_ground_station_pointer[1]][transmitter_pointer] = dict()
            data_receiver[Earth_ground_station_pointer[1]][transmitter_pointer]['Observation time at receiver'] = list()
            data_receiver[Earth_ground_station_pointer[1]][transmitter_pointer]['Elevation at receiver'] = list()
            data_receiver[Earth_ground_station_pointer[1]][transmitter_pointer]['Azimuth at receiver'] = list()

            # Iterate along transmitting times
            for transmitter_time_pointer in range(0,len(data_transmitter[transmitter_pointer]['Time at transmitter'])):
                # Compute azimuth, elevation angles and range for transmitter
                bool_transmitter = estimation.compute_target_angles_and_range(bodies,('Earth',transmitter_pointer),'Mars',
                [data_transmitter[transmitter_pointer]['Time at transmitter'][transmitter_time_pointer]],True)

                # Compute observation time, azimuth, elevation angles and range for reflector
                reflector_time = data_transmitter[transmitter_pointer]['Time at transmitter'][transmitter_time_pointer]+bool_transmitter[list(bool_transmitter.keys())[0]][2]/constants.SPEED_OF_LIGHT_LONG
                bool_reflector = estimation.compute_target_angles_and_range(bodies,('Mars',reflector_name),'Earth',[reflector_time],False)

                # Compute observation time, azimuth, elevation angles and range for receiver
                receiver_time = reflector_time+bool_transmitter[list(bool_transmitter.keys())[0]][2]/constants.SPEED_OF_LIGHT_LONG
                bool_receiver = estimation.compute_target_angles_and_range(bodies,('Earth',Earth_ground_station_pointer[1]),'Mars',[receiver_time],False)

                # Viability setting
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

    plt.figure(figsize=(10,7))
    plt.rcParams.update({'font.size': 16})
    colors = [plt.cm.jet(i) for i in np.linspace(0, 1, len(zones.keys()))]
    plt.gca().set_prop_cycle(plt.cycler('color', colors))
    for zones_pointer in list(zones.keys()):
        one_ground_station_per_location = 0
        for transmitter_pointer in zones[zones_pointer]:
            if one_ground_station_per_location==0:
                plt.scatter((np.array(data_transmitter[transmitter_pointer]['Time at reflector'])-observation_start_epoch)/constants.JULIAN_DAY,
                    np.rad2deg(data_transmitter[transmitter_pointer]['Earth elevation']),s=25,label=zones_pointer,color=colors[list(zones.keys()).index(zones_pointer)])
            else:
                plt.scatter((np.array(data_transmitter[transmitter_pointer]['Time at reflector'])-observation_start_epoch)/constants.JULIAN_DAY,
                    np.rad2deg(data_transmitter[transmitter_pointer]['Earth elevation']),s=25,color=colors[list(zones.keys()).index(zones_pointer)])
            one_ground_station_per_location+=1
    plt.ylabel('Earth elevation as seen by '+reflector_name+' [deg]')
    plt.xlabel('Time after landing [Earth days]')
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.legend(loc='best')
    plt.grid()
    plt.savefig(output_folder_path+'/Earth_elevation_seen_'+reflector_name+'_vs_time.pdf',bbox_inches='tight')
    plt.show()
    plt.close('all')

    plt.figure(figsize=(10,7))
    plt.rcParams.update({'font.size': 16})
    colors = [plt.cm.jet(i) for i in np.linspace(0, 1, len(zones.keys()))]
    plt.gca().set_prop_cycle(plt.cycler('color', colors))
    for zones_pointer in list(zones.keys()):
        one_ground_station_per_location = 0
        for transmitter_pointer in zones[zones_pointer]:
            if one_ground_station_per_location==0:
                plt.scatter((np.array(data_transmitter[transmitter_pointer]['Time at reflector'])-observation_start_epoch)/constants.JULIAN_DAY,
                    np.mod(-np.rad2deg(data_transmitter[transmitter_pointer]['Earth azimuth'])+90,360),s=25,label=zones_pointer,color=colors[list(zones.keys()).index(zones_pointer)])
            else:
                plt.scatter((np.array(data_transmitter[transmitter_pointer]['Time at reflector'])-observation_start_epoch)/constants.JULIAN_DAY,
                    np.mod(-np.rad2deg(data_transmitter[transmitter_pointer]['Earth azimuth'])+90,360),s=25,color=colors[list(zones.keys()).index(zones_pointer)])
            one_ground_station_per_location+=1
    plt.ylabel('Earth azimuth as seen by '+reflector_name+' [deg]')
    plt.xlabel('Time after landing [Earth days]')
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.legend(loc='best')
    plt.grid()
    plt.savefig(output_folder_path+'/Earth_azimuth_seen_'+reflector_name+'_vs_time.pdf',bbox_inches='tight')
    plt.show()
    plt.close('all')

    plt.figure(figsize=(10,7))
    plt.rcParams.update({'font.size': 16})
    colors = [plt.cm.jet(i) for i in np.linspace(0, 1, len(zones.keys()))]
    plt.gca().set_prop_cycle(plt.cycler('color', colors))
    for zones_pointer in list(zones.keys()):
        one_ground_station_per_location = 0
        for transmitter_pointer in zones[zones_pointer]:
            if one_ground_station_per_location==0:
                plt.scatter(np.mod(-np.rad2deg(data_transmitter[transmitter_pointer]['Earth azimuth'])+90,360),
                    np.rad2deg(data_transmitter[transmitter_pointer]['Earth elevation']),s=25,label=zones_pointer,color=colors[list(zones.keys()).index(zones_pointer)])
            else:
                plt.scatter(np.mod(-np.rad2deg(data_transmitter[transmitter_pointer]['Earth azimuth'])+90,360),
                    np.rad2deg(data_transmitter[transmitter_pointer]['Earth elevation']),s=25,color=colors[list(zones.keys()).index(zones_pointer)])
            one_ground_station_per_location+=1
    plt.ylabel('Earth elevation as seen by '+reflector_name+' [deg]')
    plt.xlabel('Earth azimuth as seen by '+reflector_name+' [deg]')
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.legend(loc='best')
    plt.grid()
    plt.savefig(output_folder_path+'/Earth_elevation_seen_'+reflector_name+'_vs_Earth_azimuth_seen_'+reflector_name+'.pdf',bbox_inches='tight')
    plt.show()
    plt.close('all')

    plt.figure(figsize=(10,7))
    plt.rcParams.update({'font.size': 16})
    colors = [plt.cm.jet(i) for i in np.linspace(0, 1, len(zones.keys()))]
    plt.gca().set_prop_cycle(plt.cycler('color', colors))
    for zones_pointer in list(zones.keys()):
        one_ground_station_per_location = 0
        for transmitter_pointer in zones[zones_pointer]:
            if one_ground_station_per_location==0:
                plt.scatter((np.array(data_transmitter[transmitter_pointer]['Time at transmitter'])-observation_start_epoch)/constants.JULIAN_DAY,
                    np.rad2deg(data_transmitter[transmitter_pointer]['Elevation at transmitter']),s=25,label=zones_pointer,color=colors[list(zones.keys()).index(zones_pointer)])
            else:
                plt.scatter((np.array(data_transmitter[transmitter_pointer]['Time at transmitter'])-observation_start_epoch)/constants.JULIAN_DAY,
                    np.rad2deg(data_transmitter[transmitter_pointer]['Elevation at transmitter']),s=25,color=colors[list(zones.keys()).index(zones_pointer)])
            one_ground_station_per_location+=1
    plt.ylabel('Elevation at DSN transmitter [deg]')
    plt.xlabel('Time after landing [Earth days]')
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.legend(loc='best')
    plt.grid()
    plt.savefig(output_folder_path+'/Transmitter_elevation_vs_time_transmission.pdf',bbox_inches='tight')
    plt.show()
    plt.close('all')

    plt.figure(figsize=(15,6))
    colors = [plt.cm.jet(i) for i in np.linspace(0, 1, len(zones.keys()))]
    plt.gca().set_prop_cycle(plt.cycler('color', colors))
    for zones_pointer in list(zones.keys()):
        one_ground_station_per_location = 0
        for transmitter_pointer in zones[zones_pointer]:
            if one_ground_station_per_location==0:
                plt.scatter((np.array(data_transmitter[transmitter_pointer]['Time at receiver'])-observation_start_epoch)/constants.JULIAN_DAY,
                    np.rad2deg(data_transmitter[transmitter_pointer]['Elevation at receiver']),s=5,label=zones_pointer,color=colors[list(zones.keys()).index(zones_pointer)])
            else:
                plt.scatter((np.array(data_transmitter[transmitter_pointer]['Time at receiver'])-observation_start_epoch)/constants.JULIAN_DAY,
                    np.rad2deg(data_transmitter[transmitter_pointer]['Elevation at receiver']),s=5,color=colors[list(zones.keys()).index(zones_pointer)])
            one_ground_station_per_location+=1
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
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(zones.keys())))))
    for zones_pointer in list(zones.keys()):
        transmitter_time_concatenate = list()
        transmitter_elevation_concatenate = list()
        for transmitter_pointer in zones[zones_pointer]:
            transmitter_time_concatenate.extend(data_transmitter[transmitter_pointer]['Time at transmitter'])
            transmitter_elevation_concatenate.extend(data_transmitter[transmitter_pointer]['Elevation at transmitter'])
        index_sort = np.argsort(transmitter_time_concatenate)
        transmitter_time_concatenate = [transmitter_time_concatenate[i] for i in index_sort]
        transmitter_elevation_concatenate = [transmitter_elevation_concatenate[i] for i in index_sort]
        plt.plot(np.arange(0,(transmitter_time_concatenate[-1]-observation_start_epoch)/constants.JULIAN_DAY,1),
            np.polyval(np.polyfit((np.array(transmitter_time_concatenate)-observation_start_epoch)/constants.JULIAN_DAY,
            np.rad2deg(transmitter_elevation_concatenate),4),
            np.arange(0,(transmitter_time_concatenate[-1]-observation_start_epoch)/constants.JULIAN_DAY,1)),label='Mean - '+zones_pointer)
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
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(zones.keys())))))
    for zones_pointer in list(zones.keys()):
        transmitter_time_concatenate = list()
        transmitter_elevation_concatenate = list()
        for transmitter_pointer in zones[zones_pointer]:
            transmitter_time_concatenate.extend(data_transmitter[transmitter_pointer]['Time at receiver'])
            transmitter_elevation_concatenate.extend(data_transmitter[transmitter_pointer]['Elevation at receiver'])
        index_sort = np.argsort(transmitter_time_concatenate)
        transmitter_time_concatenate = [transmitter_time_concatenate[i] for i in index_sort]
        transmitter_elevation_concatenate = [transmitter_elevation_concatenate[i] for i in index_sort]
        plt.plot(np.arange(0,(transmitter_time_concatenate[-1]-observation_start_epoch)/constants.JULIAN_DAY,1),
            np.polyval(np.polyfit((np.array(transmitter_time_concatenate)-observation_start_epoch)/constants.JULIAN_DAY,
            np.rad2deg(transmitter_elevation_concatenate),4),
            np.arange(0,(transmitter_time_concatenate[-1]-observation_start_epoch)/constants.JULIAN_DAY,1)),label='Mean - '+zones_pointer)
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
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(Earth_ground_station_list)-len(transmitter_names)))))
    for Earth_ground_station_pointer in Earth_ground_station_list:
        if Earth_ground_station_pointer[1] in transmitter_names:
            continue
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
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(Earth_ground_station_list)-len(transmitter_names)))))
    for Earth_ground_station_pointer in Earth_ground_station_list:
        if Earth_ground_station_pointer[1] in transmitter_names:
            continue
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
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(Earth_ground_station_list)-len(transmitter_names)))))
    station_names = list()
    receiver_observation_number = list()
    for Earth_ground_station_pointer in Earth_ground_station_list:
        if Earth_ground_station_pointer[1] in transmitter_names:
            continue
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
    plot1 = plt.barh(np.arange(len(Earth_ground_station_list)-len(transmitter_names)),receiver_observation_number)
    plt.yticks(np.arange(len(station_names)),labels=station_names)
    plt.xlabel('Number of Observations')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.bar_label(plot1)
    plt.grid()
    plt.savefig(output_folder_path+'/Receivers_hbar.pdf',bbox_inches='tight')
    plt.show()
    plt.close('all')

print("--- %s seconds ---" % (time.time() - run_time))