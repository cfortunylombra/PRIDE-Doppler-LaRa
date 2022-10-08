""" 
Description: Plot correlation matrices

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
    import matplotlib.pyplot as plt
    from multiprocessing import Pool
    import scipy.interpolate
    import scipy.sparse
    import scipy.sparse.linalg
    from tudatpy.kernel import constants
    from tudatpy.kernel.interface import spice_interface
    from tudatpy.kernel.astro import element_conversion
    from tudatpy.kernel.numerical_simulation import environment_setup
    from tudatpy.kernel.numerical_simulation.estimation_setup import observation

    CPU_par = 4

    # Booleans to understand whether we want to simulate together RISE and LaRa missions, or separetely
    RISE_boolean = True
    LaRa_boolean = True
    
    if LaRa_boolean:
        # PRIDE stations boolean
        PRIDE_boolean = True
        
        # Boolean for removing PRIDE stations only after the POD
        remove_PRIDE_weight_boolean = False

        # Define fixed correlation between PRIDE and DSN stations
        correlation = 0
    else:
        # PRIDE stations boolean
        PRIDE_boolean = False
        
        # Boolean for removing PRIDE stations only after the POD 
        remove_PRIDE_weight_boolean = False
        
        # Define fixed correlation between PRIDE and DSN stations
        correlation = 0

    # Evaluation step 
    step_eval = 1

    # Output folder
    if LaRa_boolean:
        output_folder_path = os.path.dirname(os.path.realpath(__file__)).replace('/src','/output/POD_RISE'+str(RISE_boolean)+'_LaRa'+str(LaRa_boolean)+'_PRIDE'+str(PRIDE_boolean)+str(remove_PRIDE_weight_boolean)+'_corr'+str(correlation))
    else:
        output_folder_path = os.path.dirname(os.path.realpath(__file__)).replace('/src','/output/POD_RISE'+str(RISE_boolean)+'_LaRa'+str(LaRa_boolean))

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

    observation_RISE_start = 596587943.0
    observation_LaRa_start = 694495488.0
    observation_end = 851524433.0

    correlation_time_eval = (observation_LaRa_start+observation_end)/2

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
            if name_ground_station=="HART15M" or name_ground_station=="HOBART12":
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

    # Create list of observation settings
    observation_settings_list = list()

    # Define link ends for RISE
    if RISE_boolean:
        for pointer_RISE_transmitter in RISE_transmitter_names:
            two_way_link_ends = dict()
            two_way_link_ends[observation.transmitter] = ("Earth",pointer_RISE_transmitter)
            two_way_link_ends[observation.reflector1] = ("Mars",RISE_reflector_name)
            two_way_link_ends[observation.receiver] = ("Earth",pointer_RISE_transmitter)

            observation_settings_list.append(two_way_link_ends)
        
    # Number of link ends for InSight mission
    RISE_link_ends_length = len(observation_settings_list)

    # Define link ends for LaRa
    if LaRa_boolean:
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
        if transmitter_pointer[1].startswith("DSS") and LaRa_boolean:
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
        if transmitter_pointer[1].startswith("DSS") and RISE_boolean:
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


    # Load files
    weights_matrix_diagonal = np.loadtxt(output_folder_path+"/weights_diagonal.dat")
    estimation_information_matrix = np.loadtxt(output_folder_path+"/estimation_information_matrix.dat")
    estimation_information_matrix_normalization = np.loadtxt(output_folder_path+"/estimation_information_matrix_normalization.dat")
    concatenated_times = np.loadtxt(output_folder_path+"/concatenated_times.dat")
    concatenated_link_ends = np.loadtxt(output_folder_path+"/concatenated_link_ends.dat")
    f = open(output_folder_path+"/concatenated_link_end_names.dat",'r')
    concatenated_link_end_names_list = []
    for line in f.readlines():
        concatenated_link_end_names_list.append(line.replace('\n',''))
    f.close()
    vector_weights = np.loadtxt(output_folder_path+"/vector_weights.dat")

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
        return estimation_information_matrix[index][:]

    # Sort concatenated times
    sort_estimation_information_matrix_dict = Pool(CPU_par)
    print("Sorting estimation information matrices")
    estimation_information_matrix_sort = sort_estimation_information_matrix_dict.map(sort_estimation_information_matrix_func,index_sort)
    sort_estimation_information_matrix_dict.close()
    sort_estimation_information_matrix_dict.join()

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

    # Initialize inverted weighting matrix
    #inv_weight_complex = (scipy.sparse.diags(1/np.array(vector_weights_sort)**2)).tocsr() #Same as correlation to 1
    inv_weight_complex = scipy.sparse.coo_matrix((0,0))
    
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
            for row_index in range(0,np.shape(estimation_information_matrix_sort_short)[0]):
                estimation_information_matrix_sort_short[row_index] = estimation_information_matrix_sort[start_index+receiver_link_ends[row_index]]
            estimation_information_matrix_sort[start_index:end_index] = estimation_information_matrix_sort_short
            concatenated_link_end_names_list_sort[start_index:end_index] = np.array(concatenated_link_end_names_list_sort[start_index:end_index])[receiver_link_ends]
            vector_weights_sort[start_index:end_index] = np.array(vector_weights_sort[start_index:end_index])[receiver_link_ends]
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
                        # When transmitter !=r receiver -> delete index
                        else:
                            indices_delete.append(start_index_block+row_index)

                # When all the PRIDE observations are taken into account
                else:
                    # Number of receivers for a transmitter
                    split_count = transmitter_count_number[0]
                    block_weight = (scipy.sparse.coo_matrix((split_count,split_count))).toarray()
                    for row_index in range(0,split_count):
                        for column_index in range(0,split_count):
                            # Diagonal element
                            if row_index == column_index:
                                block_weight[row_index][column_index] = vector_weights_sort[start_index_block+row_index]**2
                            # Non-diagonal element
                            else:
                                block_weight[row_index][column_index] = correlation*vector_weights_sort[start_index_block+row_index]*vector_weights_sort[start_index_block+column_index]
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
                            for row_index in range(0,split_count):
                                for column_index in range(0,split_count):
                                    # Diagonal element
                                    if row_index == column_index:
                                        block_weight[row_index][column_index] = vector_weights_sort[start_index_block+row_index]**2
                                    # Non-diagonal element
                                    else:
                                        block_weight[row_index][column_index] = correlation*vector_weights_sort[start_index_block+row_index]*vector_weights_sort[start_index_block+column_index]
                            inv_weight_complex = scipy.sparse.block_diag((inv_weight_complex,np.linalg.inv(block_weight)))

                    # Remove all the observations from other transmitters
                    else:
                        indices_delete.extend(range(start_index_block,start_index_block+split_count))

                    # Increasing the start_index_block and total_split
                    start_index_block += split_count
                    total_split+=split_count
    
    # Convert the inverted weighting matrix to Compressed Sparse Row format
    inv_weight_complex_total = (inv_weight_complex).tocsr()

    estimation_information_matrix_sort = np.loadtxt(output_folder_path+"/estimation_information_matrix_sort.dat")
    concatenated_times_sort = np.loadtxt(output_folder_path+"/concatenated_times_sort.dat")
    concatenated_times_no_duplicated = np.loadtxt(output_folder_path+"/concatenated_times_sort_no_duplicated.dat")
    concatenated_times_index = [int(i) for i in np.loadtxt(output_folder_path+"/concatenated_times_sort_index.dat")]
    concatenated_times_count = [int(i) for i in np.loadtxt(output_folder_path+"/concatenated_times_sort_count.dat")]
    concatenated_link_ends_sort = np.loadtxt(output_folder_path+"/concatenated_link_ends_sort.dat")
    f = open(output_folder_path+"/concatenated_link_end_names_sort.dat",'r')
    concatenated_link_end_names_list_sort = []
    for line in f.readlines():
        concatenated_link_end_names_list_sort.append(line.replace('\n',''))
    f.close()
    vector_weights_sort = np.loadtxt(output_folder_path+"/vector_weights_sort.dat")

    # A priori
    if RISE_boolean and LaRa_boolean:
        add_par = 3
    else:
        add_par = 0 
    apriori_vector = np.zeros(39+add_par)
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

    # Step evaluation
    arange_eval = np.arange(0,len(concatenated_times_no_duplicated),step_eval) 
    time_eval = [concatenated_times_no_duplicated[i] for i in arange_eval]

    arange_eval = [np.abs(correlation_time_eval-np.array(time_eval)).argmin()]
    time_eval = [time_eval[arange_eval[0]]]

    # Define a priori covariance
    inverse_a_priori_covariance = np.diag(1/apriori_vector**2)

    # Normalized inverse a priori covariance
    norm_inverse_a_priori_covariance = np.diag(inverse_a_priori_covariance.diagonal()/(estimation_information_matrix_normalization**2))

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
                   (estimation_information_matrix_normalization[i]*estimation_information_matrix_normalization[j])
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

    if RISE_boolean and LaRa_boolean and correlation_time_eval>observation_LaRa_start:
        index_start_LaRa = np.abs(observation_LaRa_start-np.array(concatenated_times_sort)).argmin()
        
        inv_norm_covariance_LaRa_value = np.linalg.inv(np.transpose(estimation_information_matrix_sort[index_start_LaRa:concatenated_times_index[arange_eval[0]]+1])@inv_weight_complex_total[index_start_LaRa:concatenated_times_index[arange_eval[0]]+1,index_start_LaRa:concatenated_times_index[arange_eval[0]]+1]@estimation_information_matrix_sort[index_start_LaRa:concatenated_times_index[arange_eval[0]]+1]\
                +norm_inverse_a_priori_covariance)

        covariance_matrix_LaRa = np.zeros(np.shape(inv_norm_covariance_LaRa_value))
        for i in range(0,np.shape(inv_norm_covariance_LaRa_value)[0]):
            for j in range(0,np.shape(inv_norm_covariance_LaRa_value)[1]):
                covariance_matrix_LaRa[i][j] = inv_norm_covariance_LaRa_value[i][j]/\
                    (estimation_information_matrix_normalization[i]*estimation_information_matrix_normalization[j])

        sigma_LaRa_value = np.sqrt(covariance_matrix_LaRa.diagonal())

        correlation_matrix_LaRa = np.zeros(np.shape(covariance_matrix_LaRa))
        for i in range(0,np.shape(covariance_matrix_LaRa)[0]):
            for j in range(0,np.shape(covariance_matrix_LaRa)[1]):
                correlation_matrix_LaRa[i][j] =covariance_matrix_LaRa[i][j]/\
                    (sigma_LaRa_value[i]*sigma_LaRa_value[j])
    
    # Correlation matrix (start from RISE)
    plt.figure(figsize=(18,18))
    plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
    plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True
    plt.imshow(np.abs(correlation_values[-1]))
    plt.colorbar()
    if LaRa_boolean==False:
        plt.xticks(range(0,len(apriori_vector)),labels=['$x$','$y$','$z$',r'$\dot{x}$',r'$\dot{y}$',r'$\dot{z}$','F',
            r'$\sigma_{FCN}$',r'$x_{{RISE}}$',r'$y_{{RISE}}$',r'$z_{{RISE}}$',
            r'$\phi^c_1$',r'$\phi^s_1$',r'$\phi^c_2$',r'$\phi^s_2$',r'$\phi^c_3$',r'$\phi^s_3$',r'$\phi^c_4$',r'$\phi^s_4$',
            r'$Xp^c_1$',r'$Xp^s_1$',r'$Yp^c_1$',r'$Yp^s_1$',r'$Xp^c_2$',r'$Xp^s_2$',r'$Yp^c_2$',r'$Yp^s_2$',
            r'$Xp^c_3$',r'$Xp^s_3$',r'$Yp^c_3$',r'$Yp^s_3$',r'$Xp^c_4$',r'$Xp^s_4$',r'$Yp^c_4$',r'$Yp^s_4$',
            r'$Xp^c_5$',r'$Xp^s_5$',r'$Yp^c_5$',r'$Yp^s_5$'])
        plt.yticks(range(0,len(apriori_vector)),labels=['$x$','$y$','$z$',r'$\dot{x}$',r'$\dot{y}$',r'$\dot{z}$','F',
            r'$\sigma_{FCN}$',r'$x_{{RISE}}$',r'$y_{{RISE}}$',r'$z_{{RISE}}$',
            r'$\phi^c_1$',r'$\phi^s_1$',r'$\phi^c_2$',r'$\phi^s_2$',r'$\phi^c_3$',r'$\phi^s_3$',r'$\phi^c_4$',r'$\phi^s_4$',
            r'$Xp^c_1$',r'$Xp^s_1$',r'$Yp^c_1$',r'$Yp^s_1$',r'$Xp^c_2$',r'$Xp^s_2$',r'$Yp^c_2$',r'$Yp^s_2$',
            r'$Xp^c_3$',r'$Xp^s_3$',r'$Yp^c_3$',r'$Yp^s_3$',r'$Xp^c_4$',r'$Xp^s_4$',r'$Yp^c_4$',r'$Yp^s_4$',
            r'$Xp^c_5$',r'$Xp^s_5$',r'$Yp^c_5$',r'$Yp^s_5$'])
    elif RISE_boolean==False:
        plt.xticks(range(0,len(apriori_vector)),labels=['$x$','$y$','$z$',r'$\dot{x}$',r'$\dot{y}$',r'$\dot{z}$','F',
            r'$\sigma_{FCN}$',r'$x_{{LaRa}}$',r'$y_{{LaRa}}$',r'$z_{{LaRa}}$',
            r'$\phi^c_1$',r'$\phi^s_1$',r'$\phi^c_2$',r'$\phi^s_2$',r'$\phi^c_3$',r'$\phi^s_3$',r'$\phi^c_4$',r'$\phi^s_4$',
            r'$Xp^c_1$',r'$Xp^s_1$',r'$Yp^c_1$',r'$Yp^s_1$',r'$Xp^c_2$',r'$Xp^s_2$',r'$Yp^c_2$',r'$Yp^s_2$',
            r'$Xp^c_3$',r'$Xp^s_3$',r'$Yp^c_3$',r'$Yp^s_3$',r'$Xp^c_4$',r'$Xp^s_4$',r'$Yp^c_4$',r'$Yp^s_4$',
            r'$Xp^c_5$',r'$Xp^s_5$',r'$Yp^c_5$',r'$Yp^s_5$'])
        plt.yticks(range(0,len(apriori_vector)),labels=['$x$','$y$','$z$',r'$\dot{x}$',r'$\dot{y}$',r'$\dot{z}$','F',
            r'$\sigma_{FCN}$',r'$x_{{LaRa}}$',r'$y_{{LaRa}}$',r'$z_{{LaRa}}$',
            r'$\phi^c_1$',r'$\phi^s_1$',r'$\phi^c_2$',r'$\phi^s_2$',r'$\phi^c_3$',r'$\phi^s_3$',r'$\phi^c_4$',r'$\phi^s_4$',
            r'$Xp^c_1$',r'$Xp^s_1$',r'$Yp^c_1$',r'$Yp^s_1$',r'$Xp^c_2$',r'$Xp^s_2$',r'$Yp^c_2$',r'$Yp^s_2$',
            r'$Xp^c_3$',r'$Xp^s_3$',r'$Yp^c_3$',r'$Yp^s_3$',r'$Xp^c_4$',r'$Xp^s_4$',r'$Yp^c_4$',r'$Yp^s_4$',
            r'$Xp^c_5$',r'$Xp^s_5$',r'$Yp^c_5$',r'$Yp^s_5$'])
    else:
        plt.xticks(range(0,len(apriori_vector)),labels=['$x$','$y$','$z$',r'$\dot{x}$',r'$\dot{y}$',r'$\dot{z}$','F',
            r'$\sigma_{FCN}$',r'$x_{{RISE}}$',r'$y_{{RISE}}$',r'$z_{{RISE}}$',r'$x_{{LaRa}}$',r'$y_{{LaRa}}$',r'$z_{{LaRa}}$',
            r'$\phi^c_1$',r'$\phi^s_1$',r'$\phi^c_2$',r'$\phi^s_2$',r'$\phi^c_3$',r'$\phi^s_3$',r'$\phi^c_4$',r'$\phi^s_4$',
            r'$Xp^c_1$',r'$Xp^s_1$',r'$Yp^c_1$',r'$Yp^s_1$',r'$Xp^c_2$',r'$Xp^s_2$',r'$Yp^c_2$',r'$Yp^s_2$',
            r'$Xp^c_3$',r'$Xp^s_3$',r'$Yp^c_3$',r'$Yp^s_3$',r'$Xp^c_4$',r'$Xp^s_4$',r'$Yp^c_4$',r'$Yp^s_4$',
            r'$Xp^c_5$',r'$Xp^s_5$',r'$Yp^c_5$',r'$Yp^s_5$'])
        plt.yticks(range(0,len(apriori_vector)),labels=['$x$','$y$','$z$',r'$\dot{x}$',r'$\dot{y}$',r'$\dot{z}$','F',
            r'$\sigma_{FCN}$',r'$x_{{RISE}}$',r'$y_{{RISE}}$',r'$z_{{RISE}}$',r'$x_{{LaRa}}$',r'$y_{{LaRa}}$',r'$z_{{LaRa}}$',
            r'$\phi^c_1$',r'$\phi^s_1$',r'$\phi^c_2$',r'$\phi^s_2$',r'$\phi^c_3$',r'$\phi^s_3$',r'$\phi^c_4$',r'$\phi^s_4$',
            r'$Xp^c_1$',r'$Xp^s_1$',r'$Yp^c_1$',r'$Yp^s_1$',r'$Xp^c_2$',r'$Xp^s_2$',r'$Yp^c_2$',r'$Yp^s_2$',
            r'$Xp^c_3$',r'$Xp^s_3$',r'$Yp^c_3$',r'$Yp^s_3$',r'$Xp^c_4$',r'$Xp^s_4$',r'$Yp^c_4$',r'$Yp^s_4$',
            r'$Xp^c_5$',r'$Xp^s_5$',r'$Yp^c_5$',r'$Yp^s_5$'])
    plt.grid()
    plt.title('Correlation Matrix - Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.savefig(output_folder_path+"/correlation_matrix_"+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0]))+".pdf",bbox_inches="tight")
    plt.show()
    plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = True
    plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = False
    plt.close('all')

    # Correlation matrix (start from LaRa)
    if LaRa_boolean and RISE_boolean and correlation_time_eval>observation_LaRa_start:
        plt.figure(figsize=(18,18))
        plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
        plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True
        plt.imshow(np.abs(correlation_matrix_LaRa))
        plt.colorbar()
        if RISE_boolean==False:
            plt.xticks(range(0,len(apriori_vector)),labels=['$x$','$y$','$z$',r'$\dot{x}$',r'$\dot{y}$',r'$\dot{z}$','F',
                r'$\sigma_{FCN}$',r'$x_{{LaRa}}$',r'$y_{{LaRa}}$',r'$z_{{LaRa}}$',
                r'$\phi^c_1$',r'$\phi^s_1$',r'$\phi^c_2$',r'$\phi^s_2$',r'$\phi^c_3$',r'$\phi^s_3$',r'$\phi^c_4$',r'$\phi^s_4$',
                r'$Xp^c_1$',r'$Xp^s_1$',r'$Yp^c_1$',r'$Yp^s_1$',r'$Xp^c_2$',r'$Xp^s_2$',r'$Yp^c_2$',r'$Yp^s_2$',
                r'$Xp^c_3$',r'$Xp^s_3$',r'$Yp^c_3$',r'$Yp^s_3$',r'$Xp^c_4$',r'$Xp^s_4$',r'$Yp^c_4$',r'$Yp^s_4$',
                r'$Xp^c_5$',r'$Xp^s_5$',r'$Yp^c_5$',r'$Yp^s_5$'])
            plt.yticks(range(0,len(apriori_vector)),labels=['$x$','$y$','$z$',r'$\dot{x}$',r'$\dot{y}$',r'$\dot{z}$','F',
                r'$\sigma_{FCN}$',r'$x_{{LaRa}}$',r'$y_{{LaRa}}$',r'$z_{{LaRa}}$',
                r'$\phi^c_1$',r'$\phi^s_1$',r'$\phi^c_2$',r'$\phi^s_2$',r'$\phi^c_3$',r'$\phi^s_3$',r'$\phi^c_4$',r'$\phi^s_4$',
                r'$Xp^c_1$',r'$Xp^s_1$',r'$Yp^c_1$',r'$Yp^s_1$',r'$Xp^c_2$',r'$Xp^s_2$',r'$Yp^c_2$',r'$Yp^s_2$',
                r'$Xp^c_3$',r'$Xp^s_3$',r'$Yp^c_3$',r'$Yp^s_3$',r'$Xp^c_4$',r'$Xp^s_4$',r'$Yp^c_4$',r'$Yp^s_4$',
                r'$Xp^c_5$',r'$Xp^s_5$',r'$Yp^c_5$',r'$Yp^s_5$'])
        else:
            plt.xticks(range(0,len(apriori_vector)),labels=['$x$','$y$','$z$',r'$\dot{x}$',r'$\dot{y}$',r'$\dot{z}$','F',
                r'$\sigma_{FCN}$',r'$x_{{RISE}}$',r'$y_{{RISE}}$',r'$z_{{RISE}}$',r'$x_{{LaRa}}$',r'$y_{{LaRa}}$',r'$z_{{LaRa}}$',
                r'$\phi^c_1$',r'$\phi^s_1$',r'$\phi^c_2$',r'$\phi^s_2$',r'$\phi^c_3$',r'$\phi^s_3$',r'$\phi^c_4$',r'$\phi^s_4$',
                r'$Xp^c_1$',r'$Xp^s_1$',r'$Yp^c_1$',r'$Yp^s_1$',r'$Xp^c_2$',r'$Xp^s_2$',r'$Yp^c_2$',r'$Yp^s_2$',
                r'$Xp^c_3$',r'$Xp^s_3$',r'$Yp^c_3$',r'$Yp^s_3$',r'$Xp^c_4$',r'$Xp^s_4$',r'$Yp^c_4$',r'$Yp^s_4$',
                r'$Xp^c_5$',r'$Xp^s_5$',r'$Yp^c_5$',r'$Yp^s_5$'])
            plt.yticks(range(0,len(apriori_vector)),labels=['$x$','$y$','$z$',r'$\dot{x}$',r'$\dot{y}$',r'$\dot{z}$','F',
                r'$\sigma_{FCN}$',r'$x_{{RISE}}$',r'$y_{{RISE}}$',r'$z_{{RISE}}$',r'$x_{{LaRa}}$',r'$y_{{LaRa}}$',r'$z_{{LaRa}}$',
                r'$\phi^c_1$',r'$\phi^s_1$',r'$\phi^c_2$',r'$\phi^s_2$',r'$\phi^c_3$',r'$\phi^s_3$',r'$\phi^c_4$',r'$\phi^s_4$',
                r'$Xp^c_1$',r'$Xp^s_1$',r'$Yp^c_1$',r'$Yp^s_1$',r'$Xp^c_2$',r'$Xp^s_2$',r'$Yp^c_2$',r'$Yp^s_2$',
                r'$Xp^c_3$',r'$Xp^s_3$',r'$Yp^c_3$',r'$Yp^s_3$',r'$Xp^c_4$',r'$Xp^s_4$',r'$Yp^c_4$',r'$Yp^s_4$',
                r'$Xp^c_5$',r'$Xp^s_5$',r'$Yp^c_5$',r'$Yp^s_5$'])
        plt.grid()
        plt.title('LaRa Correlation Matrix - Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
        plt.savefig(output_folder_path+"/correlation_matrix_LaRa_"+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0]))+".pdf",bbox_inches="tight")
        plt.show()
        plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = True
        plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = False
        plt.close('all')


print("--- %s seconds ---" % (time.time() - run_time))