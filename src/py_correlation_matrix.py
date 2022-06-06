"""
Description: Plot v7v

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
    import matplotlib.pyplot as plt
    from multiprocessing import Pool
    import scipy.interpolate
    import scipy.sparse
    import scipy.sparse.linalg

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

    # Load files
    weights_matrix_diagonal = np.loadtxt(output_folder_path+"/weights_diagonal.dat")
    estimation_information_matrix = np.loadtxt(output_folder_path+"/estimation_information_matrix.dat")
    estimation_information_matrix_normalization = np.loadtxt(output_folder_path+"/estimation_information_matrix_normalization.dat")
    concatenated_times = np.loadtxt(output_folder_path+"/concatenated_times.dat")
    concatenated_link_ends = np.loadtxt(output_folder_path+"/concatenated_link_ends.dat")
    concatenated_link_end_names_list = np.loadtxt(output_folder_path+"/concatenated_link_end_names.dat")
    doppler_residuals = np.loadtxt(output_folder_path+"/doppler_residuals.dat")
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
            residuals_sort_short = [None]*count_time_obs
            for row_index in range(0,np.shape(estimation_information_matrix_sort_short)[0]):
                estimation_information_matrix_sort_short[row_index] = estimation_information_matrix_sort[start_index+receiver_link_ends[row_index]]
                residuals_sort_short[row_index] = residuals_sort[start_index+receiver_link_ends[row_index]]
            estimation_information_matrix_sort[start_index:end_index] = estimation_information_matrix_sort_short
            residuals_sort[start_index:end_index] = residuals_sort_short
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
    
    # Delete terms
    for index_delete in indices_delete[::-1]:
        estimation_information_matrix_sort.pop(index_delete)
        residuals_sort.pop(index_delete)
        concatenated_times_sort.pop(index_delete)
        concatenated_link_ends_sort.pop(index_delete)
        concatenated_link_end_names_list_sort.pop(index_delete)
        vector_weights_sort.pop(index_delete)

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
    if PRIDE_boolean==False or correlation==0:
        partial_cov = np.transpose(estimation_information_matrix_sort)@(inv_weight_complex_total.sqrt())



    estimation_information_matrix_sort = np.loadtxt(output_folder_path+"/estimation_information_matrix_sort.dat")
    concatenated_times_sort = np.loadtxt(output_folder_path+"/concatenated_times_sort.dat")
    concatenated_times_no_duplicated = np.loadtxt(output_folder_path+"/concatenated_times_sort_no_duplicated.dat")
    concatenated_times_index = np.loadtxt(output_folder_path+"/concatenated_times_sort_index.dat")
    concatenated_times_count = np.loadtxt(output_folder_path+"/concatenated_times_sort_count.dat")
    concatenated_link_ends_sort = np.loadtxt(output_folder_path+"/concatenated_link_ends_sort.dat")
    concatenated_link_end_names_list_sort = np.loadtxt(output_folder_path+"/concatenated_link_end_names_sort.dat")
    residuals_sort = np.loadtxt(output_folder_path+"/doppler_residuals_sort.dat")
    vector_weights_sort = np.loadtxt(output_folder_path+"/vector_weights_sort.dat")

    # Normalized inverse a priori covariance
    norm_inverse_a_priori_covariance = np.diag(inverse_a_priori_covariance.diagonal()/(estimation_information_matrix_normalization**2))

print("--- %s seconds ---" % (time.time() - run_time))