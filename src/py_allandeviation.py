"""
Description: Compute Allan Deviation

Author: C. Fortuny-Lombraña
"""

import time

run_time = time.time()
if __name__=="__main__":
    ########################################################################################################################
    ################################################## IMPORT PACKAGES #####################################################
    ########################################################################################################################

    import os
    import json
    import copy
    import matplotlib.pyplot as plt
    import allantools
    import numpy as np
    import numpy as np
    from scipy import stats
    from scipy.stats import norm
    from collections import Counter
    from scipy import stats

    integration_time = 60 #s taken from ED045 folder

    base_frequency = 8400.5*10**6 #Hz

    output_folder_path = os.path.dirname(os.path.realpath(__file__)).replace('/src','/output/InSight')
    os.makedirs(output_folder_path,exist_ok=True)

    ########################################################################################################################
    ################################################## FUNCTIONS ###########################################################
    ########################################################################################################################

    # Define a function for interquartile range (SNR) 
    def IQR_y(value,max=75,min=25):
        upper_quartile = np.percentile( np.array(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y_axis_fdets]), max)
        lower_quartile = np.percentile( np.array(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y_axis_fdets]), min)
        IQR = (upper_quartile - lower_quartile)
        quartileSet = lower_quartile - value*IQR
        
        for keys_pointer in list(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer].keys())[:6]:
            data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][keys_pointer] =  \
                list(np.array(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][keys_pointer])[
                np.where(np.array(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y_axis_fdets]) >= quartileSet)])
        pass
    
    # Define a function for interquartile range (Doppler noise) 
    def IQR_y2(value,max=75,min=25):
        upper_quartile = np.percentile( np.array(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets]), max)
        lower_quartile = np.percentile( np.array(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets]), min)
        IQR = (upper_quartile - lower_quartile)
        quartileSet = (lower_quartile - value*IQR, upper_quartile + value*IQR)

        for keys_pointer in list(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer].keys())[:6]:
            data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][keys_pointer] =  \
                list(np.array(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][keys_pointer])[
                np.where((np.array(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets]) >= quartileSet[0]) &\
                    (np.array(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets]) <= quartileSet[1]))])
        pass

    # Define a function for standard score (SNR) 
    def zscore_y(value):
        z_score = stats.zscore(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y_axis_fdets])
        z_score_index = np.where(np.abs(z_score)<=value)[0]

        for keys_pointer in list(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer].keys())[:6]:
            data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][keys_pointer] = [
                data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][keys_pointer][i] for i in z_score_index]
        pass
    
    # Define a function for standard score (Doppler noise)
    def zscore_y2(value):
        z_score = stats.zscore(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets])
        z_score_index = np.where(np.abs(z_score)<=value)[0]

        for keys_pointer in list(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer].keys())[:6]:
            data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][keys_pointer] = [
                data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][keys_pointer][i] for i in z_score_index]
        pass

    # Define a function for plotting Gaussians
    def gaussian_plot(dict_path,name):
        plt.figure()
        x_sorted = sorted(dict_path)
        plt.plot(x_sorted, norm.pdf(x_sorted, np.mean(x_sorted), np.std(x_sorted)),'b--')
        plt.title(fdets_station_pointer+' station at '+fdets_station_starttime_pointer)
        plt.grid()
        plt.xlabel(name)
        plt.show()
        plt.close('all')
        pass


    ########################################################################################################################
    ################################################## OPEN JSON FILE ######################################################
    ########################################################################################################################

    # Open the Fdets JSON file
    with open(output_folder_path+'/Fdets_data.json', 'r') as fp:
        data_fdets = json.load(fp)

    ########################################################################################################################
    ################################################## CONSTRAINTS #########################################################
    ########################################################################################################################

    # Define list with all ground stations
    ground_station_list = list(data_fdets.keys())

    # Create a dictionary for the constraints
    constraints_dict = {}
    
    # Iterate along all the ground stations, and state that constraints are not needed with False statement
    for station_pointer in ground_station_list:
        constraints_dict[station_pointer]={}
        for time_stamp_pointer in list(data_fdets[station_pointer].keys()):
            constraints_dict[station_pointer][time_stamp_pointer]={}
            constraints_dict[station_pointer][time_stamp_pointer]['IQR_bool'] = False
            constraints_dict[station_pointer][time_stamp_pointer]['z_bool'] = False

    # The False statements are changed to True when the constraints are applied
    constraints_dict['Bd']['2020-02-22 01:54:16']['IQR_bool'] = True
    constraints_dict['Bd']['2020-02-22 01:54:16']['IQR_y'] = 1.25
    constraints_dict['Bd']['2020-02-22 01:54:16']['IQR_y2'] = 1.25

    constraints_dict['Cd']['2020-02-22 01:53:25']['IQR_bool'] = True
    constraints_dict['Cd']['2020-02-22 01:53:25']['IQR_y'] = 1.25
    constraints_dict['Cd']['2020-02-22 01:53:25']['IQR_y2'] = 1.25

    constraints_dict['Hh']['2020-02-22 01:54:15']['IQR_bool'] = True
    constraints_dict['Hh']['2020-02-22 01:54:15']['IQR_y'] = 1.25
    constraints_dict['Hh']['2020-02-22 01:54:15']['IQR_y2'] = 1.25

    constraints_dict['T6']['2020-02-22 01:35:35']['IQR_bool'] = True
    constraints_dict['T6']['2020-02-22 01:35:35']['IQR_y'] = 1.25
    constraints_dict['T6']['2020-02-22 01:35:35']['IQR_y2'] = 1.25

    constraints_dict['T6']['2020-02-22 01:54:15']['IQR_bool'] = True
    constraints_dict['T6']['2020-02-22 01:54:15']['IQR_y'] = 1.25
    constraints_dict['T6']['2020-02-22 01:54:15']['IQR_y2'] = 1.25

    #OLDS
    #constraints_dict['Hh']['2020-02-22 01:37:57']['IQR_bool'] = True
    #constraints_dict['Hh']['2020-02-22 01:37:57']['IQR_y'] = 2
    #constraints_dict['Hh']['2020-02-22 01:37:57']['IQR_y2'] = 2

    #constraints_dict['Ww']['2020-02-22 01:30:10']['z_bool'] = True
    #constraints_dict['Ww']['2020-02-22 01:30:10']['z_y'] = 3
    #constraints_dict['Ww']['2020-02-22 01:30:10']['z_y2'] = 3

    constraints_dict['Ef']['2020-10-21 03:35:05']['IQR_bool'] = True
    constraints_dict['Ef']['2020-10-21 03:35:05']['IQR_y'] = 1.5
    constraints_dict['Ef']['2020-10-21 03:35:05']['IQR_y2'] = 1.5

    constraints_dict['Ir']['2020-05-29 08:40:05']['IQR_bool'] = True
    constraints_dict['Ir']['2020-05-29 08:40:05']['IQR_y'] = 1.5
    constraints_dict['Ir']['2020-05-29 08:40:05']['IQR_y2'] = 1.5

    constraints_dict['Mc']['2020-05-30 03:25:05']['IQR_bool'] = True
    constraints_dict['Mc']['2020-05-30 03:25:05']['IQR_y'] = 1.25
    constraints_dict['Mc']['2020-05-30 03:25:05']['IQR_y2'] = 1.25

    constraints_dict['Mc']['2020-05-30 03:35:15']['IQR_bool'] = True
    constraints_dict['Mc']['2020-05-30 03:35:15']['IQR_y'] = 1.5
    constraints_dict['Mc']['2020-05-30 03:35:15']['IQR_y2'] = 1.5

    constraints_dict['Mc']['2020-10-22 04:03:55']['z_bool'] = True
    constraints_dict['Mc']['2020-10-22 04:03:55']['z_y'] = 3
    constraints_dict['Mc']['2020-10-22 04:03:55']['z_y2'] = 3

    constraints_dict['Wz']['2020-05-30 03:35:15']['IQR_bool'] = True
    constraints_dict['Wz']['2020-05-30 03:35:15']['IQR_y'] = 1.5
    constraints_dict['Wz']['2020-05-30 03:35:15']['IQR_y2'] = 1.5

    constraints_dict['Wz']['2020-10-22 04:04:25']['IQR_bool'] = True
    constraints_dict['Wz']['2020-10-22 04:04:25']['IQR_y'] = 0
    constraints_dict['Wz']['2020-10-22 04:04:25']['IQR_y2'] = 0

    constraints_dict['O6']['2020-10-21 03:35:05']['IQR_bool'] = True
    constraints_dict['O6']['2020-10-21 03:35:05']['IQR_y'] = 1.5
    constraints_dict['O6']['2020-10-21 03:35:05']['IQR_y2'] = 1.5

    constraints_dict['Wb']['2020-10-21 03:35:16']['IQR_bool'] = True
    constraints_dict['Wb']['2020-10-21 03:35:16']['IQR_y'] = 1.5
    constraints_dict['Wb']['2020-10-21 03:35:16']['IQR_y2'] = 1.5

    ########################################################################################################################
    ################################################## PLOTS ###############################################################
    ########################################################################################################################

    # Choose x- and y-axis (Comment out the useful x_axis and y_axis), we can use as well 'Doppler noise [Hz]', 'Spectral max', 'Freq. detection [Hz]'
    x_axis_fdets = 'Time(UTC) [s]'
    x2_axis_fdets = 'Doppler noise [Hz]'
    y_axis_fdets = 'Signal-to-Noise ratio [dB]'
    y2_axis_fdets = 'Doppler noise [Hz]'
    y3_axis_fdets = 'Freq. detection [Hz]'

    # Original files without any automatic cut
    boolean_fdets = True
    if boolean_fdets:
        # Iterate along the ground stations inside the Fdets JSON file
        for fdets_station_pointer in data_fdets.keys():
            # Iterate along the first time stamp for each ground station
            for fdets_station_starttime_pointer in data_fdets[fdets_station_pointer].keys():
                # Print station name & starting time
                print(fdets_station_pointer+' station at '+fdets_station_starttime_pointer)

                # Create path for saving the figures
                output_figures_path = output_folder_path+'/Fdets/'+fdets_station_pointer+'/'+fdets_station_starttime_pointer.replace(':','-')
                os.makedirs(output_figures_path,exist_ok=True)

                # Plot Guassian plots
                #gaussian_plot(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y_axis_fdets],y_axis_fdets)

                #gaussian_plot(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets],y2_axis_fdets)

                # First plot: y_axis_fdets vs x_axis_fdets
                plt.figure()
                plt.plot(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][x_axis_fdets],\
                    data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y_axis_fdets],'bo-')
                plt.xlabel(x_axis_fdets)
                plt.ylabel(y_axis_fdets)
                plt.title(fdets_station_pointer+' station at '+fdets_station_starttime_pointer)
                plt.grid()
                plt.savefig(output_figures_path+'/SNR_vs_time_initial.pdf',bbox_inches="tight")
                plt.show()
                plt.close('all')

                # Second plot: y_axis_fdets vs x2_axis_fdets
                plt.figure()
                plt.plot(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][x2_axis_fdets],\
                    data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y_axis_fdets],'rx')
                plt.xlabel(x2_axis_fdets)
                plt.ylabel(y_axis_fdets)
                plt.title(fdets_station_pointer+' station at '+fdets_station_starttime_pointer)
                plt.grid()
                plt.savefig(output_figures_path+'/SNR_vs_Doppler_noise_initial.pdf',bbox_inches="tight")
                plt.show()
                plt.close('all')

                # Third plot: y2_axis_fdets vs x_axis_fdets
                plt.figure()
                plt.plot(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][x_axis_fdets],\
                    data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets],'ko-')
                plt.xlabel(x_axis_fdets)
                plt.ylabel(y2_axis_fdets)
                plt.title(fdets_station_pointer+' station at '+fdets_station_starttime_pointer)
                plt.grid()
                plt.savefig(output_figures_path+'/Doppler_noise_vs_time_initial.pdf',bbox_inches="tight")
                plt.show()
                plt.close('all')

                # Fourth plot: y3_axis_fdets vs x_axis_fdets
                plt.figure()
                plt.plot(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][x_axis_fdets],\
                    data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y3_axis_fdets],'go-')
                plt.xlabel(x_axis_fdets)
                plt.ylabel(y3_axis_fdets)
                plt.title(fdets_station_pointer+' station at '+fdets_station_starttime_pointer)
                plt.grid()
                plt.savefig(output_figures_path+'/Freq_vs_time_initial.pdf',bbox_inches="tight")
                plt.show()
                plt.close('all')

                # Fifth plot: Allan deviation using the y2_axis_fdets
                plt.figure()
                rate_fdets = 1/Counter(np.diff(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][x_axis_fdets])).most_common(1)[0][0]
                taus2, adevs, errors, ns = allantools.mdev(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets]/(np.array(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y3_axis_fdets])+base_frequency),\
                    rate = rate_fdets, data_type = 'freq',taus='decade')
                #y = np.random.normal(np.mean(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets]),np.std(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets]),
                #        size=len(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets]))
                #taus2_white, adevs_white, errors_white, ns_white = allantools.mdev(y, rate=rate_fdets, data_type="freq",taus='decade')
                adevs_white = list()
                adevs_white.append(adevs[0])
                for i in range(1,len(taus2)):
                    adevs_white.append(10**(-1/2*(np.log10(taus2[i])-np.log10(taus2[i-1]))+np.log10(adevs_white[i-1]))) 
                plot1 = plt.errorbar(taus2,adevs,yerr=errors,ecolor='green')
                plot2 = plt.errorbar(taus2, adevs_white,ecolor='red')
                plt.xscale('log')
                plt.yscale('log')
                plt.xlabel('Tau [s]')
                plt.ylabel('Allan Deviation')
                plt.legend([plot1,plot2],['Real Measurements','White Noise'])
                plt.title(fdets_station_pointer+' station at '+fdets_station_starttime_pointer)
                plt.grid()
                plt.axis('equal')
                plt.savefig(output_figures_path+'/allan_deviation.pdf',bbox_inches="tight")
                plt.show()
                plt.close('all')
                
                # Compute the correlation, mean slope and the sigma at integration time
                slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(taus2), np.log10(adevs))
                slope_white, intercept_white, r_value_white, p_value_white, std_err_white = stats.linregress(np.log10(taus2), np.log10(adevs_white))
                correlation = np.corrcoef(adevs,adevs_white)[0][1]
                print("Correlation: ",correlation)
                print("Slope Allan Deviation: ",slope)
                print("Slope White Allan Deviation: ",slope_white)

    boolean_fdets = True
    if boolean_fdets:
        # Make a deep copy of the data_fdets dictionary
        data_fdets_updated = copy.deepcopy(data_fdets)
        # Iterate along the ground stations inside the Fdets JSON file
        for fdets_station_pointer in data_fdets.keys():
            # Iterate along the first time stamp for each ground station
            for fdets_station_starttime_pointer in data_fdets[fdets_station_pointer].keys():
                # Print station name & starting time
                print(fdets_station_pointer+' station at '+fdets_station_starttime_pointer)

                # Remove the list values for making list of lists
                types_datafdets = [type(k) for k in data_fdets[fdets_station_pointer][fdets_station_starttime_pointer].values()]
                for types_pointer in range(0,len(types_datafdets)):
                    if types_datafdets[types_pointer] == list:
                        key = list(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer].keys())[types_pointer]
                        del data_fdets_updated[fdets_station_pointer][fdets_station_starttime_pointer][key] 
                        data_fdets_updated[fdets_station_pointer][fdets_station_starttime_pointer][key] = list()

                # Plot Guassian plots
                #gaussian_plot(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y_axis_fdets],y_axis_fdets)

                #gaussian_plot(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets],y2_axis_fdets)

                # Statistical constraints: IQR & z_score
                if constraints_dict[fdets_station_pointer][fdets_station_starttime_pointer]['IQR_bool']:
                    IQR_y(constraints_dict[fdets_station_pointer][fdets_station_starttime_pointer]['IQR_y'])
                    IQR_y2(constraints_dict[fdets_station_pointer][fdets_station_starttime_pointer]['IQR_y2'])

                if constraints_dict[fdets_station_pointer][fdets_station_starttime_pointer]['z_bool']:
                    zscore_y(constraints_dict[fdets_station_pointer][fdets_station_starttime_pointer]['z_y'])
                    zscore_y2(constraints_dict[fdets_station_pointer][fdets_station_starttime_pointer]['z_y2'])

                # Identify time jumps
                time_diff = np.round(np.diff(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][x_axis_fdets]),decimals=10)
                time_jumps = np.where(time_diff!=Counter(time_diff).most_common(1)[0][0])[0]

                # Save the indices deleted from data_fdets
                delete_indices = []

                # Create empty lists for saving the results of the Allan deviation
                taus_total = list()
                adevs_total = list()
                errors_total = list()
                ns_total = list()
                slope_total = list()
                slope_white_total = list()
                correlation_total = list()
                sigma_total = list()
                white_boolean = list()

                # Operation to identify the start index and end index
                jump_label = -1
                for jump_pointer in range(0,len(time_jumps)+1):
                    jump_label += 1
                    if jump_pointer == 0:
                        start_index = 0
                    else:
                        start_index = time_jumps[jump_pointer-1]+1
                    
                    if jump_pointer == len(time_jumps):
                        end_index = len(time_diff)+1
                    else: 
                        end_index = time_jumps[jump_pointer]

                    if end_index-start_index < 10:
                        print("Warning:",end_index-start_index)
                        if jump_pointer == len(time_jumps): 
                            delete_indices.extend(list(range(start_index,end_index)))
                        else:
                            delete_indices.extend(list(range(start_index,end_index+1)))
                        print(delete_indices)
                        jump_label -= 1
                        continue
                    
                    # Create path for saving the figures
                    output_figures_path = output_folder_path+'/Fdets/'+fdets_station_pointer+'/'+fdets_station_starttime_pointer.replace(':','-')+'/Part_'+str(jump_label)
                    os.makedirs(output_figures_path,exist_ok=True)

                    # First plot: y_axis_fdets vs x_axis_fdets
                    plt.figure()
                    plt.rcParams.update({'font.size': 16})
                    plt.plot(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][x_axis_fdets][start_index:end_index],\
                        data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y_axis_fdets][start_index:end_index],'bo-')
                    plt.xlabel(x_axis_fdets)
                    plt.ylabel(y_axis_fdets)
                    #plt.title(fdets_station_pointer+' station at '+fdets_station_starttime_pointer+' - Part: '+str(jump_label))
                    plt.grid()
                    plt.savefig(output_figures_path+'/SNR_vs_time_'+fdets_station_pointer+'.png',bbox_inches="tight")
                    plt.show()
                    plt.close('all')
                    
                    # Second plot: y_axis_fdets vs x2_axis_fdets
                    plt.figure()
                    plt.plot(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][x2_axis_fdets][start_index:end_index],\
                        data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y_axis_fdets][start_index:end_index],'rx')
                    plt.xlabel(x2_axis_fdets)
                    plt.ylabel(y_axis_fdets)
                    plt.title(fdets_station_pointer+' station at '+fdets_station_starttime_pointer+' - Part: '+str(jump_label))
                    plt.grid()
                    plt.savefig(output_figures_path+'/SNR_vs_Doppler_noise.pdf',bbox_inches="tight")
                    plt.show()
                    plt.close('all')
                
                    # Third plot: y2_axis_fdets vs x_axis_fdets
                    plt.figure()
                    plt.plot(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][x_axis_fdets][start_index:end_index],\
                        data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets][start_index:end_index],'ko-')
                    plt.xlabel(x_axis_fdets)
                    plt.ylabel(y2_axis_fdets)
                    #plt.title(fdets_station_pointer+' station at '+fdets_station_starttime_pointer+' - Part: '+str(jump_label))
                    plt.grid()
                    plt.savefig(output_figures_path+'/Doppler_noise_vs_time_'+fdets_station_pointer+'.png',bbox_inches="tight")
                    plt.show()
                    plt.close('all')

                    # Fourth plot: y3_axis_fdets vs x_axis_fdets
                    plt.figure()
                    plt.plot(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][x_axis_fdets][start_index:end_index],\
                        data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y3_axis_fdets][start_index:end_index],'ko-')
                    plt.xlabel(x_axis_fdets)
                    plt.ylabel(y3_axis_fdets)
                    #plt.title(fdets_station_pointer+' station at '+fdets_station_starttime_pointer+' - Part: '+str(jump_label))
                    plt.grid()
                    plt.savefig(output_figures_path+'/Freq_vs_time.pdf',bbox_inches="tight")
                    plt.show()
                    plt.close('all')
                    
                    # Fifth plot: Allan deviation using the y2_axis_fdets
                    plt.figure()
                    rate_fdets = 1/Counter(np.diff(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][x_axis_fdets][start_index:end_index])).most_common(1)[0][0]
                    taus2, adevs, errors, ns = allantools.mdev(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets][start_index:end_index]/(np.array(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y3_axis_fdets][start_index:end_index])+base_frequency),\
                        rate = rate_fdets, data_type = 'freq',taus='decade')
                    #y = np.random.normal(np.mean(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets]),np.std(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets]),
                    #        size=len(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets]))
                    #taus2_white, adevs_white, errors_white, ns_white = allantools.mdev(y, rate=rate_fdets, data_type="freq",taus='decade')
                    adevs_white = list()
                    adevs_white.append(adevs[0])
                    for i in range(1,len(taus2)):
                        adevs_white.append(10**(-1/2*(np.log10(taus2[i])-np.log10(taus2[i-1]))+np.log10(adevs_white[i-1])))
                    plot1 = plt.errorbar(taus2,adevs,yerr=errors,ecolor='green')
                    plot2 = plt.errorbar(taus2,adevs_white,ecolor='red')
                    plt.xscale('log')
                    plt.yscale('log')
                    plt.xlabel('Tau [s]')
                    plt.ylabel('Allan Deviation')
                    plt.legend([plot1,plot2],['Observed','White Noise'])
                    #plt.title(fdets_station_pointer+' station at '+fdets_station_starttime_pointer+' - Part: '+str(jump_label))
                    plt.grid()
                    plt.axis('equal')
                    plt.savefig(output_figures_path+'/allan_deviation_'+fdets_station_pointer+'.png',bbox_inches="tight")
                    plt.show()
                    plt.close('all')

                    # Compute the correlation, mean slope and the sigma at integration time
                    slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(taus2), np.log10(adevs))
                    slope_white, intercept_white, r_value_white, p_value_white, std_err_white = stats.linregress(np.log10(taus2), np.log10(adevs_white))
                    correlation = np.corrcoef(adevs,adevs_white)[0][1]
                    print("Correlation: ",correlation)
                    print("Slope Allan Deviation: ",slope)
                    print("Slope White Allan Deviation: ",slope_white)
                    if taus2[-1] < integration_time:
                        print("Sigma: 0")
                    else:
                        print("Sigma: ",allantools.mdev(np.array(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets][start_index:end_index])/(np.array(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y3_axis_fdets][start_index:end_index])+base_frequency),\
                        rate = rate_fdets, data_type = 'freq',taus=integration_time)[1][0])
                     
                    #slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(taus2), np.log10(adevs))
                    #y = np.random.normal(np.mean(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets][start_index:end_index]),np.std(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets][start_index:end_index]),
                    #    size=len(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets][start_index:end_index]))
                    #taus2_white, adevs_white, errors_white, ns_white = allantools.mdev(y, rate=rate_fdets, data_type="freq",taus='decade')
                    #slope_white, intercept_white, r_value_white, p_value_white, std_err_white = stats.linregress(np.log10(taus2_white), np.log10(adevs_white))
                    #slopes_white.append(slope_white)
                    #correlation = np.corrcoef(adevs,adevs_white)[0][1]
                    #correlations.append(correlation)
                    #print("Correlation: ",np.mean(correlations))
                    #print("Slope Allan Deviation: ",slope)
                    #print("Slope White Allan Deviation: ",np.mean(slopes_white))
                    #plot2 = plt.errorbar(taus2_white, adevs_white,yerr=errors_white,ecolor='red')
                    #plt.xscale('log')
                    #plt.yscale('log')
                    #plt.xlabel('Tau [s]')
                    #plt.ylabel('Allan Deviation')
                    #plt.legend([plot1,plot2],['Real-data','Simulated White Noise'])
                    #plt.title(fdets_station_pointer+' station at '+fdets_station_starttime_pointer+' - Part: '+str(jump_label))
                    #plt.grid()
                    #plt.axis('equal')
                    #plt.savefig(output_figures_path+'/allan_deviation.pdf',bbox_inches="tight")
                    #plt.show()
                    #plt.close('all')
                        
                    # Appending Allan deviation results
                    taus_total.append(taus2.tolist())
                    adevs_total.append(adevs.tolist())
                    errors_total.append(errors.tolist())
                    ns_total.append(ns.tolist())
                    slope_total.append([slope])
                    slope_white_total.append([slope_white])
                    correlation_total.append([correlation])
                    if taus2[-1] < integration_time:
                        sigma_total.append([0])
                    else:
                        sigma_total.append([allantools.mdev(np.array(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets][start_index:end_index])/(np.array(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y3_axis_fdets][start_index:end_index])+base_frequency),\
                        rate = rate_fdets, data_type = 'freq',taus=integration_time)[1][0]])
                    if np.mean(correlation)>0.6 and np.abs(slope_white-slope) <0.75 and slope<0:
                        white_boolean.append([True])
                    else:
                        white_boolean.append([False])
                    print(white_boolean[-1])
                    
                    # Saving the lists of lists in the new dictionary
                    for types_pointer in range(0,len(types_datafdets)):
                        if types_datafdets[types_pointer] == list:
                            key = list(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer].keys())[types_pointer]
                            data_fdets_updated[fdets_station_pointer][fdets_station_starttime_pointer][key].append(
                                data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][key][start_index:end_index]
                            )
                
                # Saving Allan deviation results
                data_fdets_updated[fdets_station_pointer][fdets_station_starttime_pointer]['Taus [s]']=taus_total
                data_fdets_updated[fdets_station_pointer][fdets_station_starttime_pointer]['Allan Variance']=adevs_total
                data_fdets_updated[fdets_station_pointer][fdets_station_starttime_pointer]['Estimated Errors']=errors_total
                data_fdets_updated[fdets_station_pointer][fdets_station_starttime_pointer]['Number of Pairs']=ns_total
                data_fdets_updated[fdets_station_pointer][fdets_station_starttime_pointer]['Slope Allan Deviation']=slope_total
                data_fdets_updated[fdets_station_pointer][fdets_station_starttime_pointer]['Slope White Allan Deviation']=slope_white_total
                data_fdets_updated[fdets_station_pointer][fdets_station_starttime_pointer]['Correlation']=correlation_total
                data_fdets_updated[fdets_station_pointer][fdets_station_starttime_pointer]['Sigma']=sigma_total
                data_fdets_updated[fdets_station_pointer][fdets_station_starttime_pointer]['White Noise']=white_boolean

        # Save the updated data in a JSON file
        with open(output_folder_path+'/Fdets_data_updated.json', 'w') as fp:
            json.dump(data_fdets_updated, fp)             

    ########################################################################################################################
    ################################################## OPEN JSON FILE ######################################################
    ########################################################################################################################               

    # Open the Phases JSON file
    with open(output_folder_path+'/Phases_data.json', 'r') as fp:
        data_phases = json.load(fp)

    ########################################################################################################################
    ################################################## PLOTS ###############################################################
    ########################################################################################################################

    # Choose x- and y-axis (Comment out the useful x_axis and y_axis)
    x_axis_phases = 'Time stamp [s]'
    y_axis_phases = 'Phase [rad]'
    y2_axis_phases = 'Phase [s]'
    y3_axis_phases = 'Frequency [Hz]'

    boolean_phases = True
    if boolean_phases:
        # Iterate along the ground stations inside the Phases JSON file
        for phases_station_pointer in data_phases.keys():
            # Iterate along the first time stamp for each ground station
            for phases_station_starttime_pointer in data_phases[phases_station_pointer].keys():
                # Print station name & starting time
                print(phases_station_pointer+' station at '+phases_station_starttime_pointer)

                # Create path for saving the figures
                output_figures_path = output_folder_path+'/Phases/'+phases_station_pointer+'/'+phases_station_starttime_pointer.replace(':','-')
                os.makedirs(output_figures_path,exist_ok=True)

                # First plot: y_axis_phases vs x_axis_phases 
                plt.figure()
                plt.plot(data_phases[phases_station_pointer][phases_station_starttime_pointer][x_axis_phases],\
                    data_phases[phases_station_pointer][phases_station_starttime_pointer][y_axis_phases],'b')
                plt.xlabel(x_axis_phases)
                plt.ylabel(y_axis_phases)
                plt.title(phases_station_pointer+' station at '+phases_station_starttime_pointer)
                plt.grid()
                plt.savefig(output_figures_path+'/Phase_rad_vs_time.pdf',bbox_inches="tight")
                plt.show()
                plt.close('all')

                # Compute the rate
                rate_phases = 1/Counter(np.diff(data_phases[phases_station_pointer][phases_station_starttime_pointer][x_axis_phases])).most_common(1)[0][0]
                
                # Convert the phase from radians to seconds
                data_phases[phases_station_pointer][phases_station_starttime_pointer][y2_axis_phases] = list(data_phases[phases_station_pointer][phases_station_starttime_pointer][y_axis_phases]/(2*np.pi*rate_phases))
                 
                # Second plot: y2_axis_phases vs x_axis_phases
                plt.figure()
                plt.plot(data_phases[phases_station_pointer][phases_station_starttime_pointer][x_axis_phases],\
                    data_phases[phases_station_pointer][phases_station_starttime_pointer][y2_axis_phases],'r')
                plt.xlabel(x_axis_phases)
                plt.ylabel(y2_axis_phases)
                plt.title(phases_station_pointer+' station at '+phases_station_starttime_pointer)
                plt.grid()
                plt.savefig(output_figures_path+'/Phase_sec_vs_time.pdf',bbox_inches="tight")
                plt.show()
                plt.close('all')

                # Check whether the conversion is correct: from radians to seconds
                #plt.figure()
                #plt.plot(data_phases[phases_station_pointer][phases_station_starttime_pointer][x_axis_phases],\
                #    allantools.phase2radians(data_phases[phases_station_pointer][phases_station_starttime_pointer]['Phase [s]'],rate_phases),'g')
                #plt.xlabel(x_axis_phases)
                #plt.ylabel(y_axis_phases)
                #plt.title(phases_station_pointer+' station at '+phases_station_starttime_pointer)
                #plt.grid()
                #plt.show()
                #plt.close('all')

                # Convert phase to frequency, NOTE that the frequency list has an element less than the phase list
                data_phases[phases_station_pointer][phases_station_starttime_pointer][y3_axis_phases] = \
                    list(allantools.phase2frequency(data_phases[phases_station_pointer][phases_station_starttime_pointer][y2_axis_phases],rate_phases))
                
                # Fourth plot: y3_axis_phases vs x_axis_phases
                plt.figure()
                plt.plot(data_phases[phases_station_pointer][phases_station_starttime_pointer][x_axis_phases][:-1],
                    data_phases[phases_station_pointer][phases_station_starttime_pointer][y3_axis_phases],'k')
                plt.xlabel(x_axis_phases)
                plt.ylabel(y3_axis_phases)
                plt.title(phases_station_pointer+' station at '+phases_station_starttime_pointer)
                plt.grid()
                plt.savefig(output_figures_path+'/frequency_vs_time.pdf',bbox_inches="tight")
                plt.show()
                plt.close('all')

                # Third plot: Allan deviation using the y3_axis_phases
                plt.figure()
                taus2, adevs, errors, ns = allantools.mdev(np.array(data_phases[phases_station_pointer][phases_station_starttime_pointer][y3_axis_phases])/base_frequency,\
                    rate = rate_phases, data_type = 'freq',taus='decade')
                #y = np.random.normal(np.mean(data_phases[phases_station_pointer][phases_station_starttime_pointer][y3_axis_phases]),np.std(data_phases[phases_station_pointer][phases_station_starttime_pointer][y3_axis_phases]),
                #    size = len(data_phases[phases_station_pointer][phases_station_starttime_pointer][y3_axis_phases]))
                #taus2_white, adevs_white, errors_white, ns_white = allantools.mdev(y, rate=rate_phases, data_type="freq",taus='decade')
                #adevs_white = 10**(-1/2*np.log(taus2[0])+np.log(adevs[0])-1/2*np.log(taus2))
                adevs_white = list()
                adevs_white.append(adevs[0])
                for i in range(1,len(taus2)):
                    adevs_white.append(10**(-1/2*(np.log10(taus2[i])-np.log10(taus2[i-1]))+np.log10(adevs_white[i-1])))        
                plot1 = plt.errorbar(taus2,adevs,yerr=errors,ecolor='green')
                plot2 =plt.errorbar(taus2, adevs_white ,ecolor='red')
                plt.xscale('log')
                plt.yscale('log')
                plt.xlabel('Tau [s]')
                plt.ylabel('Allan Deviation')
                plt.legend([plot1,plot2],['Real Measurements','White Noise'])
                plt.title(phases_station_pointer+' station at '+phases_station_starttime_pointer)
                plt.grid()
                plt.axis('equal')
                plt.savefig(output_figures_path+'/allan_deviation.pdf',bbox_inches="tight")
                plt.show()
                plt.close('all')

                # Compute the correlation, mean slope and the sigma at integration time
                slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(taus2), np.log10(adevs))
                slope_white, intercept_white, r_value_white, p_value_white, std_err_white = stats.linregress(np.log10(taus2), np.log10(adevs_white))
                correlation = np.corrcoef(adevs,adevs_white)[0][1]
                print("Correlation: ",correlation)
                print("Slope Allan Deviation: ",slope)
                print("Slope White Allan Deviation: ",slope_white)
                if taus2[-1] < integration_time:
                    print("Sigma: 0")
                else:
                    print("Sigma: ",allantools.mdev(np.array(data_phases[phases_station_pointer][phases_station_starttime_pointer][y3_axis_phases])/base_frequency,\
                    rate = rate_phases, data_type = 'freq',taus=integration_time)[1][0])

                # Saving Allan deviation results
                data_phases[phases_station_pointer][phases_station_starttime_pointer]['Correlation'] = correlation
                data_phases[phases_station_pointer][phases_station_starttime_pointer]['Slope Allan Deviation'] = slope
                data_phases[phases_station_pointer][phases_station_starttime_pointer]['Slope White Allan Deviation'] = slope_white
                if taus2[-1] < integration_time:
                    data_phases[phases_station_pointer][phases_station_starttime_pointer]['Sigma'] = 0
                else:
                    data_phases[phases_station_pointer][phases_station_starttime_pointer]['Sigma'] = allantools.mdev(np.array(data_phases[phases_station_pointer][phases_station_starttime_pointer][y3_axis_phases])/base_frequency,\
                    rate = rate_phases, data_type = 'freq',taus=integration_time)[1][0]

        # Save the updated data in a JSON file
        with open(output_folder_path+'/Phases_data_updated.json', 'w') as fp:
            json.dump(data_phases, fp) 

print("--- %s seconds ---" % (time.time() - run_time))