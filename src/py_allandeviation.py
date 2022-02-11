"""
Description: Allan Deviation

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
    import matplotlib.pyplot as plt
    import allantools
    import numpy as np
    import numpy as np
    from scipy import stats
    from scipy.stats import norm
    from collections import Counter

    output_folder_path = os.path.dirname(os.path.realpath(__file__)).replace('/src','/output/InSight')
    os.makedirs(output_folder_path,exist_ok=True)

    ########################################################################################################################
    ################################################## PLOTS ###############################################################
    ########################################################################################################################

    # Open the Fdets JSON file
    with open(output_folder_path+'/Fdets_data.json', 'r') as fp:
        data_fdets = json.load(fp)

    # Choose x- and y-axis (Comment out the useful x_axis and y_axis), we can use as well 'Doppler noise [Hz]', 'Spectral max', 'Freq. detection [Hz]'
    x_axis_fdets = 'Time(UTC) [s]'
    x2_axis_fdets = 'Doppler noise [Hz]'
    y_axis_fdets = 'Signal-to-Noise ratio [dB]'
    y2_axis_fdets = 'Doppler noise [Hz]'

    boolean_fdets = True
    if boolean_fdets:
        # Iterate along the ground stations inside the Fdets JSON file
        for fdets_station_pointer in ['Wb']:#data_fdets.keys():
            # Iterate along the first time stamp for each ground station
            for fdets_station_starttime_pointer in ['2020-10-22 04:04:46']:#data_fdets[fdets_station_pointer].keys():
                print(fdets_station_pointer+' station at '+fdets_station_starttime_pointer)

                x_sorted = sorted(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y_axis_fdets])
                plt.plot(x_sorted, norm.pdf(x_sorted, np.mean(x_sorted), np.std(x_sorted)),'b--')
                plt.title(fdets_station_pointer+' station at '+fdets_station_starttime_pointer)
                plt.grid()
                plt.xlabel(y_axis_fdets)
                plt.show()

                x_sorted = sorted(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets])
                plt.plot(x_sorted, norm.pdf(x_sorted, np.mean(x_sorted), np.std(x_sorted)),'k--')
                plt.title(fdets_station_pointer+' station at '+fdets_station_starttime_pointer)
                plt.grid()
                plt.xlabel(y2_axis_fdets)
                plt.show()

                # First plot: y_axis_fdets vs x_axis_fdets
                plt.plot(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][x_axis_fdets],\
                    data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y_axis_fdets],'bo-')
                plt.xlabel(x_axis_fdets)
                plt.ylabel(y_axis_fdets)
                plt.title(fdets_station_pointer+' station at '+fdets_station_starttime_pointer)
                plt.grid()
                plt.show()

                # Second plot: y_axis_fdets vs x2_axis_fdets
                plt.plot(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][x2_axis_fdets],\
                    data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y_axis_fdets],'rx')
                plt.xlabel(x2_axis_fdets)
                plt.ylabel(y_axis_fdets)
                plt.title(fdets_station_pointer+' station at '+fdets_station_starttime_pointer)
                plt.grid()
                plt.show()

                # Third plot: y2_axis_fdets vs x_axis_fdets
                plt.plot(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][x_axis_fdets],\
                    data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets],'ko-')
                plt.xlabel(x_axis_fdets)
                plt.ylabel(y2_axis_fdets)
                plt.title(fdets_station_pointer+' station at '+fdets_station_starttime_pointer)
                plt.grid()
                plt.show()

                # Allan deviation using the y2_axis_fdets
                rate_fdets = 1/Counter(np.diff(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][x_axis_fdets])).most_common(1)[0][0]
                taus2, adevs, errors, ns = allantools.mdev(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets],\
                    rate = rate_fdets, data_type = 'freq',taus='decade')
                y = allantools.noise.white(len(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets]),b0=np.std(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets]),fs=rate_fdets)
                taus2_white, adevs_white, errors_white, ns_white = allantools.mdev(y, rate=rate_fdets, data_type="freq",taus='decade')
                plot1 = plt.errorbar(taus2,adevs,yerr=errors,ecolor='green')
                plot2 = plt.errorbar(taus2_white, adevs_white,yerr=errors_white,ecolor='red')
                plt.xscale('log')
                plt.yscale('log')
                plt.xlabel('Tau [s]')
                plt.ylabel('Allan Deviation')
                plt.legend([plot1,plot2],['Real-data','Simulated White Noise'])
                plt.title(fdets_station_pointer+' station at '+fdets_station_starttime_pointer)
                plt.grid()
                plt.axis('equal')
                plt.show()
                

    boolean_fdets1 = True
    if boolean_fdets1:
        # Iterate along the ground stations inside the Fdets JSON file
        for fdets_station_pointer in ['Wb']:#data_fdets.keys():
            # Iterate along the first time stamp for each ground station
            for fdets_station_starttime_pointer in ['2020-10-22 04:04:46']:#data_fdets[fdets_station_pointer].keys():
                print(fdets_station_pointer+' station at '+fdets_station_starttime_pointer)
                x_sorted = sorted(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y_axis_fdets])
                plt.plot(x_sorted, norm.pdf(x_sorted, np.mean(x_sorted), np.std(x_sorted)),'b--')
                plt.title(fdets_station_pointer+' station at '+fdets_station_starttime_pointer)
                plt.grid()
                plt.show()

                x_sorted = sorted(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets])
                plt.plot(x_sorted, norm.pdf(x_sorted, np.mean(x_sorted), np.std(x_sorted)),'k--')
                plt.title(fdets_station_pointer+' station at '+fdets_station_starttime_pointer)
                plt.grid()
                plt.show()

                time_diff = np.round(np.diff(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][x_axis_fdets]),decimals=10)

                print(len(x_sorted))
                '''
                upper_quartile = np.percentile( np.array(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y_axis_fdets]), 75)
                lower_quartile = np.percentile( np.array(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y_axis_fdets]), 25)
                IQR = (upper_quartile - lower_quartile)
                quartileSet = lower_quartile - 1.5*IQR
                
                data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets] =  \
                    np.array(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets])[
                    np.where(np.array(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y_axis_fdets]) >= quartileSet)]

                data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][x_axis_fdets] =  \
                    np.array(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][x_axis_fdets])[
                    np.where(np.array(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y_axis_fdets]) >= quartileSet)]

                data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y_axis_fdets] =  \
                    np.array(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y_axis_fdets])[
                    np.where(np.array(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y_axis_fdets]) >= quartileSet)]
                
                
                upper_quartile = np.percentile( np.array(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets]), 75)
                lower_quartile = np.percentile( np.array(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets]), 25)
                IQR = (upper_quartile - lower_quartile)
                quartileSet = (lower_quartile - 1.5*IQR, upper_quartile + 1.5*IQR)
                
                data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y_axis_fdets] =  \
                    np.array(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y_axis_fdets])[
                    np.where((np.array(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets]) >= quartileSet[0]) &\
                        (np.array(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets]) <= quartileSet[1]))]

                data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][x_axis_fdets] =  \
                    np.array(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][x_axis_fdets])[
                    np.where((np.array(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets]) >= quartileSet[0]) &\
                        (np.array(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets]) <= quartileSet[1]))]

                data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets] =  \
                    np.array(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets])[
                    np.where((np.array(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets]) >= quartileSet[0]) &\
                        (np.array(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets]) <= quartileSet[1]))]
                '''
                '''
                z_score = stats.zscore(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y_axis_fdets])
                z_score_index = np.where(np.abs(z_score)<=3)[0]

                data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y_axis_fdets] = [
                    data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y_axis_fdets][i] for i in z_score_index]

                data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets] = [
                    data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets][i] for i in z_score_index]

                data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][x_axis_fdets] = [
                    data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][x_axis_fdets][i] for i in z_score_index]

                z_score = stats.zscore(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets])
                z_score_index = np.where(np.abs(z_score)<=3)[0]

                data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y_axis_fdets] = [
                    data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y_axis_fdets][i] for i in z_score_index]

                data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets] = [
                    data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets][i] for i in z_score_index]

                data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][x_axis_fdets] = [
                    data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][x_axis_fdets][i] for i in z_score_index]
                '''
                x_sorted = sorted(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y_axis_fdets])
                plt.plot(x_sorted, norm.pdf(x_sorted, np.mean(x_sorted), np.std(x_sorted)),'b--')
                plt.title(fdets_station_pointer+' station at '+fdets_station_starttime_pointer)
                plt.grid()
                plt.show()

                x_sorted = sorted(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets])
                plt.plot(x_sorted, norm.pdf(x_sorted, np.mean(x_sorted), np.std(x_sorted)),'k--')
                plt.title(fdets_station_pointer+' station at '+fdets_station_starttime_pointer)
                plt.grid()
                plt.show()

                time_diff = np.round(np.diff(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][x_axis_fdets]),decimals=10)
                time_jumps = np.where(time_diff!=Counter(time_diff).most_common(1)[0][0])[0]

                print(time_diff)
                print(len(x_sorted))

                for jump_pointer in range(0,len(time_jumps)+1):
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
                        continue

                    # First plot: y_axis_fdets vs x_axis_fdets
                    plt.plot(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][x_axis_fdets][start_index:end_index],\
                        data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y_axis_fdets][start_index:end_index],'bo-')
                    plt.xlabel(x_axis_fdets)
                    plt.ylabel(y_axis_fdets)
                    plt.title(fdets_station_pointer+' station at '+fdets_station_starttime_pointer+' - Part: '+str(jump_pointer))
                    plt.grid()
                    plt.show()

                    # Second plot: y_axis_fdets vs x2_axis_fdets
                    plt.plot(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][x2_axis_fdets][start_index:end_index],\
                        data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y_axis_fdets][start_index:end_index],'rx')
                    plt.xlabel(x2_axis_fdets)
                    plt.ylabel(y_axis_fdets)
                    plt.title(fdets_station_pointer+' station at '+fdets_station_starttime_pointer+' - Part: '+str(jump_pointer))
                    plt.grid()
                    plt.show()

                    # Third plot: y2_axis_fdets vs x_axis_fdets
                    plt.plot(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][x_axis_fdets][start_index:end_index],\
                        data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets][start_index:end_index],'ko-')
                    plt.xlabel(x_axis_fdets)
                    plt.ylabel(y2_axis_fdets)
                    plt.title(fdets_station_pointer+' station at '+fdets_station_starttime_pointer+' - Part: '+str(jump_pointer))
                    plt.grid()
                    plt.show()

                    # Allan deviation using the y2_axis_fdets
                    rate_fdets = 1/Counter(np.diff(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][x_axis_fdets][start_index:end_index])).most_common(1)[0][0]
                    taus2, adevs, errors, ns = allantools.mdev(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets][start_index:end_index],\
                        rate = rate_fdets, data_type = 'freq',taus='decade')
                    y = allantools.noise.white(len(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets][start_index:end_index]),b0=np.std(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y2_axis_fdets][start_index:end_index]),fs=rate_fdets)
                    taus2_white, adevs_white, errors_white, ns_white = allantools.mdev(y, rate=rate_fdets, data_type="freq",taus='decade')
                    plot1 = plt.errorbar(taus2,adevs,yerr=errors,ecolor='green')
                    plot2 = plt.errorbar(taus2_white, adevs_white,yerr=errors_white,ecolor='red')
                    plt.xscale('log')
                    plt.yscale('log')
                    plt.xlabel('Tau [s]')
                    plt.ylabel('Allan Deviation')
                    plt.legend([plot1,plot2],['Real-data','Simulated White Noise'])
                    plt.title(fdets_station_pointer+' station at '+fdets_station_starttime_pointer+' - Part: '+str(jump_pointer))
                    plt.grid()
                    plt.axis('equal')
                    plt.show()
                
                    

    # Open the Phases JSON file
    with open(output_folder_path+'/Phases_data.json', 'r') as fp:
        data_phases = json.load(fp)

    # Choose x- and y-axis (Comment out the useful x_axis and y_axis)
    x_axis_phases = 'Time stamp [s]'
    y_axis_phases = 'Phase [rad]'

    boolean_phases = False
    if boolean_phases:
        # Iterate along the ground stations inside the Phases JSON file
        for phases_station_pointer in data_phases.keys():
            # Iterate along the first time stamp for each ground station
            for phases_station_starttime_pointer in data_phases[phases_station_pointer].keys():
                # First plot: y_axis_phases vs x_axis_phases 
                plt.plot(data_phases[phases_station_pointer][phases_station_starttime_pointer][x_axis_phases],\
                    data_phases[phases_station_pointer][phases_station_starttime_pointer][y_axis_phases],'b')
                plt.xlabel(x_axis_phases)
                plt.ylabel(y_axis_phases)
                plt.title(phases_station_pointer+' station at '+phases_station_starttime_pointer)
                plt.grid()
                plt.show()

                rate_phases = 1/np.diff(data_phases[phases_station_pointer][phases_station_starttime_pointer][x_axis_phases])[0]
                freq = allantools.phase2frequency(data_phases[phases_station_pointer][phases_station_starttime_pointer][y_axis_phases],rate_phases)

                #plt.plot(freq,\
                #    data_phases[phases_station_pointer][phases_station_starttime_pointer][y_axis_phases],'k')
                #plt.xlabel(x_axis_phases)
                #plt.ylabel('Freq [Hz')
                #plt.title(phases_station_pointer+' station at '+phases_station_starttime_pointer)
                #plt.grid()
                #plt.show()

                # Allan deviation using x_axis_phases
                taus2, adevs, errors, ns = allantools.adev(freq,rate = rate_phases, data_type = 'freq')
                plt.loglog(taus2, adevs,'go-')
                plt.xlabel('Tau')
                plt.ylabel('Allan Deviation')
                plt.title(phases_station_pointer+' station at '+phases_station_starttime_pointer)
                plt.grid()
                plt.axis('equal')
                plt.show()

print("--- %s seconds ---" % (time.time() - run_time))