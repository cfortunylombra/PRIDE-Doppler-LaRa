"""
Description: Save the minimum Allan Deviation for InSight observation span

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
    from scipy import stats
    from scipy.stats import norm
    from collections import Counter
    import statistics

    output_folder_path = os.path.dirname(os.path.realpath(__file__)).replace('/src','/output/POD_RISETrue_LaRaTrue_PRIDEcomplexTrueFalse')
    os.makedirs(output_folder_path,exist_ok=True)

    correlation_list = np.loadtxt(output_folder_path+"/correlation_list.dat")
    weather_list = np.loadtxt(output_folder_path+"/weather_list.dat")
    indicator_Doppler_list = np.loadtxt(output_folder_path+"/indicator_Doppler_list.dat")
    indicator_weather_list = np.loadtxt(output_folder_path+"/indicator_weather_list.dat")
    Doppler_total_list = np.loadtxt(output_folder_path+"/Doppler_total_list.dat")
    weather_total_list = np.loadtxt(output_folder_path+"/weather_total_list.dat")

    # Histogram correlation
    
    plt.figure(figsize=(8, 5))
    plt.rcParams.update({'font.size': 14})
    plt.hist(correlation_list, bins = 25)
    average_correlation = np.mean(correlation_list)
    plt.axvline(x=average_correlation, color='k', linestyle='--')
    plt.axvline(x=statistics.median(correlation_list), color='r', linestyle='--')
    plt.ylabel('Count [-]')
    plt.xlabel(r'$\rho$ [-]')
    plt.grid()
    plt.savefig(output_folder_path+"/correlation.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(8, 5))
    plt.rcParams.update({'font.size': 14})
    plt.hist(weather_total_list, bins = 25)
    plt.ylabel('Count [-]')
    plt.xlabel(r'$\rho$ [-]')
    plt.grid()
    plt.savefig(output_folder_path+"/weather.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(8, 5))
    plt.rcParams.update({'font.size': 14})
    plt.hist(Doppler_total_list, bins = 25)
    plt.ylabel('Count [-]')
    plt.xlabel(r'$\rho$ [-]')
    plt.grid()
    plt.savefig(output_folder_path+"/Doppler.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    print('Mean = '+str(average_correlation))
    print('Median = '+str(statistics.median(correlation_list)))

    fig=plt.figure(figsize=(8,5))
    plt.rcParams.update({'font.size': 14})
    ax=fig.add_subplot(111, label="1")
    ax2=fig.add_subplot(111, label="2", frame_on=False)
    
    ax.hist(correlation_list, bins = 25, color="C0")
    ax.set_xlabel(r"$\rho_{ij}$ [-]", color="C0")
    ax.set_ylabel("Count [-]", color="C0")
    ax.set_xlim([0, 1])
    ax.tick_params(axis='x', colors="C0")
    ax.tick_params(axis='y', colors="C0")

    ax2.hist(weather_list, bins = 25, color="C1")
    ax2.xaxis.tick_top()
    ax2.yaxis.tick_right()
    ax2.set_xlabel(r'$w_{weather}$ [-]', color="C1") 
    ax2.set_ylabel('Count [-]', color="C1")       
    ax2.xaxis.set_label_position('top') 
    ax2.yaxis.set_label_position('right') 
    ax2.set_xlim([0, 1])
    ax2.tick_params(axis='x', colors="C1")
    ax2.tick_params(axis='y', colors="C1")

    plt.savefig(output_folder_path+"/correlation_combined.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(8, 5))
    plt.rcParams.update({'font.size': 14})
    plt.hist(Doppler_total_list, bins = 20,alpha=0.75, label="Doppler")
    plt.hist(weather_total_list, bins = 50,alpha=0.75, label="weather")
    plt.ylabel('Count [-]')
    plt.xlabel(r'$w \cdot \rho$ [-]')
    plt.legend()
    plt.grid()
    plt.savefig(output_folder_path+"/Doppler_weather.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    fig=plt.figure(figsize=(8,5))
    plt.rcParams.update({'font.size': 14})
    ax=fig.add_subplot(111, label="1")
    ax2=fig.add_subplot(111, label="2", frame_on=False)
    
    ax.hist(Doppler_total_list, bins = 25, color="C0")
    ax.set_xlabel(r"$w_{Doppler}$ $\cdot$ $\rho_{Doppler}$ [-]", color="C0")
    ax.set_ylabel("Count [-]", color="C0")
    ax.set_xlim([0, 1])
    ax.tick_params(axis='x', colors="C0")
    ax.tick_params(axis='y', colors="C0")

    ax2.hist(weather_total_list, bins = 25, color="C1")
    ax2.xaxis.tick_top()
    ax2.yaxis.tick_right()
    ax2.set_xlabel(r'$w_{weather}$ $\cdot$ $\rho_{weather}$ [-]', color="C1") 
    ax2.set_ylabel('Count [-]', color="C1")       
    ax2.xaxis.set_label_position('top') 
    ax2.yaxis.set_label_position('right') 
    ax2.set_xlim([0, 1])
    ax2.tick_params(axis='x', colors="C1")
    ax2.tick_params(axis='y', colors="C1")

    plt.savefig(output_folder_path+"/Doppler_weather.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    with open(output_folder_path+'/correlations.json', 'r') as fp:
        correlation_dict = json.load(fp)

    LaRa_stations = ['DSS 43','DSS 63','DSS 14','MEDICINA','WETTZELL','ONSALA60','EFLSBERG','WRT0','YEBES40M','TIANMA65','CEDUNA','BADARY','HARTRAO']
    correlation_station_matrix = np.zeros((len(LaRa_stations),len(LaRa_stations)))
    for i_station in LaRa_stations:
        for j_station in LaRa_stations:
            print(i_station,j_station)
            if i_station!=j_station and (len(correlation_dict[i_station][j_station])+len(correlation_dict[j_station][i_station]))!=0:
                print(LaRa_stations.index(i_station),LaRa_stations.index(j_station))
                correlation_station_matrix[LaRa_stations.index(i_station)][LaRa_stations.index(j_station)] =sum(correlation_dict[i_station][j_station]+correlation_dict[j_station][i_station])/(len(correlation_dict[i_station][j_station])+len(correlation_dict[j_station][i_station]))
            else:
                correlation_station_matrix[LaRa_stations.index(i_station)][LaRa_stations.index(j_station)] = np.nan

    plt.figure(figsize=(7,7))
    #plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
    #plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True
    plt.imshow(correlation_station_matrix)
    plt.colorbar(label=r"$\rho$ [-]",shrink=0.75,pad=0.02)
    plt.clim(0,1)
    plt.xticks(range(0,len(LaRa_stations)),labels=LaRa_stations,rotation=90)
    plt.yticks(range(0,len(LaRa_stations)),labels=LaRa_stations)
    plt.savefig(output_folder_path+"/correlations_stations.pdf",bbox_inches="tight")
    plt.show()
    #plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = True
    #plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = False
    plt.close('all')

    LaRa_mop = np.loadtxt(output_folder_path+"/LaRa_mop_correlation.dat")
    RISE_mop = np.loadtxt(output_folder_path+"/RISE_mop_correlation.dat")
    final_mop = np.loadtxt(output_folder_path+"/final_mop_correlation.dat")

    plt.figure(figsize=(7.5,7.5))
    plt.rcParams.update({'font.size': 12})
    #plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
    #plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True
    plt.imshow(np.abs(LaRa_mop))
    plt.colorbar(label=r"$\rho$ [-]",shrink=0.75,pad=0.02)
    plt.xticks(range(0,len(LaRa_mop[0])),labels=['','$y$','',r'$\dot{x}$',r'',r'$\dot{z}$','',
            r'$\sigma_{FCN}$',r'',r'$y_{{RISE}}$',r'',r'$x_{{LaRa}}$',r'',r'$z_{{LaRa}}$',
            r'',r'$\phi^s_1$',r'',r'$\phi^s_2$',r'',r'$\phi^s_3$',r'',r'$\phi^s_4$',
            r'',r'${X_p}^s_1$',r'',r'${Y_p}^s_1$',r'',r'${X_p}^s_2$',r'',r'${Y_p}^s_2$',
            r'',r'${X_p}^s_3$',r'',r'${Y_p}^s_3$',r'',r'${X_p}^s_4$',r'',r'${Y_p}^s_4$',
            r'',r'${X_p}^s_5$',r'',r'${Y_p}^s_5$'],rotation=90)
    plt.yticks(range(0,len(LaRa_mop[0])),labels=['$x$','','$z$',r'',r'$\dot{y}$',r'','F',
            r'',r'$x_{{RISE}}$',r'',r'$z_{{RISE}}$',r'',r'$y_{{LaRa}}$',r'',
            r'$\phi^c_1$',r'',r'$\phi^c_2$',r'',r'$\phi^c_3$',r'',r'$\phi^c_4$',r'',
            r'${X_p}^c_1$',r'',r'${Y_p}^c_1$',r'',r'${X_p}^c_2$',r'',r'${Y_p}^c_2$',r'',
            r'${X_p}^c_3$',r'',r'${Y_p}^c_3$',r'',r'${X_p}^c_4$',r'',r'${Y_p}^c_4$',r'',
            r'${X_p}^c_5$',r'',r'${Y_p}^c_5$',r''])
    plt.savefig(output_folder_path+"/LaRa_small_correlation_abs.pdf",bbox_inches="tight")
    plt.show()
    #plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = True
    #plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = False
    plt.close('all')

    plt.figure(figsize=(7.5,7.5))
    plt.rcParams.update({'font.size': 12})
    #plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
    #plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True
    plt.imshow(LaRa_mop)
    plt.colorbar(label=r"$\rho$ [-]",shrink=0.75,pad=0.02)
    plt.clim(-1,1)
    plt.xticks(range(0,len(RISE_mop[0])),labels=['0','','','','','5','','','','','10','','','','','15','','','','','20',
        '','','','','25','','','','','30','','','','','35','','','','','40',''])
    plt.yticks(range(0,len(RISE_mop[0])),labels=['0','','','','','5','','','','','10','','','','','15','','','','','20',
        '','','','','25','','','','','30','','','','','35','','','','','40',''])
    plt.ylabel("Index Parameters")
    plt.xlabel("Index Parameters")
    plt.savefig(output_folder_path+"/LaRa_small_correlation.pdf",bbox_inches="tight")
    plt.show()
    #plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = True
    #plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = False
    plt.close('all')

    plt.figure(figsize=(7.5,7.5))
    plt.rcParams.update({'font.size': 12})
    plt.imshow(np.abs(RISE_mop))
    plt.colorbar(label=r"$\rho$ [-]",shrink=0.75,pad=0.02)
    plt.xticks(range(0,len(RISE_mop[0])),labels=['','$y$','',r'$\dot{x}$',r'',r'$\dot{z}$','',
            r'$\sigma_{FCN}$',r'',r'$y_{{RISE}}$',r'',r'$x_{{LaRa}}$',r'',r'$z_{{LaRa}}$',
            r'',r'$\phi^s_1$',r'',r'$\phi^s_2$',r'',r'$\phi^s_3$',r'',r'$\phi^s_4$',
            r'',r'${X_p}^s_1$',r'',r'${Y_p}^s_1$',r'',r'${X_p}^s_2$',r'',r'${Y_p}^s_2$',
            r'',r'${X_p}^s_3$',r'',r'${Y_p}^s_3$',r'',r'${X_p}^s_4$',r'',r'${Y_p}^s_4$',
            r'',r'${X_p}^s_5$',r'',r'${Y_p}^s_5$'],rotation=90)
    plt.yticks(range(0,len(RISE_mop[0])),labels=['$x$','','$z$',r'',r'$\dot{y}$',r'','F',
            r'',r'$x_{{RISE}}$',r'',r'$z_{{RISE}}$',r'',r'$y_{{LaRa}}$',r'',
            r'$\phi^c_1$',r'',r'$\phi^c_2$',r'',r'$\phi^c_3$',r'',r'$\phi^c_4$',r'',
            r'${X_p}^c_1$',r'',r'${Y_p}^c_1$',r'',r'${X_p}^c_2$',r'',r'${Y_p}^c_2$',r'',
            r'${X_p}^c_3$',r'',r'${Y_p}^c_3$',r'',r'${X_p}^c_4$',r'',r'${Y_p}^c_4$',r'',
            r'${X_p}^c_5$',r'',r'${Y_p}^c_5$',r''])
    plt.savefig(output_folder_path+"/RISE_small_correlation_abs.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')


    plt.figure(figsize=(7.5,7.5))
    plt.rcParams.update({'font.size': 12})
    plt.imshow(RISE_mop)
    plt.colorbar(label=r"$\rho$ [-]",shrink=0.75,pad=0.02)
    plt.clim(-1,1)
    plt.xticks(range(0,len(RISE_mop[0])),labels=['0','','','','','5','','','','','10','','','','','15','','','','','20',
        '','','','','25','','','','','30','','','','','35','','','','','40',''])
    plt.yticks(range(0,len(RISE_mop[0])),labels=['0','','','','','5','','','','','10','','','','','15','','','','','20',
        '','','','','25','','','','','30','','','','','35','','','','','40',''])
    plt.ylabel("Index Parameters")
    plt.xlabel("Index Parameters")
    #plt.xticks(range(0,len(RISE_mop[0])),labels=['','$y$','',r'$\dot{x}$',r'',r'$\dot{z}$','',
    #        r'$\sigma_{FCN}$',r'',r'$y_{{RISE}}$',r'',r'$x_{{LaRa}}$',r'',r'$z_{{LaRa}}$',
    #        r'',r'$\phi^s_1$',r'',r'$\phi^s_2$',r'',r'$\phi^s_3$',r'',r'$\phi^s_4$',
    #        r'',r'${X_p}^s_1$',r'',r'${Y_p}^s_1$',r'',r'${X_p}^s_2$',r'',r'${Y_p}^s_2$',
    #        r'',r'${X_p}^s_3$',r'',r'${Y_p}^s_3$',r'',r'${X_p}^s_4$',r'',r'${Y_p}^s_4$',
    #        r'',r'${X_p}^s_5$',r'',r'${Y_p}^s_5$'],rotation=90)
    #plt.yticks(range(0,len(RISE_mop[0])),labels=['$x$','','$z$',r'',r'$\dot{y}$',r'','F',
    #        r'',r'$x_{{RISE}}$',r'',r'$z_{{RISE}}$',r'',r'$y_{{LaRa}}$',r'',
    #        r'$\phi^c_1$',r'',r'$\phi^c_2$',r'',r'$\phi^c_3$',r'',r'$\phi^c_4$',r'',
    #        r'${X_p}^c_1$',r'',r'${Y_p}^c_1$',r'',r'${X_p}^c_2$',r'',r'${Y_p}^c_2$',r'',
    #        r'${X_p}^c_3$',r'',r'${Y_p}^c_3$',r'',r'${X_p}^c_4$',r'',r'${Y_p}^c_4$',r'',
    #        r'${X_p}^c_5$',r'',r'${Y_p}^c_5$',r''])
    plt.savefig(output_folder_path+"/RISE_small_correlation.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(7.5,7.5))
    plt.rcParams.update({'font.size': 12})
    plt.imshow(np.abs(final_mop))
    plt.colorbar(label=r"$\rho$ [-]",shrink=0.75,pad=0.02)
    plt.xticks(range(0,len(final_mop[0])),labels=['','$y$','',r'$\dot{x}$',r'',r'$\dot{z}$','',
            r'$\sigma_{FCN}$',r'',r'$y_{{RISE}}$',r'',r'$x_{{LaRa}}$',r'',r'$z_{{LaRa}}$',
            r'',r'$\phi^s_1$',r'',r'$\phi^s_2$',r'',r'$\phi^s_3$',r'',r'$\phi^s_4$',
            r'',r'${X_p}^s_1$',r'',r'${Y_p}^s_1$',r'',r'${X_p}^s_2$',r'',r'${Y_p}^s_2$',
            r'',r'${X_p}^s_3$',r'',r'${Y_p}^s_3$',r'',r'${X_p}^s_4$',r'',r'${Y_p}^s_4$',
            r'',r'${X_p}^s_5$',r'',r'${Y_p}^s_5$'],rotation=90)
    plt.yticks(range(0,len(final_mop[0])),labels=['$x$','','$z$',r'',r'$\dot{y}$',r'','F',
            r'',r'$x_{{RISE}}$',r'',r'$z_{{RISE}}$',r'',r'$y_{{LaRa}}$',r'',
            r'$\phi^c_1$',r'',r'$\phi^c_2$',r'',r'$\phi^c_3$',r'',r'$\phi^c_4$',r'',
            r'${X_p}^c_1$',r'',r'${Y_p}^c_1$',r'',r'${X_p}^c_2$',r'',r'${Y_p}^c_2$',r'',
            r'${X_p}^c_3$',r'',r'${Y_p}^c_3$',r'',r'${X_p}^c_4$',r'',r'${Y_p}^c_4$',r'',
            r'${X_p}^c_5$',r'',r'${Y_p}^c_5$',r''])
    plt.savefig(output_folder_path+"/final_small_correlation_abs.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(7.5,7.5))
    plt.rcParams.update({'font.size': 12})
    plt.imshow(final_mop)
    plt.colorbar(label=r"$\rho$ [-]",shrink=0.75,pad=0.02)
    plt.clim(-1,1)
    plt.xticks(range(0,len(RISE_mop[0])),labels=['0','','','','','5','','','','','10','','','','','15','','','','','20',
        '','','','','25','','','','','30','','','','','35','','','','','40',''])
    plt.yticks(range(0,len(RISE_mop[0])),labels=['0','','','','','5','','','','','10','','','','','15','','','','','20',
        '','','','','25','','','','','30','','','','','35','','','','','40',''])
    plt.ylabel("Index Parameters")
    plt.xlabel("Index Parameters")
    plt.savefig(output_folder_path+"/final_small_correlation.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    LaRa_mop_noinit = np.loadtxt(output_folder_path+"/LaRa_mop_correlation_noinit.dat")
    RISE_mop_noinit = np.loadtxt(output_folder_path+"/RISE_mop_correlation_noinit.dat")
    final_mop_noinit = np.loadtxt(output_folder_path+"/final_mop_correlation_noinit.dat")

    plt.figure(figsize=(7.5,7.5))
    plt.rcParams.update({'font.size': 12})
    #plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
    #plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True
    plt.imshow(LaRa_mop_noinit)
    plt.colorbar(label=r"$\rho$ [-]",shrink=0.75,pad=0.02)
    plt.clim(-1,1)
    plt.xticks(range(0,len(LaRa_mop_noinit[0])),labels=['',
            r'$\sigma_{FCN}$',r'',r'$y_{{RISE}}$',r'',r'$x_{{LaRa}}$',r'',r'$z_{{LaRa}}$',
            r'',r'$\phi^s_1$',r'',r'$\phi^s_2$',r'',r'$\phi^s_3$',r'',r'$\phi^s_4$',
            r'',r'${X_p}^s_1$',r'',r'${Y_p}^s_1$',r'',r'${X_p}^s_2$',r'',r'${Y_p}^s_2$',
            r'',r'${X_p}^s_3$',r'',r'${Y_p}^s_3$',r'',r'${X_p}^s_4$',r'',r'${Y_p}^s_4$',
            r'',r'${X_p}^s_5$',r'',r'${Y_p}^s_5$'],rotation=90)
    plt.yticks(range(0,len(LaRa_mop_noinit[0])),labels=['F',
            r'',r'$x_{{RISE}}$',r'',r'$z_{{RISE}}$',r'',r'$y_{{LaRa}}$',r'',
            r'$\phi^c_1$',r'',r'$\phi^c_2$',r'',r'$\phi^c_3$',r'',r'$\phi^c_4$',r'',
            r'${X_p}^c_1$',r'',r'${Y_p}^c_1$',r'',r'${X_p}^c_2$',r'',r'${Y_p}^c_2$',r'',
            r'${X_p}^c_3$',r'',r'${Y_p}^c_3$',r'',r'${X_p}^c_4$',r'',r'${Y_p}^c_4$',r'',
            r'${X_p}^c_5$',r'',r'${Y_p}^c_5$',r''])
    plt.savefig(output_folder_path+"/LaRa_small_correlation_noinit.pdf",bbox_inches="tight")
    plt.show()
    #plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = True
    #plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = False
    plt.close('all')

    plt.figure(figsize=(7.5,7.5))
    plt.rcParams.update({'font.size': 12})
    #plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
    #plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True
    colormap = plt.cm.get_cmap('viridis')
    sm = plt.cm.ScalarMappable(cmap=colormap)
    sm.set_clim(vmin=0, vmax=1)
    plt.colorbar(sm,label=r"$\rho$ [-]",shrink=0.75,pad=0.02)
    plt.imshow(np.abs(LaRa_mop_noinit),cmap=colormap)
    plt.xticks(range(0,len(LaRa_mop_noinit[0])),labels=['',
            r'$\sigma_{FCN}$',r'',r'$y_{{RISE}}$',r'',r'$x_{{LaRa}}$',r'',r'$z_{{LaRa}}$',
            r'',r'$\phi^s_1$',r'',r'$\phi^s_2$',r'',r'$\phi^s_3$',r'',r'$\phi^s_4$',
            r'',r'${X_p}^s_1$',r'',r'${Y_p}^s_1$',r'',r'${X_p}^s_2$',r'',r'${Y_p}^s_2$',
            r'',r'${X_p}^s_3$',r'',r'${Y_p}^s_3$',r'',r'${X_p}^s_4$',r'',r'${Y_p}^s_4$',
            r'',r'${X_p}^s_5$',r'',r'${Y_p}^s_5$'],rotation=90)
    plt.yticks(range(0,len(LaRa_mop_noinit[0])),labels=['F',
            r'',r'$x_{{RISE}}$',r'',r'$z_{{RISE}}$',r'',r'$y_{{LaRa}}$',r'',
            r'$\phi^c_1$',r'',r'$\phi^c_2$',r'',r'$\phi^c_3$',r'',r'$\phi^c_4$',r'',
            r'${X_p}^c_1$',r'',r'${Y_p}^c_1$',r'',r'${X_p}^c_2$',r'',r'${Y_p}^c_2$',r'',
            r'${X_p}^c_3$',r'',r'${Y_p}^c_3$',r'',r'${X_p}^c_4$',r'',r'${Y_p}^c_4$',r'',
            r'${X_p}^c_5$',r'',r'${Y_p}^c_5$',r''])
    plt.savefig(output_folder_path+"/LaRa_small_correlation_noinit_abs.pdf",bbox_inches="tight")
    plt.show()
    #plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = True
    #plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = False
    plt.close('all')

    plt.figure(figsize=(7.5,7.5))
    plt.rcParams.update({'font.size': 12})
    #plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
    #plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True
    plt.imshow(RISE_mop_noinit)
    plt.colorbar(label=r"$\rho$ [-]",shrink=0.75,pad=0.02)
    plt.clim(-1,1)
    plt.xticks(range(0,len(RISE_mop_noinit[0])),labels=['',
            r'$\sigma_{FCN}$',r'',r'$y_{{RISE}}$',r'',r'$x_{{LaRa}}$',r'',r'$z_{{LaRa}}$',
            r'',r'$\phi^s_1$',r'',r'$\phi^s_2$',r'',r'$\phi^s_3$',r'',r'$\phi^s_4$',
            r'',r'${X_p}^s_1$',r'',r'${Y_p}^s_1$',r'',r'${X_p}^s_2$',r'',r'${Y_p}^s_2$',
            r'',r'${X_p}^s_3$',r'',r'${Y_p}^s_3$',r'',r'${X_p}^s_4$',r'',r'${Y_p}^s_4$',
            r'',r'${X_p}^s_5$',r'',r'${Y_p}^s_5$'],rotation=90)
    plt.yticks(range(0,len(RISE_mop_noinit[0])),labels=['F',
            r'',r'$x_{{RISE}}$',r'',r'$z_{{RISE}}$',r'',r'$y_{{LaRa}}$',r'',
            r'$\phi^c_1$',r'',r'$\phi^c_2$',r'',r'$\phi^c_3$',r'',r'$\phi^c_4$',r'',
            r'${X_p}^c_1$',r'',r'${Y_p}^c_1$',r'',r'${X_p}^c_2$',r'',r'${Y_p}^c_2$',r'',
            r'${X_p}^c_3$',r'',r'${Y_p}^c_3$',r'',r'${X_p}^c_4$',r'',r'${Y_p}^c_4$',r'',
            r'${X_p}^c_5$',r'',r'${Y_p}^c_5$',r''])
    plt.savefig(output_folder_path+"/RISE_small_correlation_noinit.pdf",bbox_inches="tight")
    plt.show()
    #plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = True
    #plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = False
    plt.close('all')

    plt.figure(figsize=(7.5,7.5))
    plt.rcParams.update({'font.size': 12})
    #plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
    #plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True
    colormap = plt.cm.get_cmap('viridis')
    sm = plt.cm.ScalarMappable(cmap=colormap)
    sm.set_clim(vmin=0, vmax=1)
    plt.colorbar(sm,label=r"$\rho$ [-]",shrink=0.75,pad=0.02)
    plt.imshow(np.abs(RISE_mop_noinit),cmap=colormap)
    plt.xticks(range(0,len(RISE_mop_noinit[0])),labels=['',
            r'$\sigma_{FCN}$',r'',r'$y_{{RISE}}$',r'',r'$x_{{LaRa}}$',r'',r'$z_{{LaRa}}$',
            r'',r'$\phi^s_1$',r'',r'$\phi^s_2$',r'',r'$\phi^s_3$',r'',r'$\phi^s_4$',
            r'',r'${X_p}^s_1$',r'',r'${Y_p}^s_1$',r'',r'${X_p}^s_2$',r'',r'${Y_p}^s_2$',
            r'',r'${X_p}^s_3$',r'',r'${Y_p}^s_3$',r'',r'${X_p}^s_4$',r'',r'${Y_p}^s_4$',
            r'',r'${X_p}^s_5$',r'',r'${Y_p}^s_5$'],rotation=90)
    plt.yticks(range(0,len(RISE_mop_noinit[0])),labels=['F',
            r'',r'$x_{{RISE}}$',r'',r'$z_{{RISE}}$',r'',r'$y_{{LaRa}}$',r'',
            r'$\phi^c_1$',r'',r'$\phi^c_2$',r'',r'$\phi^c_3$',r'',r'$\phi^c_4$',r'',
            r'${X_p}^c_1$',r'',r'${Y_p}^c_1$',r'',r'${X_p}^c_2$',r'',r'${Y_p}^c_2$',r'',
            r'${X_p}^c_3$',r'',r'${Y_p}^c_3$',r'',r'${X_p}^c_4$',r'',r'${Y_p}^c_4$',r'',
            r'${X_p}^c_5$',r'',r'${Y_p}^c_5$',r''])
    plt.savefig(output_folder_path+"/RISE_small_correlation_noinit_abs.pdf",bbox_inches="tight")
    plt.show()
    #plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = True
    #plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = False
    plt.close('all')

    plt.figure(figsize=(7.5,7.5))
    plt.rcParams.update({'font.size': 12})
    #plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
    #plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True
    colormap = plt.cm.get_cmap('viridis') #turbo
    sm = plt.cm.ScalarMappable(cmap=colormap)
    sm.set_clim(vmin=-1, vmax=1)
    plt.colorbar(sm,label=r"$\rho$ [-]",shrink=0.75,pad=0.02)
    plt.imshow(final_mop_noinit,cmap=colormap)
    plt.xticks(range(0,len(final_mop_noinit[0])),labels=['',
            r'$\sigma_{FCN}$',r'',r'$y_{{RISE}}$',r'',r'$x_{{LaRa}}$',r'',r'$z_{{LaRa}}$',
            r'',r'$\phi^s_1$',r'',r'$\phi^s_2$',r'',r'$\phi^s_3$',r'',r'$\phi^s_4$',
            r'',r'${X_p}^s_1$',r'',r'${Y_p}^s_1$',r'',r'${X_p}^s_2$',r'',r'${Y_p}^s_2$',
            r'',r'${X_p}^s_3$',r'',r'${Y_p}^s_3$',r'',r'${X_p}^s_4$',r'',r'${Y_p}^s_4$',
            r'',r'${X_p}^s_5$',r'',r'${Y_p}^s_5$'],rotation=90)
    plt.yticks(range(0,len(final_mop_noinit[0])),labels=['F',
            r'',r'$x_{{RISE}}$',r'',r'$z_{{RISE}}$',r'',r'$y_{{LaRa}}$',r'',
            r'$\phi^c_1$',r'',r'$\phi^c_2$',r'',r'$\phi^c_3$',r'',r'$\phi^c_4$',r'',
            r'${X_p}^c_1$',r'',r'${Y_p}^c_1$',r'',r'${X_p}^c_2$',r'',r'${Y_p}^c_2$',r'',
            r'${X_p}^c_3$',r'',r'${Y_p}^c_3$',r'',r'${X_p}^c_4$',r'',r'${Y_p}^c_4$',r'',
            r'${X_p}^c_5$',r'',r'${Y_p}^c_5$',r''])
    plt.savefig(output_folder_path+"/final_small_correlation_noinit.pdf",bbox_inches="tight")
    plt.show()
    #plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = True
    #plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = False
    plt.close('all')

    plt.figure(figsize=(7.5,7.5))
    plt.rcParams.update({'font.size': 12})
    #plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
    #plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True
    colormap = plt.cm.get_cmap('viridis')
    sm = plt.cm.ScalarMappable(cmap=colormap)
    sm.set_clim(vmin=0, vmax=1)
    plt.colorbar(sm,label=r"$\rho$ [-]",shrink=0.75,pad=0.02)
    plt.imshow(np.abs(final_mop_noinit),cmap=colormap)
    plt.xticks(range(0,len(final_mop_noinit[0])),labels=['',
            r'$\sigma_{FCN}$',r'',r'$y_{{RISE}}$',r'',r'$x_{{LaRa}}$',r'',r'$z_{{LaRa}}$',
            r'',r'$\phi^s_1$',r'',r'$\phi^s_2$',r'',r'$\phi^s_3$',r'',r'$\phi^s_4$',
            r'',r'${X_p}^s_1$',r'',r'${Y_p}^s_1$',r'',r'${X_p}^s_2$',r'',r'${Y_p}^s_2$',
            r'',r'${X_p}^s_3$',r'',r'${Y_p}^s_3$',r'',r'${X_p}^s_4$',r'',r'${Y_p}^s_4$',
            r'',r'${X_p}^s_5$',r'',r'${Y_p}^s_5$'],rotation=90)
    plt.yticks(range(0,len(final_mop_noinit[0])),labels=['F',
            r'',r'$x_{{RISE}}$',r'',r'$z_{{RISE}}$',r'',r'$y_{{LaRa}}$',r'',
            r'$\phi^c_1$',r'',r'$\phi^c_2$',r'',r'$\phi^c_3$',r'',r'$\phi^c_4$',r'',
            r'${X_p}^c_1$',r'',r'${Y_p}^c_1$',r'',r'${X_p}^c_2$',r'',r'${Y_p}^c_2$',r'',
            r'${X_p}^c_3$',r'',r'${Y_p}^c_3$',r'',r'${X_p}^c_4$',r'',r'${Y_p}^c_4$',r'',
            r'${X_p}^c_5$',r'',r'${Y_p}^c_5$',r''])
    plt.savefig(output_folder_path+"/final_small_correlation_noinit_abs.pdf",bbox_inches="tight")
    plt.show()
    #plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = True
    #plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = False
    plt.close('all')

    final_formalerrors_RISE = np.loadtxt(output_folder_path+"/final_formalerrors_plot_RISE.dat")
    final_formalerrors_corrvar_noinit = np.loadtxt(output_folder_path+"/final_formalerrors_plot_corrvar_noinit.dat")
    final_formalerrors_corrvar = np.loadtxt(output_folder_path+"/final_formalerrors_plot_corrvar.dat")
    final_formalerrors_099corr = np.loadtxt(output_folder_path+"/final_formalerrors_plot_095corr.dat")
    final_formalerrors_04677corr = np.loadtxt(output_folder_path+"/final_formalerrors_plot_04677corr.dat")
    final_formalerrors_mediancorr = np.loadtxt(output_folder_path+"/final_formalerrors_plot_mediancorr.dat")
    final_formalerrors_0corr = np.loadtxt(output_folder_path+"/final_formalerrors_plot_0corr.dat")
    final_formalerrors_noPRIDEcorr = np.loadtxt(output_folder_path+"/final_formalerrors_plot_noPRIDEcorr.dat")

    plt.figure(figsize=(3,3))
    #plt.scatter(np.arange(0,3),final_formalerrors_RISE[0:3],label="RISE")
    plt.scatter(np.arange(0,3),final_formalerrors_corrvar[0:3],label=r"RISE+LaRa $\rho$=variable")
    plt.scatter(np.arange(0,3),final_formalerrors_099corr[0:3],label=r"RISE+LaRa $\rho$=0.99")
    plt.scatter(np.arange(0,3),final_formalerrors_0corr[0:3],label=r"RISE+LaRa $\rho$=0")
    #plt.scatter(np.arange(0,3),final_formalerrors_04677corr[0:3],label=r"RISE+LaRa $\rho$=0.4677")
    plt.ylabel(r"1-$\sigma$ x,y,z [m]")
    plt.yscale('log')
    plt.xticks(range(0,3),labels=["x","y","z"])
    plt.grid()
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(output_folder_path+"/formal_errors_position_plot.pdf",bbox_inches="tight")
    plt.show()  
    plt.close('all')

    plt.figure(figsize=(3,3))
    #plt.scatter(np.arange(0,3),final_formalerrors_RISE[3:6],label="RISE")
    plt.scatter(np.arange(0,3),final_formalerrors_corrvar[3:6],label=r"RISE+LaRa $\rho$=variable")
    plt.scatter(np.arange(0,3),final_formalerrors_099corr[3:6],label=r"RISE+LaRa $\rho$=0.99")
    plt.scatter(np.arange(0,3),final_formalerrors_0corr[3:6],label=r"RISE+LaRa $\rho$=0")
    #plt.scatter(np.arange(0,3),final_formalerrors_04677corr[3:6],label=r"RISE+LaRa $\rho$=0.4677")
    plt.ylabel(r"1-$\sigma$ $\dot{x}$,$\dot{y}$,$\dot{z}$ [m/s]")
    plt.yscale('log')
    plt.xticks(range(0,3),labels=[r"$\dot{x}$",r"$\dot{y}$",r"$\dot{z}$"])
    plt.grid()
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(output_folder_path+"/formal_errors_velocity_plot.pdf",bbox_inches="tight")
    plt.show()  
    plt.close('all')

    plt.figure(figsize=(3,3))
    #plt.scatter(np.arange(0,1),final_formalerrors_RISE[6],label="RISE")
    plt.scatter(np.arange(0,1),final_formalerrors_corrvar[6],label=r"RISE+LaRa $\rho$=variable")
    plt.scatter(np.arange(0,1),final_formalerrors_099corr[6],label=r"RISE+LaRa $\rho$=0.99")
    plt.scatter(np.arange(0,1),final_formalerrors_0corr[6],label=r"RISE+LaRa $\rho$=0")
    #plt.scatter(np.arange(0,1),final_formalerrors_04677corr[6],label=r"RISE+LaRa $\rho$=0.4677")
    #plt.scatter(np.arange(0,1),final_formalerrors_corrvar_noinit[0],label=r"RISE+LaRa $\rho$=var no init")
    plt.ylabel(r"1-$\sigma$ $F$ [-]")
    plt.yscale('log')
    plt.xticks(range(0,1),labels=[r"F"])
    plt.grid()
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(output_folder_path+"/formal_errors_corefactor_plot.pdf",bbox_inches="tight")
    plt.show()  
    plt.close('all')

    i = 6
    print("F",100*np.abs(final_formalerrors_corrvar[i]-final_formalerrors_099corr[i])/final_formalerrors_099corr[i])

    plt.figure(figsize=(3,3))
    #plt.scatter(np.arange(0,1),final_formalerrors_RISE[7],label="RISE")
    plt.scatter(np.arange(0,1),final_formalerrors_corrvar[7],label=r"RISE+LaRa $\rho$=variable")
    plt.scatter(np.arange(0,1),final_formalerrors_099corr[7],label=r"RISE+LaRa $\rho$=0.99")
    plt.scatter(np.arange(0,1),final_formalerrors_0corr[7],label=r"RISE+LaRa $\rho$=0")
    #plt.scatter(np.arange(0,1),final_formalerrors_04677corr[7],label=r"RISE+LaRa $\rho$=0.4677")
    #plt.scatter(np.arange(0,1),final_formalerrors_corrvar_noinit[1],label=r"RISE+LaRa $\rho$=var no init")
    plt.ylabel(r"1-$\sigma$ $\sigma_{FCN}$ [rad/s]")
    plt.yscale('log')
    plt.xticks(range(0,1),labels=[r"$\sigma_{FCN}$"])
    plt.grid()
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(output_folder_path+"/formal_errors_fcnrate_plot.pdf",bbox_inches="tight")
    plt.show()  
    plt.close('all')  

    i = 7
    print("FCN",100*np.abs(final_formalerrors_corrvar[i]-final_formalerrors_099corr[i])/final_formalerrors_099corr[i])

    print("nutation",(100*np.abs(final_formalerrors_corrvar[i]-final_formalerrors_noPRIDEcorr[i])/final_formalerrors_noPRIDEcorr[i]+100*np.abs(final_formalerrors_corrvar[i]-final_formalerrors_noPRIDEcorr[i])/final_formalerrors_noPRIDEcorr[i])/2)

    plt.figure(figsize=(3,3))
    #plt.scatter(np.arange(0,3),final_formalerrors_RISE[8:11],label="RISE")
    plt.scatter(np.arange(0,6),final_formalerrors_corrvar[8:14],label=r"RISE+LaRa $\rho$=variable")
    plt.scatter(np.arange(0,6),final_formalerrors_099corr[8:14],label=r"RISE+LaRa $\rho$=0.99")
    plt.scatter(np.arange(0,6),final_formalerrors_0corr[8:14],label=r"RISE+LaRa $\rho$=0")
    #plt.scatter(np.arange(0,6),final_formalerrors_04677corr[8:14],label=r"RISE+LaRa $\rho$=0.4677")
    #plt.scatter(np.arange(0,6),final_formalerrors_corrvar_noinit[2:8],label=r"RISE+LaRa $\rho$=var no init")
    plt.ylabel(r"1-$\sigma$ x,y,z [rad/s]")
    plt.yscale('log')
    plt.xticks(range(0,6),labels=[r"$x_{RISE}$",r"$y_{RISE}$",r"$z_{RISE}$",r"$x_{LaRa}$",r"$y_{LaRa}$",r"$z_{LaRa}$"])
    plt.grid()
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(output_folder_path+"/formal_errors_landers_plot.pdf",bbox_inches="tight")
    plt.show()  
    plt.close('all')

    mas =np.pi/(180.0*1000.0*3600.0)

    plt.figure(figsize=(3,3))
    #plt.scatter(np.arange(0,8),final_formalerrors_RISE[11:19]/mas,label="RISE")
    plt.scatter(np.arange(0,8),final_formalerrors_corrvar[14:22]/mas,label=r"RISE+LaRa $\rho$=variable")
    plt.scatter(np.arange(0,8),final_formalerrors_099corr[14:22]/mas,label=r"RISE+LaRa $\rho$=0.99")
    plt.scatter(np.arange(0,8),final_formalerrors_0corr[14:22]/mas,label=r"RISE+LaRa $\rho$=0")
    #plt.scatter(np.arange(0,8),final_formalerrors_04677corr[14:22]/mas,label=r"RISE+LaRa $\rho$=0.4677")
    #plt.scatter(np.arange(0,8),final_formalerrors_corrvar_noinit[8:16]/mas,label=r"RISE+LaRa $\rho$=var no init")
    plt.ylabel(r"1-$\sigma$ $\phi$ [map]")
    plt.yscale('log')
    plt.xticks(range(0,8),labels=[r"$\phi^c_1$",r"$\phi^s_1$",r"$\phi^c_2$",r"$\phi^s_2$",r"$\phi^c_3$",r"$\phi^s_3$",r"$\phi^c_4$",r"$\phi^s_4$"])
    plt.grid()
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(output_folder_path+"/formal_errors_spinvariations_plot.pdf",bbox_inches="tight")
    plt.show()  
    plt.close('all')

    phi_sum = 0
    for i in range(14,22):
        print(100*np.abs(final_formalerrors_corrvar[i]-final_formalerrors_noPRIDEcorr[i])/final_formalerrors_noPRIDEcorr[i])
        phi_sum += 100*np.abs(final_formalerrors_corrvar[i]-final_formalerrors_noPRIDEcorr[i])/final_formalerrors_noPRIDEcorr[i]/len(range(14,22))
    print("phi",phi_sum)

    plt.figure(figsize=(9,3))
    #plt.scatter(np.arange(0,20),final_formalerrors_RISE[19:]/mas,label="RISE")
    plt.scatter(np.arange(0,20),final_formalerrors_corrvar[22:]/mas,label=r"RISE+LaRa $\rho$=variable")
    plt.scatter(np.arange(0,20),final_formalerrors_099corr[22:]/mas,label=r"RISE+LaRa $\rho$=0.99")
    plt.scatter(np.arange(0,20),final_formalerrors_0corr[22:]/mas,label=r"RISE+LaRa $\rho$=0")
    #plt.scatter(np.arange(0,20),final_formalerrors_04677corr[22:]/mas,label=r"RISE+LaRa $\rho$=0.4677")
    #plt.scatter(np.arange(0,20),final_formalerrors_corrvar_noinit[16:]/mas,label=r"RISE+LaRa $\rho$=var no init")
    plt.ylabel(r"1-$\sigma$ $X_p$,$Y_p$ [map]")
    plt.yscale('log')
    plt.xticks(range(0,20),labels=[r'$Xp^c_1$',r'$Xp^s_1$',r'$Yp^c_1$',r'$Yp^s_1$',r'$Xp^c_2$',r'$Xp^s_2$',r'$Yp^c_2$',r'$Yp^s_2$',
                r'$Xp^c_3$',r'$Xp^s_3$',r'$Yp^c_3$',r'$Yp^s_3$',r'$Xp^c_4$',r'$Xp^s_4$',r'$Yp^c_4$',r'$Yp^s_4$',
                r'$Xp^c_5$',r'$Xp^s_5$',r'$Yp^c_5$',r'$Yp^s_5$'])
    plt.grid()
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(output_folder_path+"/formal_errors_polarmotion_plot.pdf",bbox_inches="tight")
    plt.show()  
    plt.close('all')

    xyp_sum = 0
    for i in range(22,42):
        print(100*np.abs(final_formalerrors_corrvar[i]-final_formalerrors_099corr[i])/final_formalerrors_099corr[i])
        xyp_sum += 100*np.abs(final_formalerrors_corrvar[i]-final_formalerrors_099corr[i])/final_formalerrors_099corr[i]/len(range(22,42))
    print("xyp",xyp_sum)

    fig, (axs1,axs2,axs3) = plt.subplots(1,3,figsize=(15,4))
    axs1.scatter(np.arange(0,1),final_formalerrors_corrvar[6],label=r"RISE+LaRa $\rho$=variable")
    axs1.scatter(np.arange(0,1),final_formalerrors_noPRIDEcorr[6],label=r"RISE+LaRa $\rho$=1")
    axs1.scatter(np.arange(0,1),final_formalerrors_0corr[6],label=r"RISE+LaRa $\rho$=0")
    axs1.set(ylabel = r"1-$\sigma$ $F$ [-]")
    axs1.set_yscale('log')
    axs1.set_xticks(range(0,1),labels=[r"F"])
    axs1.grid()

    axs2.scatter(np.arange(0,1),final_formalerrors_corrvar[7],label=r"RISE+LaRa $\rho$=variable")
    axs2.scatter(np.arange(0,1),final_formalerrors_noPRIDEcorr[7],label=r"RISE+LaRa $\rho$=1",color='red')
    axs2.scatter(np.arange(0,1),final_formalerrors_0corr[7],label=r"RISE+LaRa $\rho$=0")
    axs2.set(ylabel = r"1-$\sigma$ $\sigma_{FCN}$ [rad/s]")
    axs2.set_yscale('log')
    axs2.set_xticks(range(0,1),labels=[r"$\sigma_{FCN}$"])
    axs2.grid()

    axs3.scatter(np.arange(0,8),final_formalerrors_corrvar[14:22]/mas,label=r"RISE+LaRa $\rho$=variable")
    axs3.scatter(np.arange(0,8),final_formalerrors_noPRIDEcorr[14:22]/mas,label=r"RISE+LaRa $\rho$=1")
    axs3.scatter(np.arange(0,8),final_formalerrors_0corr[14:22]/mas,label=r"RISE+LaRa $\rho$=0")
    axs3.set(ylabel = r"1-$\sigma$ $\phi$ [map]")
    axs3.set_yscale('log')
    axs3.set_xticks(range(0,8),labels=[r"$\phi^c_1$",r"$\phi^s_1$",r"$\phi^c_2$",r"$\phi^s_2$",r"$\phi^c_3$",r"$\phi^s_3$",r"$\phi^c_4$",r"$\phi^s_4$"])
    axs3.grid()
    axs3.legend()

    fig.tight_layout()

    plt.rcParams.update({'font.size': 9})

    fig, (axs1) = plt.subplots(1,figsize=(2,1.7))
    axs1.scatter(np.arange(0,1),final_formalerrors_corrvar[6],label=r"RISE+LaRa $\rho$=variable",color='blue')
    axs1.scatter(np.arange(0,1),final_formalerrors_noPRIDEcorr[6],label=r"RISE+LaRa $\rho$=1",color='red')
    axs1.scatter(np.arange(0,1),final_formalerrors_0corr[6],label=r"RISE+LaRa $\rho$=0",color='green')
    #axs1.scatter(np.arange(0,1),final_formalerrors_04677corr[6],label=r"RISE+LaRa $\rho$=11",color='black')
    axs1.set(ylabel = r"1-$\sigma$ $F$ [-]")
    axs1.set_xticks(range(0,1),labels=[r"F"])
    axs1.grid()
    fig.savefig(output_folder_path+"/formal_errors_corefactor_plot.pdf",bbox_inches="tight")

    fig, (axs1) = plt.subplots(1,figsize=(2,1.7))
    axs1.scatter(np.arange(0,1),final_formalerrors_corrvar[7],label=r"RISE+LaRa $\rho$=variable",color='blue')
    axs1.scatter(np.arange(0,1),final_formalerrors_noPRIDEcorr[7],label=r"RISE+LaRa $\rho$=1",color='red')
    axs1.scatter(np.arange(0,1),final_formalerrors_0corr[7],label=r"RISE+LaRa $\rho$=0",color='green')
    #axs1.scatter(np.arange(0,1),final_formalerrors_04677corr[7],label=r"RISE+LaRa $\rho$=11",color='black')
    axs1.set(ylabel = r"1-$\sigma$ $\sigma_{FCN}$ [rad/s]")
    axs1.set_xticks(range(0,1),labels=[r"$\sigma_{FCN}$"])
    axs1.grid()
    fig.savefig(output_folder_path+"/formal_errors_fcnrate_plot.pdf",bbox_inches="tight")

    fig, (axs1) = plt.subplots(1,figsize=(8,1.7))
    axs1.scatter(np.arange(0,8),final_formalerrors_noPRIDEcorr[14:22]/mas,label=r"RISE+LaRa $\rho$=1",color='red')
    axs1.scatter(np.arange(0,8),final_formalerrors_0corr[14:22]/mas,label=r"RISE+LaRa $\rho$=0",color='green')
    axs1.scatter(np.arange(0,8),final_formalerrors_corrvar[14:22]/mas,label=r"RISE+LaRa $\rho$=variable",color='blue')
    #axs1.scatter(np.arange(0,8),final_formalerrors_04677corr[14:22]/mas,label=r"RISE+LaRa $\rho$=11",color='black')
    axs1.set(ylabel = r"1-$\sigma$ $\phi$ [mas]")
    axs1.set_xticks(range(0,8),labels=[r"$\phi^c_1$",r"$\phi^s_1$",r"$\phi^c_2$",r"$\phi^s_2$",r"$\phi^c_3$",r"$\phi^s_3$",r"$\phi^c_4$",r"$\phi^s_4$"])
    axs1.grid()
    axs1.legend()
    fig.savefig(output_folder_path+"/formal_errors_spinvariations_plot.pdf",bbox_inches="tight")

    fig, (axs1) = plt.subplots(1,figsize=(12,1.7))
    axs1.scatter(np.arange(0,20),final_formalerrors_corrvar[22:]/mas,label=r"RISE+LaRa $\rho$=variable",color='blue')
    axs1.scatter(np.arange(0,20),final_formalerrors_noPRIDEcorr[22:]/mas,label=r"RISE+LaRa $\rho$=1",color='red')
    axs1.scatter(np.arange(0,20),final_formalerrors_0corr[22:]/mas,label=r"RISE+LaRa $\rho$=0",color='green')
    #axs1.scatter(np.arange(0,20),final_formalerrors_04677corr[22:]/mas,label=r"RISE+LaRa $\rho$=11",color='black')
    axs1.set(ylabel = r"1-$\sigma$ $X_p$,$Y_p$ [mas]")
    axs1.set_xticks(range(0,20),labels=[r'$Xp^c_1$',r'$Xp^s_1$',r'$Yp^c_1$',r'$Yp^s_1$',r'$Xp^c_2$',r'$Xp^s_2$',r'$Yp^c_2$',r'$Yp^s_2$',
                r'$Xp^c_3$',r'$Xp^s_3$',r'$Yp^c_3$',r'$Yp^s_3$',r'$Xp^c_C$',r'$Xp^s_C$',r'$Yp^c_C$',r'$Yp^s_C$',
                r'$Xp^c_4$',r'$Xp^s_4$',r'$Yp^c_4$',r'$Yp^s_4$'])
    axs1.grid()
    #axs1.legend()
    fig.savefig(output_folder_path+"/formal_errors_polarmotion_plot.pdf",bbox_inches="tight")

    RISE_sum = 0
    list_sum = list()
    list_sum.extend(range(6,8))
    list_sum.extend(range(14,42))
    print(list_sum)
    for i in list_sum:
        if i>8:
                print(100*np.abs(final_formalerrors_corrvar[i]-final_formalerrors_RISE[i-3])/final_formalerrors_RISE[i-3])
                RISE_sum += 100*np.abs(final_formalerrors_corrvar[i]-final_formalerrors_RISE[i-3])/final_formalerrors_RISE[i-3]/len(list_sum)
        else:
                print(100*np.abs(final_formalerrors_corrvar[i]-final_formalerrors_RISE[i])/final_formalerrors_RISE[i])
                RISE_sum += 100*np.abs(final_formalerrors_corrvar[i]-final_formalerrors_RISE[i])/final_formalerrors_RISE[i]/len(list_sum)

    print("RISE",RISE_sum)

print("--- %s seconds ---" % (time.time() - run_time))