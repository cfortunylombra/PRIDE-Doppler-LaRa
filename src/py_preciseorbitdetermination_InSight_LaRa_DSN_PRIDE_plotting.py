"""
Description: Environment Setup for the Precise Orbit Determination (RISE and LaRa with DSN-PRIDE)

Author: C. Fortuny-Lombraña
"""
import time
run_time = time.time()
if __name__=="__main__":

    import sys
    sys.path.insert(0, "/home/cfortunylombra/tudat-bundle/cmake-build-release-wsl/tudatpy/")
    import os
    import datetime
    import numpy as np
    import matplotlib.pyplot as plt
    from tudatpy.kernel import constants

    ########################################################################################################################
    ################################################## FILES ###############################################################
    ########################################################################################################################

    benchmark_folder = os.path.dirname(os.path.realpath(__file__)).replace('/src','/output/POD_RISEFalse_LaRaTrue_PRIDETrue_corr0')
    main_folder = os.path.dirname(os.path.realpath(__file__)).replace('/src','/output/POD_RISEFalse_LaRaTrue_PRIDETrue_corr0.95')

    # Booleans to understand whether we want to simulate together RISE and LaRa missions, or separetely
    RISE_boolean = False
    LaRa_boolean = True

    # Number of estimated parameters
    if RISE_boolean and LaRa_boolean:
        parameters_total_number = 42
    else:
        parameters_total_number = 39
    
    step_eval = 500

    concatenated_times_no_duplicated = np.loadtxt(benchmark_folder+'/concatenated_times_sort_no_duplicated.dat')

    arange_eval = np.arange(0,len(concatenated_times_no_duplicated),step_eval) 
    time_eval = [concatenated_times_no_duplicated[i] for i in arange_eval]

    # 1-sigma position as a function of time
    plt.figure(figsize=(15, 6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 3))))
    x_values = np.loadtxt(main_folder+'/xposition_plot.dat')
    y_values = np.loadtxt(main_folder+'/yposition_plot.dat') 
    z_values = np.loadtxt(main_folder+'/zposition_plot.dat')
    x_bench_values = np.loadtxt(benchmark_folder+'/xposition_plot.dat')
    y_bench_values = np.loadtxt(benchmark_folder+'/yposition_plot.dat')
    z_bench_values = np.loadtxt(benchmark_folder+'/zposition_plot.dat')

    plt.plot((time_eval-concatenated_times_no_duplicated[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        x_values,'-o',label='$x$')
    plt.plot((time_eval-concatenated_times_no_duplicated[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        y_values,'-o',label='$y$')
    plt.plot((time_eval-concatenated_times_no_duplicated[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        z_values,'-o',label='$z$')

    plt.plot((time_eval-concatenated_times_no_duplicated[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        x_bench_values,'--',label='$x_b$')
    plt.plot((time_eval-concatenated_times_no_duplicated[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        y_bench_values,'--',label='$y_b$')
    plt.plot((time_eval-concatenated_times_no_duplicated[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        z_bench_values,'--',label='$z_b$')
    #if RISE_boolean and LaRa_boolean:
    #    plt.axvline(x=(LaRa_observation_times_list[0]-observation_start_epoch)/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ x,y,z [m]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=concatenated_times_no_duplicated[0])))
    plt.grid()
    plt.legend()
    #plt.savefig(output_folder_path+"/positionvalues_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')
'''
    # 1-sigma position as a velocity of time
    plt.figure(figsize=(15, 6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(parameter_perturbation[3:6])))))
    xdot_values = list()
    ydot_values = list()
    zdot_values = list()
    for time_index in range(0,len(time_eval)):
        xdot_values.append(sigma_values[time_index][3])
        ydot_values.append(sigma_values[time_index][4])
        zdot_values.append(sigma_values[time_index][5])
    np.savetxt(output_folder_path+"/xdotvelocity_plot.dat",xdot_values,fmt='%.15e')
    np.savetxt(output_folder_path+"/ydotvelocity_plot.dat",ydot_values,fmt='%.15e')
    np.savetxt(output_folder_path+"/zdotvelocity_plot.dat",zdot_values,fmt='%.15e')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        xdot_values,'-o',label=r'$\dot{x}$')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        ydot_values,'-o',label=r'$\dot{y}$')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        zdot_values,'-o',label=r'$\dot{z}$')
    plt.ylabel(r'1-$\sigma$ $\dot{x}$,$\dot{y}$,$\dot{z}$ [m/s]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.grid()
    plt.legend()
    plt.savefig(output_folder_path+"/velocityvalues_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # 1-sigma F as a function of time
    plt.figure(figsize=(15, 6))
    F_values = list()
    for time_index in range(0,len(time_eval)):
        F_values.append(sigma_values[time_index][6])
    np.savetxt(output_folder_path+"/corefactor_plot.dat",F_values,fmt='%.15e')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        F_values,'-o')
    plt.ylabel(r'1-$\sigma$ F [-]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.grid()
    plt.savefig(output_folder_path+"/Fvalues_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # 1-sigma sigma_FCN as a function of time
    plt.figure(figsize=(15, 6))
    sigma_FCN_values = list()
    for time_index in range(0,len(time_eval)):
        sigma_FCN_values.append(sigma_values[time_index][7])
    np.savetxt(output_folder_path+"/sigmaFCN_plot.dat",sigma_FCN_values,fmt='%.15e')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        sigma_FCN_values,'-o')
    plt.ylabel(r'1-$\sigma$ $\sigma_{FCN}$ [rad/s]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.grid()
    plt.savefig(output_folder_path+"/sigmaFCNvalues_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # 1-sigma x_RISE,y_RISE,z_RISE,x_LaRa,y_LaRa,z_LaRa as a function of time
    plt.figure(figsize=(15, 6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(parameter_perturbation[8:11+add_par])))))
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
            xRISElander_values.append(sigma_values[time_index][8])
            yRISElander_values.append(sigma_values[time_index][9])
            zRISElander_values.append(sigma_values[time_index][10])
        
        if len(LaRa_observation_times_list)!=0:
            xLaRalander_values.append(sigma_values[time_index][8+add_par])
            yLaRalander_values.append(sigma_values[time_index][9+add_par])
            zLaRalander_values.append(sigma_values[time_index][10+add_par])

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
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(parameter_perturbation[11+add_par:19+add_par])))))
    cos1spin_values = list()
    sin1spin_values = list()
    cos2spin_values = list()
    sin2spin_values = list()
    cos3spin_values = list()
    sin3spin_values = list()
    cos4spin_values = list()
    sin4spin_values = list()
    for time_index in range(0,len(time_eval)):
        cos1spin_values.append(sigma_values[time_index][11+add_par])
        sin1spin_values.append(sigma_values[time_index][12+add_par])
        cos2spin_values.append(sigma_values[time_index][13+add_par])
        sin2spin_values.append(sigma_values[time_index][14+add_par])
        cos3spin_values.append(sigma_values[time_index][15+add_par])
        sin3spin_values.append(sigma_values[time_index][16+add_par])
        cos4spin_values.append(sigma_values[time_index][17+add_par])
        sin4spin_values.append(sigma_values[time_index][18+add_par])
    np.savetxt(output_folder_path+"/cos1spin_plot.dat",cos1spin_values,fmt='%.15e')
    np.savetxt(output_folder_path+"/sin1spin_plot.dat",sin1spin_values,fmt='%.15e')
    np.savetxt(output_folder_path+"/cos2spin_plot.dat",cos2spin_values,fmt='%.15e')
    np.savetxt(output_folder_path+"/sin2spin_plot.dat",sin2spin_values,fmt='%.15e')
    np.savetxt(output_folder_path+"/cos3spin_plot.dat",cos3spin_values,fmt='%.15e')
    np.savetxt(output_folder_path+"/sin3spin_plot.dat",sin3spin_values,fmt='%.15e')
    np.savetxt(output_folder_path+"/cos4spin_plot.dat",cos4spin_values,fmt='%.15e')
    np.savetxt(output_folder_path+"/sin4spin_plot.dat",sin4spin_values,fmt='%.15e')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(cos1spin_values)/mas,'-o',label=r'$\psi^c_1$')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(sin1spin_values)/mas,'-o',label=r'$\psi^s_1$')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(cos2spin_values)/mas,'-o',label=r'$\psi^c_2$')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(sin2spin_values)/mas,'-o',label=r'$\psi^s_2$')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(cos3spin_values)/mas,'-o',label=r'$\psi^c_3$')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(sin3spin_values)/mas,'-o',label=r'$\psi^s_3$')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(cos4spin_values)/mas,'-o',label=r'$\psi^c_4$')
    plt.plot((time_eval-observation_start_epoch*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(sin4spin_values)/mas,'-o',label=r'$\psi^s_4$')
    plt.ylabel(r'1-$\sigma$ $\psi$ [mas]')
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
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(parameter_perturbation[19+add_par:23+add_par])))))
    xpcos1_values = list()
    xpsin1_values = list()
    ypcos1_values = list()
    ypsin1_values = list()
    for time_index in range(0,len(time_eval)):
        xpcos1_values.append(sigma_values[time_index][19+add_par])
        xpsin1_values.append(sigma_values[time_index][20+add_par])
        ypcos1_values.append(sigma_values[time_index][21+add_par])
        ypsin1_values.append(sigma_values[time_index][22+add_par])
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
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(parameter_perturbation[23+add_par:27+add_par])))))
    xpcos2_values = list()
    xpsin2_values = list()
    ypcos2_values = list()
    ypsin2_values = list()
    for time_index in range(0,len(time_eval)):
        xpcos2_values.append(sigma_values[time_index][23+add_par])
        xpsin2_values.append(sigma_values[time_index][24+add_par])
        ypcos2_values.append(sigma_values[time_index][25+add_par])
        ypsin2_values.append(sigma_values[time_index][26+add_par])
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
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(parameter_perturbation[27+add_par:31+add_par])))))
    xpcos3_values = list()
    xpsin3_values = list()
    ypcos3_values = list()
    ypsin3_values = list()
    for time_index in range(0,len(time_eval)):
        xpcos3_values.append(sigma_values[time_index][27+add_par])
        xpsin3_values.append(sigma_values[time_index][28+add_par])
        ypcos3_values.append(sigma_values[time_index][29+add_par])
        ypsin3_values.append(sigma_values[time_index][30+add_par])
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
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(parameter_perturbation[31+add_par:35+add_par])))))
    xpcos4_values = list()
    xpsin4_values = list()
    ypcos4_values = list()
    ypsin4_values = list()
    for time_index in range(0,len(time_eval)):
        xpcos4_values.append(sigma_values[time_index][31+add_par])
        xpsin4_values.append(sigma_values[time_index][32+add_par])
        ypcos4_values.append(sigma_values[time_index][33+add_par])
        ypsin4_values.append(sigma_values[time_index][34+add_par]) 
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
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(parameter_perturbation[35+add_par:])))))
    xpcos5_values = list()
    xpsin5_values = list()
    ypcos5_values = list()
    ypsin5_values = list()
    for time_index in range(0,len(time_eval)):
        xpcos5_values.append(sigma_values[time_index][35+add_par])
        xpsin5_values.append(sigma_values[time_index][36+add_par])
        ypcos5_values.append(sigma_values[time_index][37+add_par])
        ypsin5_values.append(sigma_values[time_index][38+add_par])
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
    plt.ylabel(r'1-$\sigma$ $Xp, Yp$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=observation_start_epoch)))
    plt.legend()
    plt.grid()
    plt.savefig(output_folder_path+"/polarmotionamp5_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    print('Is there any duplicated total time value? :',any(list(concatenated_times).count(x) > 1 for x in list(concatenated_times)))

print("--- %s seconds ---" % (time.time() - run_time))
'''