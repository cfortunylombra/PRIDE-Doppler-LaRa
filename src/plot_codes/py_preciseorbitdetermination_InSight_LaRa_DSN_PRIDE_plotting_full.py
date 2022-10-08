"""
Description: Plot comparison between data

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

    benchmark_folder = os.path.dirname(os.path.realpath(__file__)).replace('/src','/output/POD_RISETrue_LaRaTrue_PRIDEFalseFalse_corr0')
    main0_folder = os.path.dirname(os.path.realpath(__file__)).replace('/src','/output/PODlat4.5_RISETrue_LaRaTrue_PRIDEFalseFalse_corr0')
    main1_folder = os.path.dirname(os.path.realpath(__file__)).replace('/src','/output/POD_RISETrue_LaRaTrue_PRIDETrueFalse_corr0')
    main2_folder = os.path.dirname(os.path.realpath(__file__)).replace('/src','/output/PODlat4.5_RISETrue_LaRaTrue_PRIDETrueFalse_corr0')
    output_folder_path = os.path.dirname(os.path.realpath(__file__)).replace('/src','/output/POD_comparison_plot_full')
    os.makedirs(output_folder_path,exist_ok=True)

    # Booleans to understand whether we want to simulate together RISE and LaRa missions, or separetely
    RISE_boolean = True
    LaRa_boolean = True

    mas =np.pi/(180.0*1000.0*3600.0) # Conversion from milli arc seconds to seconds 
    
    time_bench_eval = np.loadtxt(benchmark_folder+'/time_plot.dat')
    time0_eval = np.loadtxt(main0_folder+'/time_plot.dat')
    time1_eval = np.loadtxt(main1_folder+'/time_plot.dat')
    time2_eval = np.loadtxt(main2_folder+'/time_plot.dat')

    label_bench_eval = r'DSN & $\rho$=0'
    label0_eval = r'DSN & $\rho$=0 & lat-LaRa = 4.5$^\circ$'
    label1_eval = r'DSN & PRIDE & $\rho$=0'
    label2_eval = r'DSN & PRIDE & $\rho$=0 & lat-LaRa = 4.5$^\circ$'

    markersize_value = 4

    marker_bench = '--x'
    marker_0 = '-o'
    marker_1 = '--s'
    marker_2 = '-v'

    # 1-sigma position as a function of time
    plt.figure(figsize=(15, 6))
    x_bench_values = np.loadtxt(benchmark_folder+'/xposition_plot.dat')
    x0_values = np.loadtxt(main0_folder+'/xposition_plot.dat')
    x1_values = np.loadtxt(main1_folder+'/xposition_plot.dat')
    x2_values = np.loadtxt(main2_folder+'/xposition_plot.dat')

    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        x_bench_values,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        x0_values,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        x1_values,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        x2_values,marker_2,label=label2_eval,markersize=markersize_value)
    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ x [m]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.grid()
    plt.legend()
    plt.savefig(output_folder_path+"/xpositionvalues_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(15, 6))
    y_bench_values = np.loadtxt(benchmark_folder+'/yposition_plot.dat')
    y0_values = np.loadtxt(main0_folder+'/yposition_plot.dat')
    y1_values = np.loadtxt(main1_folder+'/yposition_plot.dat')
    y2_values = np.loadtxt(main2_folder+'/yposition_plot.dat')

    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        y_bench_values,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        y0_values,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        y1_values,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        y2_values,marker_2,label=label2_eval,markersize=markersize_value)
    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ y [m]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.grid()
    plt.legend()
    plt.savefig(output_folder_path+"/ypositionvalues_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(15, 6))
    z_bench_values = np.loadtxt(benchmark_folder+'/zposition_plot.dat')
    z0_values = np.loadtxt(main0_folder+'/zposition_plot.dat')
    z1_values = np.loadtxt(main1_folder+'/zposition_plot.dat')
    z2_values = np.loadtxt(main2_folder+'/zposition_plot.dat')

    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        z_bench_values,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        z0_values,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        z1_values,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        z2_values,marker_2,label=label2_eval,markersize=markersize_value)
    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ z [m]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.grid()
    plt.legend()
    plt.savefig(output_folder_path+"/ypositionvalues_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # 1-sigma position as a velocity of time
    plt.figure(figsize=(15, 6))
    xdot_bench_values = np.loadtxt(benchmark_folder+'/xdotvelocity_plot.dat')
    xdot0_values = np.loadtxt(main0_folder+'/xdotvelocity_plot.dat')
    xdot1_values = np.loadtxt(main1_folder+'/xdotvelocity_plot.dat')
    xdot2_values = np.loadtxt(main2_folder+'/xdotvelocity_plot.dat')

    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        xdot_bench_values,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        xdot0_values,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        xdot1_values,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        xdot2_values,marker_2,label=label2_eval,markersize=markersize_value)
    
    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $\dot{x}$ [m/s]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.grid()
    plt.legend()
    plt.savefig(output_folder_path+"/xvelocityvalues_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(15, 6))
    ydot_bench_values = np.loadtxt(benchmark_folder+'/ydotvelocity_plot.dat')
    ydot0_values = np.loadtxt(main0_folder+'/ydotvelocity_plot.dat')
    ydot1_values = np.loadtxt(main1_folder+'/ydotvelocity_plot.dat')
    ydot2_values = np.loadtxt(main2_folder+'/ydotvelocity_plot.dat')

    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        ydot_bench_values,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        ydot0_values,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        ydot1_values,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        ydot2_values,marker_2,label=label2_eval,markersize=markersize_value)
    
    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $\dot{y}$ [m/s]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.grid()
    plt.legend()
    plt.savefig(output_folder_path+"/yvelocityvalues_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(15, 6))
    zdot_bench_values = np.loadtxt(benchmark_folder+'/zdotvelocity_plot.dat')
    zdot0_values = np.loadtxt(main0_folder+'/zdotvelocity_plot.dat')
    zdot1_values = np.loadtxt(main1_folder+'/zdotvelocity_plot.dat')
    zdot2_values = np.loadtxt(main2_folder+'/zdotvelocity_plot.dat')

    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        zdot_bench_values,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        zdot0_values,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        zdot1_values,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        zdot2_values,marker_2,label=label2_eval,markersize=markersize_value)
    
    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $\dot{z}$ [m/s]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.grid()
    plt.legend()
    plt.savefig(output_folder_path+"/zvelocityvalues_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # 1-sigma F as a function of time
    plt.figure(figsize=(15, 6))
    F_bench_values = np.loadtxt(benchmark_folder+'/corefactor_plot.dat')
    F0_values = np.loadtxt(main0_folder+'/corefactor_plot.dat')
    F1_values = np.loadtxt(main1_folder+'/corefactor_plot.dat')
    F2_values = np.loadtxt(main2_folder+'/corefactor_plot.dat')
    
    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        F_bench_values,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        F0_values,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        F1_values,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        F2_values,marker_2,label=label2_eval,markersize=markersize_value)
    
    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ F [-]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.grid()
    plt.legend()
    plt.savefig(output_folder_path+"/Fvalues_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # 1-sigma sigma_FCN as a function of time
    plt.figure(figsize=(15, 6))
    sigma_FCN_bench_values = np.loadtxt(benchmark_folder+'/sigmaFCN_plot.dat')
    sigma_FCN0_values = np.loadtxt(main0_folder+'/sigmaFCN_plot.dat')
    sigma_FCN1_values = np.loadtxt(main1_folder+'/sigmaFCN_plot.dat')
    sigma_FCN2_values = np.loadtxt(main2_folder+'/sigmaFCN_plot.dat')
    
    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        sigma_FCN_bench_values,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        sigma_FCN0_values,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        sigma_FCN1_values,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        sigma_FCN2_values,marker_2,label=label2_eval,markersize=markersize_value)

    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $\sigma_{FCN}$ [rad/s]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.grid()
    plt.legend()
    plt.savefig(output_folder_path+"/sigmaFCNvalues_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # 1-sigma x_RISE,y_RISE,z_RISE,x_LaRa,y_LaRa,z_LaRa as a function of time
    if RISE_boolean:
        plt.figure(figsize=(15, 6))
        xRISElander_bench_values = np.loadtxt(benchmark_folder+'/xRISE_plot.dat')
        xRISElander0_values = np.loadtxt(main0_folder+'/xRISE_plot.dat')
        xRISElander1_values = np.loadtxt(main1_folder+'/xRISE_plot.dat')
        xRISElander2_values = np.loadtxt(main2_folder+'/xRISE_plot.dat')

        plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
            xRISElander_bench_values,marker_bench,label=label_bench_eval,markersize=markersize_value)
        plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
            xRISElander0_values,marker_0,label=label0_eval,markersize=markersize_value)
        plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
            xRISElander1_values,marker_1,label=label1_eval,markersize=markersize_value)
        plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
            xRISElander2_values,marker_2,label=label2_eval,markersize=markersize_value)

        plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
        plt.ylabel(r'1-$\sigma$ $x_{RISE}$ [m]')
        plt.xlabel('Time [days]')
        plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
        plt.legend()
        plt.grid()
        plt.savefig(output_folder_path+"/xRISElander_time.pdf",bbox_inches="tight")
        plt.show()
        plt.close('all')
    
        plt.figure(figsize=(15, 6))
        yRISElander_bench_values = np.loadtxt(benchmark_folder+'/yRISE_plot.dat')
        yRISElander0_values = np.loadtxt(main0_folder+'/yRISE_plot.dat')
        yRISElander1_values = np.loadtxt(main1_folder+'/yRISE_plot.dat')
        yRISElander2_values = np.loadtxt(main2_folder+'/yRISE_plot.dat')

        plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
            yRISElander_bench_values,marker_bench,label=label_bench_eval,markersize=markersize_value)
        plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
            yRISElander0_values,marker_0,label=label0_eval,markersize=markersize_value)
        plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
            yRISElander1_values,marker_1,label=label1_eval,markersize=markersize_value)
        plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
            yRISElander2_values,marker_2,label=label2_eval,markersize=markersize_value)

        plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
        plt.ylabel(r'1-$\sigma$ $y_{RISE}$ [m]')
        plt.xlabel('Time [days]')
        plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
        plt.legend()
        plt.grid()
        plt.savefig(output_folder_path+"/yRISElander_time.pdf",bbox_inches="tight")
        plt.show()
        plt.close('all')

        plt.figure(figsize=(15, 6))
        zRISElander_bench_values = np.loadtxt(benchmark_folder+'/zRISE_plot.dat')
        zRISElander0_values = np.loadtxt(main0_folder+'/zRISE_plot.dat')
        zRISElander1_values = np.loadtxt(main1_folder+'/zRISE_plot.dat')
        zRISElander2_values = np.loadtxt(main2_folder+'/zRISE_plot.dat')

        plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
            zRISElander_bench_values,marker_bench,label=label_bench_eval,markersize=markersize_value)
        plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
            zRISElander0_values,marker_0,label=label0_eval,markersize=markersize_value)
        plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
            zRISElander1_values,marker_1,label=label1_eval,markersize=markersize_value)
        plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
            zRISElander2_values,marker_2,label=label2_eval,markersize=markersize_value)

        plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
        plt.ylabel(r'1-$\sigma$ $z_{RISE}$ [m]')
        plt.xlabel('Time [days]')
        plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
        plt.legend()
        plt.grid()
        plt.savefig(output_folder_path+"/zRISElander_time.pdf",bbox_inches="tight")
        plt.show()
        plt.close('all')
    
    if LaRa_boolean:
        plt.figure(figsize=(15, 6))
        xLaRalander_bench_values = np.loadtxt(benchmark_folder+'/xLaRa_plot.dat')
        xLaRalander0_values = np.loadtxt(main0_folder+'/xLaRa_plot.dat')
        xLaRalander1_values = np.loadtxt(main1_folder+'/xLaRa_plot.dat')
        xLaRalander2_values = np.loadtxt(main2_folder+'/xLaRa_plot.dat')

        plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
            xLaRalander_bench_values,marker_bench,label=label_bench_eval,markersize=markersize_value)
        plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
            xLaRalander0_values,marker_0,label=label0_eval,markersize=markersize_value)
        plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
            xLaRalander1_values,marker_1,label=label1_eval,markersize=markersize_value)
        plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
            xLaRalander2_values,marker_2,label=label2_eval,markersize=markersize_value)

        plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
        plt.ylabel(r'1-$\sigma$ ${x_{LaRa}}$ [m]')
        plt.xlabel('Time [days]')
        plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
        plt.legend()
        plt.grid()
        plt.savefig(output_folder_path+"/xLaRalander_time.pdf",bbox_inches="tight")
        plt.show()
        plt.close('all')
    
        plt.figure(figsize=(15, 6))
        yLaRalander_bench_values = np.loadtxt(benchmark_folder+'/yLaRa_plot.dat')
        yLaRalander0_values = np.loadtxt(main0_folder+'/yLaRa_plot.dat')
        yLaRalander1_values = np.loadtxt(main1_folder+'/yLaRa_plot.dat')
        yLaRalander2_values = np.loadtxt(main2_folder+'/yLaRa_plot.dat')

        plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
            yLaRalander_bench_values,marker_bench,label=label_bench_eval,markersize=markersize_value)
        plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
            yLaRalander0_values,marker_0,label=label0_eval,markersize=markersize_value)
        plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
            yLaRalander1_values,marker_1,label=label1_eval,markersize=markersize_value)
        plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
            yLaRalander2_values,marker_2,label=label2_eval,markersize=markersize_value)

        plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
        plt.ylabel(r'1-$\sigma$ ${y_{LaRa}}$ [m]')
        plt.xlabel('Time [days]')
        plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
        plt.legend()
        plt.grid()
        plt.savefig(output_folder_path+"/yLaRalander_time.pdf",bbox_inches="tight")
        plt.show()
        plt.close('all')

        plt.figure(figsize=(15, 6))
        zLaRalander_bench_values = np.loadtxt(benchmark_folder+'/zLaRa_plot.dat')
        zLaRalander0_values = np.loadtxt(main0_folder+'/zLaRa_plot.dat')
        zLaRalander1_values = np.loadtxt(main1_folder+'/zLaRa_plot.dat')
        zLaRalander2_values = np.loadtxt(main2_folder+'/zLaRa_plot.dat')

        plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
            zLaRalander_bench_values,marker_bench,label=label_bench_eval,markersize=markersize_value)
        plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
            zLaRalander0_values,marker_0,label=label0_eval,markersize=markersize_value)
        plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
            zLaRalander1_values,marker_1,label=label1_eval,markersize=markersize_value)
        plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
            zLaRalander2_values,marker_2,label=label2_eval,markersize=markersize_value)

        plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
        plt.ylabel(r'1-$\sigma$ ${z_{LaRa}}$ [m]')
        plt.xlabel('Time [days]')
        plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
        plt.legend()
        plt.grid()
        plt.savefig(output_folder_path+"/zLaRalander_time.pdf",bbox_inches="tight")
        plt.show()
        plt.close('all')

    # 1-sigma spin variations as a function of time
    plt.figure(figsize=(15, 6))
    cos1spin_bench_values = np.loadtxt(benchmark_folder+'/cos1spin_plot.dat')
    cos1spin0_values = np.loadtxt(main0_folder+'/cos1spin_plot.dat')
    cos1spin1_values = np.loadtxt(main1_folder+'/cos1spin_plot.dat')
    cos1spin2_values = np.loadtxt(main2_folder+'/cos1spin_plot.dat')
    
    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(cos1spin_bench_values)/mas,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        np.array(cos1spin0_values)/mas,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        np.array(cos1spin1_values)/mas,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        np.array(cos1spin2_values)/mas,marker_2,label=label2_eval,markersize=markersize_value)

    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $\phi^c_1$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.legend(ncol=1)
    plt.grid()
    plt.savefig(output_folder_path+"/psispincos1_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(15, 6))
    sin1spin_bench_values = np.loadtxt(benchmark_folder+'/sin1spin_plot.dat')
    sin1spin0_values = np.loadtxt(main0_folder+'/sin1spin_plot.dat')
    sin1spin1_values = np.loadtxt(main1_folder+'/sin1spin_plot.dat')
    sin1spin2_values = np.loadtxt(main2_folder+'/sin1spin_plot.dat')
    
    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(sin1spin_bench_values)/mas,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        np.array(sin1spin0_values)/mas,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        np.array(sin1spin1_values)/mas,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        np.array(sin1spin2_values)/mas,marker_2,label=label2_eval,markersize=markersize_value)

    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $\phi^s_1$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.legend(ncol=1)
    plt.grid()
    plt.savefig(output_folder_path+"/psispinsin1_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(15, 6))
    cos2spin_bench_values = np.loadtxt(benchmark_folder+'/cos2spin_plot.dat')
    cos2spin0_values = np.loadtxt(main0_folder+'/cos2spin_plot.dat')
    cos2spin1_values = np.loadtxt(main1_folder+'/cos2spin_plot.dat')
    cos2spin2_values = np.loadtxt(main2_folder+'/cos2spin_plot.dat')
    
    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(cos2spin_bench_values)/mas,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        np.array(cos2spin0_values)/mas,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        np.array(cos2spin1_values)/mas,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        np.array(cos2spin2_values)/mas,marker_2,label=label2_eval,markersize=markersize_value)

    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $\phi^c_2$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.legend(ncol=1)
    plt.grid()
    plt.savefig(output_folder_path+"/psispincos2_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(15, 6))
    sin2spin_bench_values = np.loadtxt(benchmark_folder+'/sin2spin_plot.dat')
    sin2spin0_values = np.loadtxt(main0_folder+'/sin2spin_plot.dat')
    sin2spin1_values = np.loadtxt(main1_folder+'/sin2spin_plot.dat')
    sin2spin2_values = np.loadtxt(main2_folder+'/sin2spin_plot.dat')
    
    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(sin2spin_bench_values)/mas,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        np.array(sin2spin0_values)/mas,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        np.array(sin2spin1_values)/mas,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        np.array(sin2spin2_values)/mas,marker_2,label=label2_eval,markersize=markersize_value)

    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $\phi^s_2$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.legend(ncol=1)
    plt.grid()
    plt.savefig(output_folder_path+"/psispinsin2_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(15, 6))
    cos3spin_bench_values = np.loadtxt(benchmark_folder+'/cos3spin_plot.dat')
    cos3spin0_values = np.loadtxt(main0_folder+'/cos3spin_plot.dat')
    cos3spin1_values = np.loadtxt(main1_folder+'/cos3spin_plot.dat')
    cos3spin2_values = np.loadtxt(main2_folder+'/cos3spin_plot.dat')
    
    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(cos3spin_bench_values)/mas,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        np.array(cos3spin0_values)/mas,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        np.array(cos3spin1_values)/mas,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        np.array(cos3spin2_values)/mas,marker_2,label=label2_eval,markersize=markersize_value)

    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $\phi^c_3$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.legend(ncol=1)
    plt.grid()
    plt.savefig(output_folder_path+"/psispincos3_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(15, 6))
    sin3spin_bench_values = np.loadtxt(benchmark_folder+'/sin3spin_plot.dat')
    sin3spin0_values = np.loadtxt(main0_folder+'/sin3spin_plot.dat')
    sin3spin1_values = np.loadtxt(main1_folder+'/sin3spin_plot.dat')
    sin3spin2_values = np.loadtxt(main2_folder+'/sin3spin_plot.dat')
    
    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(sin3spin_bench_values)/mas,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        np.array(sin3spin0_values)/mas,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        np.array(sin3spin1_values)/mas,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        np.array(sin3spin2_values)/mas,marker_2,label=label2_eval,markersize=markersize_value)

    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $\phi^s_3$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.legend(ncol=1)
    plt.grid()
    plt.savefig(output_folder_path+"/psispinsin3_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(15, 6))
    cos4spin_bench_values = np.loadtxt(benchmark_folder+'/cos4spin_plot.dat')
    cos4spin0_values = np.loadtxt(main0_folder+'/cos4spin_plot.dat')
    cos4spin1_values = np.loadtxt(main1_folder+'/cos4spin_plot.dat')
    cos4spin2_values = np.loadtxt(main2_folder+'/cos4spin_plot.dat')
    
    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(cos4spin_bench_values)/mas,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        np.array(cos4spin0_values)/mas,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        np.array(cos4spin1_values)/mas,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        np.array(cos4spin2_values)/mas,marker_2,label=label2_eval,markersize=markersize_value)

    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $\phi^c_4$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.legend(ncol=1)
    plt.grid()
    plt.savefig(output_folder_path+"/psispincos4_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(15, 6))
    sin4spin_bench_values = np.loadtxt(benchmark_folder+'/sin4spin_plot.dat')
    sin4spin0_values = np.loadtxt(main0_folder+'/sin4spin_plot.dat')
    sin4spin1_values = np.loadtxt(main1_folder+'/sin4spin_plot.dat')
    sin4spin2_values = np.loadtxt(main2_folder+'/sin4spin_plot.dat')
    
    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(sin4spin_bench_values)/mas,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        np.array(sin4spin0_values)/mas,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        np.array(sin4spin1_values)/mas,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        np.array(sin4spin2_values)/mas,marker_2,label=label2_eval,markersize=markersize_value)

    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $\phi^s_4$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.legend(ncol=1)
    plt.grid()
    plt.savefig(output_folder_path+"/psispinsin4_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # 1-sigma polar motion (order 1) as a function of time
    plt.figure(figsize=(15, 6))
    xpcos1_bench_values = np.loadtxt(benchmark_folder+'/xpcos1_plot.dat')
    xpcos10_values = np.loadtxt(main0_folder+'/xpcos1_plot.dat')
    xpcos11_values = np.loadtxt(main1_folder+'/xpcos1_plot.dat')
    xpcos12_values = np.loadtxt(main2_folder+'/xpcos1_plot.dat')
    
    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(xpcos1_bench_values)/mas,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        np.array(xpcos10_values)/mas,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        np.array(xpcos11_values)/mas,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        np.array(xpcos12_values)/mas,marker_2,label=label2_eval,markersize=markersize_value)

    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $Xp^c_1$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.legend(ncol=1)
    plt.grid()
    plt.savefig(output_folder_path+"/xpcos_polarmotionamp1_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(15, 6))
    xpsin1_bench_values = np.loadtxt(benchmark_folder+'/xpsin1_plot.dat')
    xpsin10_values = np.loadtxt(main0_folder+'/xpsin1_plot.dat')
    xpsin11_values = np.loadtxt(main1_folder+'/xpsin1_plot.dat')
    xpsin12_values = np.loadtxt(main2_folder+'/xpsin1_plot.dat')
    
    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(xpsin1_bench_values)/mas,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        np.array(xpsin10_values)/mas,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        np.array(xpsin11_values)/mas,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        np.array(xpsin12_values)/mas,marker_2,label=label2_eval,markersize=markersize_value)

    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $Xp^s_1$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.legend(ncol=1)
    plt.grid()
    plt.savefig(output_folder_path+"/xpsin_polarmotionamp1_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(15, 6))
    ypcos1_bench_values = np.loadtxt(benchmark_folder+'/ypcos1_plot.dat')
    ypcos10_values = np.loadtxt(main0_folder+'/ypcos1_plot.dat')
    ypcos11_values = np.loadtxt(main1_folder+'/ypcos1_plot.dat')
    ypcos12_values = np.loadtxt(main2_folder+'/ypcos1_plot.dat')
    
    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(ypcos1_bench_values)/mas,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        np.array(ypcos10_values)/mas,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        np.array(ypcos11_values)/mas,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        np.array(ypcos12_values)/mas,marker_2,label=label2_eval,markersize=markersize_value)

    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $Yp^c_1$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.legend(ncol=1)
    plt.grid()
    plt.savefig(output_folder_path+"/ypcos_polarmotionamp1_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(15, 6))
    ypsin1_bench_values = np.loadtxt(benchmark_folder+'/ypsin1_plot.dat')
    ypsin10_values = np.loadtxt(main0_folder+'/ypsin1_plot.dat')
    ypsin11_values = np.loadtxt(main1_folder+'/ypsin1_plot.dat')
    ypsin12_values = np.loadtxt(main2_folder+'/ypsin1_plot.dat')
    
    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(ypsin1_bench_values)/mas,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        np.array(ypsin10_values)/mas,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        np.array(ypsin11_values)/mas,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        np.array(ypsin12_values)/mas,marker_2,label=label2_eval,markersize=markersize_value)

    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $Yp^s_1$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.legend(ncol=1)
    plt.grid()
    plt.savefig(output_folder_path+"/ypsin_polarmotionamp1_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # 1-sigma polar motion (order 2) as a function of time
    plt.figure(figsize=(15, 6))
    xpcos2_bench_values = np.loadtxt(benchmark_folder+'/xpcos2_plot.dat')
    xpcos20_values = np.loadtxt(main0_folder+'/xpcos2_plot.dat')
    xpcos21_values = np.loadtxt(main1_folder+'/xpcos2_plot.dat')
    xpcos22_values = np.loadtxt(main2_folder+'/xpcos2_plot.dat')
    
    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(xpcos2_bench_values)/mas,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        np.array(xpcos20_values)/mas,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        np.array(xpcos21_values)/mas,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        np.array(xpcos22_values)/mas,marker_2,label=label2_eval,markersize=markersize_value)

    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $Xp^c_2$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.legend(ncol=1)
    plt.grid()
    plt.savefig(output_folder_path+"/xpcos_polarmotionamp2_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(15, 6))
    xpsin2_bench_values = np.loadtxt(benchmark_folder+'/xpsin2_plot.dat')
    xpsin20_values = np.loadtxt(main0_folder+'/xpsin2_plot.dat')
    xpsin21_values = np.loadtxt(main1_folder+'/xpsin2_plot.dat')
    xpsin22_values = np.loadtxt(main2_folder+'/xpsin2_plot.dat')
    
    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(xpsin2_bench_values)/mas,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        np.array(xpsin20_values)/mas,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        np.array(xpsin21_values)/mas,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        np.array(xpsin22_values)/mas,marker_2,label=label2_eval,markersize=markersize_value)

    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $Xp^s_2$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.legend(ncol=1)
    plt.grid()
    plt.savefig(output_folder_path+"/xpsin_polarmotionamp2_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(15, 6))
    ypcos2_bench_values = np.loadtxt(benchmark_folder+'/ypcos2_plot.dat')
    ypcos20_values = np.loadtxt(main0_folder+'/ypcos2_plot.dat')
    ypcos21_values = np.loadtxt(main1_folder+'/ypcos2_plot.dat')
    ypcos22_values = np.loadtxt(main2_folder+'/ypcos2_plot.dat')
    
    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(ypcos2_bench_values)/mas,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        np.array(ypcos20_values)/mas,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        np.array(ypcos21_values)/mas,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        np.array(ypcos22_values)/mas,marker_2,label=label2_eval,markersize=markersize_value)

    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $Yp^c_2$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.legend(ncol=1)
    plt.grid()
    plt.savefig(output_folder_path+"/ypcos_polarmotionamp2_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(15, 6))
    ypsin2_bench_values = np.loadtxt(benchmark_folder+'/ypsin2_plot.dat')
    ypsin20_values = np.loadtxt(main0_folder+'/ypsin2_plot.dat')
    ypsin21_values = np.loadtxt(main1_folder+'/ypsin2_plot.dat')
    ypsin22_values = np.loadtxt(main2_folder+'/ypsin2_plot.dat')
    
    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(ypsin2_bench_values)/mas,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        np.array(ypsin20_values)/mas,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        np.array(ypsin21_values)/mas,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        np.array(ypsin22_values)/mas,marker_2,label=label2_eval,markersize=markersize_value)

    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $Yp^s_2$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.legend(ncol=1)
    plt.grid()
    plt.savefig(output_folder_path+"/ypsin_polarmotionamp2_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # 1-sigma polar motion (order 3) as a function of time
    plt.figure(figsize=(15, 6))
    xpcos3_bench_values = np.loadtxt(benchmark_folder+'/xpcos3_plot.dat')
    xpcos30_values = np.loadtxt(main0_folder+'/xpcos3_plot.dat')
    xpcos31_values = np.loadtxt(main1_folder+'/xpcos3_plot.dat')
    xpcos32_values = np.loadtxt(main2_folder+'/xpcos3_plot.dat')
    
    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(xpcos3_bench_values)/mas,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        np.array(xpcos30_values)/mas,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        np.array(xpcos31_values)/mas,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        np.array(xpcos32_values)/mas,marker_2,label=label2_eval,markersize=markersize_value)

    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $Xp^c_3$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.legend(ncol=1)
    plt.grid()
    plt.savefig(output_folder_path+"/xpcos_polarmotionamp3_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(15, 6))
    xpsin3_bench_values = np.loadtxt(benchmark_folder+'/xpsin3_plot.dat')
    xpsin30_values = np.loadtxt(main0_folder+'/xpsin3_plot.dat')
    xpsin31_values = np.loadtxt(main1_folder+'/xpsin3_plot.dat')
    xpsin32_values = np.loadtxt(main2_folder+'/xpsin3_plot.dat')
    
    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(xpsin3_bench_values)/mas,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        np.array(xpsin30_values)/mas,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        np.array(xpsin31_values)/mas,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        np.array(xpsin32_values)/mas,marker_2,label=label2_eval,markersize=markersize_value)

    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $Xp^s_3$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.legend(ncol=1)
    plt.grid()
    plt.savefig(output_folder_path+"/xpsin_polarmotionamp3_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(15, 6))
    ypcos3_bench_values = np.loadtxt(benchmark_folder+'/ypcos3_plot.dat')
    ypcos30_values = np.loadtxt(main0_folder+'/ypcos3_plot.dat')
    ypcos31_values = np.loadtxt(main1_folder+'/ypcos3_plot.dat')
    ypcos32_values = np.loadtxt(main2_folder+'/ypcos3_plot.dat')
    
    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(ypcos3_bench_values)/mas,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        np.array(ypcos30_values)/mas,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        np.array(ypcos31_values)/mas,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        np.array(ypcos32_values)/mas,marker_2,label=label2_eval,markersize=markersize_value)

    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $Yp^c_3$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.legend(ncol=1)
    plt.grid()
    plt.savefig(output_folder_path+"/ypcos_polarmotionamp3_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(15, 6))
    ypsin3_bench_values = np.loadtxt(benchmark_folder+'/ypsin3_plot.dat')
    ypsin30_values = np.loadtxt(main0_folder+'/ypsin3_plot.dat')
    ypsin31_values = np.loadtxt(main1_folder+'/ypsin3_plot.dat')
    ypsin32_values = np.loadtxt(main2_folder+'/ypsin3_plot.dat')
    
    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(ypsin3_bench_values)/mas,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        np.array(ypsin30_values)/mas,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        np.array(ypsin31_values)/mas,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        np.array(ypsin32_values)/mas,marker_2,label=label2_eval,markersize=markersize_value)

    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $Yp^s_3$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.legend(ncol=1)
    plt.grid()
    plt.savefig(output_folder_path+"/ypsin_polarmotionamp3_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # 1-sigma polar motion (order 4) as a function of time
    plt.figure(figsize=(15, 6))
    xpcos4_bench_values = np.loadtxt(benchmark_folder+'/xpcos4_plot.dat')
    xpcos40_values = np.loadtxt(main0_folder+'/xpcos4_plot.dat')
    xpcos41_values = np.loadtxt(main1_folder+'/xpcos4_plot.dat')
    xpcos42_values = np.loadtxt(main2_folder+'/xpcos4_plot.dat')
    
    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(xpcos4_bench_values)/mas,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        np.array(xpcos40_values)/mas,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        np.array(xpcos41_values)/mas,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        np.array(xpcos42_values)/mas,marker_2,label=label2_eval,markersize=markersize_value)

    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $Xp^c_4$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.legend(ncol=1)
    plt.grid()
    plt.savefig(output_folder_path+"/xpcos_polarmotionamp4_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(15, 6))
    xpsin4_bench_values = np.loadtxt(benchmark_folder+'/xpsin4_plot.dat')
    xpsin40_values = np.loadtxt(main0_folder+'/xpsin4_plot.dat')
    xpsin41_values = np.loadtxt(main1_folder+'/xpsin4_plot.dat')
    xpsin42_values = np.loadtxt(main2_folder+'/xpsin4_plot.dat')
    
    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(xpsin4_bench_values)/mas,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        np.array(xpsin40_values)/mas,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        np.array(xpsin41_values)/mas,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        np.array(xpsin42_values)/mas,marker_2,label=label2_eval,markersize=markersize_value)

    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $Xp^s_4$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.legend(ncol=1)
    plt.grid()
    plt.savefig(output_folder_path+"/xpsin_polarmotionamp4_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(15, 6))
    ypcos4_bench_values = np.loadtxt(benchmark_folder+'/ypcos4_plot.dat')
    ypcos40_values = np.loadtxt(main0_folder+'/ypcos4_plot.dat')
    ypcos41_values = np.loadtxt(main1_folder+'/ypcos4_plot.dat')
    ypcos42_values = np.loadtxt(main2_folder+'/ypcos4_plot.dat')
    
    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(ypcos4_bench_values)/mas,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        np.array(ypcos40_values)/mas,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        np.array(ypcos41_values)/mas,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        np.array(ypcos42_values)/mas,marker_2,label=label2_eval,markersize=markersize_value)

    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $Yp^c_4$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.legend(ncol=1)
    plt.grid()
    plt.savefig(output_folder_path+"/ypcos_polarmotionamp4_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(15, 6))
    ypsin4_bench_values = np.loadtxt(benchmark_folder+'/ypsin4_plot.dat')
    ypsin40_values = np.loadtxt(main0_folder+'/ypsin4_plot.dat')
    ypsin41_values = np.loadtxt(main1_folder+'/ypsin4_plot.dat')
    ypsin42_values = np.loadtxt(main2_folder+'/ypsin4_plot.dat')
    
    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(ypsin4_bench_values)/mas,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        np.array(ypsin40_values)/mas,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        np.array(ypsin41_values)/mas,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        np.array(ypsin42_values)/mas,marker_2,label=label2_eval,markersize=markersize_value)

    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $Yp^s_4$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.legend(ncol=1)
    plt.grid()
    plt.savefig(output_folder_path+"/ypsin_polarmotionamp4_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # 1-sigma polar motion (order 5) as a function of time
    plt.figure(figsize=(15, 6))
    xpcos5_bench_values = np.loadtxt(benchmark_folder+'/xpcos5_plot.dat')
    xpcos50_values = np.loadtxt(main0_folder+'/xpcos5_plot.dat')
    xpcos51_values = np.loadtxt(main1_folder+'/xpcos5_plot.dat')
    xpcos52_values = np.loadtxt(main2_folder+'/xpcos5_plot.dat')
    
    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(xpcos5_bench_values)/mas,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        np.array(xpcos50_values)/mas,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        np.array(xpcos51_values)/mas,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        np.array(xpcos52_values)/mas,marker_2,label=label2_eval,markersize=markersize_value)

    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $Xp^c_5$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.legend(ncol=1)
    plt.grid()
    plt.savefig(output_folder_path+"/xpcos_polarmotionamp5_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(15, 6))
    xpsin5_bench_values = np.loadtxt(benchmark_folder+'/xpsin5_plot.dat')
    xpsin50_values = np.loadtxt(main0_folder+'/xpsin5_plot.dat')
    xpsin51_values = np.loadtxt(main1_folder+'/xpsin5_plot.dat')
    xpsin52_values = np.loadtxt(main2_folder+'/xpsin5_plot.dat')
    
    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(xpsin5_bench_values)/mas,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        np.array(xpsin50_values)/mas,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        np.array(xpsin51_values)/mas,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        np.array(xpsin52_values)/mas,marker_2,label=label2_eval,markersize=markersize_value)

    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $Xp^s_5$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.legend(ncol=1)
    plt.grid()
    plt.savefig(output_folder_path+"/xpsin_polarmotionamp5_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(15, 6))
    ypcos5_bench_values = np.loadtxt(benchmark_folder+'/ypcos5_plot.dat')
    ypcos50_values = np.loadtxt(main0_folder+'/ypcos5_plot.dat')
    ypcos51_values = np.loadtxt(main1_folder+'/ypcos5_plot.dat')
    ypcos52_values = np.loadtxt(main2_folder+'/ypcos5_plot.dat')
    
    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(ypcos5_bench_values)/mas,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        np.array(ypcos50_values)/mas,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        np.array(ypcos51_values)/mas,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        np.array(ypcos52_values)/mas,marker_2,label=label2_eval,markersize=markersize_value)

    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $Yp^s_5$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.legend(ncol=1)
    plt.grid()
    plt.savefig(output_folder_path+"/ypcos_polarmotionamp5_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(15, 6))
    ypsin5_bench_values = np.loadtxt(benchmark_folder+'/ypsin5_plot.dat')
    ypsin50_values = np.loadtxt(main0_folder+'/ypsin5_plot.dat')
    ypsin51_values = np.loadtxt(main1_folder+'/ypsin5_plot.dat')
    ypsin52_values = np.loadtxt(main2_folder+'/ypsin5_plot.dat')
    
    plt.plot((time_bench_eval-time0_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(ypsin5_bench_values)/mas,marker_bench,label=label_bench_eval,markersize=markersize_value)
    plt.plot((time0_eval-time0_eval[0]*np.ones(len(time0_eval)))/constants.JULIAN_DAY,
        np.array(ypsin50_values)/mas,marker_0,label=label0_eval,markersize=markersize_value)
    plt.plot((time1_eval-time0_eval[0]*np.ones(len(time1_eval)))/constants.JULIAN_DAY,
        np.array(ypsin51_values)/mas,marker_1,label=label1_eval,markersize=markersize_value)
    plt.plot((time2_eval-time0_eval[0]*np.ones(len(time2_eval)))/constants.JULIAN_DAY,
        np.array(ypsin52_values)/mas,marker_2,label=label2_eval,markersize=markersize_value)

    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $Yp^s_5$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.legend(ncol=1)
    plt.grid()
    plt.savefig(output_folder_path+"/ypsin_polarmotionamp5_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # Observations plot
    plt.figure(figsize=(15,4))
    time_concatenated_bench_eval = np.loadtxt(benchmark_folder+'/concatenated_times_sort.dat')
    time_concatenated0_eval = np.loadtxt(main0_folder+'/concatenated_times_sort.dat')
    time_concatenated1_eval = np.loadtxt(main1_folder+'/concatenated_times_sort.dat')
    time_concatenated2_eval = np.loadtxt(main2_folder+'/concatenated_times_sort.dat')
    data = [(time_concatenated_bench_eval-time0_eval[0]*np.ones(len(time_concatenated_bench_eval)))/constants.JULIAN_DAY,\
        (time_concatenated0_eval-time0_eval[0]*np.ones(len(time_concatenated0_eval)))/constants.JULIAN_DAY,\
            (time_concatenated1_eval-time0_eval[0]*np.ones(len(time_concatenated1_eval)))/constants.JULIAN_DAY,\
                (time_concatenated2_eval-time0_eval[0]*np.ones(len(time_concatenated2_eval)))/constants.JULIAN_DAY]
    plt.hist(data,bins=int((max(time_bench_eval)-time0_eval[0])/constants.JULIAN_DAY)//84,histtype='bar',label=[label_bench_eval,label0_eval,label1_eval,label2_eval])
    plt.axvline(x=(6.944957280000000e+08-time0_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'N$^\circ$ of Obs per 12 Weeks')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time0_eval[0])))
    plt.grid()
    plt.legend()
    plt.savefig(output_folder_path+"/obs_bars.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all') 

print("--- %s seconds ---" % (time.time() - run_time))
