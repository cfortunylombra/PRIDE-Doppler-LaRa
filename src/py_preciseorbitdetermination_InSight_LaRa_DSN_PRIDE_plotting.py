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

    benchmark_folder = os.path.dirname(os.path.realpath(__file__)).replace('/src','/output/POD_RISEFalse_LaRaTrue_PRIDEFalseFalse_corr0')
    main_folder = os.path.dirname(os.path.realpath(__file__)).replace('/src','/output/POD_RISEFalse_LaRaTrue_PRIDETrueTrue_corr0')
    output_folder_path = os.path.dirname(os.path.realpath(__file__)).replace('/src','/output/POD_comparison_plot')
    os.makedirs(output_folder_path,exist_ok=True)

    # Booleans to understand whether we want to simulate together RISE and LaRa missions, or separetely
    RISE_boolean = False
    LaRa_boolean = True

    mas =np.pi/(180.0*1000.0*3600.0) # Conversion from milli arc seconds to seconds 
    
    time_eval = np.loadtxt(main_folder+'/time_plot.dat')
    time_bench_eval = np.loadtxt(benchmark_folder+'/time_plot.dat')

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

    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        x_values,'-o',label='$x$')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        y_values,'-o',label='$y$')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        z_values,'-o',label='$z$')

    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        x_bench_values,'--',label='$x_b$')
    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        y_bench_values,'--',label='$y_b$')
    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        z_bench_values,'--',label='$z_b$')
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ x,y,z [m]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    plt.legend()
    plt.savefig(output_folder_path+"/positionvalues_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # 1-sigma position as a velocity of time
    plt.figure(figsize=(15, 6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 3))))
    xdot_values = np.loadtxt(main_folder+'/xdotvelocity_plot.dat')
    ydot_values = np.loadtxt(main_folder+'/ydotvelocity_plot.dat')
    zdot_values = np.loadtxt(main_folder+'/zdotvelocity_plot.dat')
    xdot_bench_values = np.loadtxt(benchmark_folder+'/xdotvelocity_plot.dat')
    ydot_bench_values = np.loadtxt(benchmark_folder+'/ydotvelocity_plot.dat')
    zdot_bench_values = np.loadtxt(benchmark_folder+'/zdotvelocity_plot.dat')
    
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        xdot_values,'-o',label=r'$\dot{x}$')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        ydot_values,'-o',label=r'$\dot{y}$')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        zdot_values,'-o',label=r'$\dot{z}$')

    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        xdot_bench_values,'--',label=r'$\dot{x}_b$')
    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        ydot_bench_values,'--',label=r'$\dot{y}_b$')
    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        zdot_bench_values,'--',label=r'$\dot{z}_b$')
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $\dot{x}$,$\dot{y}$,$\dot{z}$ [m/s]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    plt.legend()
    plt.savefig(output_folder_path+"/velocityvalues_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # 1-sigma F as a function of time
    plt.figure(figsize=(15, 6))
    F_values = np.loadtxt(main_folder+'/corefactor_plot.dat')
    F_bench_values = np.loadtxt(benchmark_folder+'/corefactor_plot.dat')
    
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        F_values,'-o',label='F')
    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        F_bench_values,'--',label=r'F$_b$')
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ F [-]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    plt.legend()
    plt.savefig(output_folder_path+"/Fvalues_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # 1-sigma sigma_FCN as a function of time
    plt.figure(figsize=(15, 6))
    sigma_FCN_values = np.loadtxt(main_folder+'/sigmaFCN_plot.dat')
    sigma_FCN_bench_values = np.loadtxt(benchmark_folder+'/sigmaFCN_plot.dat')
    
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        sigma_FCN_values,'-o',label=r'$\sigma_{FCN}$')
    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        sigma_FCN_bench_values,'--',label=r'${\sigma_{FCN}}_b$')

    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $\sigma_{FCN}$ [rad/s]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    plt.legend()
    plt.savefig(output_folder_path+"/sigmaFCNvalues_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # 1-sigma x_RISE,y_RISE,z_RISE,x_LaRa,y_LaRa,z_LaRa as a function of time
    plt.figure(figsize=(15, 6))
    colormap = plt.cm.gist_ncar
    if RISE_boolean and LaRa_boolean:
        plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 6))))
    else:
        plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 3))))
    if RISE_boolean:
        xRISElander_values = np.loadtxt(main_folder+'/xRISE_plot.dat')
        yRISElander_values = np.loadtxt(main_folder+'/yRISE_plot.dat')
        zRISElander_values = np.loadtxt(main_folder+'/zRISE_plot.dat')

        plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
            xRISElander_values,'-o',label='$x_{RISE}$')
        plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
            yRISElander_values,'-o',label='$y_{RISE}$')
        plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
            zRISElander_values,'-o',label='$z_{RISE}$')
    
    if LaRa_boolean:
        xLaRalander_values = np.loadtxt(main_folder+'/xLaRa_plot.dat')
        yLaRalander_values = np.loadtxt(main_folder+'/yLaRa_plot.dat')
        zLaRalander_values = np.loadtxt(main_folder+'/zLaRa_plot.dat')

        plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
            xLaRalander_values,'-o',label='$x_{LaRa}$')
        plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
            yLaRalander_values,'-o',label='$y_{LaRa}$')
        plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
            zLaRalander_values,'-o',label='$z_{LaRa}$')

    if RISE_boolean:
        xRISElander_bench_values = np.loadtxt(benchmark_folder+'/xRISE_plot.dat')
        yRISElander_bench_values = np.loadtxt(benchmark_folder+'/yRISE_plot.dat')
        zRISElander_bench_values = np.loadtxt(benchmark_folder+'/zRISE_plot.dat')

        plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
            xRISElander_bench_values,'--',label='${x_{RISE}}_b$')
        plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
            xRISElander_bench_values,'--',label='${y_{RISE}}_b$')
        plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
            xRISElander_bench_values,'--',label='${z_{RISE}}_b$')
    
    if LaRa_boolean:
        xLaRalander_bench_values = np.loadtxt(benchmark_folder+'/xLaRa_plot.dat')
        yLaRalander_bench_values = np.loadtxt(benchmark_folder+'/yLaRa_plot.dat')
        zLaRalander_bench_values = np.loadtxt(benchmark_folder+'/zLaRa_plot.dat')

        plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
            xLaRalander_bench_values,'--',label='${x_{LaRa}}_b$')
        plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
            yLaRalander_bench_values,'--',label='${y_{LaRa}}_b$')
        plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
            zLaRalander_bench_values,'--',label='${z_{LaRa}}_b$')

    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ x,y,z [m]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.legend()
    plt.grid()
    plt.savefig(output_folder_path+"/xyzlander_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # 1-sigma spin variations as a function of time
    plt.figure(figsize=(15, 6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 8))))
    cos1spin_values = np.loadtxt(main_folder+'/cos1spin_plot.dat')
    sin1spin_values = np.loadtxt(main_folder+'/sin1spin_plot.dat')
    cos2spin_values = np.loadtxt(main_folder+'/cos2spin_plot.dat')
    sin2spin_values = np.loadtxt(main_folder+'/sin2spin_plot.dat')
    cos3spin_values = np.loadtxt(main_folder+'/cos3spin_plot.dat')
    sin3spin_values = np.loadtxt(main_folder+'/sin3spin_plot.dat')
    cos4spin_values = np.loadtxt(main_folder+'/cos4spin_plot.dat')
    sin4spin_values = np.loadtxt(main_folder+'/sin4spin_plot.dat')

    cos1spin_bench_values = np.loadtxt(benchmark_folder+'/cos1spin_plot.dat')
    sin1spin_bench_values = np.loadtxt(benchmark_folder+'/sin1spin_plot.dat')
    cos2spin_bench_values = np.loadtxt(benchmark_folder+'/cos2spin_plot.dat')
    sin2spin_bench_values = np.loadtxt(benchmark_folder+'/sin2spin_plot.dat')
    cos3spin_bench_values = np.loadtxt(benchmark_folder+'/cos3spin_plot.dat')
    sin3spin_bench_values = np.loadtxt(benchmark_folder+'/sin3spin_plot.dat')
    cos4spin_bench_values = np.loadtxt(benchmark_folder+'/cos4spin_plot.dat')
    sin4spin_bench_values = np.loadtxt(benchmark_folder+'/sin4spin_plot.dat')
    
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(cos1spin_values)/mas,'-o',label=r'$\psi^c_1$')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(sin1spin_values)/mas,'-o',label=r'$\psi^s_1$')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(cos2spin_values)/mas,'-o',label=r'$\psi^c_2$')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(sin2spin_values)/mas,'-o',label=r'$\psi^s_2$')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(cos3spin_values)/mas,'-o',label=r'$\psi^c_3$')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(sin3spin_values)/mas,'-o',label=r'$\psi^s_3$')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(cos4spin_values)/mas,'-o',label=r'$\psi^c_4$')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(sin4spin_values)/mas,'-o',label=r'$\psi^s_4$')

    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(cos1spin_bench_values)/mas,'--',label=r'${\psi^c_1}_b$')
    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(sin1spin_bench_values)/mas,'--',label=r'${\psi^s_1}_b$')
    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(cos2spin_bench_values)/mas,'--',label=r'${\psi^c_2}_b$')
    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(sin2spin_bench_values)/mas,'--',label=r'${\psi^s_2}_b$')
    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(cos3spin_bench_values)/mas,'--',label=r'${\psi^c_3}_b$')
    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(sin3spin_bench_values)/mas,'--',label=r'${\psi^s_3}_b$')
    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(cos4spin_bench_values)/mas,'--',label=r'${\psi^c_4}_b$')
    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(sin4spin_bench_values)/mas,'--',label=r'${\psi^s_4}_b$')
    
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $\psi$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.legend(ncol=2)
    plt.grid()
    plt.savefig(output_folder_path+"/psispin_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # 1-sigma polar motion (order 1) as a function of time
    plt.figure(figsize=(15, 6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 4))))
    xpcos1_values = np.loadtxt(main_folder+'/xpcos1_plot.dat')
    xpsin1_values = np.loadtxt(main_folder+'/xpsin1_plot.dat')
    ypcos1_values = np.loadtxt(main_folder+'/ypcos1_plot.dat')
    ypsin1_values = np.loadtxt(main_folder+'/ypsin1_plot.dat')

    xpcos1_bench_values = np.loadtxt(benchmark_folder+'/xpcos1_plot.dat')
    xpsin1_bench_values = np.loadtxt(benchmark_folder+'/xpsin1_plot.dat')
    ypcos1_bench_values = np.loadtxt(benchmark_folder+'/ypcos1_plot.dat')
    ypsin1_bench_values = np.loadtxt(benchmark_folder+'/ypsin1_plot.dat')
    
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(xpcos1_values)/mas,'-o',label=r'$Xp^c_1$')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(xpsin1_values)/mas,'-o',label=r'$Xp^s_1$')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(ypcos1_values)/mas,'-o',label=r'$Yp^c_1$')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(ypsin1_values)/mas,'-o',label=r'$Yp^s_1$')

    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(xpcos1_bench_values)/mas,'--',label=r'${Xp^c_1}_b$')
    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(xpsin1_bench_values)/mas,'--',label=r'${Xp^s_1}_b$')
    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(ypcos1_bench_values)/mas,'--',label=r'${Yp^c_1}_b$')
    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(ypsin1_bench_values)/mas,'--',label=r'${Yp^s_1}_b$')

    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $Xp, Yp$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.legend(ncol=2)
    plt.grid()
    plt.savefig(output_folder_path+"/polarmotionamp1_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # 1-sigma polar motion (order 2) as a function of time
    plt.figure(figsize=(15, 6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 4))))
    xpcos2_values = np.loadtxt(main_folder+'/xpcos2_plot.dat')
    xpsin2_values = np.loadtxt(main_folder+'/xpsin2_plot.dat')
    ypcos2_values = np.loadtxt(main_folder+'/ypcos2_plot.dat')
    ypsin2_values = np.loadtxt(main_folder+'/ypsin2_plot.dat')

    xpcos2_bench_values = np.loadtxt(benchmark_folder+'/xpcos2_plot.dat')
    xpsin2_bench_values = np.loadtxt(benchmark_folder+'/xpsin2_plot.dat')
    ypcos2_bench_values = np.loadtxt(benchmark_folder+'/ypcos2_plot.dat')
    ypsin2_bench_values = np.loadtxt(benchmark_folder+'/ypsin2_plot.dat')
    
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(xpcos2_values)/mas,'-o',label=r'$Xp^c_2$')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(xpsin2_values)/mas,'-o',label=r'$Xp^s_2$')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(ypcos2_values)/mas,'-o',label=r'$Yp^c_2$')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(ypsin2_values)/mas,'-o',label=r'$Yp^s_2$')

    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(xpcos2_bench_values)/mas,'--',label=r'${Xp^c_2}_b$')
    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(xpsin2_bench_values)/mas,'--',label=r'${Xp^s_2}_b$')
    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(ypcos2_bench_values)/mas,'--',label=r'${Yp^c_2}_b$')
    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(ypsin2_bench_values)/mas,'--',label=r'${Yp^s_2}_b$')

    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $Xp, Yp$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.legend(ncol=2)
    plt.grid()
    plt.savefig(output_folder_path+"/polarmotionamp2_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # 1-sigma polar motion (order 3) as a function of time
    plt.figure(figsize=(15, 6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 4))))
    xpcos3_values = np.loadtxt(main_folder+'/xpcos3_plot.dat')
    xpsin3_values = np.loadtxt(main_folder+'/xpsin3_plot.dat')
    ypcos3_values = np.loadtxt(main_folder+'/ypcos3_plot.dat')
    ypsin3_values = np.loadtxt(main_folder+'/ypsin3_plot.dat')

    xpcos3_bench_values = np.loadtxt(benchmark_folder+'/xpcos3_plot.dat')
    xpsin3_bench_values = np.loadtxt(benchmark_folder+'/xpsin3_plot.dat')
    ypcos3_bench_values = np.loadtxt(benchmark_folder+'/ypcos3_plot.dat')
    ypsin3_bench_values = np.loadtxt(benchmark_folder+'/ypsin3_plot.dat')
    
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(xpcos3_values)/mas,'-o',label=r'$Xp^c_3$')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(xpsin3_values)/mas,'-o',label=r'$Xp^s_3$')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(ypcos3_values)/mas,'-o',label=r'$Yp^c_3$')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(ypsin3_values)/mas,'-o',label=r'$Yp^s_3$')

    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(xpcos3_bench_values)/mas,'--',label=r'${Xp^c_3}_b$')
    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(xpsin3_bench_values)/mas,'--',label=r'${Xp^s_3}_b$')
    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(ypcos3_bench_values)/mas,'--',label=r'${Yp^c_3}_b$')
    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(ypsin3_bench_values)/mas,'--',label=r'${Yp^s_3}_b$')

    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $Xp, Yp$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.legend(ncol=2)
    plt.grid()
    plt.savefig(output_folder_path+"/polarmotionamp3_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # 1-sigma polar motion (order 4) as a function of time
    plt.figure(figsize=(15, 6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 4))))
    xpcos4_values = np.loadtxt(main_folder+'/xpcos4_plot.dat')
    xpsin4_values = np.loadtxt(main_folder+'/xpsin4_plot.dat')
    ypcos4_values = np.loadtxt(main_folder+'/ypcos4_plot.dat')
    ypsin4_values = np.loadtxt(main_folder+'/ypsin4_plot.dat')

    xpcos4_bench_values = np.loadtxt(benchmark_folder+'/xpcos4_plot.dat')
    xpsin4_bench_values = np.loadtxt(benchmark_folder+'/xpsin4_plot.dat')
    ypcos4_bench_values = np.loadtxt(benchmark_folder+'/ypcos4_plot.dat')
    ypsin4_bench_values = np.loadtxt(benchmark_folder+'/ypsin4_plot.dat')
    
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(xpcos4_values)/mas,'-o',label=r'$Xp^c_4$')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(xpsin4_values)/mas,'-o',label=r'$Xp^s_4$')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(ypcos4_values)/mas,'-o',label=r'$Yp^c_4$')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(ypsin4_values)/mas,'-o',label=r'$Yp^s_4$')

    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(xpcos4_bench_values)/mas,'--',label=r'${Xp^c_4}_b$')
    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(xpsin4_bench_values)/mas,'--',label=r'${Xp^s_4}_b$')
    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(ypcos4_bench_values)/mas,'--',label=r'${Yp^c_4}_b$')
    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(ypsin4_bench_values)/mas,'--',label=r'${Yp^s_4}_b$')

    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $Xp, Yp$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.legend(ncol=2)
    plt.grid()
    plt.savefig(output_folder_path+"/polarmotionamp4_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

    # 1-sigma polar motion (order 5) as a function of time
    plt.figure(figsize=(15, 6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 4))))
    xpcos5_values = np.loadtxt(main_folder+'/xpcos5_plot.dat')
    xpsin5_values = np.loadtxt(main_folder+'/xpsin5_plot.dat')
    ypcos5_values = np.loadtxt(main_folder+'/ypcos5_plot.dat')
    ypsin5_values = np.loadtxt(main_folder+'/ypsin5_plot.dat')

    xpcos5_bench_values = np.loadtxt(benchmark_folder+'/xpcos5_plot.dat')
    xpsin5_bench_values = np.loadtxt(benchmark_folder+'/xpsin5_plot.dat')
    ypcos5_bench_values = np.loadtxt(benchmark_folder+'/ypcos5_plot.dat')
    ypsin5_bench_values = np.loadtxt(benchmark_folder+'/ypsin5_plot.dat')
    
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(xpcos5_values)/mas,'-o',label=r'$Xp^c_5$')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(xpsin5_values)/mas,'-o',label=r'$Xp^s_5$')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(ypcos5_values)/mas,'-o',label=r'$Yp^c_5$')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY,
        np.array(ypsin5_values)/mas,'-o',label=r'$Yp^s_5$')

    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(xpcos5_bench_values)/mas,'--',label=r'${Xp^c_5}_b$')
    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(xpsin5_bench_values)/mas,'--',label=r'${Xp^s_5}_b$')
    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(ypcos5_bench_values)/mas,'--',label=r'${Yp^c_5}_b$')
    plt.plot((time_bench_eval-time_eval[0]*np.ones(len(time_bench_eval)))/constants.JULIAN_DAY,
        np.array(ypsin5_bench_values)/mas,'--',label=r'${Yp^s_5}_b$')

    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY, color='k', linestyle='--',label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $Xp, Yp$ [mas]')
    plt.xlabel('Time [days]')
    plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.legend(ncol=2)
    plt.grid()
    plt.savefig(output_folder_path+"/polarmotionamp5_time.pdf",bbox_inches="tight")
    plt.show()
    plt.close('all')

print("--- %s seconds ---" % (time.time() - run_time))
