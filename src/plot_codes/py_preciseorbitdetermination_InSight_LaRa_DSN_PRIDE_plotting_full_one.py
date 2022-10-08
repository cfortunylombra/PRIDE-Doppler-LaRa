"""
Description: Plot data

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

    main_folder = os.path.dirname(os.path.realpath(__file__)).replace('/src','/output/POD_RISETrue_LaRaTrue_PRIDEcomplexTrueFalse')
    output_folder_path = os.path.dirname(os.path.realpath(__file__)).replace('/src','/output/POD_comparison_plot_one')
    os.makedirs(output_folder_path,exist_ok=True)

    # Booleans to understand whether we want to simulate together RISE and LaRa missions, or separetely
    RISE_boolean = True
    LaRa_boolean = True

    mas =np.pi/(180.0*1000.0*3600.0) # Conversion from milli arc seconds to seconds 
    
    time_eval = np.loadtxt(main_folder+'/time_plot.dat')

    label_eval = r'- DSN + PRIDE & variable $\rho$'

    plt.figure(figsize=(4,4))
    plt.rcParams.update({'font.size': 18})
    F_values = np.loadtxt(main_folder+'/corefactor_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        F_values,'ko-',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        F_values[0],'ro',label="A priori",markersize=8)
    plt.plot((time_eval[23600]-time_eval[0])/constants.JULIAN_DAY/7,
        F_values[23600],'gs',label="only RISE",markersize=8)
    plt.plot((time_eval[-1]-time_eval[0])/constants.JULIAN_DAY/7,
        F_values[-1],'s',color='orange',label="RISE + LaRa\n with PRIDE",markersize=8)
    #plt.legend(loc=1, prop={'size': 15})
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'Uncertainty of $F$ [-]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/Fvalue_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(3,3))
    plt.rcParams.update({'font.size': 18})
    #plt.gca().invert_yaxis()
    fcnrate_values = np.loadtxt(main_folder+'/sigmaFCN_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        np.rad2deg(fcnrate_values)*constants.JULIAN_DAY,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        np.rad2deg(fcnrate_values[0])*constants.JULIAN_DAY,'s',markersize=8)
    print(np.rad2deg(fcnrate_values[0])*constants.JULIAN_DAY)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $\sigma_{FCN}$ [deg/day]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/fcnratevalue_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    mas =np.pi/(180.0*1000.0*3600.0)

    plt.figure(figsize=(4,4))
    plt.rcParams.update({'font.size': 18})
    cos1spin_values = np.loadtxt(main_folder+'/cos1spin_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        cos1spin_values/mas,'ko-',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        cos1spin_values[0]/mas,'ro',label="A priori",markersize=8)
    plt.plot((time_eval[23600]-time_eval[0])/constants.JULIAN_DAY/7,
        cos1spin_values[23600]/mas,'gs',label="only RISE",markersize=8)
    plt.plot((time_eval[-1]-time_eval[0])/constants.JULIAN_DAY/7,
        cos1spin_values[-1]/mas,'s',color='orange',label="RISE + LaRa\n with PRIDE",markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'Uncertainty of $\phi_1^c$ [mas]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/cos1spinvalue_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(3,3))
    plt.rcParams.update({'font.size': 18})
    sin1spin_values = np.loadtxt(main_folder+'/sin1spin_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        sin1spin_values/mas,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        sin1spin_values[0]/mas,'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $\phi_1^s$ [mas]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/sin1spinvalue_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(3,3))
    plt.rcParams.update({'font.size': 18})
    cos2spin_values = np.loadtxt(main_folder+'/cos2spin_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        cos2spin_values/mas,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        cos2spin_values[0]/mas,'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $\phi_2^c$ [mas]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/cos2spinvalue_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(3,3))
    plt.rcParams.update({'font.size': 18})
    sin2spin_values = np.loadtxt(main_folder+'/sin2spin_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        sin2spin_values/mas,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        sin2spin_values[0]/mas,'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $\phi_2^s$ [mas]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/sin2spinvalue_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(3,3))
    plt.rcParams.update({'font.size': 18})
    cos3spin_values = np.loadtxt(main_folder+'/cos3spin_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        cos3spin_values/mas,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        cos3spin_values[0]/mas,'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $\phi_3^c$ [mas]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/cos3spinvalue_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(3,3))
    plt.rcParams.update({'font.size': 18})
    sin3spin_values = np.loadtxt(main_folder+'/sin3spin_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        sin3spin_values/mas,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        sin3spin_values[0]/mas,'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $\phi_3^s$ [mas]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/sin3spinvalue_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(3,3))
    plt.rcParams.update({'font.size': 18})
    cos4spin_values = np.loadtxt(main_folder+'/cos4spin_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        cos4spin_values/mas,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        cos4spin_values[0]/mas,'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $\phi_4^c$ [mas]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/cos4spinvalue_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(3,3))
    plt.rcParams.update({'font.size': 18})
    sin4spin_values = np.loadtxt(main_folder+'/sin4spin_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        sin4spin_values/mas,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        sin4spin_values[0]/mas,'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $\phi_4^s$ [mas]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/sin4spinvalue_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(4,4))
    plt.rcParams.update({'font.size': 18})
    xpcos1_values = np.loadtxt(main_folder+'/xpcos1_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        xpcos1_values/mas,'ko-',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        xpcos1_values[0]/mas,'ro',label="A priori",markersize=8)
    plt.plot((time_eval[23600]-time_eval[0])/constants.JULIAN_DAY/7,
        xpcos1_values[23600]/mas,'gs',label="only RISE",markersize=8)
    plt.plot((time_eval[-1]-time_eval[0])/constants.JULIAN_DAY/7,
        xpcos1_values[-1]/mas,'s',color='orange',label="RISE + LaRa\n with PRIDE",markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.legend(loc=1, prop={'size': 12})
    plt.ylabel(r'Uncertainty of ${X_p}_1^c$ [mas]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/xpcos1value_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(3,3))
    plt.rcParams.update({'font.size': 18})
    xpsin1_values = np.loadtxt(main_folder+'/xpsin1_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        xpsin1_values/mas,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        xpsin1_values[0]/mas,'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ ${X_p}_1^s$ [mas]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/xpsin1value_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(3,3))
    plt.rcParams.update({'font.size': 18})
    ypcos1_values = np.loadtxt(main_folder+'/ypcos1_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        ypcos1_values/mas,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        ypcos1_values[0]/mas,'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ ${Y_p}_1^c$ [mas]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/ypcos1value_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(3,3))
    plt.rcParams.update({'font.size': 18})
    ypsin1_values = np.loadtxt(main_folder+'/ypsin1_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        ypsin1_values/mas,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        ypsin1_values[0]/mas,'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ ${Y_p}_1^s$ [mas]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/ypsin1value_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(3,3))
    plt.rcParams.update({'font.size': 18})
    xpcos2_values = np.loadtxt(main_folder+'/xpcos2_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        xpcos2_values/mas,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        xpcos2_values[0]/mas,'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ ${X_p}_2^c$ [mas]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/xpcos2value_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(3,3))
    plt.rcParams.update({'font.size': 18})
    xpsin2_values = np.loadtxt(main_folder+'/xpsin2_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        xpsin2_values/mas,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        xpsin2_values[0]/mas,'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ ${X_p}_2^s$ [mas]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/xpsin2value_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(3,3))
    plt.rcParams.update({'font.size': 18})
    ypcos2_values = np.loadtxt(main_folder+'/ypcos2_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        ypcos2_values/mas,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        ypcos2_values[0]/mas,'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ ${Y_p}_2^c$ [mas]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/ypcos2value_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(3,3))
    plt.rcParams.update({'font.size': 18})
    ypsin2_values = np.loadtxt(main_folder+'/ypsin2_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        ypsin2_values/mas,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        ypsin2_values[0]/mas,'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ ${Y_p}_2^s$ [mas]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/ypsin2value_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(3,3))
    plt.rcParams.update({'font.size': 18})
    xpcos3_values = np.loadtxt(main_folder+'/xpcos3_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        xpcos3_values/mas,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        xpcos3_values[0]/mas,'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ ${X_p}_3^c$ [mas]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/xpcos3value_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(3,3))
    plt.rcParams.update({'font.size': 18})
    xpsin3_values = np.loadtxt(main_folder+'/xpsin3_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        xpsin3_values/mas,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        xpsin3_values[0]/mas,'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ ${X_p}_3^s$ [mas]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/xpsin3value_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(3,3))
    plt.rcParams.update({'font.size': 18})
    ypcos3_values = np.loadtxt(main_folder+'/ypcos3_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        ypcos3_values/mas,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        ypcos3_values[0]/mas,'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ ${Y_p}_3^c$ [mas]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/ypcos3value_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(3,3))
    plt.rcParams.update({'font.size': 18})
    ypsin3_values = np.loadtxt(main_folder+'/ypsin3_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        ypsin3_values/mas,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        ypsin3_values[0]/mas,'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ ${Y_p}_3^s$ [mas]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/ypsin3value_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(3,3))
    plt.rcParams.update({'font.size': 18})
    xpcos4_values = np.loadtxt(main_folder+'/xpcos4_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        xpcos4_values/mas,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        xpcos4_values[0]/mas,'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ ${X_p}_C^c$ [mas]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/xpcos4value_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(3,3))
    plt.rcParams.update({'font.size': 18})
    xpsin4_values = np.loadtxt(main_folder+'/xpsin4_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        xpsin4_values/mas,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        xpsin4_values[0]/mas,'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ ${X_p}_C^s$ [mas]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/xpsin4value_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(3,3))
    plt.rcParams.update({'font.size': 18})
    ypcos4_values = np.loadtxt(main_folder+'/ypcos4_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        ypcos4_values/mas,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        ypcos4_values[0]/mas,'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ ${Y_p}_C^c$ [mas]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/ypcos4value_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(3,3))
    plt.rcParams.update({'font.size': 18})
    ypsin4_values = np.loadtxt(main_folder+'/ypsin4_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        ypsin4_values/mas,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        ypsin4_values[0]/mas,'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ ${Y_p}_C^s$ [mas]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/ypsin4value_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(3,3))
    plt.rcParams.update({'font.size': 18})
    xpcos5_values = np.loadtxt(main_folder+'/xpcos5_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        xpcos5_values/mas,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        xpcos5_values[0]/mas,'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ ${X_p}_4^c$ [mas]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/xpcos5value_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(3,3))
    plt.rcParams.update({'font.size': 18})
    xpsin5_values = np.loadtxt(main_folder+'/xpsin5_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        xpsin5_values/mas,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        xpsin5_values[0]/mas,'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ ${X_p}_4^s$ [mas]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/xpsin5value_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(3,3))
    plt.rcParams.update({'font.size': 18})
    ypcos5_values = np.loadtxt(main_folder+'/ypcos5_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        ypcos5_values/mas,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        ypcos5_values[0]/mas,'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ ${Y_p}_4^c$ [mas]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/ypcos5value_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(3,3))
    plt.rcParams.update({'font.size': 18})
    ypsin5_values = np.loadtxt(main_folder+'/ypsin5_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        ypsin5_values/mas,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        ypsin5_values[0]/mas,'s',markersize=8,label='A priori')
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of ExoMars mission')
    plt.ylabel(r'1-$\sigma$ ${Y_p}_4^s$ [mas]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    plt.legend()
    plt.savefig(output_folder_path+"/ypsin5value_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')  

    plt.figure(figsize=(5,5))
    plt.rcParams.update({'font.size': 18})
    x_values = np.loadtxt(main_folder+'/xposition_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        x_values,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        x_values[0],'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $x$ [m]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/xvalue_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')  

    plt.figure(figsize=(5,5))
    plt.rcParams.update({'font.size': 18})
    y_values = np.loadtxt(main_folder+'/yposition_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        y_values,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        y_values[0],'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $y$ [m]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/yvalue_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(5,5))
    plt.rcParams.update({'font.size': 18})
    z_values = np.loadtxt(main_folder+'/zposition_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        z_values,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        z_values[0],'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $z$ [m]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/zvalue_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all') 

    plt.figure(figsize=(5,5))
    plt.rcParams.update({'font.size': 18})
    xdot_values = np.loadtxt(main_folder+'/xdotvelocity_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        xdot_values,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        xdot_values[0],'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $\dot{x}$ [m/s]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/xdotvalue_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all') 

    plt.figure(figsize=(5,5))
    plt.rcParams.update({'font.size': 18})
    ydot_values = np.loadtxt(main_folder+'/ydotvelocity_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        ydot_values,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        ydot_values[0],'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $\dot{y}$ [m/s]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/ydotvalue_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(5,5))
    plt.rcParams.update({'font.size': 18})
    zdot_values = np.loadtxt(main_folder+'/zdotvelocity_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        zdot_values,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        zdot_values[0],'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $\dot{z}$ [m/s]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/zdotvalue_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(5,5))
    plt.rcParams.update({'font.size': 18})
    xRISE_values = np.loadtxt(main_folder+'/xRISE_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        xRISE_values,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        xRISE_values[0],'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $x_{RISE}$ [m]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/xRISEvalue_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(5,5))
    plt.rcParams.update({'font.size': 18})
    yRISE_values = np.loadtxt(main_folder+'/yRISE_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        yRISE_values,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        yRISE_values[0],'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $y_{RISE}$ [m]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/yRISEvalue_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(5,5))
    plt.rcParams.update({'font.size': 18})
    zRISE_values = np.loadtxt(main_folder+'/zRISE_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        zRISE_values,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        zRISE_values[0],'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $z_{RISE}$ [m]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/zRISEvalue_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(5,5))
    plt.rcParams.update({'font.size': 18})
    xLaRa_values = np.loadtxt(main_folder+'/xLaRa_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        xLaRa_values,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        xLaRa_values[0],'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $x_{LaRa}$ [m]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/xLaRavalue_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(5,5))
    plt.rcParams.update({'font.size': 18})
    yLaRa_values = np.loadtxt(main_folder+'/yLaRa_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        yLaRa_values,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        yLaRa_values[0],'s',markersize=8)
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $y_{LaRa}$ [m]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    #plt.legend()
    plt.savefig(output_folder_path+"/yLaRavalue_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

    plt.figure(figsize=(5,5))
    plt.rcParams.update({'font.size': 18})
    zLaRa_values = np.loadtxt(main_folder+'/zLaRa_plot.dat')
    plt.plot((time_eval-time_eval[0]*np.ones(len(time_eval)))/constants.JULIAN_DAY/7,
        zLaRa_values,'-o',markersize=3)
    plt.plot((time_eval[0]-time_eval[0])/constants.JULIAN_DAY/7,
        zLaRa_values[0],'s',markersize=8,label='A priori')
    plt.axvline(x=(6.944957280000000e+08-time_eval[0])/constants.JULIAN_DAY/7, color='k', linestyle='--')#,label='Start of LaRa mission')
    plt.ylabel(r'1-$\sigma$ $z_{LaRa}$ [m]',labelpad=0)
    plt.xlabel('Weeks',labelpad=0)
    #plt.title('Start Date: '+str(datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=time_eval[0])))
    plt.grid()
    plt.legend()
    plt.savefig(output_folder_path+"/zLaRavalue_time.png",bbox_inches="tight")
    plt.show()
    plt.close('all')

print("--- %s seconds ---" % (time.time() - run_time))
