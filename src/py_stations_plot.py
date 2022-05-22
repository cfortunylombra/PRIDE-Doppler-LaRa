"""
Description: Plot DSN and EVN station

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
    from mpl_toolkits.basemap import Basemap # conda install basemap
    from tudatpy.kernel.astro import element_conversion

    np.set_printoptions(suppress=True,precision=15)

    ########################################################################################################################
    ################################################## CREATE GROUND STATIONS AND LANDER ###################################
    ########################################################################################################################
    
    # Earth-based transmitters
    transmitters_dict = dict() #Taken from JPL web site
    transmitters_dict['Canberra']=np.array([-4460894.9170,2682361.5070,-3674748.1517]) # DSS 43
    transmitters_dict['Madrid']=np.array([4849092.5175,-360180.3480,4115109.2506]) # DSS 63
    transmitters_dict['Goldstone']=np.array([-2353621.4197,-4641341.4717,3677052.3178]) # DSS 14

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
    
    # Function for drawing the shadedrelief
    def draw_map(m, scale=1):
        # draw a shaded-relief image
        m.shadedrelief(scale=scale)

    # Plot the antennas in a map
    cm = plt.get_cmap('gist_rainbow')
    fig = plt.figure(figsize=(8,6), edgecolor='w')
    m = Basemap(projection='cyl', resolution = 'h',
        llcrnrlat = -90, urcrnrlat = 90,
        llcrnrlon= -180, urcrnrlon = 180)

    # Add scatter points (DSN transmitters)
    j = 0
    for i in transmitters_dict.keys():
        # Conversion from cartesian to spherical state
        cartesian_state = np.zeros(6)
        cartesian_state[0] = transmitters_dict[i][0]
        cartesian_state[1] = transmitters_dict[i][1]
        cartesian_state[2] = transmitters_dict[i][2]
        cartesian_state[3] = 0
        cartesian_state[4] = 0
        cartesian_state[5] = 0
        spherical_state = element_conversion.cartesian_to_spherical( cartesian_state )
        m.scatter(np.rad2deg(spherical_state[2]), np.rad2deg(spherical_state[1]), s = 20, color = cm(j*50), label=i, edgecolors='black')
        j+=1

    # Add scatter points (radio telescopes)
    j = 0
    markers = ['x','s','^']
    for i in radio_telescopes_dict.keys():
        # Conversion from cartesian to spherical state
        cartesian_state = np.zeros(6)
        cartesian_state[0] = radio_telescopes_dict[i][0]
        cartesian_state[1] = radio_telescopes_dict[i][1]
        cartesian_state[2] = radio_telescopes_dict[i][2]
        cartesian_state[3] = 0
        cartesian_state[4] = 0
        cartesian_state[5] = 0
        spherical_state = element_conversion.cartesian_to_spherical( cartesian_state )
        m.scatter(np.rad2deg(spherical_state[2]), np.rad2deg(spherical_state[1]), s = 10, color = cm(j*20), marker=markers[j%len(markers)],label=i)
        j+=1
        
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    draw_map(m)

    fig.savefig('output/DSN-EVN-stations.pdf',bbox_inches="tight")

print("--- %s seconds ---" % (time.time() - run_time))