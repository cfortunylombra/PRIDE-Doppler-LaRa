"""
Description: Metric

Author: C. Fortuny-Lombraña
"""

import time

run_time = time.time()
if __name__=="__main__":
    ########################################################################################################################
    ################################################## IMPORT PACKAGES #####################################################
    ########################################################################################################################

    import numpy as np

    coordinates = dict()
    coordinates['DSS 43']=np.array([-4460894.9170,2682361.5070,-3674748.1517])
    coordinates['DSS 63']=np.array([4849092.5175,-360180.3480,4115109.2506])
    coordinates['DSS 14']=np.array([-2353621.4197,-4641341.4717,3677052.3178])
    coordinates['BADARY'] = np.array([-838201.2618,3865751.5589,4987670.8708])
    coordinates['CEDUNA'] = np.array([-3753442.7457,3912709.7530,-3348067.6095])
    coordinates['HARTRAO'] = np.array([5085442.7721,2668263.9300,-2768696.6299])
    coordinates['TIANMA65'] = np.array([-2826708.8081,4679236.9722,3274667.4495])
    coordinates['EFLSBERG'] = np.array([4033947.1525,486990.8961,4900431.0604])
    coordinates['IRBENE'] = np.array([3183649.341,1276902.985,5359264.715])
    coordinates['YEBES40M'] = np.array([4848761.7579,-261484.0570,4123085.1343])
    coordinates['MEDICINA'] = np.array([4461369.5682,919597.2489,4449559.4702])
    coordinates['WETTZELL'] = np.array([4075539.5173,931735.6497,4801629.6028])
    coordinates['ONSALA60'] = np.array([3370605.7035,711917.8146,5349830.9852])
    coordinates['WRT0'] = np.array([3828767.1338,442446.1588,5064921.5700])

    def rho_distance(station1,station2):
        distance = np.sqrt((coordinates[station1][0] - coordinates[station2][0])**2\
            + (coordinates[station1][1] - coordinates[station2][1])**2\
                + (coordinates[station1][2] - coordinates[station2][2])**2)
        print(distance)
        return np.exp(-distance/(2*413000))

    print(rho_distance('DSS 63','MEDICINA'))

    SNR = dict()
    SNR['DSS 43'] = 20.9
    SNR['DSS 63'] = 20.9
    SNR['DSS 14'] = 20.9
    SNR['BADARY'] = 21.188963221944505
    SNR['CEDUNA'] = 21.68991848541254
    SNR['HARTRAO'] = 20.870687469362966
    SNR['TIANMA65'] = 27.388973475002572
    SNR['EFLSBERG'] = 39.80319950457843
    SNR['IRBENE'] = 17.31395371999131
    SNR['YEBES40M'] = 29.539625649012528
    SNR['MEDICINA'] = 30.227433564275287
    SNR['WETTZELL'] = 29.379914821805926
    SNR['ONSALA60'] = 24.947058655336654
    SNR['WRT0'] = 21.777191716353105

    min_SNR = np.inf
    max_SNR = -np.inf

    for SNR1 in SNR:
        for SNR2 in SNR:
            if SNR1!=SNR2:
                x = (SNR[SNR1]-20)+(SNR[SNR2]-20)
                if x>max_SNR:
                    max_SNR = x
                if x<min_SNR:
                    min_SNR = x
    
    print("min",min_SNR)
    print("max",max_SNR)
        
    def rho_SNR(station1,station2):
        x = ((SNR[station1]-20)+(SNR[station2]-20))
        return 1-(x-min_SNR)/(max_SNR-min_SNR)

    print("SNR")
    print(rho_SNR('DSS 63','DSS 43'))
    print(rho_SNR('EFLSBERG','YEBES40M'))
    print(rho_SNR('WETTZELL','YEBES40M'))
    print(rho_SNR('IRBENE','EFLSBERG'))
    print(rho_SNR('IRBENE','YEBES40M'))
    print(rho_SNR('IRBENE','HARTRAO'))

print("--- %s seconds ---" % (time.time() - run_time))