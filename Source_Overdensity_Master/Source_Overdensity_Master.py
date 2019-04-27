import pandas as pd
import os
from os import system
import sys
#from statistics import mean
import numpy as np
path_Modules="/Volumes/xray/anthony/Research_Git"
sys.path.append(os.path.abspath(path_Modules))
from Galaxy_Name_Reducer import Galaxy_Name_Reducer
dir = os.path.dirname(__file__)
path=os.path.realpath('../')
sys.path.append(os.path.abspath(path))
from Source_Overdensity_Plot import Source_Overdensity_Plot
from D25_Excess_Test import D25_Excess_Test
def Source_Overdensity_Master(Gname_L,Source_Overdensity_Plot_Run_B=False):
    Tidal_Stream_Data_Path="/Volumes/xray/anthony/Research_Git/Data_Processing/Source_Overdensity/Tidal_Stream_Data/Galaxy_Sample_Info.csv"
    Tidal_Stream_Data=pd.read_csv(Tidal_Stream_Data_Path)
    #print "Tidal_Stream_Data : \n", Tidal_Stream_Data
    Gname_Raw_DF=Tidal_Stream_Data["Gname"]
    #print "Gname_Raw_DF : ", Gname_Raw_DF
    Gname_Raw_L=list(Gname_Raw_DF)
    #print "Gname_Raw_L : ", Gname_Raw_L
    Gname_TS_L=[]
    for Gname_Raw in Gname_Raw_L:
        Gname_M=Gname_Raw.split("(")[0]
        #print "Gname_M : ", Gname_M
        Gname_TS_M_Reduced=Galaxy_Name_Reducer.Galaxy_Name_Reducer(Gname_M)
        #print "Gname_TS_M_Reduced:", Gname_TS_M_Reduced
        Gname_TS_Reduced=Gname_TS_M_Reduced.replace("M","MESSIER")
        #print "Gname_TS_Reduced:", Gname_TS_Reduced
        #Gname_TS_Reduced_2=Galaxy_Name_Reducer.Galaxy_Name_Reducer(Gname_TS_Reduced)
        #print "Gname_TS_Reduced_2:", Gname_TS_Reduced_2
        Gname_TS_L.append(Gname_TS_Reduced)
    #print "Gname_TS_L : ", Gname_TS_L
    Tidal_Stream_YN_DF=Tidal_Stream_Data["Tidal Streams"]
    #print "Tidal_Stream_B_DF : ", Tidal_Stream_B_DF
    Tidal_Stream_YN_L=list(Tidal_Stream_YN_DF)
    #print "Tidal_Stream_YN_L : ", Tidal_Stream_YN_L
    Tidal_Stream_B_L=[]
    for YN in Tidal_Stream_YN_L:
        Cur_Tidal_Stream_B=(YN=="Y")
        Tidal_Stream_B_L.append(Cur_Tidal_Stream_B)
    #print "Tidal_Stream_B_L : ", Tidal_Stream_B_L
    Galaxy_TS_Data_L=[]
    for i in range(0,len(Gname_TS_L)):
        Cur_Gname_TS=Gname_TS_L[i]
        Cur_TS_B=Tidal_Stream_B_L[i]
        Cur_Galaxy_TS_Data_L=[Cur_Gname_TS,Cur_TS_B]
        Galaxy_TS_Data_L.append(Cur_Galaxy_TS_Data_L)
    #print "Galaxy_TS_Data_L : ", Galaxy_TS_Data_L
    D25_Excess_B_HL=[]
    for Gname in Gname_L:
        Gname_Modifed=Galaxy_Name_Reducer.Galaxy_Name_Reducer(Gname)
        if(Source_Overdensity_Plot_Run_B):
            Source_Overdensity_Plot.Source_Histogram_Plot(Gname,Outpath="/Volumes/xray/anthony/Research_Git/Data_Processing/Source_Overdensity/Source_Overdensity_Plot/Source_Overdensity_Histograms/")
        D25_Excess_B=D25_Excess_Test.D25_Excess_Test(Gname)
        Cur_D25_Excess_B_L=[Gname_Modifed,D25_Excess_B]
        D25_Excess_B_HL.append(Cur_D25_Excess_B_L)
    print "D25_Excess_B_HL : ", D25_Excess_B_HL
    Num_Excess=0.0
    Num_Galaxies=len(D25_Excess_B_HL)
    for D25_Excess_B_L in D25_Excess_B_HL:
        D25_Excess_Bool=D25_Excess_B_L[1]
        if(D25_Excess_Bool):
            Num_Excess=Num_Excess+1.0
    Excess_Ratio=Num_Excess/Num_Galaxies
    print "Excess_Ratio : ", Excess_Ratio
    Y_TS_L=[]
    N_TS_L=[]
    for Tidal_Stream_Data in Galaxy_TS_Data_L:
        Gname_TS=Tidal_Stream_Data[0]
        print "Gname_TS:", Gname_TS
        TS_B=Tidal_Stream_Data[1]
        #print "TS_B:",TS_B
        for D25_Excess_Data in D25_Excess_B_HL:
            Gname_Excess=D25_Excess_Data[0]
            print "Gname_Excess:",Gname_Excess
            Excess_B=D25_Excess_Data[1]
            #print "Excess_B:",Excess_B
            if(Gname_TS==Gname_Excess):
                print "Name Match"
                if(TS_B==True):
                    if(Excess_B==True):
                        Y_TS_L.append(1.0)
                    if(Excess_B==False):
                        Y_TS_L.append(0.0)
                if(TS_B==False):
                    if(Excess_B==True):
                        N_TS_L.append(1.0)
                    if(Excess_B==False):
                        N_TS_L.append(0.0)
    print "Excess_Ratio : ", Excess_Ratio
    print "Y_TS_L : ", Y_TS_L
    print "N_TS_L : ", N_TS_L
    print "len(Y_TS_L) : ", len(Y_TS_L)
    print "len(N_TS_L) : ", len(N_TS_L)
    Y_TS_A=np.array(Y_TS_L)
    #print "Y_TS_A : ", Y_TS_A
    N_TS_A=np.array(N_TS_L)
    #print "N_TS_A : ", N_TS_A
    Y_Avg=np.mean(Y_TS_A)
    #print "Y_Avg : ", Y_Avg
    N_Avg=np.mean(N_TS_A)
    #print "N_Avg : ", N_Avg
    return [Excess_Ratio,Y_Avg,N_Avg]
#Source_Overdensity_Master(["NGC_5813"])
#Source_Overdensity_Master(["NGC_253","NGC_5813"])
#print Source_Overdensity_Master(["NGC_3631"])
#print Source_Overdensity_Master(['NGC 4278', 'NGC 5204', 'NGC 2841', 'NGC 3877', 'MESSIER 106', 'NGC 5194', 'MESSIER 104', 'MESSIER 105', 'MESSIER 101', 'NGC 5054', 'NGC 5813', 'MESSIER 108', 'MESSIER 066', 'MESSIER 061', 'MESSIER 063', 'MESSIER 086', 'MESSIER 084', 'MESSIER 083', 'MESSIER 082', 'MESSIER 081', 'NGC 0278', 'MESSIER 088', 'NGC 1042', 'NGC 6744', 'IC 5332', 'NGC 3585', 'NGC 4478', 'NGC 7507', 'NGC 1637', 'NGC 4473', 'NGC 1365', 'MESSIER 074', 'NGC 4476', 'NGC 4570', 'NGC 5576', 'NGC 4321', 'NGC 5474', 'NGC 7090', 'MESSIER 094', 'MESSIER 095', 'NGC 4494', 'NGC 0247', 'NGC 4490', 'IC 1613', 'NGC 4477', 'NGC 4365', 'NGC 2787', 'NGC 3557', 'IC 5267', 'NGC 4388', 'NGC 3923', 'NGC 4945', 'NGC 891', 'NGC 1300', 'UGC 05340', 'NGC 3631', 'UGCA 166', 'NGC 4314', 'NGC 4550', 'Holmberg IX                   ', 'NGC 4559', 'IC 1459', 'NGC 1399', 'NGC 4039', 'NGC 4038', 'NGC 1316', 'NGC 1097', 'NGC 0383', 'NGC 2681', 'NGC 5018', 'NGC 5253', 'NGC 4631', 'NGC 4308', 'MESSIER 060', 'NGC 4742', 'NGC 1672', 'NGC 5846', 'NGC 4725', 'NGC 2403', 'NGC 3507', 'MESSIER 087', 'NGC 0891', 'NGC 4698', 'NGC 3384', 'NGC 6946', 'NGC 1291:[LFF2012] 084', 'NGC 3115', 'NGC 1332', 'NGC 1700', 'NGC 5584', 'NGC 4527', 'NGC 7552', 'NGC 2997', 'NGC 4449', 'MESSIER 049', 'NGC 3198', 'NGC 0855', 'NGC 7793', 'NGC 0119', 'NGC 2865', 'Circinus Galaxy               ', 'MESSIER 059', 'NGC 7331', 'NGC 1427', 'NGC 3628', 'NGC 3608', 'NGC 0055', 'NGC 4457', 'NGC 4214', 'NGC 4459', 'NGC 3521','NGC 4565', 'NGC 1313', 'NGC 0253'],Source_Overdensity_Plot_Run_B=True)
#print Source_Overdensity_Master(['NGC 4278', 'NGC 2841', 'NGC 3877', 'MESSIER 106', 'NGC 5194', 'MESSIER 104', 'MESSIER 105', 'MESSIER 101', 'NGC 5054', 'NGC 5813', 'MESSIER 108', 'MESSIER 066', 'MESSIER 061', 'MESSIER 063', 'MESSIER 086', 'MESSIER 084', 'MESSIER 083', 'MESSIER 082', 'MESSIER 081', 'NGC 0278', 'MESSIER 088', 'NGC 6744', 'IC 5332', 'NGC 3585', 'NGC 7507', 'NGC 1637', 'NGC 4473', 'NGC 1365', 'MESSIER 074', 'NGC 4476', 'NGC 4570', 'NGC 5576', 'NGC 4321', 'NGC 5474', 'NGC 7090', 'MESSIER 094', 'MESSIER 095', 'NGC 4494', 'NGC 4490', 'IC 1613', 'NGC 4477', 'NGC 4365', 'NGC 2787', 'NGC 3557', 'IC 5267', 'NGC 4388', 'NGC 3923', 'NGC 4945', 'NGC 891', 'NGC 1300', 'UGC 05340', 'NGC 3631', 'UGCA 166', 'NGC 4314', 'NGC 4550', 'Holmberg IX                   ', 'NGC 4559', 'NGC 1399', 'NGC 4039', 'NGC 4038', 'NGC 1316', 'NGC 1097', 'NGC 0383', 'NGC 2681', 'NGC 5018', 'NGC 5253', 'NGC 4631', 'MESSIER 060', 'NGC 4742', 'NGC 1672', 'NGC 5846', 'NGC 4725', 'NGC 2403', 'NGC 3507', 'MESSIER 087', 'NGC 0891', 'NGC 3384', 'NGC 6946', 'NGC 1291:[LFF2012] 084', 'NGC 3115', 'NGC 1332', 'NGC 1700', 'NGC 5584', 'NGC 7552', 'NGC 2997', 'NGC 4449', 'MESSIER 049', 'NGC 3198', 'NGC 0855', 'NGC 7793', 'NGC 0119', 'NGC 2865', 'MESSIER 059', 'NGC 1427', 'NGC 3628', 'NGC 3608', 'NGC 0055', 'NGC 4457', 'NGC 4214', 'NGC 4459', 'NGC 3521', 'NGC 4565', 'NGC 1313', 'NGC 0253'],Source_Overdensity_Plot_Run_B=True)
#NGC_1365
#Source_Overdensity_Master(["NGC_1365"])
#print Source_Overdensity_Master(['NGC 4278', 'NGC 2841', 'NGC 3877', 'MESSIER 106', 'NGC 5194', 'MESSIER 104', 'MESSIER 105', 'MESSIER 101', 'NGC 5054', 'NGC 5813', 'MESSIER 108', 'MESSIER 066', 'MESSIER 061', 'MESSIER 063', 'MESSIER 086', 'MESSIER 084', 'MESSIER 083', 'MESSIER 082', 'MESSIER 081', 'NGC 0278', 'MESSIER 088', 'NGC 6744', 'IC 5332', 'NGC 3585', 'NGC 7507', 'NGC 1637', 'NGC 4473', 'NGC 1365', 'MESSIER 074', 'NGC 4476', 'NGC 4570', 'NGC 5576', 'NGC 4321', 'NGC 5474', 'NGC 7090', 'MESSIER 094', 'MESSIER 095', 'NGC 4494', 'NGC 4490', 'IC 1613', 'NGC 4477', 'NGC 4365', 'NGC 2787', 'NGC 3557', 'IC 5267', 'NGC 4388', 'NGC 3923', 'NGC 4945', 'NGC 891', 'NGC 1300', 'UGC 05340', 'NGC 3631', 'UGCA 166', 'NGC 4314', 'NGC 4550', 'Holmberg IX                   ', 'NGC 4559', 'NGC 1399', 'NGC 4039', 'NGC 4038', 'NGC 1316', 'NGC 1097', 'NGC 0383', 'NGC 2681', 'NGC 5018', 'NGC 5253', 'NGC 4631', 'MESSIER 060', 'NGC 4742', 'NGC 1672', 'NGC 5846', 'NGC 4725', 'NGC 2403', 'NGC 3507', 'MESSIER 087', 'NGC 0891', 'NGC 3384', 'NGC 6946', 'NGC 1291:[LFF2012] 084', 'NGC 3115', 'NGC 1332', 'NGC 1700', 'NGC 5584', 'NGC 7552', 'NGC 2997', 'NGC 4449', 'MESSIER 049', 'NGC 3198', 'NGC 0855', 'NGC 7793', 'NGC 0119', 'NGC 2865', 'MESSIER 059', 'NGC 1427', 'NGC 3628', 'NGC 3608', 'NGC 0055', 'NGC 4457', 'NGC 4214', 'NGC 4459', 'NGC 3521', 'NGC 4565', 'NGC 1313', 'NGC 0253'])

#print Source_Overdensity_Master(['NGC 4278', 'NGC 2841', 'NGC 3877', 'NGC 5194', 'NGC 5054', 'NGC 5813', 'MESSIER 108', 'MESSIER 066', 'MESSIER 061', 'MESSIER 063', 'MESSIER 086', 'MESSIER 084', 'MESSIER 083', 'MESSIER 082', 'NGC 0278', 'MESSIER 088', 'NGC 3585', 'NGC 7507', 'NGC 1637', 'NGC 4473', 'NGC 1365', 'MESSIER 074', 'NGC 4570', 'NGC 5576', 'NGC 4321', 'NGC 5474', 'NGC 7090', 'MESSIER 094', 'MESSIER 095', 'NGC 4494', 'IC 1613', 'NGC 4477', 'NGC 4365', 'NGC 2787', 'NGC 3557', 'IC 5267', 'NGC 4388', 'NGC 3923', 'NGC 891', 'NGC 1300', 'UGC 05340', 'NGC 3631', 'UGCA 166', 'NGC 4314', 'NGC 4550', 'Holmberg IX                   ', 'NGC 4559', 'NGC 1399', 'NGC 1316', 'NGC 1097', 'NGC 2681', 'NGC 5018', 'NGC 5253', 'NGC 4631', 'MESSIER 060', 'NGC 4742', 'NGC 1672', 'NGC 5846', 'NGC 4725', 'NGC 3507', 'MESSIER 087', 'NGC 0891', 'NGC 3384', 'NGC 6946', 'NGC 1291:[LFF2012] 084', 'NGC 3115', 'NGC 1332', 'NGC 1700', 'NGC 5584', 'NGC 7552', 'NGC 2997', 'NGC 4449', 'MESSIER 049', 'NGC 3198', 'NGC 0855', 'NGC 7793', 'NGC 0119', 'NGC 2865', 'MESSIER 059', 'NGC 1427', 'NGC 3628', 'NGC 3608', 'NGC 0055', 'NGC 4457', 'NGC 4214', 'NGC 4459', 'NGC 3521', 'NGC 4565', 'NGC 1313'])

#print Source_Overdensity_Master(['NGC 4278', 'NGC 2841', 'NGC 3877', 'NGC 5194', 'NGC 5054', 'NGC 5813', 'MESSIER 108', 'MESSIER 066', 'MESSIER 061', 'MESSIER 063', 'MESSIER 086', 'MESSIER 084', 'MESSIER 083', 'MESSIER 082', 'NGC 0278', 'MESSIER 088', 'NGC 3585', 'NGC 7507', 'NGC 1637', 'NGC 4473', 'NGC 1365', 'MESSIER 074', 'NGC 4570', 'NGC 5576', 'NGC 4321', 'NGC 5474', 'NGC 7090', 'MESSIER 094', 'MESSIER 095', 'NGC 4494', 'IC 1613', 'NGC 4477', 'NGC 4365', 'NGC 2787', 'NGC 3557', 'IC 5267', 'NGC 4388', 'NGC 3923', 'NGC 891', 'NGC 1300', 'UGC 05340', 'NGC 3631', 'UGCA 166', 'NGC 4314', 'NGC 4550', 'Holmberg IX                   ', 'NGC 4559', 'NGC 1399', 'NGC 1316', 'NGC 1097', 'NGC 2681', 'NGC 5018', 'NGC 5253', 'NGC 4631', 'MESSIER 060', 'NGC 4742', 'NGC 1672', 'NGC 5846', 'NGC 4725', 'NGC 3507', 'MESSIER 087', 'NGC 0891', 'NGC 3384', 'NGC 6946', 'NGC 1291:[LFF2012] 084', 'NGC 3115', 'NGC 1332', 'NGC 1700', 'NGC 5584', 'NGC 7552', 'NGC 2997', 'NGC 4449', 'MESSIER 049', 'NGC 3198', 'NGC 0855', 'NGC 7793', 'NGC 0119', 'NGC 2865', 'MESSIER 059', 'NGC 1427', 'NGC 3628', 'NGC 3608', 'NGC 0055', 'NGC 4457', 'NGC 4214', 'NGC 4459', 'NGC 3521', 'NGC 4565', 'NGC 1313'],Source_Overdensity_Plot_Run_B=True)

#print Source_Overdensity_Master(['NGC 2841', 'NGC 3877', 'NGC 5054', 'NGC 5813', 'MESSIER 108', 'MESSIER 066', 'MESSIER 061', 'MESSIER 063', 'MESSIER 086', 'MESSIER 084', 'MESSIER 083', 'MESSIER 082', 'NGC 0278', 'MESSIER 088', 'NGC 3585', 'NGC 7507', 'NGC 1637', 'NGC 4473', 'NGC 1365', 'MESSIER 074', 'NGC 4570', 'NGC 4321', 'NGC 5474', 'NGC 7090', 'MESSIER 094', 'MESSIER 095', 'NGC 4494', 'IC 1613', 'NGC 4477', 'NGC 2787', 'IC 5267', 'NGC 3923', 'NGC 891', 'NGC 1300', 'UGC 05340', 'NGC 3631', 'UGCA 166', 'NGC 4314', 'NGC 4559', 'NGC 2681', 'NGC 5018', 'NGC 5253', 'NGC 4742', 'NGC 1672', 'NGC 4725', 'NGC 0891', 'NGC 6946', 'NGC 1291:[LFF2012] 084', 'NGC 3115', 'NGC 1332', 'NGC 1700', 'NGC 5584', 'NGC 7552', 'NGC 2997', 'NGC 4449', 'MESSIER 049', 'NGC 3198', 'NGC 0855', 'NGC 7793', 'NGC 0119', 'NGC 2865', 'MESSIER 059', 'NGC 1427', 'NGC 3628', 'NGC 0055', 'NGC 4457', 'NGC 4214', 'NGC 4459', 'NGC 3521'])

print Source_Overdensity_Master(['NGC 2841', 'NGC 3877', 'NGC 5054', 'NGC 5813', 'MESSIER 108', 'MESSIER 066', 'MESSIER 061', 'MESSIER 063', 'MESSIER 086', 'MESSIER 084', 'MESSIER 083', 'MESSIER 082', 'NGC 0278', 'MESSIER 088', 'NGC 3585', 'NGC 7507', 'NGC 1637', 'NGC 4473', 'NGC 1365', 'MESSIER 074', 'NGC 4570', 'NGC 4321', 'NGC 5474', 'NGC 7090', 'MESSIER 094', 'MESSIER 095', 'NGC 4494', 'IC 1613', 'NGC 4477', 'NGC 2787', 'IC 5267', 'NGC 3923', 'NGC 891', 'NGC 1300', 'UGC 05340', 'NGC 3631', 'NGC 4314', 'NGC 4559', 'NGC 2681', 'NGC 5018', 'NGC 5253', 'NGC 4742', 'NGC 1672', 'NGC 4725', 'NGC 0891', 'NGC 6946', 'NGC 1291:[LFF2012] 084', 'NGC 3115', 'NGC 1332', 'NGC 1700', 'NGC 5584', 'NGC 7552', 'NGC 2997', 'NGC 4449', 'MESSIER 049', 'NGC 3198', 'NGC 0855', 'NGC 7793', 'NGC 2865', 'MESSIER 059', 'NGC 1427', 'NGC 3628', 'NGC 0055', 'NGC 4457', 'NGC 4214', 'NGC 4459', 'NGC 3521']) # Real One
