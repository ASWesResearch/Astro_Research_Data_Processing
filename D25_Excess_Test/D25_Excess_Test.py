import math
import os
from os import system
import sys
path_Modules="/Volumes/xray/anthony/Research_Git"
sys.path.append(os.path.abspath(path_Modules))
from Galaxy_Name_Reducer import Galaxy_Name_Reducer
from File_Query_Code import File_Query_Code_5
def D25_Excess_Test(Gname,Overdensity_Ratio_Threshold=3.0):
    Gname_Modifed=Galaxy_Name_Reducer.Galaxy_Name_Reducer(Gname)
    #print "D25 Excess Gname: ",Gname
    Evt2_File_H_L_Max_Exp=File_Query_Code_5.File_Query(Gname,"evt2",Exp_Max_B=True)
    print "Evt2_File_H_L_Max_Exp : ", Evt2_File_H_L_Max_Exp   
    Evt2_File_L_Max_Exp=Evt2_File_H_L_Max_Exp[0]
    Obs_ID_Max=Evt2_File_L_Max_Exp[0]
    #/Volumes/xray/anthony/Research_Git/Data_Processing/Source_Overdensity/Source_Overdensity_Plot/Source_Overdensity_Histograms/NGC_5813/Histograms/12952
    Overdensity_Ratio_Data_Fpath="/Volumes/xray/anthony/Research_Git/Data_Processing/Source_Overdensity/Source_Overdensity_Plot/Source_Overdensity_Histograms/"+Gname_Modifed+"/Histograms/"+str(Obs_ID_Max)+"/"+Gname_Modifed+"_"+str(Obs_ID_Max)+"_Overdensity_Ratios.txt"
    file=open(Overdensity_Ratio_Data_Fpath,"r")
    data=file.read()
    #print data
    data_L=data.split("|")
    D25_Str=data_L[0]
    D25_Str_Reduced=D25_Str.split(":")[1]
    D25=float(D25_Str_Reduced)
    #print "D25 : ", D25
    Histogram_Data_Str=data_L[1]
    Histogram_Data_Str_L=Histogram_Data_Str.split("\n")
    #print "Histogram_Data_Str_L Before: ", Histogram_Data_Str_L
    Histogram_Data_Str_L.pop(0)
    Histogram_Data_Str_L.pop(len(Histogram_Data_Str_L)-1)
    #print "Histogram_Data_Str_L After: ", Histogram_Data_Str_L
    Overdensity_Ratio_L=[]
    Num_Sources_L=[]
    for Histogram_Data_Str in Histogram_Data_Str_L:
        Cur_Histogram_Data_L=Histogram_Data_Str.split(",")
        #print "Cur_Histogram_Data_L : ", Cur_Histogram_Data_L
        Cur_Overdensity_Ratio_Str=Cur_Histogram_Data_L[0]
        #print "Cur_Overdensity_Ratio_Str : ", Cur_Overdensity_Ratio_Str
        Cur_Overdensity_Ratio=float(Cur_Overdensity_Ratio_Str)
        Overdensity_Ratio_L.append(Cur_Overdensity_Ratio)
        Cur_Num_Sources_Str=Cur_Histogram_Data_L[1]
        Cur_Num_Sources=float(Cur_Num_Sources_Str)
        Num_Sources_L.append(Cur_Num_Sources)
    D25_Rounded_Up=math.ceil(D25)
    #print "D25_Rounded_Up : ", D25_Rounded_Up
    i=int(D25_Rounded_Up) #Only checks bins outside D25 for source overdensities
    if(i>8):
        return False
    while(i<len(Overdensity_Ratio_L)-1):
        Overdensity_Ratio=Overdensity_Ratio_L[i]
        if(Overdensity_Ratio>Overdensity_Ratio_Threshold):
            return True
        i=i+1
    return False
#print D25_Excess_Test("NGC_5813")
#print D25_Excess_Test("NGC_5813",Overdensity_Ratio_Threshold=20.0)
#print D25_Excess_Test("NGC_253")
