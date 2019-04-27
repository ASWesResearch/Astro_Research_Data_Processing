import astropy.io.ascii as ascii
import tarfile
import numpy as np
import math
import matplotlib.pyplot as plt
import os
from os import system
from astropy.io import fits
from astroquery.ned import Ned
import sys
import pandas as pd
import math
from astropy.io import ascii
import matplotlib.pyplot as plt
#path_GR=os.path.realpath('../') #Path needs to be updated
#path_GR=os.path.realpath('../../../') #Path needs to be updated
#"/Volumes/xray/anthony/Research_Git"
path_GR="/Volumes/xray/anthony/Research_Git"
#print "path_GR=",path_GR
sys.path.append(os.path.abspath(path_GR))
from Galaxy_Name_Reducer import Galaxy_Name_Reducer
from D25_Finder import D25_Finder
from File_Query_Code import File_Query_Code_5
path_BG_Source_Calc="/Volumes/xray/anthony/Research_Git/Data_Processing/Source_Overdensity/"
path_BG_Source_Calc_Reduced=os.path.realpath(path_BG_Source_Calc)
sys.path.append(os.path.abspath(path_BG_Source_Calc_Reduced))
from Background_Sources_Calc import Background_Sources_Calc
def Source_Histogram_Sum_Plot(Gname_L,Fname_Key="",Fileout_B=True,Outpath=False,Bins=100):
    """
    Gname: str- The name of the galaxy in the form NGC #
    Data: array- a table of data containg the coordinates of each object

    This fucntion takes the inputs for the name of the galaxy and returns a histrogram that plots the number of objects per bin by the
    area enclosed by the cirlce centered on the center of the galaxy that inculdes the objects in each bin in
    square degrees divided by the visible area of the galaxy in square degrees.
    This function plots the visible Major axis of the galaxy area enclosed by a circle that uses the Major axis as
    the diameter of the cirlce divided by itself to give 1 on histogram.
    This function uses astroquary to get data directly from NED

    #THIS IS THE CURRENT RUNNING VERSION OF THIS CODE
    """
    Failed_Galaxy_List=[]
    dis_Total_L=[]
    #BG_Source_D25_Total_L=[0,0,0,0,0,0,0,0,0,0]
    BG_Source_D25_Total_L=[0,0,0,0,0,0,0,0,0,0,0] #Cheap Fix to binning problem will fix later
    #BG_Source_D25_Total_A=np.zeros(10)
    for Gname in Gname_L:
        try:
            Gname_Modifed=Galaxy_Name_Reducer.Galaxy_Name_Reducer(Gname)
            dir = os.path.dirname(__file__)
            #Evt2_File_H_L=File_Query_Code_5.File_Query(Gname,"evt2")
            Evt2_File_H_L=File_Query_Code_5.File_Query(Gname,"evt2",Exp_Max_B=True)
            #print "Evt2_File_H_L : ", Evt2_File_H_L
            if(Evt2_File_H_L==False):
                print "Invalid Galaxy"
                return
            Galaxy_Obs_ID_L=[]
            for Evt2_File_L in Evt2_File_H_L:
                Cur_Galaxy_Obs_ID=Evt2_File_L[0]
                Galaxy_Obs_ID_L.append(Cur_Galaxy_Obs_ID)
            for Obs_ID in Galaxy_Obs_ID_L:
                G_Data= Ned.query_object(Gname) #G_Data:-astropy.table.table.Table, Galaxy_Data, The Galaxy Data Table queried from NED
                #print "G_Data : \n",G_Data
                #print type(G_Data)
                #/Volumes/xray/anthony/Research_Git/Master_Code/Master_Output/MESSIER_066/Coords_Lists/9548
                #MESSIER_066_ObsID_9548_Coords.csv
                Coords_Path="/Volumes/xray/anthony/Research_Git/Master_Code/Master_Output/"+Gname_Modifed+"/Coords_Lists/"+str(Obs_ID)+"/"+Gname_Modifed+"_ObsID_"+str(Obs_ID)+"_Coords.csv"
                Coords_Table=pd.read_csv(Coords_Path)
                #print "Coords_Table: \n",Coords_Table
                #Phys_X,Phys_Y,Chip_X,Chip_Y,Chip_ID,RA,DEC,Offaxis_Angle
                #Dec_Match_L,RA_Match_L
                RA_Match_Deg_DF=Coords_Table["RA"] #In decimal degrees
                RA_Match_Deg_L=list(RA_Match_Deg_DF)
                RA_Match_DF=RA_Match_Deg_DF*60.0 #In units of arcmins
                #print "RA_Match_DF :\n", RA_Match_DF
                RA_Match_L=list(RA_Match_DF)
                #print "RA_Match_L :\n", RA_Match_L
                Dec_Match_Deg_DF=Coords_Table["DEC"] #In decimal degrees
                Dec_Match_Deg_L=list(Dec_Match_Deg_DF)
                Dec_Match_DF=Dec_Match_Deg_DF*60.0 #In units of arcmins
                Dec_Match_L=list(Dec_Match_DF)
                #print "Dec_Match_L :\n", Dec_Match_L
                #return "STOPPING PROGRAM"
                #print "Area Matching_Index_List : ", Matching_Index_List
                D25_S_Maj_Deg=D25_Finder.D25_Finder(Gname)
                D25_S_Maj_Arcmin=D25_S_Maj_Deg*60.0
                area_T=((D25_S_Maj_Deg)**2)*math.pi
                try:
                    raGC=float(G_Data['RA(deg)'])
                    decGC=float(G_Data['DEC(deg)'])
                except:
                    raGC=float(G_Data['RA'])
                    decGC=float(G_Data['DEC'])
                #area_A=[((((((decGC-dec)**2)+((raGC-ra)**2)))*(math.pi))/area_T) for dec,ra in zip(decA,raA)]



                ##area_A=[((((((decGC-dec)**2)+((raGC-ra)**2)))*(math.pi))/area_T) for dec,ra in zip(Dec_Match_Deg_L,RA_Match_Deg_L)] #REAL ONE
                #disA=[math.sqrt(((decGC-dec)**2)+((raGC-ra)**2)) for dec,ra in zip(decA,raA)] #REAL ONE
                #print area_A
                #area_max=max(area_A)
                #print area_max
                #plt.vlines(0.00001,0,10,color='red')
                #plt.vlines(1,0,10,color='red')
                #plt.hist(area_A)
                ##Hist_A=plt.hist(area_A)
                ##Bin_Hight_A=Hist_A[0]
                ##Bin_Hight_Max=max(Bin_Hight_A)
                #print Bin_Hight_Max
                ##plt.vlines(1,0,Bin_Hight_Max,color='red')
                ##plt.xlabel('A')
                ##plt.ylabel('N')
                #Hist_Max=max(area_A)
                #print "Hist_Max ", Hist_Max
                Flux_90_Path="/Volumes/xray/anthony/Research_Git/Master_Code/Master_Output/"+Gname_Modifed+"/Flux_90_Files/"+str(Obs_ID)+"/"+Gname_Modifed+"_ObsID_"+str(Obs_ID)+"_Flux_90.txt"
                BG_Source_HL=Background_Sources_Calc.Big_Background_Source_Calc(Flux_90_Path)
                #print "BG_Source_HL : ", BG_Source_HL
                """
                BG_Source_Sum_Area_L=BG_Source_HL[4]
                BG_Source_Sum_Frac_L=[]
                for Sum_Area in BG_Source_Sum_Area_L:
                    Sum_Frac=Sum_Area/area_T
                    BG_Source_Sum_Frac_L.append(Sum_Frac)
                """
                #BG_Source_Sum_N_L=BG_Source_HL[5]
                BG_Source_L=BG_Source_HL[3]
                """
                BG_Source_L_Reduced=BG_Source_L[0:len(BG_Source_L)-1]
                #BG_Source_A_Reduced=np.array(BG_Source_L_Reduced)
                #np.delete(BG_Source_A, -1) #Deletes right most value so that bins are [0,1),[1,)...[n-1,n] were the value of any bin [i,i+1] is the value of the left edge of the bin (ie. the ith value, so if the number of expected BG source at an offaxis of 0' is 2 souces and at 1' 3 sources then the value for the bin [0,1) is 2 not 3. This may be reversed to be 3 later)
                """
                #print "BG_Source_Sum_Frac_L : ", BG_Source_Sum_Frac_L
                #print "BG_Source_L : ", BG_Source_L
                print "D25_S_Maj_Arcmin : ", D25_S_Maj_Arcmin
                D25_S_Maj_Arcmin_Whole=round(D25_S_Maj_Arcmin)
                print "D25_S_Maj_Arcmin_Whole : ", D25_S_Maj_Arcmin_Whole
                D25_S_Maj_Arcmin_Whole_Int=int(D25_S_Maj_Arcmin_Whole)
                print "D25_S_Maj_Arcmin_Whole_Int : ", D25_S_Maj_Arcmin_Whole_Int
                Num_BG_Bins=10.0/D25_S_Maj_Arcmin
                #print "Num_BG_Bins : ", Num_BG_Bins
                Num_BG_Bins_Whole=round(Num_BG_Bins)
                #print "Num_BG_Bins_Whole : ", Num_BG_Bins_Whole
                Num_BG_Bins_Whole_Int=int(Num_BG_Bins_Whole)
                #print "Num_BG_Bins_Whole_Int : ", Num_BG_Bins_Whole_Int
                Remainder_BG_Bins=10.0%Num_BG_Bins_Whole
                #print "Remainder_BG_Bins : ",Remainder_BG_Bins
                #Subdivided_A_Len=10*Num_BG_Bins_Whole_Int
                Subdivided_A_Len=11*Num_BG_Bins_Whole_Int
                #Subdivided_A_Len=10*D25_S_Maj_Arcmin_Whole_Int
                Subdivided_A=np.zeros(Subdivided_A_Len)
                BG_Source_SD_HL=[]
                #for BG_Source in BG_Source_L_Reduced:
                for BG_Source in BG_Source_L:
                    BG_Source_SD=BG_Source/Num_BG_Bins_Whole
                    #BG_Source_SD=BG_Source/Num_BG_Bins_Whole
                    Cur_BG_Source_SD_A=np.empty(Num_BG_Bins_Whole_Int)
                    Cur_BG_Source_SD_A.fill(BG_Source_SD)
                    Cur_BG_Source_SD_L=list(Cur_BG_Source_SD_A)
                    BG_Source_SD_HL.append(Cur_BG_Source_SD_L)
                #print "BG_Source_L : ", BG_Source_L
                #print "len(BG_Source_L) : ", len(BG_Source_L)
                #print "BG_Source_SD_HL : ", BG_Source_SD_HL
                #print "len(BG_Source_SD_HL) : ", len(BG_Source_SD_HL)
                i=0
                for BG_Source_SD_L in BG_Source_SD_HL:
                    for BG_Source in BG_Source_SD_L:
                        Subdivided_A[i]=BG_Source
                        i=i+1
                #print "Subdivided_A : ", Subdivided_A
                #print "len(Subdivided_A) : ", len(Subdivided_A)
                Subdivided_L=list(Subdivided_A)
                #print "Subdivided_L : ", Subdivided_L
                #[l[i:i + n] for i in xrange(0, len(l), n)]
                #Subdivided_HL=[Subdivided_L[j:j + Num_BG_Bins_Whole_Int] for j in xrange(0, len(Subdivided_L), Num_BG_Bins_Whole_Int)]
                #Subdivided_HL=[Subdivided_L[j:j + D25_S_Maj_Arcmin_Whole_Int] for j in xrange(0, len(Subdivided_L), D25_S_Maj_Arcmin_Whole_Int)]
                #Subdivided_HL=[Subdivided_L[j:j + 10] for j in xrange(0, len(Subdivided_L), 10)]
                Subdivided_HL=[Subdivided_L[j:j + 11] for j in xrange(0, len(Subdivided_L), 11)]
                #print "Subdivided_HL : ", Subdivided_HL
                #print "len(Subdivided_HL) : ", len(Subdivided_HL)
                Subdivided_D25_A=np.array(Subdivided_HL)
                D25_BG_Source_A=np.sum(Subdivided_D25_A, axis=1)
                #print "D25_BG_Source_A : ", D25_BG_Source_A
                D25_BG_Source_L=list(D25_BG_Source_A)
                #print "D25_BG_Source_L : ", D25_BG_Source_L
                #print "len(D25_BG_Source_L) : ", len(D25_BG_Source_L)
                #for k in range(0,len(D25_BG_Source_L)):
                #for k in range(1,len(D25_BG_Source_L)): # 1 instead of zero in order to get plt.step to work correctly
                #for k in range(0,len(D25_BG_Source_L)):
                for k in range(0,len(D25_BG_Source_L)):
                    #BG_Source_D25_Total_L[k]=BG_Source_D25_Total_L[k]+D25_BG_Source_L[k]
                    BG_Source_D25_Total_L[k+1]=BG_Source_D25_Total_L[k+1]+D25_BG_Source_L[k]
                """
                plt.step(BG_Source_Sum_Frac_L, BG_Source_L)
                plt.plot()
                #plt.savefig('Test1.png')
                #path_2=os.path.realpath('../Master_Code/Master_Output/') #Goes to Master_Output folder, which will hold all the data calculated for the galaxies including the histogram pictures
                path_2=os.path.realpath('./Source_Overdensity_Histograms/')
                if(Outpath!=False):
                    path_2=Outpath
                #print "Path_2=", path_2
                #Gname_Modifed=Galaxy_Name_Reducer.Galaxy_Name_Reducer(Gname)
                #print Gname_Modifed
                path_3=path_2+'/'+Gname_Modifed+'/'
                directory = os.path.dirname(path_3)
                if not os.path.exists(directory):
                    os.makedirs(directory)
                #os.chdir(path_3) #Goes to Current Galaxies Folder
                path_Hist=path_3+'Histograms/'
                directory_Hist=os.path.dirname(path_Hist)
                if not os.path.exists(directory_Hist):
                    os.makedirs(directory_Hist)
                print "path_Hist=",path_Hist
                #os.chdir(path_Hist)
                path_Obs_ID=path_Hist+str(Obs_ID)+'/'
                directory_Obs_ID=os.path.dirname(path_Obs_ID)
                if not os.path.exists(directory_Obs_ID):
                    os.makedirs(directory_Obs_ID)
                plt.savefig(path_Obs_ID+Gname_Modifed+'_'+str(Obs_ID)+'_Frac.png') #Saves angluar histogram figure
                #system('pwd')
                #path_4=os.path.realpath('../../../../Histogram_Code/')
                #print "Path_4=",path_4
                #os.chdir(path_4) #Goes back to where this code (the histogram code) is being run, ie. Desktop/GitHub
                plt.close()
                """
                #plt.show()
                #For Area Histrograms
                raGC_Arcmin=raGC*60.0
                decGC_Arcmin=decGC*60.0
                disA=[math.sqrt(((decGC_Arcmin-dec)**2)+((raGC_Arcmin-ra)**2)) for dec,ra in zip(Dec_Match_L,RA_Match_L)] #REAL ONE in units of arcmins
                #print "disA : ", disA
                #D25_S_Maj_Arcmin
                dis_D25_A=disA/D25_S_Maj_Arcmin
                #print "dis_D25_A : ", dis_D25_A
                dis_D25_L=list(dis_D25_A)
                for dis_D25 in dis_D25_L:
                    Dist_Arcmin=dis_D25*D25_S_Maj_Arcmin
                    if(Dist_Arcmin>10.0): #Throws out all sources outside the reasonable Chandra FOV (Offaxis > 10')
                        continue
                    dis_Total_L.append(dis_D25)
                #dis_max=np.max(disA)
                #N_Bins_Float=dis_max/2.0
                #N_Bins=int(round(N_Bins_Float))
                #print "dis_max : ", dis_max
                #print area_A
                #area_max=max(area_A)
                #print area_max
                #plt.vlines(0.00001,0,10,color='red')
                #plt.vlines(1,0,10,color='red')
                #plt.hist(area_A)
                #Hist_A=plt.hist(disA,[0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0])
                #Hist_A=plt.hist(disA,[0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
                #Hist_A=plt.hist(disA,[0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
                ##Hist_A=plt.hist(disA,range=(0.0,10.0))
        except:
            Failed_Galaxy_List.append(Gname)
            print "Failed_Galaxy_List : ", Failed_Galaxy_List
    dis_Total_A=np.array(dis_Total_L)
    #Hist_A=plt.hist(disA,range=(0.0,10.0/D25_S_Maj_Arcmin))
    #Hist_A=plt.hist(dis_Total_A,range=(0.0,20.0))
    Hist_A=plt.hist(dis_Total_A,bins=Bins,range=(0.0,10.0))
    Bin_Hight_A=Hist_A[0]
    Bin_A=Hist_A[1]
    print "Bin_A : ", Bin_A
    Bin_Hight_Max=max(Bin_Hight_A)
    #Hist_A=plt.hist(disA,N_Bins)
    #print "Angular Hist_A : ", Hist_A
    #print "Hist_A : ", Hist_A
    print "BG_Source_D25_Total_L : ", BG_Source_D25_Total_L
    BG_Source_D25_Total_A=np.array(BG_Source_D25_Total_L)
    BG_Source_D25_Total_A=BG_Source_D25_Total_A*(10.0/Bins)
    BG_Source_D25_Total_L=list(BG_Source_D25_Total_A)
    BG_Source_Sigificance_A=3.0*BG_Source_D25_Total_A
    """
    Bin_Hight_A=Hist_A[0]
    Bin_Hight_Max=max(Bin_Hight_A)
    #print "Bin_Hight_Max : ", Bin_Hight_Max
    #/Volumes/xray/anthony/Research_Git/Master_Code/Master_Output/NGC_3631/Flux_90_Files/3951/NGC_3631_ObsID_3951_Flux_90.txt #Example Path
    #Flux_90_Path="/Volumes/xray/anthony/Research_Git/Master_Code/Master_Output/"+Gname_Modifed+"/Flux_90_Files/"+str(Obs_ID)+"/"+Gname_Modifed+"_ObsID_"+str(Obs_ID)+"_Flux_90.txt"
    #BG_Source_HL=Background_Sources_Calc.Big_Background_Source_Calc(Flux_90_Path)
    #print "BG_Source_HL : ", BG_Source_HL
    BG_Bins=BG_Source_HL[2]
    BG_Source_L=BG_Source_HL[3]
    #Hist_BG_A=plt.hist(BG_Source_L,[0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0],alpha=0.5)
    #Bar_BG_A=plt.bar([0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0], BG_Source_L,color="red",alpha=0.5)
    #plt.plot([0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0], BG_Source_L)
    #plt.step([0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0], BG_Source_L)
    plt.step([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], BG_Source_L)
    ##plt.step([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]/D25_S_Maj_Arcmin, BG_Source_L)
    #BG_Bin_Hight=
    #print Bin_Hight_Max
    #plt.vlines(D25_S_Maj_Deg,0,Bin_Hight_Max,color='red') #Plots red line at D25
    plt.vlines(D25_S_Maj_Arcmin,0,Bin_Hight_Max,color='red') #Plots red line at D25
    """
    plt.step([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], BG_Source_D25_Total_L)
    plt.step([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], BG_Source_Sigificance_A,color="orange")
    plt.step([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], BG_Source_Sigificance_A*5.0,color="yellow")
    plt.vlines(1,0,Bin_Hight_Max,color='red') #Plots red line at D25
    #plt.xlabel('R (arcmin)')
    plt.xlabel('R (D25)')
    plt.ylabel('N')
    #print D25_S_Maj_Deg
    #Hist_Max=max(area_A)
    #print "Hist_Max ", Hist_Max
    plt.plot()
    #plt.savefig('Test2.png')
    #path_2=os.path.realpath('../Master_Code/Master_Output/') #Goes to Master_Output folder, which will hold all the data calculated for the galaxies including the histogram pictures, Quoted Out because path_2 is already defined
    #print "Path_2=", path_2
    """
    Gname_Modifed=Galaxy_Name_Reducer.Galaxy_Name_Reducer(Gname)
    #print Gname_Modifed
    path_3=path_2+'/'+Gname_Modifed+'/'
    directory = os.path.dirname(path_3)
    if not os.path.exists(directory):
        os.makedirs(directory)
    #os.chdir(path_3) #Goes to Current Galaxies Folder
    path_Hist=path_3+'Histograms/'
    directory_Hist=os.path.dirname(path_Hist)
    if not os.path.exists(directory_Hist):
        os.makedirs(directory_Hist)
    #print "path_Hist=",path_Hist
    #os.chdir(path_Hist)
    path_Obs_ID=path_Hist+str(Obs_ID)+'/'
    print "path_Obs_ID : ",path_Obs_ID
    directory_Obs_ID=os.path.dirname(path_Obs_ID)
    if not os.path.exists(directory_Obs_ID):
        os.makedirs(directory_Obs_ID)
    plt.savefig(path_Obs_ID + Gname_Modifed+'_'+str(Obs_ID)+'_Ang.png') #Saves angluar histogram figure
    """
    Path_Sum="/Volumes/xray/anthony/Research_Git/Data_Processing/Source_Overdensity/Source_Overdensity_Sum_Plot/Sum_Histograms/"
    """
    directory = os.path.dirname(path_3)
    if not os.path.exists(directory):
        os.makedirs(directory)
    """
    plt.savefig(Path_Sum+'Sum_D25_'+Fname_Key+'_Histogram.png') #Saves angluar histogram figure
    #BG_Source_A=np.array(BG_Source_L)
    #system('pwd')
    #path_4=os.path.realpath('../../../../Histogram_Code/')
    #print "Path_4=",path_4
    #os.chdir(path_4) #Goes back to where this code (the histogram code) is being run, ie. Desktop/GitHub
    plt.close()
    """
    if(Fileout_B):
        O_R_fname=path_Obs_ID + Gname_Modifed+'_'+str(Obs_ID)+'_Overdensity_Ratios.txt'
        file=open(O_R_fname,"w")
        print "O_R_fname:",O_R_fname
        D25_Line="D25:"+str(D25_S_Maj_Arcmin)+"|\n"
        file.write(D25_Line)
        BG_Source_A=np.array(BG_Source_L)
        Bin_Hight_L=list(Bin_Hight_A)
        #print "Bin_Hight_A : ", Bin_Hight_A
        #print "BG_Source_A: ", BG_Source_A
        #BG_Source_L_Reduced=BG_Source_L[0:len(BG_Source_L)-2]
        BG_Source_L_Reduced=BG_Source_L[0:len(BG_Source_L)-1]
        BG_Source_A_Reduced=np.array(BG_Source_L_Reduced)
        np.delete(BG_Source_A, -1) #Deletes right most value so that bins are [0,1),[1,)...[n-1,n] were the value of any bin [i,i+1] is the value of the left edge of the bin (ie. the ith value, so if the number of expected BG source at an offaxis of 0' is 2 souces and at 1' 3 sources then the value for the bin [0,1) is 2 not 3. This may be reversed to be 3 later)
        #print "BG_Source_A_Reduced: ", BG_Source_A_Reduced
        Overdensity_Ratio_A=(Bin_Hight_A/BG_Source_A_Reduced) #Somthing Might be wrong here
        Overdensity_Ratio_L=list(Overdensity_Ratio_A)
        #for i in range(0,len(Overdensity_Ratio_L)-1):
        for i in range(0,len(Overdensity_Ratio_L)):
            Num_Sources=Bin_Hight_L[i]
            Overdensity_Ratio=Overdensity_Ratio_L[i]
            Cur_Line=str(Overdensity_Ratio)+","+str(Num_Sources)+"\n"
            file.write(Cur_Line)
    """
    if(Fileout_B):
        #O_R_fname=path_Obs_ID + Gname_Modifed+'_'+str(Obs_ID)+'_Overdensity_Ratios.txt'
        #file=open(Path_Sum+'Sum_D25_'+Fname_Key+'_Data.txt',"w")
        file=open(Path_Sum+'Sum_D25_'+Fname_Key+'_Data.csv',"w")
        #print "O_R_fname:",O_R_fname
        #D25_Line="D25:"+str(D25_S_Maj_Arcmin)+"|\n"
        #file.write(D25_Line)
        #BG_Source_A=np.array(BG_Source_L)
        Bin_Hight_L=list(Bin_Hight_A)
        #print "Bin_Hight_A : ", Bin_Hight_A
        #print "BG_Source_A: ", BG_Source_A
        Bin_L=list(Bin_A)
        #print "Bin_L : ", Bin_L
        #BG_Source_L_Reduced=BG_Source_L[0:len(BG_Source_L)-2]
        BG_Source_L_Reduced=BG_Source_L[0:len(BG_Source_L)-1]
        #BG_Source_A_Reduced=np.array(BG_Source_L_Reduced)
        #np.delete(BG_Source_A, -1) #Deletes right most value so that bins are [0,1),[1,)...[n-1,n] were the value of any bin [i,i+1] is the value of the left edge of the bin (ie. the ith value, so if the number of expected BG source at an offaxis of 0' is 2 souces and at 1' 3 sources then the value for the bin [0,1) is 2 not 3. This may be reversed to be 3 later)
        #print "BG_Source_A_Reduced: ", BG_Source_A_Reduced
        #Overdensity_Ratio_A=(Bin_Hight_A/BG_Source_A_Reduced) #Somthing Might be wrong here
        #Overdensity_Ratio_L=list(Overdensity_Ratio_A)
        #for i in range(0,len(Overdensity_Ratio_L)-1):
        Header_Line_1="BG_Source_D25_Total_L:"+str(BG_Source_D25_Total_L)+"\n"+"Diluted:"+str((10.0/Bins))+"|\n"
        file.write(Header_Line_1)
        Header_Line_2="Bins,Num_Sources\n"
        file.write(Header_Line_2)
        for i in range(0,len(Bin_Hight_L)):
            Bin=Bin_L[i]
            Num_Sources=Bin_Hight_L[i]
            #BG_Source=BG_Source_L_Reduced[i]
            #Overdensity_Ratio=Overdensity_Ratio_L[i]
            #Cur_Line=str(Overdensity_Ratio)+","+str(Num_Sources)+"\n"
            #Cur_Line=str(Bin)+","+str(Num_Sources)+","+str(BG_Source)"\n"
            Cur_Line=str(Bin)+","+str(Num_Sources)+"\n"
            file.write(Cur_Line)
#Source_Histogram_Sum_Plot(["NGC 3631"],Fileout_B=True,Outpath=False)
#Source_Histogram_Sum_Plot(['NGC 4278'])
#Source_Histogram_Sum_Plot(['NGC 4278', 'NGC 2841', 'NGC 3877'])
#Source_Histogram_Sum_Plot(['NGC 4278', 'NGC 2841', 'NGC 3877', 'NGC 5194', 'NGC 5054', 'NGC 5813', 'MESSIER 108', 'MESSIER 066', 'MESSIER 061', 'MESSIER 063', 'MESSIER 086', 'MESSIER 084', 'MESSIER 083', 'MESSIER 082', 'NGC 0278', 'MESSIER 088', 'NGC 3585', 'NGC 7507', 'NGC 1637', 'NGC 4473', 'NGC 1365', 'MESSIER 074', 'NGC 4570', 'NGC 5576', 'NGC 4321', 'NGC 5474', 'NGC 7090', 'MESSIER 094', 'MESSIER 095', 'NGC 4494', 'IC 1613', 'NGC 4477', 'NGC 4365', 'NGC 2787', 'NGC 3557', 'IC 5267', 'NGC 4388', 'NGC 3923', 'NGC 891', 'NGC 1300', 'UGC 05340', 'NGC 3631', 'UGCA 166', 'NGC 4314', 'NGC 4550', 'Holmberg IX                   ', 'NGC 4559', 'NGC 1399', 'NGC 1316', 'NGC 1097', 'NGC 2681', 'NGC 5018', 'NGC 5253', 'NGC 4631', 'MESSIER 060', 'NGC 4742', 'NGC 1672', 'NGC 5846', 'NGC 4725', 'NGC 3507', 'MESSIER 087', 'NGC 0891', 'NGC 3384', 'NGC 6946', 'NGC 1291:[LFF2012] 084', 'NGC 3115', 'NGC 1332', 'NGC 1700', 'NGC 5584', 'NGC 7552', 'NGC 2997', 'NGC 4449', 'MESSIER 049', 'NGC 3198', 'NGC 0855', 'NGC 7793', 'NGC 0119', 'NGC 2865', 'MESSIER 059', 'NGC 1427', 'NGC 3628', 'NGC 3608', 'NGC 0055', 'NGC 4457', 'NGC 4214', 'NGC 4459', 'NGC 3521', 'NGC 4565', 'NGC 1313'])

#Tests:
#Source_Histogram_Sum_Plot(["NGC 3631"],"Test")
#Source_Histogram_Sum_Plot(['NGC 4278', 'NGC 2841', 'NGC 3877'],"Test")
#Source_Histogram_Sum_Plot(['NGC 2841', 'NGC 3877', 'NGC 5054', 'NGC 5813', 'MESSIER 108', 'MESSIER 066', 'MESSIER 061', 'MESSIER 063', 'MESSIER 086', 'MESSIER 084', 'MESSIER 083', 'MESSIER 082', 'NGC 0278', 'MESSIER 088', 'NGC 3585', 'NGC 7507', 'NGC 1637', 'NGC 4473', 'NGC 1365', 'MESSIER 074', 'NGC 4570', 'NGC 4321', 'NGC 5474', 'NGC 7090', 'MESSIER 094', 'MESSIER 095', 'NGC 4494', 'IC 1613', 'NGC 4477', 'NGC 2787', 'IC 5267', 'NGC 3923', 'NGC 891', 'NGC 1300', 'UGC 05340', 'NGC 3631', 'UGCA 166', 'NGC 4314', 'NGC 4559', 'NGC 2681', 'NGC 5018', 'NGC 5253', 'NGC 4742', 'NGC 1672', 'NGC 4725', 'NGC 0891', 'NGC 6946', 'NGC 1291:[LFF2012] 084', 'NGC 3115', 'NGC 1332', 'NGC 1700', 'NGC 5584', 'NGC 7552', 'NGC 2997', 'NGC 4449', 'MESSIER 049', 'NGC 3198', 'NGC 0855', 'NGC 7793', 'NGC 0119', 'NGC 2865', 'MESSIER 059', 'NGC 1427', 'NGC 3628', 'NGC 0055', 'NGC 4457', 'NGC 4214', 'NGC 4459', 'NGC 3521'],"Test_All") #All
#Source_Histogram_Sum_Plot(['NGC 2841', 'NGC 3877', 'NGC 5054', 'NGC 5813', 'MESSIER 108', 'MESSIER 066', 'MESSIER 061', 'MESSIER 063', 'MESSIER 086', 'MESSIER 084', 'MESSIER 083', 'MESSIER 082', 'NGC 0278', 'MESSIER 088', 'NGC 3585', 'NGC 7507', 'NGC 1637', 'NGC 4473', 'NGC 1365', 'MESSIER 074', 'NGC 4570', 'NGC 4321', 'NGC 5474', 'NGC 7090', 'MESSIER 094', 'MESSIER 095', 'NGC 4494', 'IC 1613', 'NGC 4477', 'NGC 2787', 'IC 5267', 'NGC 3923', 'NGC 891', 'NGC 1300', 'UGC 05340', 'NGC 3631', 'UGCA 166', 'NGC 4314', 'NGC 4559', 'NGC 2681', 'NGC 5018', 'NGC 5253', 'NGC 4742', 'NGC 1672', 'NGC 4725', 'NGC 0891', 'NGC 6946', 'NGC 1291:[LFF2012] 084', 'NGC 3115', 'NGC 1332', 'NGC 1700', 'NGC 5584', 'NGC 7552', 'NGC 2997', 'NGC 4449', 'MESSIER 049', 'NGC 3198', 'NGC 0855', 'NGC 7793', 'NGC 0119', 'NGC 2865', 'MESSIER 059', 'NGC 1427', 'NGC 3628', 'NGC 0055', 'NGC 4457', 'NGC 4214', 'NGC 4459', 'NGC 3521'],"Test_All_10_Bins",Bins=10) #All
#Source_Histogram_Sum_Plot(["NGC 3631"],"Test_10_Bins",Bins=10)
#Source_Histogram_Sum_Plot(['NGC 2841', 'NGC 3877', 'NGC 5054', 'NGC 5813', 'MESSIER 108', 'MESSIER 066', 'MESSIER 061', 'MESSIER 063', 'MESSIER 086', 'MESSIER 084', 'MESSIER 083', 'MESSIER 082', 'NGC 0278', 'MESSIER 088', 'NGC 3585', 'NGC 7507', 'NGC 1637', 'NGC 4473', 'NGC 1365', 'MESSIER 074', 'NGC 4570', 'NGC 4321', 'NGC 5474', 'NGC 7090', 'MESSIER 094', 'MESSIER 095', 'NGC 4494', 'IC 1613', 'NGC 4477', 'NGC 2787', 'IC 5267', 'NGC 3923', 'NGC 891', 'NGC 1300', 'UGC 05340', 'NGC 3631', 'UGCA 166', 'NGC 4314', 'NGC 4559', 'NGC 2681', 'NGC 5018', 'NGC 5253', 'NGC 4742', 'NGC 1672', 'NGC 4725', 'NGC 0891', 'NGC 6946', 'NGC 1291:[LFF2012] 084', 'NGC 3115', 'NGC 1332', 'NGC 1700', 'NGC 5584', 'NGC 7552', 'NGC 2997', 'NGC 4449', 'MESSIER 049', 'NGC 3198', 'NGC 0855', 'NGC 7793', 'NGC 0119', 'NGC 2865', 'MESSIER 059', 'NGC 1427', 'NGC 3628', 'NGC 0055', 'NGC 4457', 'NGC 4214', 'NGC 4459', 'NGC 3521'],"Test_All_100_Bins") #All

# Real Ones:
Source_Histogram_Sum_Plot(['NGC 2841', 'NGC 3877', 'NGC 5054', 'NGC 5813', 'MESSIER 108', 'MESSIER 066', 'MESSIER 061', 'MESSIER 063', 'MESSIER 086', 'MESSIER 084', 'MESSIER 083', 'MESSIER 082', 'NGC 0278', 'MESSIER 088', 'NGC 3585', 'NGC 7507', 'NGC 1637', 'NGC 4473', 'NGC 1365', 'MESSIER 074', 'NGC 4570', 'NGC 4321', 'NGC 5474', 'NGC 7090', 'MESSIER 094', 'MESSIER 095', 'NGC 4494', 'IC 1613', 'NGC 4477', 'NGC 2787', 'IC 5267', 'NGC 3923', 'NGC 891', 'NGC 1300', 'UGC 05340', 'NGC 3631', 'NGC 4314', 'NGC 4559', 'NGC 2681', 'NGC 5018', 'NGC 5253', 'NGC 4742', 'NGC 1672', 'NGC 4725', 'NGC 0891', 'NGC 6946', 'NGC 1291:[LFF2012] 084', 'NGC 3115', 'NGC 1332', 'NGC 1700', 'NGC 5584', 'NGC 7552', 'NGC 2997', 'NGC 4449', 'MESSIER 049', 'NGC 3198', 'NGC 0855', 'NGC 7793', 'NGC 2865', 'MESSIER 059', 'NGC 1427', 'NGC 3628', 'NGC 0055', 'NGC 4457', 'NGC 4214', 'NGC 4459', 'NGC 3521'],"All") #All
Source_Histogram_Sum_Plot(['NGC 2841', 'NGC 3877', 'NGC 5054', 'MESSIER 108', 'MESSIER 066', 'MESSIER 061', 'MESSIER 063', 'MESSIER 083', 'NGC 0278', 'MESSIER 088', 'NGC 1637', 'NGC 1365', 'MESSIER 074', 'NGC 4321', 'NGC 5474', 'NGC 7090', 'MESSIER 094', 'MESSIER 095', 'NGC 4477', 'NGC 2787', 'IC 5267', 'NGC 891', 'NGC 1300', 'NGC 3631', 'NGC 4314', 'NGC 4559', 'NGC 2681', 'NGC 5018', 'NGC 1672', 'NGC 4725', 'NGC 0891', 'NGC 6946', 'NGC 1291:[LFF2012] 084', 'NGC 3115', 'NGC 1332', 'NGC 5584', 'NGC 7552', 'NGC 2997', 'NGC 3198', 'NGC 0855', 'NGC 7793', 'NGC 3628', 'NGC 4457', 'NGC 4459', 'NGC 3521'],"Spirals") #Just Spirals
Source_Histogram_Sum_Plot(['NGC 5813', 'MESSIER 086', 'MESSIER 084', 'NGC 3585', 'NGC 7507', 'NGC 4473', 'NGC 4570', 'NGC 4494', 'NGC 3923', 'NGC 4742', 'NGC 1700', 'MESSIER 049', 'NGC 2865', 'MESSIER 059', 'NGC 1427'],"Ellipticals") #Just Ellipticals
Source_Histogram_Sum_Plot(['MESSIER 082', 'IC 1613', 'UGC 05340', 'NGC 5253', 'NGC 4449', 'NGC 0055', 'NGC 4214'],"Irregulars") #Just Irregulars
#Source_Histogram_Sum_Plot(['NGC 4278', 'NGC 5813', 'MESSIER 086', 'MESSIER 084', 'NGC 3585', 'NGC 7507', 'NGC 4473', 'NGC 4570', 'NGC 5576', 'NGC 4494', 'NGC 4365', 'NGC 3557', 'NGC 3923', 'NGC 4550', 'NGC 1399', 'NGC 1316', 'MESSIER 060', 'NGC 4742', 'NGC 5846', 'MESSIER 087', 'NGC 1700', 'MESSIER 049', 'NGC 2865', 'MESSIER 059', 'NGC 1427', 'NGC 3608'],"Ellipticals_BG_Galaxies") #Just Ellipticals with Background Galaxies
