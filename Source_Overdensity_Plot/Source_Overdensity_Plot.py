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
def Source_Histogram_Plot(Gname,Fileout_B=True,Outpath=False):
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
    Gname_Modifed=Galaxy_Name_Reducer.Galaxy_Name_Reducer(Gname)
    import math
    from astropy.io import ascii
    import matplotlib.pyplot as plt
    #system('pwd')
    #system('cd ~/Desktop/SQL_Standard_File/')
    #import os
    dir = os.path.dirname(__file__)
    #filename= os.path.join(dir, '~','Desktop','SQL_Standard_File',)
    #filepath=os.path.abspath("~/Desktop/SQL_Standard_File")
    #print "Filepath =",filepath
    #path= os.path.join(dir,'~','Desktop','SQL_Standard_File',)
    #path=os.path.realpath('~/Desktop/SQL_Standard_File/SQL_Sandard_File.csv')
    #path=os.path.realpath('../SQL_Standard_File/SQL_Standard_File.csv')
    #path=os.path.realpath('../SQL_Standard_File/Source_Flux_Table.csv') #To Do: Update this path
    #path=os.path.realpath('../../../SQL_Standard_File/Source_Flux_Table.csv') #To Do: Update this path
    #path="/Volumes/xray/anthony/Research_Git/SQL_Standard_File/Source_Flux_Table.csv" #THIS CSC FILE DOES NOT HAVE ALL THE SOURCES ! ! ! IT IS WRONG ! ! !
    #print "Path=",path
    #os.chdir(path)
    #os.chdir('~/Desktop/SQL_Standard_File/')
    #system('cd ~/Desktop/Big_Object_Regions/')
    #system('cd ../SQL_Standard_File/')
    #system('pwd')
    #system('ls')
    #data = ascii.read("SQL_Sandard_File.csv") #data:-astropy.table.table.Table, data, The data from the SQL_Standard_File
    #data = ascii.read(filename) #data:-astropy.table.table.Table, data, The data from the SQL_Standard_File
    #data = ascii.read(filepath) #data:-astropy.table.table.Table, data, The data from the SQL_Standard_File
    Evt2_File_H_L=File_Query_Code_5.File_Query(Gname,"evt2")
    #print "Evt2_File_H_L : ", Evt2_File_H_L
    if(Evt2_File_H_L==False):
        print "Invalid Galaxy"
        return
    Galaxy_Obs_ID_L=[]
    for Evt2_File_L in Evt2_File_H_L:
        Cur_Galaxy_Obs_ID=Evt2_File_L[0]
        Galaxy_Obs_ID_L.append(Cur_Galaxy_Obs_ID)
    """
    data = ascii.read(path) #data:-astropy.table.table.Table, data, The data from the SQL_Standard_File
    #data2=open("SQL_Sandard_File.csv","r")
    #print data2
    #system('cd ~/Desktop/galaxies/out')
    RA_A=data['sourceRA'] #RA_A:-astropy.table.column.Column, Right_Ascension_Array, The array containing all Right Ascensions in the SQL Standard File
    #print type(RA_A)
    RA_L=list(RA_A) #RA_L:-list, Right_Ascension_List, The list containing all Right Ascensions in the SQL Standard File
    #print RA_L
    RA_A_Arcmin=RA_A*60.0
    Dec_A=data['sourceDec'] #Dec_A:-astropy.table.column.Column, Declination_Array, The array containing all Declinations in the SQL Standard File
    Dec_L=list(Dec_A) #Dec_L:-List, Declination_List, The list containing all Declinations in the SQL Standard File
    Dec_A_Arcmin=Dec_A*60.0
    #print Dec_L
    #Obs_ID_A=data["obsid"] #Obs_ID_A:-astropy.table.column.Column, Observation_Idenification_Array, The array containing all Observation IDs in the SQL_Standard_File (not indexable)
    #print type(Obs_ID_A)
    #Obs_ID_L=list(Obs_ID_A) #Obs_ID_L:-List, Observation_Idenification_List, The list containing all Observation IDs in the SQL_Standard_File (So it is indexable)
    Obs_ID_A=data["OBSID"] #Obs_ID_A:-astropy.table.column.Column, Observation_Idenification_Array, The array containing all Observation IDs in the SQL_Standard_File (not indexable)
    #print type(Obs_ID_A)
    Obs_ID_L=list(Obs_ID_A) #Obs_ID_L:-List, Observation_Idenification_List, The list containing all Observation IDs in the SQL_Standard_File (So it is indexable)
    #print "Obs_ID_L ", Obs_ID_L
    """
    for Obs_ID in Galaxy_Obs_ID_L:
        #print type(Obs_ID_L)
        #print Obs_ID_A
        #FGname_A=data["foundName"]
        #FGname_L=list(FGname_A)
        #print FGname_A
        #QGname_A=data["queriedName"] #QGname_A:-Obs_ID_A:-astropy.table.column.Column, Query_Galaxy_Name_Array, The array containing all Query Galaxy Names in the SQL_Standard_File (not indexable)
        """
        QGname_A=data["resolvedObject"] #QGname_A:-Obs_ID_A:-astropy.table.column.Column, Query_Galaxy_Name_Array, The array containing all Query Galaxy Names in the SQL_Standard_File (not indexable)
        QGname_L=list(QGname_A) #QGname_L:-List, Query_Galaxy_Name_Array, The list containing all Query Galaxy Names in the SQL_Standard_File (So it is indexable)
        #print type(QGname_A)
        #print QGname_A
        Matching_Index_List=[] #Matching_Index_List:-List, Matching_Index_List, The list of all indexes (ref. QGname_L) that corresepond to the input Galaxy Name, All arrays are of equal lenth, and "ith" value of an array is the correseponding value for any other arrays "ith" value, so for example Obs_ID_L[228]=794 and the Galaxy in the Observation is QGname_L[228]="NGC 891", Note both lists have the same index
        for i in range(0,len(QGname_L)): # i:-int, i, the "ith" index of QGname_L
            #print "i ", i
            QGname=QGname_L[i] #QGname:-string, Query_Galaxy_Name, The current test Galaxy Name, if this Galaxy name equals the input Galaxy Name (Gname) then this Matching_Index, i (ref. QGname_L) will be appended to the Matching_Index_List
            #QGname_Reduced=QGname.replace(" ", "")
            #print "QGname ", QGname
            #print "QGname_Reduced ", QGname_Reduced
            if(Gname==QGname): #Checks to see if the current test Galaxy Name is the same as the input Galaxy Name, if so it appends the current index (ref. QGname_L) to the Matching_Index_List
                #print "i ", i
                Matching_Index_List.append(i) #Appends the current index (ref. QGname_L) to the Matching_Index_List
        """
        """
        Matching_Index_List=[] #Matching_Index_List:-List, Matching_Index_List, The list of all indexes (ref. QGname_L) that corresepond to the input Galaxy Name, All arrays are of equal lenth, and "ith" value of an array is the correseponding value for any other arrays "ith" value, so for example Obs_ID_L[228]=794 and the Galaxy in the Observation is QGname_L[228]="NGC 891", Note both lists have the same index
        #Matching_Obs_ID_List=[] #Matching_Obs_ID_List:-List, Matching_Obs_ID_List, The list of all ob
        for i in range(0,len(Obs_ID_L)): # i:-int, i, the "ith" index of QGname_L
            #print "i ", i
            QObs_ID=Obs_ID_L[i] #QGname:-string, Query_Galaxy_Name, The current test Galaxy Name, if this Galaxy name equals the input Galaxy Name (Gname) then this Matching_Index, i (ref. QGname_L) will be appended to the Matching_Index_List
            #QGname_Reduced=QGname.replace(" ", "")
            #print "QGname ", QGname
            #print "QGname_Reduced ", QGname_Reduced
            if(Obs_ID==QObs_ID): #Checks to see if the current test Galaxy Name is the same as the input Galaxy Name, if so it appends the current index (ref. QGname_L) to the Matching_Index_List
                #print "i ", i
                Matching_Index_List.append(i) #Appends the current index (ref. QGname_L) to the Matching_Index_List
        RA_Match_L=[] #RA_Match_L:-List, Right_Ascension_Match_List, The list of all source RA's for the input Galaxy Name in decimal degrees
        Dec_Match_L=[] #Dec_Match_L:-List, Declination_Match_List, The list of all source Dec's for the input Galaxy Name in decimal degrees
        for Cur_Matching_Index in Matching_Index_List: #Cur_Matching_Index:-int, Current_Matching_Index, The current index (ref. QGname_L) in the list of matching indexes for the current input Galaxy Name (Matching_Index_List)
            Cur_Match_RA=RA_L[Cur_Matching_Index] #Cur_Match_RA:-numpy.float64, Current_Match_Right_Ascension, The RA of the current source in decimal degrees
            #print type(Cur_Match_RA)
            Cur_Match_Dec=Dec_L[Cur_Matching_Index] #Cur_Match_Dec:-numpy.float64, Current_Match_Declination, The Dec of the current source in decimal degrees
            RA_Match_L.append(Cur_Match_RA) #RA_Match_L:-list, Right_Ascension_Match_List, The list of all source RA's for the input Galaxy Name in decimal degrees
            Dec_Match_L.append(Cur_Match_Dec) #Dec_Match_L:-list, Declination_Match_List, The list of all source Dec's for the input Galaxy Name in decimal degrees
        """

        #print RA_Match_L
        #print len(RA_Match_L)
        #print Dec_Match_L
        #print len(Dec_Match_L)
        #decA=Data['dec']
        #raA=Data['ra']
        #Maj=Maj/3600
        #S_Maj=Maj/2
        #area_T=((S_Maj)**2)*math.pi
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
        """
        Dia_Table = Ned.get_table(Gname, table='diameters') #Dia_Table:-astropy.table.table.Table, Diameter_Table, The Data table queried from NED that contains the infomation about the Major Axis of the input Galaxy Name
        #print type(Dia_Table)
        #print G_Data
        #print Dia_Table
        #print Dia_Table.colnames
        #print Dia_Table.meta
        #print Dia_Table.columns
        Dia_Table_Feq=Dia_Table['Frequency targeted'] #Dia_Table_Feq:-astropy.table.column.MaskedColumn, Diameter_Table_Fequency, The Array containing all named frequencies of light that are being used for the Major Axis Measurement
        #print Dia_Table['NED Frequency']
        #print Dia_Table_Feq
        #print type(Dia_Table_Feq)
        Dia_Table_Feq_L=list(Dia_Table_Feq) #Dia_Table_Feq_L:-List, Diameter_Table_Fequency_List, The list containing all named frequencies of light that are being used for the Major Axis Measurement
        #print Dia_Table_Feq_L
        Dia_Table_Num=Dia_Table['No.'] #Dia_Table_Num:-astropy.table.column.MaskedColumn, Diameter_Table_Number, The number Ned assigns to
        #print Dia_Table_Num
        #print type(Dia_Table_Num)
        Dia_Table_Num_L=list(Dia_Table_Num)
        #print Dia_Table_Num_L
        for i in range(0,len(Dia_Table_Feq_L)-1): #There is a bug here with index matching, The matched index isn't that same index for the major axis
            Cur_Feq=Dia_Table_Feq_L[i]
            #print Cur_Feq
            if(Cur_Feq=="RC3 D_25, R_25 (blue)"):
                Match_inx=i
                Match_Feq=Dia_Table_Feq_L[Match_inx]
                Match_Num=Dia_Table_Num_L[Match_inx]
                #Match_Num
                #print "Match_Feq ", Match_Feq
                #print "Match_inx ", Match_inx
                #print "Match_Num ", Match_Num
        #Dia_Table_Maj=Dia_Table['Major Axis']
        Dia_Table_Maj=Dia_Table['NED Major Axis']
        #print Dia_Table_Maj
        Dia_Table_Maj_L=list(Dia_Table_Maj)
        #print Dia_Table_Maj_L
        Dia_Table_Maj_Units=Dia_Table['Major Axis Unit']
        #print Dia_Table_Maj_Units
        Dia_Table_Maj_Units_L=list(Dia_Table_Maj_Units)
        #print Dia_Table_Maj_Units_L
        #print "i ", i
        D25_Maj=Dia_Table_Maj_L[Match_inx]
        #print "D25_Maj ", D25_Maj
        D25_Units=Dia_Table_Maj_Units[Match_inx]
        #print "D25_Units ", D25_Units
        #print type(Dia_Table)
        #print Dia_Table.info()
        #Dia_Table_2=Dia_Table[6]
        #print Dia_Table_2
        #Maj=Dia_Table_2[18]
        #print "Maj, ! ! !", Maj
        D25_S_Maj=D25_Maj/2.0
        D25_S_Maj_Deg=D25_S_Maj/3600.0
        """
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
        area_A=[((((((decGC-dec)**2)+((raGC-ra)**2)))*(math.pi))/area_T) for dec,ra in zip(Dec_Match_Deg_L,RA_Match_Deg_L)] #REAL ONE
        #disA=[math.sqrt(((decGC-dec)**2)+((raGC-ra)**2)) for dec,ra in zip(decA,raA)] #REAL ONE
        #print area_A
        #area_max=max(area_A)
        #print area_max
        #plt.vlines(0.00001,0,10,color='red')
        #plt.vlines(1,0,10,color='red')
        #plt.hist(area_A)
        Hist_A=plt.hist(area_A)
        Bin_Hight_A=Hist_A[0]
        Bin_Hight_Max=max(Bin_Hight_A)
        #print Bin_Hight_Max
        plt.vlines(1,0,Bin_Hight_Max,color='red')
        plt.xlabel('A')
        plt.ylabel('N')
        #Hist_Max=max(area_A)
        #print "Hist_Max ", Hist_Max
        Flux_90_Path="/Volumes/xray/anthony/Research_Git/Master_Code/Master_Output/"+Gname_Modifed+"/Flux_90_Files/"+str(Obs_ID)+"/"+Gname_Modifed+"_ObsID_"+str(Obs_ID)+"_Flux_90.txt"
        BG_Source_HL=Background_Sources_Calc.Big_Background_Source_Calc(Flux_90_Path)
        #print "BG_Source_HL : ", BG_Source_HL
        BG_Source_Sum_Area_L=BG_Source_HL[4]
        BG_Source_Sum_Frac_L=[]
        for Sum_Area in BG_Source_Sum_Area_L:
            Sum_Frac=Sum_Area/area_T
            BG_Source_Sum_Frac_L.append(Sum_Frac)
        #BG_Source_Sum_N_L=BG_Source_HL[5]
        BG_Source_L=BG_Source_HL[3]
        #print "BG_Source_Sum_Frac_L : ", BG_Source_Sum_Frac_L
        #print "BG_Source_L : ", BG_Source_L
        plt.step(BG_Source_Sum_Frac_L, BG_Source_L)
        plt.plot()
        #plt.savefig('Test1.png')
        #path_2=os.path.realpath('../Master_Code/Master_Output/') #Goes to Master_Output folder, which will hold all the data calculated for the galaxies including the histogram pictures
        path_2=os.path.realpath('./Source_Overdensity_Histograms/')
        if(Outpath!=False):
            path_2=Outpath
        #print "Path_2=", path_2
        """
        path_Hist=path_2+'/Histograms/'
        directory_Hist=os.path.dirname(path_Hist)
        if not os.path.exists(directory_Hist):
            os.makedirs(directory_Hist)
        print "path_Hist=",path_Hist
        os.chdir(path_Hist)
        """
        #system('mkdir '+Gname) #Creates Current Galaxy's Folder, Folder Named after Galaxy, Note: will have to remove space from "NGC #" to change to "NGC_#", I Don't know if this works
        """
        Gname_L=Gname.split(" ")
        print "Gname_L: ", Gname_L
        if(len(Gname_L)>1):
            Gname_Modifed=Gname_L[0]+"_"+Gname_L[1] #Adds underscore to remove space from "NGC #" to change to "NGC_#" if there is a space in the name
        else:
            Gname_Modifed=Gname # Does nothing if the galaxy name has no space, ie. NGC#, For example NGC253 instead of NGC 253 or NGC_253
        """
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
        #plt.show()
        #For Area Histrograms
        raGC_Arcmin=raGC*60.0
        decGC_Arcmin=decGC*60.0
        '''
        RA_A_Arcmin=RA_A*60.0
        RA_L_Arcmin=list(RA_A_Arcmin)
        Dec_A_Arcmin=Dec_A*60.0
        Dec_L_Arcmin=list(Dec_A_Arcmin)
        #print "Distance Matching_Index_List : ", Matching_Index_List
        RA_Match_L=[]
        Dec_Match_L=[]
        for Cur_Matching_Index in Matching_Index_List: #Cur_Matching_Index:-int, Current_Matching_Index, The current index (ref. QGname_L) in the list of matching indexes for the current input Galaxy Name (Matching_Index_List)
            """
            Cur_Match_RA=RA_L[Cur_Matching_Index] #Cur_Match_RA:-numpy.float64, Current_Match_Right_Ascension, The RA of the current source in decimal degrees
            #print type(Cur_Match_RA)
            Cur_Match_Dec=Dec_L[Cur_Matching_Index] #Cur_Match_Dec:-numpy.float64, Current_Match_Declination, The Dec of the current source in decimal degrees
            """
            Cur_Match_RA=RA_L_Arcmin[Cur_Matching_Index] #Cur_Match_RA:-numpy.float64, Current_Match_Right_Ascension, The RA of the current source in decimal degrees
            #print type(Cur_Match_RA)
            Cur_Match_Dec=Dec_L_Arcmin[Cur_Matching_Index] #Cur_Match_Dec:-numpy.float64, Current_Match_Declination, The Dec of the current source in decimal degrees
            RA_Match_L.append(Cur_Match_RA) #RA_Match_L:-list, Right_Ascension_Match_List, The list of all source RA's for the input Galaxy Name in decimal degrees
            Dec_Match_L.append(Cur_Match_Dec) #Dec_Match_L:-list, Declination_Match_List, The list of all source Dec's for the input Galaxy Name in decimal degrees
            #print "Dec_Match_L : ", Dec_Match_L
            #print "decGC_Arcmin : ", decGC_Arcmin
        '''
        disA=[math.sqrt(((decGC_Arcmin-dec)**2)+((raGC_Arcmin-ra)**2)) for dec,ra in zip(Dec_Match_L,RA_Match_L)] #REAL ONE in units of arcmins
        #print "disA : ", disA
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
        Hist_A=plt.hist(disA,range=(0.0,10.0))
        #Hist_A=plt.hist(disA,N_Bins)
        #print "Angular Hist_A : ", Hist_A
        #print "Hist_A : ", Hist_A
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
        BG_Source_Plot_L=[0.0]
        for i in range(0,len(BG_Source_L)-1):
            BG_Source_Plot_L.append(BG_Source_L[i])
        print "BG_Source_L : ", BG_Source_L
        print "BG_Source_Plot_L : ", BG_Source_Plot_L
        #BG_Source_D25_Total_L=list(BG_Source_D25_Total_A)
        BG_Source_Plot_A=np.array(BG_Source_Plot_L)
        BG_Source_Sigificance_Plot_A=3.0*BG_Source_Plot_A
        #plt.step([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], BG_Source_L)
        plt.step([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], BG_Source_Plot_L)
        plt.step([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], BG_Source_Sigificance_Plot_A,color="orange")
        plt.step([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], BG_Source_Sigificance_Plot_A*3.0,color="gold")
        #plt.step([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], BG_Source_Sigificance_Plot_A*3.0,color="goldenrod")
        #plt.step([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], BG_Source_Sigificance_Plot_A*5.0,color="gold")

        #BG_Bin_Hight=
        #print Bin_Hight_Max
        #plt.vlines(D25_S_Maj_Deg,0,Bin_Hight_Max,color='red') #Plots red line at D25
        plt.vlines(D25_S_Maj_Arcmin,0,Bin_Hight_Max,color='red') #Plots red line at D25
        plt.xlabel('R (arcmin)')
        plt.ylabel('N')
        #print D25_S_Maj_Deg
        #Hist_Max=max(area_A)
        #print "Hist_Max ", Hist_Max
        plt.plot()
        #plt.savefig('Test2.png')
        #path_2=os.path.realpath('../Master_Code/Master_Output/') #Goes to Master_Output folder, which will hold all the data calculated for the galaxies including the histogram pictures, Quoted Out because path_2 is already defined
        #print "Path_2=", path_2
        """
        path_Hist=path_2+'/Histograms/'
        directory_Hist=os.path.dirname(path_Hist)
        if not os.path.exists(directory_Hist):
            os.makedirs(directory_Hist)
        print "path_Hist=",path_Hist
        os.chdir(path_Hist)
        """
        #system('mkdir '+Gname) #Creates Current Galaxy's Folder, Folder Named after Galaxy, Note: will have to remove space from "NGC #" to change to "NGC_#", I Don't know if this works
        """
        Gname_L=Gname.split(" ")
        print "Gname_L: ", Gname_L
        if(len(Gname_L)>1):
            Gname_Modifed=Gname_L[0]+"_"+Gname_L[1] #Adds underscore to remove space from "NGC #" to change to "NGC_#" if there is a space in the name
        else:
            Gname_Modifed=Gname # Does nothing if the galaxy name has no space, ie. NGC#, For example NGC253 instead of NGC 253 or NGC_253
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
        #BG_Source_A=np.array(BG_Source_L)
        #system('pwd')
        #path_4=os.path.realpath('../../../../Histogram_Code/')
        #print "Path_4=",path_4
        #os.chdir(path_4) #Goes back to where this code (the histogram code) is being run, ie. Desktop/GitHub
        plt.close()
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
        #print "Done"
        #plt.show()

#Area_GC_R_N_F_2('NGC4258')
#Area_GC_R_N('NGC4258')
#plt.savefig(Gname_Modifed+'_Ang.png') #Saves angluar histogram figure
"""
def Driver_Code(Gname):
    Area_GC_R_N_F_2(Gname)
    #Area_GC_R_N(Gname)
"""
#Driver_Code('NGC4258')
#Driver_Code('NGC 4649')
#"NGC 891"
#Driver_Code('NGC 891')
#Driver_Code('NGC_253')
#Driver_Code('NGC_3631')
#Driver_Code('NGC_5813')
#Driver_Code('MESSIER 083')
#Driver_Code('MESSIER 084')
#Driver_Code('NGC_6946')
#NGC_1365
#Source_Histogram_Plot('NGC_1365')
#Source_Histogram_Plot('NGC_3631')
#Source_Histogram_Plot('MESSIER_066')
#Outpath="/Volumes/xray/anthony/Research_Git/Data_Processing/Source_Overdensity/Source_Overdensity_Plot/Source_Overdensity_Histograms_Test/"
#Source_Histogram_Plot('NGC_3631',Outpath="/Volumes/xray/anthony/Research_Git/Data_Processing/Source_Overdensity/Source_Overdensity_Plot/Source_Overdensity_Histograms_Test/")
#Source_Histogram_Plot('NGC_3628',Outpath="/Volumes/xray/anthony/Research_Git/Data_Processing/Source_Overdensity/Source_Overdensity_Plot/Source_Overdensity_Histograms_Test/")
