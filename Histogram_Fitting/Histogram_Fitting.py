import numpy as np
import pandas as pd
import sys
import ast
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import chisquare
if sys.version_info[0] < 3:
    from StringIO import StringIO
else:
    from io import StringIO
df = pd.read_csv(StringIO("csv string..."))
def Exponential(x,A,C):
    return A*np.exp(C*x)

def Histogram_Fitting(path,Label="",Normalize_B=False,Color_L=["navy","cyan"],Points=False):
    #path="/Volumes/xray/anthony/Research_Git/Data_Processing/Source_Overdensity/Source_Overdensity_Sum_Plot/Sum_Histograms/"
    file = open(path,"r")
    Data_Str=file.read()
    #print "Data_Str : \n", Data_Str
    Data_Str_L=Data_Str.split("|")
    Data_Str_Header=Data_Str_L[0]
    #print "Data_Str_Header : ", Data_Str_Header
    Data_Str_Header_L=Data_Str_Header.split("\n")
    #print "Data_Str_Header_L : ", Data_Str_Header_L
    BG_Source_D25_Total_L_Str=Data_Str_Header_L[0].split(":")[1]
    #print "BG_Source_D25_Total_L_Str : ", BG_Source_D25_Total_L_Str
    #ast.literal_eval('["A","B" ,"C" ," D"]')
    BG_Source_D25_Total_L=ast.literal_eval(BG_Source_D25_Total_L_Str)
    #BG_Source_D25_Total_L=list(BG_Source_D25_Total_L_Str)
    #print "BG_Source_D25_Total_L : ", BG_Source_D25_Total_L
    #print "type(BG_Source_D25_Total_L) : ", type(BG_Source_D25_Total_L)
    BG_Source_D25_Total_A=np.array(BG_Source_D25_Total_L)
    BG_Source_Sigificance_A=3.0*BG_Source_D25_Total_A
    BG_Source_Bin_Dilution_Str=Data_Str_Header_L[1].split(":")[1]
    #print "BG_Source_Bin_Dilution_Str : ", BG_Source_Bin_Dilution_Str
    BG_Source_Bin_Dilution=float(BG_Source_Bin_Dilution_Str)
    #print "BG_Source_Bin_Dilution : ", BG_Source_Bin_Dilution
    Data_Str_Reduced=Data_Str_L[1]
    #print "Data_Str_Reduced : ", Data_Str_Reduced
    Data=pd.read_csv(StringIO(Data_Str_Reduced))
    #print "Data :\n",Data
    Bins_DF=Data["Bins"]
    #print "Bins_DF : ", Bins_DF
    Num_Sources_DF=Data["Num_Sources"]
    #print "Num_Sources_DF : ", Num_Sources_DF
    Bin_Hight_Max=np.max(Num_Sources_DF)
    #"""
    #plt.plot(Bins_DF,Num_Sources_DF)
    #plt.step([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], BG_Source_D25_Total_L)
    #plt.step([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], BG_Source_Sigificance_A,color="orange",label=Label)
    """
    if(Normalize_B==False):
        plt.step([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], BG_Source_Sigificance_A,label=Label+"_BG",color=Color_L[1])
    else:
        #BG_Source_Sigificance_Max=np.max(BG_Source_Sigificance_A)
        #BG_Source_Sigificance_A_Normalized=BG_Source_Sigificance_A/BG_Source_Sigificance_Max
        BG_Source_Sigificance_A_Normalized=BG_Source_Sigificance_A/Bin_Hight_Max
        plt.step([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], BG_Source_Sigificance_A_Normalized,label=Label+"_BG",color=Color_L[1])
    """
    #plt.step([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], BG_Source_Sigificance_A*2.0,color="grey")
    Bin_Hight_Max=np.max(Num_Sources_DF)
    if(Normalize_B==False):
        plt.vlines(1,0,Bin_Hight_Max,color='red') #Plots red line at D25
    else:
        plt.vlines(1,0,1,color='red') #Plots red line at D25
    #plt.show()
    #"""
    popt, pcov = curve_fit(Exponential, Bins_DF, Num_Sources_DF)
    print "popt : \n", popt
    #print "pcov : \n", pcov
    if(Normalize_B==False):
        plt.plot(Bins_DF, Exponential(Bins_DF, *popt),label=Label,color=Color_L[0])
        plt.ylabel('N')
    else:
        plt.plot(Bins_DF, Exponential(Bins_DF, *(1.0,popt[1])),label=Label,color=Color_L[0])
        plt.ylabel('N (Normalized)')
    if(Points):
        plt.plot(Bins_DF,Num_Sources_DF,".",label=Label+"_Data")
    if(Normalize_B==False):
        plt.step([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], BG_Source_Sigificance_A,label=Label+"_BG",color=Color_L[1])
    else:
        #BG_Source_Sigificance_Max=np.max(BG_Source_Sigificance_A)
        #BG_Source_Sigificance_A_Normalized=BG_Source_Sigificance_A/BG_Source_Sigificance_Max
        BG_Source_Sigificance_A_Normalized=BG_Source_Sigificance_A/Bin_Hight_Max
        plt.step([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], BG_Source_Sigificance_A_Normalized,label=Label+"_BG",color=Color_L[1])
    plt.xlabel('R (D25)')
    #plt.ylabel('N')
    plt.legend(loc='best')
    #plt.show()
    """
    Nexp = f(t, *popt)
    r = N - Nexp
    chisq = np.sum((r/stdev)**2)
    """
    """
    perr = np.sqrt(np.diag(pcov))
    stdev=np.std(Num_Sources_DF)
    print "stdev : ", stdev
    Nexp = Exponential(Bins_DF, *popt)
    #print "Nexp : \n", Nexp
    r = Num_Sources_DF - Nexp
    #print "r : ", r
    chisq = np.sum((r/stdev)**2)
    print "chisq : ", chisq
    """
    Nexp = Exponential(Bins_DF, *popt)
    chisq_T=chisquare( Num_Sources_DF, f_exp=Nexp)
    print "chisq_T : ", chisq_T
    chisq=chisq_T[0]
    print "chisq : ", chisq
    chisq_Reduced=chisq/(len(Num_Sources_DF)-len(popt))
    print "chisq_Reduced : ", chisq_Reduced
    perr = np.sqrt(np.diag(pcov))
    return popt,perr,chisq,chisq_Reduced
#Histogram_Fitting("/Volumes/xray/anthony/Research_Git/Data_Processing/Source_Overdensity/Source_Overdensity_Sum_Plot/Sum_Histograms/Sum_D25_All_Data.csv","All")

def Histogram_Fitting_Master():
    #plt.figure()
    Histogram_Fitting("/Volumes/xray/anthony/Research_Git/Data_Processing/Source_Overdensity/Source_Overdensity_Sum_Plot/Sum_Histograms/Sum_D25_All_Data.csv","All",Points=True)
    #plt.savefig("All.png")
    plt.close()
    Histogram_Fitting("/Volumes/xray/anthony/Research_Git/Data_Processing/Source_Overdensity/Source_Overdensity_Sum_Plot/Sum_Histograms/Sum_D25_Spirals_Data.csv","Spirals",Points=True)
    #plt.savefig("Spirals.png")
    plt.close()
    Histogram_Fitting("/Volumes/xray/anthony/Research_Git/Data_Processing/Source_Overdensity/Source_Overdensity_Sum_Plot/Sum_Histograms/Sum_D25_Ellipticals_Data.csv","Ellipticals",Points=True)
    #plt.savefig("Ellipticals.png")
    """
    plt.close()
    Histogram_Fitting("/Volumes/xray/anthony/Research_Git/Data_Processing/Source_Overdensity/Source_Overdensity_Sum_Plot/Sum_Histograms/Sum_D25_Spirals_Data.csv","Spirals")
    Histogram_Fitting("/Volumes/xray/anthony/Research_Git/Data_Processing/Source_Overdensity/Source_Overdensity_Sum_Plot/Sum_Histograms/Sum_D25_Ellipticals_Data.csv","Ellipticals",Color_L=["green","lime"])
    #plt.savefig("Spirals_vs_Ellipticals.png")
    plt.close()
    Histogram_Fitting("/Volumes/xray/anthony/Research_Git/Data_Processing/Source_Overdensity/Source_Overdensity_Sum_Plot/Sum_Histograms/Sum_D25_Spirals_Data.csv","Spirals",Normalize_B=True)
    Histogram_Fitting("/Volumes/xray/anthony/Research_Git/Data_Processing/Source_Overdensity/Source_Overdensity_Sum_Plot/Sum_Histograms/Sum_D25_Ellipticals_Data.csv","Ellipticals",Normalize_B=True,Color_L=["green","lime"])
    #plt.savefig("Spirals_vs_Ellipticals_Normalized.png")
    plt.close()
    #plt.show()
    """

Histogram_Fitting_Master()
