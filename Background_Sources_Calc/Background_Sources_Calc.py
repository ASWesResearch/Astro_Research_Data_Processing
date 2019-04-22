import numpy as np
def Background_Source_Calc(Flux,Hardness_Str="S"):
    if(Hardness_Str=="S"):
        #Soft Equation
        N=370.0*((Flux/(2.0E-15))**(-0.85)) #In sources per deg^2
    if(Hardness_Str=="H"):
        #Hard Equation
        N=1200.0*((Flux/(2.0E-15))**(-1.0)) #In sources per deg^2
    return N

def Big_Background_Source_Calc(fpath,r_step=1.0):
    file=open(fpath,"r")
    data=file.read()
    #print data
    data_L=data.split("|")
    Background_Str=data_L[0]
    Background_Str_Reduced=Background_Str.split(":")[1]
    Background=float(Background_Str_Reduced)
    Fluxes_Str=data_L[1]
    #print "Fluxes_Str : ",Fluxes_Str
    Fluxes_L=Fluxes_Str.split("\n")
    #print "Fluxes_L : ", Fluxes_L
    #Fluxes_L.pop(0)
    #Fluxes_L.pop(len(Fluxes_L)-1)
    Fluxes_L_Reduced=[]
    for Flux_Str in Fluxes_L:
        if(Flux_Str==''):
            continue
        Flux=float(Flux_Str)
        Fluxes_L_Reduced.append(Flux)
    #print "Fluxes_L_Reduced : ", Fluxes_L_Reduced
    BG_Source_Density_L=[]
    for Flux in Fluxes_L_Reduced:
        Source_Density=Background_Source_Calc(Flux)
        BG_Source_Density_L.append(Source_Density)
    #print "BG_Source_Density_L : ", BG_Source_Density_L
    Cur_r=0.0
    r_step_deg=(r_step/60.0)
    BG_Source_L=[]
    r_List=[]
    Sum_Area_L=[]
    Sum_N_L=[]
    for BG_Source_Density in BG_Source_Density_L:
    #for i in range(0,len(BG_Source_Density_L)-2): # 2 instead of one because the number of BG soruces per bin needed not the number of BG sources per edge
    #for i in range(0,len(BG_Source_Density_L)-1): # 2 instead of one because the number of BG soruces per bin needed not the number of BG sources per edge
        #BG_Source_Density=BG_Source_Density_L[i]
        #Cur_r=Cur_r+r_step_deg
        r_List.append(Cur_r*60.0)
        Cur_Annulus_Area=np.pi*((2.0*r_step_deg*Cur_r)+(r_step_deg**2.0))
        #print "Cur_r : ", Cur_r
        #print "Cur_Annulus_Area : ", Cur_Annulus_Area
        Sum_Area=np.pi*(Cur_r**2.0)
        Cur_N=BG_Source_Density*Cur_Annulus_Area
        Cur_Sum_N=BG_Source_Density*Sum_Area
        #print "Sum_Area : ", Sum_Area
        Sum_Area_L.append(Sum_Area)
        #print "Cur_Sum_N : ", Cur_Sum_N
        Sum_N_L.append(Cur_Sum_N)
        BG_Source_L.append(Cur_N)
        Cur_r=Cur_r+r_step_deg
    Source_Density_HL=[Fluxes_L_Reduced,BG_Source_Density_L,r_List,BG_Source_L,Sum_Area_L,Sum_N_L]
    return Source_Density_HL

#print Big_Background_Source_Calc("/Volumes/xray/anthony/Research_Git/Data_Processing/Source_Overdensity/Background_Sources_Calc/NGC_253_ObsID_13830_Flux_90.txt")
#print Big_Background_Source_Calc("/Volumes/xray/anthony/Research_Git/Master_Code/Master_Output/NGC_3631/Flux_90_Files/3951/NGC_3631_ObsID_3951_Flux_90.txt")
