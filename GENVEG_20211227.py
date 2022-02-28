#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 27 09:00:43 2021

@author: rdel1cmc
"""

  #############################################################################################################################
  #Model:		  Python Veg Model 
  #Purpose:	Simulates vegetation growth and mortality
  #Last Modified by:	Carrillo
  #Date: 		27 Dec 2021
  #############################################################################################################################
  
  #############################################################################################################################
  #IMPORT DATA and initiate variable values
  #############################################################################################################################
  

import os
import pandas as pd #abbreviates pandas command to pd...shortens code
import numpy as np
import sqlite3
import sys
import math
import array as arr
import csv


os.chdir(os.path.dirname(__file__))

Species = "corn"
Light = pd.read_csv("Model_Data/Light.csv", delimiter = ',')
Plants = pd.read_csv("Model_Data/Plants.csv", delimiter = ',')

# pull species data with corresponding parameter list and returns a dictionary for Species Type
def setPlantParameters(Plants_PD_, Species_Type_):
  input_dict = {}
  for i in range(len(Plants_PD_[["parm"]].values.tolist())):
    input_dict[Plants_PD_[["parm"]].values.tolist()[i][0]] = Plants_PD_[[Species_Type_]].values.tolist()[i][0]
  return input_dict

PlantParameters = setPlantParameters(Plants, Species)
print(PlantParameters)



print('pMax = ' + str(PlantParameters['pMax']))
print('pMax from R Code = 0.0372')

beginGrow = (PlantParameters['begingrow'])
print(beginGrow)
tillerDensity = 3
maxTillerHeight = 80
maxTillerWeight = 6
rWeightTillerHeight = 13.333
maxRootLength = 60
maxDensity = 2322
maxBiomass = 2322
minSize = 0
RSratio = 0.95
FracDM_LVG = 0.5   #amount of carbon that goes into the leaves
FracDM_STG = 0.2   #amount of carbon that goes into the stems
FracDM_RTG = 0.3   #amount of carbon that goes into the roots


kmLVG_prime = 0.03  #set from Teh 2006 Table 7.1 Leaves base maintenance respiration coefficient at 25 °C
kmSTG_prime = 0.015 #set from Teh 2006 Table 7.1 Stems base maintenance respiration coefficient at 25 °C
kmRTG_prime = 0.015 #set from Teh 2006 Table 7.1 Roots base maintenance respiration coefficient at 25 °C
glucoseReqLVG = 1.436 #set from Teh 2006 Table 7.4
glucoseReqSTG = 1.513 #set from Teh 2006 Table 7.4
glucoseReqRTG = 1.444 #set from Teh 2006 Table 7.4
TCA = 0.231           #total cross sectional area


##################################################

#Set length of simulation (Day of Year)
startsim = 1   #start of simulation
endsim = 365   #end of simulation
day = pd.Series(list(range(startsim, endsim+1)))  #Length of simulation, creates a vector from 1 to 365 given values of startsim and endsim defined

dt_balance = 1    #timestep
endGrow = 320   #date that vegetation stops growing due to day of first freeze (changes based on latitude)
leafDieOff = 163  #date that plant starts to self-shade
##################################################


##################################################
#Initialize Variables

#Initialize veg parameter variables and vectors
twlvd = 0
twstd = 0
twrtd = 0
twlvg = 0
twstg = 0
twrtg = 0
fgross = pd.Series([0,0,0])
print(fgross)
dtga = pd.Series([0,0,0])
print(dtga)
wgaus = pd.Series([0.2778, 0.4444, 0.2778])
print(wgaus)
xgauss = pd.Series([0.1127, 0.5, 0.8873])
print(xgauss)
Hi = 40
risk = 1

#Initalizes Vectors for daily outputs of critical variables 

dailyleafg = pd.Series(range(startsim,endsim+1))
dailyrootg = pd.Series(range(startsim,endsim+1))
dailystemg = pd.Series(range(startsim,endsim+1))
dailyleafd = pd.Series(range(startsim,endsim+1))
dailyrootd = pd.Series(range(startsim,endsim+1))
dailystemd = pd.Series(range(startsim,endsim+1))
dailyAGbiomass = pd.Series(range(startsim,endsim+1))
dailyrespiration = pd.Series(range(startsim,endsim+1))
dailyplantage = pd.Series(range(startsim,endsim+1))
dailyrespMaint = pd.Series(range(startsim,endsim+1))
dailyphoto = pd.Series(range(startsim,endsim+1))
dailyTOTbiomass = pd.Series(range(startsim,endsim+1))
dailyRTbiomass = pd.Series(range(startsim,endsim+1))

##################################################
#FUNCTIONS CALLED IN MODEL!!!!!!
##################################################

#FUNCTION: logistic shape
#Creates non-linear response curves for effect
#INPUTS: XY Paird for curve shape, input value for curve
#Output: value of effect, given specific curve and input


def logistic(X1, Y1, X2, Y2, x):
    G = math.log((Y1)/(1-Y1))
    C = math.log((Y2)/(1-Y2))
    B = ((G-C)/(X1-X2))
    A = (G-(B*X1))
    
    Z = math.exp(A+(B*x))

    S = (Z/(1+Z))
    
    return S

print(logistic(1, 0.5, 0, 0.5, 0.25))

#FUNCTION: PAR
#Calculates daylength, amount of photosynthetically availble light, Determines amount of plant growth per day
#INPUTS: Day of year and latitude
#Output: PAR vlues at 3 different tiems of dat for Gaussian Integration

#twlvg = 10
#k = -0.02
#pMax = 0.25
#day = 258
#lat = 30

def PAR(day, lat):
    
    degree_to_rad = 0.017453292  #required to convert degrees to radians
    #rad_to_degree = (-math.asin((math.sin(23.45*degree_to_rad))*(math.cos(2*math.pi*(day+10)/365))))
  
  
  
    tmpvec = []
    declination = (-math.asin((math.sin(23.45*degree_to_rad))*(math.cos(2*math.pi*(day+10)/365))))
    
    #Intermediate variables
    sinld = ((math.sin(lat*degree_to_rad))*(math.sin(declination)))   #radians
    cosld = math.cos(lat*degree_to_rad)*math.cos(declination)  #radians
    
    aob = (sinld/cosld)  #radians
    
    temp1 = math.asin(aob)
    
    daylength = 12 * (1 + 2 * temp1/math.pi)   #calculates daylength based on declination and latitude
    
    dsinB = 3600 * (daylength * sinld + 24 * cosld * math.sqrt(1 - aob * aob)/math.pi)
    dsinBE = 3600 * (daylength * (sinld + 0.4 * (sinld * sinld + cosld * cosld * 0.5)) + 12 * cosld * (2 + 3 * 0.4 * sinld) * math.sqrt(1 - aob * aob)/math.pi)
    
    sc = 1370 * (1 + 0.033 * math.cos(2 * math.pi * day/365))  #solar constant
    
    dso = sc * dsinB  #Daily solar radiation
    
    for hr in range(0,3):
        hour1 = 12 + (daylength * 0.5 * (xgauss[hr]))  #calculates hour in which photosynthesis is applied
        print(hour1)
        
        sinb_tmp = sinld + cosld * math.cos(2 * math.pi * (hour1 + 12)/24)
        print(sinb_tmp)
        
        sinB = max(pd.Series([0, sinb_tmp]))  #calculates sin of solar elevation, max functions prevents values less than 0
        print(sinB)
        
        PAR1 = 0.5 * dso * sinB * (1 + 0.4 * sinB) / dsinBE  #dso can be replaced with values from FAO chart
        print(PAR1)
        
        PAR1 = PAR1 * (868/208.32)  #convert to correct units
        print(PAR1)
        
        tmpvec.append(PAR1)  #output of function is vector of 3 values that represents time of day
        
        
    return tmpvec  #returns a vector of light values in MicroEinsteins

print(PAR(258, 30))


##################################################
#RUN SIMULATION
##################################################
    
#this loops through the simulation time via for loop function

for j in range(0,len(day)):   #length function gets or sets the length of a vector or other objects.
    print('day: {} '.format(day[j]))    
    if day[j] == beginGrow:
        twlvg = 0.25
        twstg = 0.1
        twrtg = 0.15
        
    ##################################################
    #Calculate Vegetation Structure Metrics Each Day
    ##################################################
    
    totLeafWeight = twlvd + twlvg #twlvd = total (t) weight(w) leaves(lv), d @ dead and g is living or green
    np.round(totLeafWeight, 4)
    print('Leaf Weight: {}'.format(totLeafWeight))
    totStemWeight = twstd + twstg #total weight of stems where d at end is dead and g is living or green
    print('Stem Weight: {}'.format(totStemWeight))
    totRootWeight = twrtd + twrtg #total weight of roots where d at end is dead and g is living or green
    print('Root Weight: {}'.format(totRootWeight))
    
    ##################################################
    #Growth and Respiration
    ##################################################
    
    
    #maintenance respiration
    if day[j] < endGrow and day[j] >= beginGrow:
        print('day: {} '.format(j)) 
        if day[j] < endGrow:
            kmLVG = kmLVG_prime * pow(2,((Light.iloc[j]['meantemp'] - 25)/10))  #repiration coefficient for lvs, temp dependence from Teh 2006
            print('kmlvg: {}'.format(kmLVG))
            kmSTG = kmSTG_prime* pow(2,((Light.iloc[j]['meantemp'] - 25)/10)) #respiration coefficient for stems, temp depencence from Teh 2006 page 134
            print('kmSTG: {}'.format(kmSTG))
            
            kmRTG = kmRTG_prime * pow(2,((Light.iloc[j]['meantemp'] - 25)/10)) #respiration coefficient for roots, temp dependence from Teh 2006 page 134
            print('kmRTG: {}'.format(kmRTG))
            
            rmPrime = (kmLVG * twlvg) + (kmSTG * twstg) + (kmRTG * twrtg)  #maintenance respiration per day from Teh 2006
            print('rmPrime: {}'.format(rmPrime))
            
            plantAge = twlvg/totLeafWeight #calculates respiration adjustment based on aboveground biomass, as plants age needs less respiration
            print('plant age: {}'.format(plantAge))
   
            if math.isnan(plantAge):
                plantAge = 0
            
            respMaint = rmPrime * plantAge  #plant age dependence from Teh 2006 page 145
            print('respMaint: {}'.format(respMaint))
      
        #if then statement stops respiration at end of growing season
        if day[j] >= endGrow:
        #if day[j] >= endGrow and day[j] < beginGrow:
            respMaint = 0
        
        #glucose requirement for growth
        glocseReg = (FracDM_LVG * glucoseReqLVG) + (FracDM_STG * glucoseReqSTG) + (FracDM_RTG * glucoseReqRTG) #from Teh 2006 page 148
        print('glocseReg: {}'.format(glocseReg))
        
        #writes results for daily respiration, plant age, and maintenance respiration
        dailyrespiration.iloc[j] = rmPrime
        print('daily respiration: {}'.format(dailyrespiration))
        dailyplantage[j] = plantAge
        print('daily plant age: {}'.format(dailyplantage))
        dailyrespMaint[j] = respMaint
        print('dailyrespMaint: {}'.format(dailyrespMaint))
        
        
        #Enter photosynthesis loop  
      
        if day[j] < endGrow:
            #dailyphoto = []
            for hr in range(0,3):  #radiation measured 3x daily, roughly correlates to morning, noon, afternoon
                parMicroE = (Light.iloc[j,hr]) * (868/208.32) #convert to correct units which is microeinsteins which is the unit measure of light and what this model is based on
                np.round_(parMicroE, decimals = 4) 
                print('day: {}'.format(day[j]))
                print('parMicroE: {}'.format(parMicroE))
                intSolarRad = parMicroE*math.exp(-(PlantParameters['k'])*twlvg)  #from Charisma instructions: tells how much of the light a plant is going to get as PAR in microeinsteins based on how many leaves are on the plant
                np.round_(intSolarRad, decimals = 4)
                print('intSolarRad: {}'.format(intSolarRad))
                intLightpH = intSolarRad/(intSolarRad+Hi) #amoung of light absorbed, per half saturaion constants from Charisma eq. 3. the monod or michaelis/menten function is adequate for describing the photosynthetic response to light
                np.round_(intLightpH, decimals = 4)
                print('intLightpH: {}'.format(intLightpH))
                photosynthesis = (PlantParameters['pMax']) * intLightpH #pMax is the maximum rate of photosynthesis, species specific
                np.round_(photosynthesis, decimals = 4) 
                print('photosynthesis: {}'.format(photosynthesis))
                fgross.iloc[hr] = photosynthesis #calculates gross assimilation of fgross(like APT) via photosynthesis at specific hour calculate growth per day at three times per day, morning, middday, and evenning and this amount is weighted based on how much light is hitting hte plant based on the latitude of your study site
                print('fgross: {}'.format(fgross[hr]))
                dtga.iloc[hr] = fgross[hr]*wgaus[hr] #weights fgross for specific time of day
                print('dtga: {}'.format(dtga[hr]))
                round(dtga[hr],5)
                
                
                
            dtgaCollapsed = sum(dtga)*twlvg  #calculates total biomass gained across plant (twlvg is amount of leaver/green matter): you feed the model total biomass and then from that we determine how much leaf mass there is and so then basically an average of how much that average leaf will produce multiplied by the number of leaves, this is assuming that all leaves are mature
            print('dtgaCollapsed: {}'.format(dtgaCollapsed))
            assimilatedCH2O = dtgaCollapsed*Light['daylength'][day[j]] #total biomass for day length
            print('Day length: {}'.format(Light['daylength'][day[j]]))
            print('assimilatedCH2): {}'.format(assimilatedCH2O))
            gphot = assimilatedCH2O*(30/44) #converts carbohydrates to glucose where photosynthesis unit is glucose and then we later convert that glucose to biomass in another section
            print('gphot: {}'.format(gphot))
            dailyphoto.iloc[j] = gphot
            print('dailyphoto: {}'.format(dailyphoto))
            
        #if then statement ends glucose generation at end of growing season 
        if day[j] >= endGrow:
            gphot = 0
        
 
            
        #############################
        #this section calculates death and mortality of the plant
        #############################
    
    
        #Mortality function
        #how much growth (lv, st, rt) are going to die each based on time of year and temperature, with assumption that there are more leaves dying as it gets hotter
        #Death rate based on julian day and temp - Best and Boyd A24-25
    
        if day[j] > leafDieOff:
            print('Day: {}'.format(day[j]))
            tempi =[-100, 19, 19.01, 30, 30.01, 40, 40.01, 100]
            print('tempi: {}'.format(tempi))
            dRi = [0.021, 0.021, 0.042, 0.042, 0.084, 0.084, 1, 1]
            print('dRi: {}'.format(dRi))
            deathFrac = np.interp(Light.iloc[j]['meantemp'],tempi, dRi) #np.interp seems to be the R equivalent of "approx" which returns a list with components x and y
                                                                 #containing n coordinates which interpolate data points according to the method and rule desired
            print('deathFrac: {}'.format(deathFrac))
            #deathFrac = int(death[2])
            
        if day[j] < leafDieOff:
            print('Day less than leaf die off: {}'.format(day[j]))
            deathFrac = 0.0
            print('deathFrac at day < leaf die off: {}'.format(deathFrac))
       
        
      
        #amount of leaves, stems, and roots that will die each day where day is morning, midday, and evening combined
        #plant death including temperature related death and probablistic death from binomial risk equation (don't have risk equation)
        #ideath = 1, risk = 1, dt_balance = 1 twlvdNEW = 1*twlvg*deathFrac^1, so this statement is actually total weight X deathFrac
        twlvdNEW = twlvg * deathFrac**(risk*dt_balance) # total weight of dead matter produced per timestep
        print('twlvdNEW: {}'.format(twlvdNEW))
        twstdNEW = twstg * deathFrac**(risk*dt_balance) #total weight of dead matter produced per timestep
        print('twstdNEW: {}'.format(twstdNEW))
        twrtdNEW = twrtg*deathFrac**(risk*dt_balance) #total weight of dead matter produced per timestep 
        print('twrtdNEW: {}'.format(twrtdNEW))
   
        #plant growth
        growthWeight = (gphot - respMaint)/glocseReg  #total weight of dry matter produced per day
        print('growth weight: {}'.format(growthWeight))
        if (twlvg < minSize):  #if plant is less than minimum size, it won't grow
            growthWeight = 0
            
        #distribute the life and death into the growing part with a floor of 0
        #FracDM variables are constants of the amount in percentage of how much growth goes towards the three plant parts
        
        twstg = twstg + FracDM_STG * growthWeight - twstdNEW
        print('twstg: {}'.format(twstg))
        twlvg = twlvg + FracDM_LVG * growthWeight - twlvdNEW
        print('twlvg: {}'.format(twlvg))
        twrtg = twrtg + FracDM_RTG * growthWeight - twrtdNEW
        print('twrtg: {}'.format(twrtg))
        
        #if plant is below minimum size, it won't grow and will lose green weight
        
        if twstg < 0.01:
            twstg = 0
            print('twstg: {}'.format(twstg))
        
        if twlvg < 0.01:
            twlvg = 0
            print('twlvg: {}'.format(twlvg))
        
        if twrtg < 0.01:
            twrtg = 0
            print('twrtg: {}'.format(twrtg))
         
        #growth when plant is too small for photosynthesis = 10% of root grows to stems and leaves
        
        reDist = 0.1 * twrtg
        print('reDist: {}'.format(reDist))
        
        if ((twlvg >= minSize) | (twlvg == 0.1)):
            reDist = 0 
            
        twstg = twstg + (reDist * (FracDM_STG/(FracDM_STG + FracDM_LVG))) #redistributes biomass to stems
        print('twstg: {}'.format(twstg))
        twlvg = twlvg + (reDist * (FracDM_LVG/(FracDM_STG + FracDM_LVG))) #redistributes biomass to leaves
        print('twlvg: {}'.format(twlvg))
        twrtg = twrtg - reDist      #bookkeeping to clear redistributed biomass from the roots
        print('twrtg: {}'.format(twrtg))
        
        #constrain the plant biomass
        aboveBiomass = twlvg + twstg
        print('above biomass: {}'.format(aboveBiomass))
        extra = aboveBiomass - maxBiomass
        print('extra: {}'.format(extra))
        
        if ((aboveBiomass < maxBiomass) | (twrtg < 0)):
            extra = 0 
            
        extraLv = extra * (FracDM_LVG/(FracDM_STG + FracDM_LVG))
        print('extraLv: {}'.format(extraLv))
        extraSt = extra * (FracDM_STG/(FracDM_STG + FracDM_LVG))
        print('extraSt: {}'.format(extraSt))
        twlvg = twlvg - extraLv
        print('twlvg: {}'.format(twlvg))
        twstg = twstg - extraSt
        print('twstg: {}'.format(twstg))
        #if root is below minimum then entire plant dies and clear out dead material to reset any potential plant age
        #problems (i.e., including dead material from previous plant)
     
        #distribute the new death into the dead part
        twstd = twstd + twstdNEW
        print('twstd: {}'.format(twstd))
        twlvd = twlvd + twlvdNEW
        print('twlvd: {}'.format(twlvd))
        twrtd = twrtd + twrtdNEW
        print('twrtd: {}'.format(twrtd))
        
        #need to update biomass variable
        biomass = twlvg + twstg
        print('biomass: {}'.format(biomass))
        fracbiomass = biomass / maxBiomass
        print('fracbiomass: {}'.format(fracbiomass))
        
        if (twstg < (minSize * FracDM_STG)):
            twstg = 0
            
        if (twlvg < (minSize * FracDM_LVG)):
            twlvg = 0
            
        if (twrtg < (minSize * FracDM_RTG)):
            twrtg = 0
            
        #need to update biomass var
        biomass = twlvg + twstg
        print('biomass: {}'.format(biomass))
        fracbiomass = biomass / maxBiomass 
        print('fracbiomass: {}'.format(fracbiomass))
        
        #distribute final biomass values into plants, stems, roots
        plant_dens = twstg / tillerDensity * maxTillerWeight 
        print('plant_dens: {}'.format(plant_dens))
        if (plant_dens > maxDensity):
            plant_dens = maxDensity
            
        till_dens = twstg / (plant_dens * tillerDensity)
        print('till_dens: {}'.format(till_dens))
      
        if math.isnan(till_dens):
            till_dens = 0
            
        if (till_dens < 0):
            till_dens = 0
          
        till_ht = till_dens * rWeightTillerHeight / 100
        print('till_ht: {}'.format(till_ht))
        
        if (till_ht > maxTillerHeight):
            till_ht = maxTillerHeight / 100
            
        root_lg = RSratio * till_ht
        print('root_lg: {}'.format(root_lg))
        
        if (root_lg > maxRootLength):
            root_lg = maxRootLength
            
        st_dia = 0
        
        if (till_dens > 0):
            st_dia = 2 * (TCA/math.pi)**(0.5/100)
            
        Density = till_dens  #density value
        print('density: {}'.format(Density))
        Height = till_ht  #height value
        print('height: {}'.format(Height))
        Roots = root_lg  #Roots value
        print('roots: {}'.format(Roots))
        Dia = st_dia   #stem diameter value
        print('dia: {}'.format(Dia))
        AGBiomass = twlvg + twlvd + twstg + twstd
        np.round_(AGBiomass, decimals = 4)
        print('AGBiomass: {}'.format(AGBiomass))
        TOTBiomass = twlvg + twlvd + twstg + twstd + twrtd + twrtg
        np.round_(TOTBiomass, decimals = 4)
        print('TOTBiomass: {}'.format(TOTBiomass))
        
    #daily results outputted to vectors
    
    dailyleafg.iloc[j] = twlvg
    dailyrootg.iloc[j] = twrtg
    dailystemg.iloc[j] = twstg
    dailyleafd.iloc[j] = twlvd
    dailyrootd.iloc[j] = twrtd
    dailystemd.iloc[j] = twstd
    dailyAGbiomass.iloc[j] = twlvg + twlvd + twstg + twstd
    dailyTOTbiomass.iloc[j] = twlvg + twlvd + twstg + twstd + twrtd + twrtg
    dailyRTbiomass.iloc[j] = twrtd + twrtg
    #dailyAGbiomass1.iloc[j] = twlvg + twlvd + twstg + twstd
    
    days = day.tolist()
    avgtemp = Light["meantemp"].tolist()
    dailyleafg_list = dailyleafg.tolist()
    dailyrootg_list = dailyrootg.tolist()
    dailystemg_list = dailystemg.tolist()
    dailyleafd_list = dailyleafd.tolist()
    dailyrootd_list = dailyrootd.tolist()
    dailystemd_list = dailystemd.tolist()
    dailyAGbiomass_list = dailyAGbiomass.tolist()
    dailyTOTbiomass_list = dailyTOTbiomass.tolist()
    dailyRTbiomass_list = dailyRTbiomass.tolist()
    

#data_needs = [days, avgtemp, dailyleafg_list, dailyrootg_list, dailystemg_list, dailyleafd_list, dailyrootd_list, dailystemd_list, dailyAGbiomass_list, dailyTOTbiomass_list]

data_needs = [day, avgtemp, dailyleafg, dailyrootg, dailystemg, dailyleafd, dailyrootd, dailystemd, dailyAGbiomass, dailyTOTbiomass]


#########OUTPUT FILE#########
'''
output_file = open("VegOutput.csv", "w")
output_file.write("Day, mean_temp, lvsg, rootsg, stemg, lvsd, rootsd, stemd, AGB, TOTB\n")

output_file.write(str(days) + "," + str(Light["meantemp"]) + "," + str(dailyleafg_list) + "," + str(dailyrootg_list) + "," + str(dailystemg_list) + "," + str(dailyleafd_list) + "," + str(dailyrootd_list) + "," + str(dailystemd_list) + "," + str(dailyAGbiomass_list) + "," + str(dailyTOTbiomass_list) + "\n")

output_file.close()
'''                


  
output_file = open("VegOutput.csv", "w")
output_file.write("Day, mean_temp, lvsg, rootsg, stemg, lvsd, rootsd, stemd, AGB, TOTB\n")  
for i in range(len(data_needs[0])):
        output_file.write(str(data_needs[0][i]) + "," + str(data_needs[1][i]) + "," + str(data_needs[2][i]) + "," + str(data_needs[3][i]) + "," + str(data_needs[4][i]) + "," + str(data_needs[5][i]) + "," + str(data_needs[6][i]) + "," + str(data_needs[7][i]) + "," + str(data_needs[8][i]) + "," + str(data_needs[9][i]) +  "\n")



output_file.close()



    

       
 
    
    
    
    














































