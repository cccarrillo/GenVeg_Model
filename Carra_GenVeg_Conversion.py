  #############################################################################################################################
  #Model:		  R Veg Model 
  #Purpose:	Simulates vegetation growth and mortality
  #Last Modified by:	Swannack/Herman/Altman
  #Date: 		13 Oct 2020
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

###################################################

#Set length of simulation (Day of Year)
startsim = 1  #start of simulation
endsim = 365  #end of simulation
day = list(range(startsim, endsim+1))


dt_balance = 1   #timestep
endGrow = 320    #date that vegetation stops growing due to day of first freeze (based on latitude)
leafDieOff = 163 #date that plant starts to self-shade
#######################################################

#Initialize variables

#initialize veg parameter variables and vectors

twlvd = 0
twstd = 0
twrtd = 0
twlvg = 0
twstg = 0
twrtg = 0
fgross = pd.Series([0,0,0])
dtga = pd.Series([0,0,0])
wgaus = pd.Series([0.2778,0.4444,0.2778])
xgauss = pd.Series([0.1127,0.5,0.8873], dtype=np.dtype("float"))
Hi = 40
risk = 1



#initialize vectors for daily outputs of critical variables
dailyleafg = pd.Series([startsim,endsim])
dailyrootg = pd.Series([startsim,endsim])
dailystemg = pd.Series([startsim,endsim])
dailyleafd = pd.Series([startsim,endsim])
dailyrootd = pd.Series([startsim,endsim])
dailystemd = pd.Series([startsim,endsim])
dailyAGbiomass =pd.Series([startsim,endsim])
dailyrespiration = pd.Series([startsim,endsim])
dailyplantage = pd.Series([startsim,endsim])
dailyrespMaint =pd.Series([startsim,endsim])
dailyphoto = pd.Series([startsim,endsim])
dailyTOTbiomass = pd.Series([startsim,endsim])
dailyRTbiomass = pd.Series([startsim,endsim])
############################################################

#FUNCTIONS CALLED IN MODEL

#Function: logistic shape

#creates non-linear response curves for effects
#inputs: XY pairs for curve shape, input value for curve
#output: value of effect, given specific curve and input
def logistic (X1, Y1, X2, Y2, x):
  G = math.log((Y1)/(1-Y1))
  C = math.log((Y2)/(1-Y2))
  B = ((G-C)/(X1-X2))
  A = (G-(B*X1))

  Z = math.exp(A+(B*x))
  S = (Z/(1+Z))
  return(S)

print(str(logistic(1,0.5,0,0.5,0.25)) + ' = 0.5 from R Code')
#Function:PAR
#Calculates daylength, amount of PAR at different times of day
#Packages required: math, numpy, pandas
#Inputs: Day of year and latitude
#Output: PAR values at 3 different times of day for Gaussian Integration



twlvg = 10
k = -0.02
pMax = 0.25
#day = 258
lat = 30

def PAR(day, lat):
  degree_to_rad = 0.017453292  #required to convert degrees to radians
  rad_to_degree = (-math.asin((math.sin(23.45*degree_to_rad))*(math.cos(2*math.pi*(day+10)/365))))

  tmpvec = pd.Series([0,0,0], dtype=np.dtype("float"))
  declination = (-math.asin(math.sin(23.45*degree_to_rad))*(math.cos(2*math.pi*(day+10)/365)))


  #intermediate variables
  sinld = ((math.sin(lat*degree_to_rad))*(math.sin(declination))) #radians
  cosld = math.cos(lat*degree_to_rad) * math.cos(declination) #radians #DOES THIS NEED TO BE ENCLOSED IN PARENTHASIS LIKE THE SINLD IS??
  aob = (sinld/cosld) #radians
  temp1 = math.asin(aob)
  global daylength
  daylength = 12*(1+2*temp1/math.pi) #calculates daylength based on declination and latitude
  dsinB = 3600*(daylength*sinld+24*cosld*math.sqrt(1-aob*aob)/math.pi)
  dsinBE = 3600 * (daylength * sinld + 0.4 * (sinld * sinld + cosld * cosld *0.5)) + 12 * cosld * (2 + 3 * 0.4 * sinld) * math.sqrt (1 - aob * aob) / math.pi
  sc = 1370 * (1 + 0.033 * math.cos(2 * math.pi * day / 365)) #solar constant
  dso = sc * dsinB #daily solar radiation

  for hr in range (1,4):
    hour1 = 12 + (daylength * 0.5 * xgauss[hr - 1]) #calculates hour in which photosynthesis is applied
    sinb_tmp = sinld + cosld * math.cos(2 * math.pi * (hour1 + 12) / 24)
    sinB = max(0, sinb_tmp) #calculates sin of solar elevation, max functions prevents values less than 0
    PAR1 = (0.5 * dso * sinB * (1 + 0.4 * sinB) / dsinBE) #NBdso can be replaced with values from the FAO chart
    PAR1 = PAR1 * (868/208.32) #converts to correct units
    #print(PAR1)
    tmpvec[hr - 1] = PAR1 #output of function is vector of 3 values that represents time of day
    #print(tmpvec)

    return tmpvec
#####################################

#Run Simulation
#This loops through the simulation time via for lop function

for j in range(1,len(day)):  #length() function gets or sets the length of a vector (list) or other objects
    if day[j] == beginGrow:  #begin grow date is rice emergence day at 10C degrees from input file
                             #beginGrow comes from the SetPlantParameters function
        twlvg = 0.25         #initial emergence
        twstg = 0.1
        twrtg = 0.15
        
        #calculate vegetation structure metrics each day
        
        totLeafWeight = twlvd + twlvg #twlvd = total(t) weight (w) leaves (lv), d @ end is dead and g is living or green
        totStemWeight = twstd + twstg #total weight of stems where d @ end is dead and g is living or green
        totRootWeight = twrtd + twrtg #total weight of roots where d @ end id dead and g is living or green
        
        #Growth and Respiration
        
        if (day[j] < endGrow) and (day[j] >= beginGrow):
            if day[j] < endGrow:
                kmLVG = kmLVG_prime * 2**((((Light["meantemp"][day[j]]-25)/10))) #respiration coefficiend for lvs, temperature dependence from Teh (2006)
                kmSTG = kmSTG_prime * 2**((((Light["meantemp"][day[j]]-25)/10))) #respiration coefficiend for lvs, temperature dependence from Teh (2006)
                kmRTG = kmRTG_prime * 2**((((Light["meantemp"][day[j]]-25)/10))) #respiration coefficiend for lvs, temperature dependence from Teh (2006)
                rmPrime = (kmLVG * twlvg) + (kmSTG * twstg) + (kmRTG * twrtg) #maintenance respiration per day from Teh (2006)
                plantAge = twlvg/totLeafWeight #calculates respiration adjustment based on aboveground biomass, as plant age needs less respiration
                plantAge.np.NA(0) #sets the variable to 0 if value if NaN
                respMaint = rmPrime * plantAge #plant age dependence from Teh (2006) page 145
                
                
        # if/then statement stops respiration at end of growing season
        if (day[j] == endGrow):
            respMaint = 0
            
            
        #glucose requirement for growth
        glocseReg = (FracDM_LVG * glucoseReqLVG) + (FracDM_STG * glucoseReqSTG) + (FracDM_RTG * glucoseReqRTG) #from Teh (2006) page 148
        
        
        #writes results for daily respiration, plant age, and maintenance respiration
        dailyrespiration[j] = rmPrime
        dailyplantage[j] = plantAge
        dailyrespMaint[j] = respMaint
        
        #Enter photosynthesis loop
        
        if (day[j] < endGrow):
            
            for hr in range(1,4):      
                parMicroE = Light[day[j],hr] * (868/208.32) #convert to correct units which is microeinsteins which is the unit measure of light and what this model is based on
                intSolarRad = parMicroE ** ((-k * twlvg)) #from Charisma instructions (tells you how much of the light a plant is goig to get as PAR in microeinstens based on number of leaves on plant)
                intLightPh = intSolarRad/(intSolarRad + Hi) #amount of light absorbed, per hald saturation constants from Charisma equation 3 The Monod or Michaelis/Menten function is adequate for describing photosynthetic response
                phytosynthesis = pMax * intLightPh  #pMax is the maximum rate of photosynthesis, species specific
                fgross[hr] = photosynthesis #calculates gross assimilation of fgross (like APT) via photosynthesis at specific hour calculate growth per day at 3X a day (morning, midday, and evening) and 
                                        #this amount is weighed based on how much light is hitting the plant based on the latitude of your study site
                dtga[hr] = fgross[hr] * wgaus[hr] #weghts fgross for specific time of day
            
            
        dtgaCollapsed = sum(dtga) * twlvg #calulates total biomass gaines across plant (twlvg is amoung of leaves/green matter).  you feed model total biomass and then from that we determine how much leaf mass there is and so then
                                          #basically an average of how much that average leaf will produce multiplied by the number of leaves, this is assuming that all leaves are mature
        assimilatedCH20 = dtgaCollapsed * Light["daylength"[day[j]]]#total biomass for day length
        
        gphot = assimilatedCH20 * (30/44) #converts carbohydrates to glucose where photosynthesis unit is glucose and then we later convert that glucose to biomass in another section
        
        dailyphoto[j] = gphot[j]
        
        
        #if/then statement ends glucose generation at end of growing season
        if (day[j] >= endGrow):
            gphot = 0
            
            
            
        # This section calculates death and mortality of the plant
        
        #mortality function
        # how much growl (lv, st, rt) are going to die each day based on time of year and temperature, with assumption that there are more leaves dying as it gets hotter
        # Death rate based on Julian day and temp - Best and Boyd A24-25
        
        if (day[j] > leafDieOff):
            tempi = pd.Series([-100, 19, 19.01, 30, 30.01, 40, 40.01, 100])
            dRi = pd.Series([0.021, 0.021, 0.042, 0.042, 0.084, 0.084, 1, 1])
            death = np.interp(tempi, dRi, Light["meantemp"][day[j]]) #the np.interp function seems to be the equivalent of the R "Approx" funciton where is returns a list with components X and Y, 
                                                                            #containing n coordinates which interpolate the given data points according to the method (and rule) desired
                                                                            
            deathFrac = int(death[2]) #int converts simply the index part of a factor into numeric integer
            
        else:
            deathFrac = 0
                
            
            #amount of leaves, stems, and roots that will die each day where day is morning, midday, and evening combined
            #plant death including temperature related death and probabilistic death from binomial risk equation (don't have risk equation),
            #ideath = 1, risk = 1, dt_balance = 1 (time step) twlvdNEW = 1 * twlvg * deathFrac(as above) ^1 * 1, so this statement is actually total weight * deathfrac
            
            twlvdNEW = twlvd * deathFrac ** (risk * dt_balance) #total weight of dead matter produced per timestep
            twstdNEW = twstg * deathFrac ** (risk * dt_balance) #total weight of dead matter produced per timestep
            twrtdNEW = twrtg * deathFrac ** (risk * dt_balance) #total weight of dead matter producted per timestep
            
            #plant growth
            growWeight = (gphot - respMaint) / glocseReg #total weight of dry matter produced per day
            
            if (twlvg < minSize):  #if plant is less thatn minimum size, it won't grow
                growthWeight = 0  
                
            #distribute the life and death into the growing part with a floor of 0 
            #FracDM variables are constants of the amount in percentage of how much growth goes towards the three plant parts
            
            twstg = twstg + FracDM_STG * growthWeight - twstdNEW
            twlvg = twlvg + FracDM_LVG * growthWeight - twlvdNEW
            twrtg = twrtg + FracDM_RTG * growthWeight - twrtdNEW
            
            #if plant is below minimum size, it won't grow and will loose green weight 
            if (twstg < 0.01):
                twstg = 0
               
               
            if (twlvg < 0.01):
                twlvg = 0
                
            if (twrtg < 0.01):
                twrtg = 0
                
                
            # growth when plant is too small for photosynthesis - 10% of root grows to stems and leaves
            reDist = 0.1 * twrtg
            
            if ((twlvg >= minSize) | (twlvg == 0.01)):
                reDist = 0
                
            twstg = twstg + (reDist * (FracDM_STG/(FracDM_STG+FracDM_LVG))) #redistributes biomass to stems
            twlvg = twlvg + (reDist * (FracDM_LVG/(FracDM_STG+FracDM_LVG))) #redistributes biomass to leaves
            twrtg = twrtg - reDist #bookkeeping to clear redistributed biomass from roots
            
            
            #constraing the plant biomass
            aboveBiomass = twlvg + twstg
            extra = aboveBiomass - maxBiomass
            
            if ((aboveBiomass < maxBiomass) | (twrtg < 0)):
                extra = 0
                
            extraLv = extra * (FracDM_LVG/(FracDM_STG + FracDM_LVG))
            extraSt = extra * (FracDM_STG/(FracDM_STG + FracDM_LVG))
            twlvg = twlvg - extraLv
            twstg = twstg - extraSt
            
            #if root is below minimum weight then entire plant dies and clear out dead material to reset any potential plant age
            #problems (i.e., including dead material from previous plants)
            
            #distribute the new death into the dead part
            twstd = twstd + twstdNEW
            twlvd = twlvd + twlvdNEW
            twrtd = twrtd + twrtdNEW
            
            #need to update biomass variable
            biomass = twlvg + twstg
            fracbiomass = biomass / maxBiomass
            
            if (twstg < (minSize * FracDM_STG)):
                twstg = 0
                
            if (twlvg < (minSize * FracDM_LVG)):
                twlvg = 0
                
            if (twrtg < (minSize * FracDM_RTG)):
                twrtg = 0
                
            #need to update biomass var
            biomass = twlvg + twstg
            fracbiomass = biomass/maxBiomass
            
            #Distribute final biomass values into plants, stems, roots
            plant_dens = twstg/tillerDensity*maxTillerWeight
            
            if (plant_dens > maaxDensity):
                plant_dens = maxDensity
            
            till_dens = twstg / (plant_dens * tillerDensity)
            till_dens[np.NaN(till_dens)] = 0
            
            if (till_dens < 0):
                till_dens = 0
                
            till_ht = till_dens * rWeightTillerHeight/100 
            
            if (till_ht > maxTillerHeight):
                till_ht = maxTillerHeight/100 
                
            root_lg = RSratio * till_ht
            
            if (root_lg > maxRootLength):
                root_lg = maxRootLength
                
            st_dia = 0
            
            if (till_dens > 0):
                st_dia = 2 * (TCA/pi)**(0.5/100)
                
            Density = till_dens #Density value
            Height = till_ht    #Height value
            Roots = root_lg     #Roots value
            Dia = st_dia        #Stem Diameter value
            AGBiomass = twlvg + twlvd + twstg + twstd
            TOTBiomass = twlvg + twlvd + twstg + twstd + twrtd + twrtg
            
            
        #Daily Results outputted to vectors
        
        dailyleafg[j] = twlvg
        dailyrootg[j] = twrtg
        dailystemg[j] = twstg
        dailyleafd[j] = twlvd
        dailyrootd[j] = twrtd
        dailystemd[j] = twstd 
        dailyAGbiomass[j] = twlvg + twlvd + twstg + twstd
        dailyTOTbiomass[j] = twlvg + twlvd + twstg + twstd + twrtd + twrtg
        dailyRTbiomass[j] = twrtd + twrtg
        #dailyAGbiomass1[j] = twlvg + twlvd + twstg + twstd
        
#CARRA"S OWN CODE WRITING!!!!!

csv_list = [day, Light["meantemp"], dailyleafg, dailyrootg, dailystemg, dailyleafd, dailyrootd, dailystemd, dailyAGbiomass, dailyTOTbiomass] 
        
#creates output file
def VegOutputCSV(filename):
    output_file_start = open(filename, "w")
    output_file_start.write("Day, Mean_Temp, lvsg, rootsg, stemg, lvsd, rootd, stemd, AGB, TOTB\n")
    output_file_start.write(csv_list)
    output_file_start.close()
    

'''
parMicroE = PAR(day, lat)


#Function: Photo1 (simplet representation)
#calculates grams of carbohydrate accumulated per day via photosynthesis
#Packages required: math, numpy, pandas
#Inputs: ParMicroE
#Output: grams of carbohydrate per day

def photo1(parMicroE):
  
  for hr in range (0,3):
    intSolarRad = parMicroE[hr] * math.exp(-k * twlvg)
    intLight = intSolarRad / (intSolarRad + Hi)
    photosyn = pMax * intLight
    fgross[hr] = photosyn
    dtga[hr] = fgross[hr] * wgaus[hr]

  dtgaCollapsed = pd.Series.sum(dtga)
  assimilatedCH2O = dtgaCollapsed * daylength #be sure to multiply by current biomass
  return(assimilatedCH2O * 30/44)

gphoto = photo1(parMicroE)

#Function: Photo2 (photosynthesis at different plant layers)
#calculates grams of carbohydrate accumulated per day via photosynthesis
#packages required: math, numpy, pandas
#inputs: ParMicroE
#outputs: grams of carbohydrate per day

def photo2(parMicroE):
  numlayer = 3 #number of layers at which photosynthesis is calculated

  for i in range(len(parMicroE)):

    return


#function: lightAtten (intended for SAV)
#Calculates how light is attenuated
#packages required: math, numpy, pandas
#Inputs: parMicroE
#outputs: grams of carbohydrate per day

def lightAtten(depth, parMicroE, current_biomass):  #i think this needs to be an underline instead of hyphen in current-biomass

  iz = parMicroE
  return 

  '''



