  #############################################################################################################################
  #Model:		  R Veg Model 
  #Purpose:	Simulates vegetation growth and mortality
  #Last Modified by:	Swannack/Herman/Altman
  #Date: 		13 Oct 2020
  #############################################################################################################################
  
  #############################################################################################################################
  #IMPORT DATA and initiate variable values
  #############################################################################################################################
  
  library(ggplot2)
  #Clear memory
  rm(list = ls())
  #getwd()
  setwd("/Users/rdel1cmc/Desktop/GenVeg_Model/R_Code/")
  
  ############################################################################################################################
  #Sets focal crop species
  Species = "corn" 
  ############################################################################################################################
  
  ############################################################################################################################
  #INPUT FILES FOR MEAN DAILY TEMPERATURE 
  Light <- read.csv("../Model_Data/Light.csv",header=TRUE)
  Plants <- read.csv("../Model_Data/Plants.csv",header=TRUE)
  ############################################################################################################################
  
  ###########################################################################################################################
  #REMOVE DURING SIIMULATIONS
  set.seed(1) ##this fixes the day (5) that a fire will occur if a fire is probable. 
  ############################################################################################################################
  
  ############################################################################################################################
  #MODEL PARAMETERS
  
  ##############################################################
  #FUNCTION: setPlantParameters
  #sets species specific plant parameters
  #Input: Input file Plants.csv required
  #output: species specific plant parameters
  
  setPlantParameters <- function(Plants, Species){
    #Locate the right entry in the database
    findSpecies <- which(colnames(Plants)==Species)
    
    #Set plant parameters from input file
    k <<- (Plants[1,findSpecies])             #k: Light attenuation coefficient, Default at 0.02 m2/g. Rice  - 0.38 from Charkraborty etal. 2018. 
    
    pMax <<- (Plants[2,findSpecies])          #Maximum rate of photosynthesis for each crop, specific daily production of the plant top at 20C degrees in the absence of resource limitation
    
    plantMature <<- (Plants[2,findSpecies])   #PlantMature: day when plant reaches maturity. 
    
    beginGrow <<- (Plants[4,findSpecies])     #day crop emerges
    
    Cs137TF <<- (Plants[5,findSpecies])       #transfer factor for Cesium
    
    Sr90TF <<- (Plants[6,findSpecies])        #transfer factor for Strontium
    
    I129TF <<- (Plants[7,findSpecies])        #transfer factor for Iodine
    
    
  } #end fxn call
  
  
  setPlantParameters(Plants, Species) #function that sets species-specific parameters
  
  print(pMax)
  
  tillerDensity <- 3
  maxTillerHeight <- 80
  maxTillerWeight <- 6
  rWeightTillerHeight <- 13.333
  maxRootLength <- 60
  maxDensity <- 129
  maxBiomass <- 2322
  minSize <- 0          
  RSratio <- 0.95
  FracDM_LVG <- 0.5       #amount of carbon that goes into leaves
  FracDM_STG <- 0.2       #amount of carbon that goes into stems
  FracDM_RTG <- 0.3       #amount of carbon that goes into roots
  
  
  kmLVG_prime <- 0.03    #set from Teh 2006 Table 7.1 Leaves base maintenance respiration coefficient at 25 °C
  kmSTG_prime <- 0.015   #set from Teh 2006 Table 7.1 Stems base maintenance respiration coefficient at 25 °C
  kmRTG_prime <- 0.015   #set from Teh 2006 Table 7.1 Roots base maintenance respiration coefficient at 25 °C
  glucoseReqLVG <- 1.436 #set from Teh 2006 Table 7.4 
  glucoseReqSTG <- 1.513 #set from Teh 2006 Table 7.4 
  glucoseReqRTG <- 1.444 #set from Teh 2006 Table 7.4 
  TCA <- 0.231          #total cross sectional area
  #############################################################
  
  #############################################################
  #Set length of simulation (DAY OF YEAR)
  startsim <- 1              #Start of simulation
  endsim <- 365              #End of simulation
  day <- c(startsim:endsim)  #Length of simulation, creates a vector from 1 to 365 given the values of startsim and endsim defined 
  
  dt_balance <- 1            #Timestep
  endGrow <- 320             #Date that vegetation stops growing due to day of first freeze (changes based on latitude)
  leafDieOff <- 163          #Date that plant starts to self-shade
  ##############################################################
  
  ##############################################################
  #INTITALIZE VARIABLES
  
  #Initialize veg parameter variables and vectors
  twlvd <- 0       
  twstd <- 0
  twrtd <- 0
  twlvg <- 0
  twstg <- 0
  twrtg <- 0
  fgross <- c(0,0,0)
  dtga <- c(0,0,0)
  wgaus <- c(0.2778,0.4444,0.2778)
  
  
  #INTIALIZES VECTORS FOR DAILY OUTPUTS OF CRITICAL VARIABLES
  dailyleafg <- c(startsim:endsim)
  dailyrootg <- c(startsim:endsim)
  dailystemg <- c(startsim:endsim)
  dailyleafd <- c(startsim:endsim)
  dailyrootd <- c(startsim:endsim)
  dailystemd <- c(startsim:endsim)
  dailyAGbiomass <-c(startsim:endsim)
  dailyrespiration <- c(startsim:endsim)
  dailyplantage <- c(startsim:endsim)
  dailyrespMaint <-c(startsim:endsim)
  dailyphoto <- c(startsim:endsim)
  dailyTOTbiomass <- c(startsim:endsim)
  dailyRTbiomass <- c(startsim:endsim)
  
  
  ##############################################################
  #FUNCTIONS CALLED IN MODEL
  
  ##############################################################
  #FUNCTION: logistic shape 
  #Creates non-linear response curves for effects
  #Inputs: XY Pairs for curve shape, input value for curve
  #Output: value of effect, given specific curve and input
  logistic <- function(X1,Y1,X2,Y2, x) {
    G <- log((Y1)/(1-Y1))
    C <- log((Y2)/(1-Y2))
    B <- ((G-C)/(X1-X2))
    A <- (G-(B * X1))
    
    Z <- exp(A+(B*x))
    S <- (Z/(1+Z))
    
    return(S)
    
  }
  print(logistic(1,0.5,0,0.5,0.25))
  
  ##############################################################
  #FUNCTION: PAR 
  #Calculates daylength, amount of photosynthetically available light, Determines amount of plant growth per day (ASTRO Procedure from )
  #Inputs: Day of year and latitude
  #output: PAR values at 3 different times of day for Gaussian integration
  
  PAR.Function <- function(day, lat){
    
    tmpvec <- c()
    declination <- (-asin ((sin(23.45 * degree.to.rad)) * (cos(2 * pi * (day + 10) / 365))))
    
    #intermediate variables
    sinld <- ((sin(lat * degree.to.rad)) * (sin(declination))) #radians
    cosld <- cos(lat * degree.to.rad) * cos(declination) #radians
    
    aob <- (sinld / cosld)  #radians
    
    temp1 <- asin(aob)
    
    daylength <- 12 * (1 + 2 * temp1 / pi) #calculates daylength based on declination and latitude
    
    
    dsinB <- 3600 * (daylength * sinld + 24 * cosld * sqrt(1 - aob * aob) / pi)
    dsinBE <- 3600 * (daylength *(sinld + 0.4 * (sinld * sinld + cosld * cosld * 0.5)) + 12 * cosld * (2 + 3 * 0.4 * sinld) * sqrt (1 - aob * aob) / pi)  
    
    sc <- 1370 * (1 + 0.033 * cos(2 * pi * day / 365)) #Solar constant
    
    
    dso <- sc * dsinB     #Daily solar radiation
    
    
    for(hr in 1:3){
      
      hour1 <- 12 + (daylength * 0.5 * xgauss[hr]) #calculates hour in which photosynthesis is applied
      print (paste("hour1r: ", hour1))
      
      sinb.tmp <-  sinld + cosld * cos(2 * pi * (hour1 + 12) / 24)
      #print (paste("sinb.tmp", sinb.tmp))  
      #sinb.tmp <-c(o,sinb.tmp1)
      
      sinB <- max(c(0,sinb.tmp)) #calculates sin of solar elevation, max functions prevents values less than 0
      #print (paste("sinB", sinB))
      
      PAR1 <- 0.5 * dso * sinB * (1 + 0.4 * sinB) / dsinBE #NB: dso can be replaced with values from the FAO chart
      #print (paste("PAR1", PAR1))
      
      PAR1 <- PAR1 * (868/208.32) #Convert to correct units
      #print (paste("hr: ", hr, " PAR1: ", PAR1, " tmpvec ", tmpvec))
      
      tmpvec[hr] <- PAR1 #output of function is vector of 3 values that represents time of day
      #print (tmpvec)
      
    }
    
    return(tmpvec) #returns a vector of light values in MicroEinsteins
    
  } #end fxn call
  
  
  
  #############################################################################################################################
  #Run Simulation
  #############################################################################################################################
  
  #This loops through the simulation time via for loop function 
  
  for (j in 1:length(day)){     #length() function gets or sets the length of a vector (list) or other objects
    if (day [j] == beginGrow){  #begin grow date is rice emergence day at 10C degrees from input file.
      twlvg <- 0.25             #initial emergence 
      twstg <- 0.1
      twrtg <- 0.15  
    }
    
    
    #############################################################################################################################
    #Calculate vegetation structure metrics each day
    #############################################################################################################################
    
    totLeafWeight <- twlvd+twlvg #twlvd=total(t) weight(w) leaves(lv), d @ end is dead & g is living or green 
    totStemWeight <- twstd+twstg #total weight of stems where d at end is dead and g is living or green 
    totRootWeight <- twrtd+twrtg #total weight of roots where d at end is dead and g is living or green
    
            
    
    
    #############################################################################################################################
    #Growth and Respiration
    #############################################################################################################################
    
    
    #maintenance respiration
    if (day[j]<endGrow && day[j] >= beginGrow){ 
      if (day[j]<endGrow){
        kmLVG <- kmLVG_prime*2^((Light$meantemp[day[j]]-25)/10) #respiration coefficient for lvs, temperature dependence from Teh (2006)
        kmSTG <- kmSTG_prime*2^((Light$meantemp[day[j]]-25)/10) #respiration coefficient for stems, temperature dependence from Teh (2006) page 134
        kmRTG <- kmRTG_prime*2^((Light$meantemp[day[j]]-25)/10) #respiration coefficient for roots, temperature dependence from Teh (2006) page 134
        
        rmPrime <- (kmLVG*twlvg)+(kmSTG*twstg)+(kmRTG*twrtg) #maintenance respiration per day from Teh (2006) 
        
        plantAge <- twlvg/totLeafWeight #calculates respiration adjustment based on aboveground biomass, as plants age needs less respiration 
        
        
        plantAge[is.na(plantAge)] <- 0  #sets the variable to 0 if value is NA
        respMaint <- rmPrime*plantAge   #plant age dependence from Teh (2006) page 145
        
      } 
      
      #if then statement stops respiration at end of growing season
      if (day[j]>= endGrow){
        respMaint <- 0
      }
      
      #glucose requirement for growth
      glocseReg <- (FracDM_LVG*glucoseReqLVG) + (FracDM_STG*glucoseReqSTG) +  (FracDM_RTG*glucoseReqRTG)  #from Teh (2006) page 148
      
      #writes results for daily respiraiton, plant age and maintenance respriation
      dailyrespiration[j] <- rmPrime
      dailyplantage[j] <- plantAge
      dailyrespMaint[j] <-respMaint
      
     #Enter photosynthesis loop
      
      if (day[j]<endGrow){
        
        for (hr in 1:3){                              #radiation measured 3X day, roughly correlates to morning, noon, afternoon  
          parMicroE <- Light[day[j],hr]*(868/208.32)  #Convert to correct units which is microeinsteins which is the unit measure of light and what this model is based on 
          intSolarRad <- parMicroE*exp(-k*twlvg)      #from Charisma instructions--> tells you how much of the light a plant is going to get as PAR in microeinsteins based
                                                      #on how many leaves are on the plant 
         intLightPh <- intSolarRad/(intSolarRad+Hi) #amount of light absorbed, per half-saturation constants from Charisma equ. 3. The Monod or Michaelis/Menten function is
          #adequate for describing the photosynthetic response to light. 
          photosynthesis <- pMax*intLightPh         #pMax is the maxiumum rate of photosynthesis, species specific       
          fgross[hr] <- photosynthesis              #calculates gross assimilation of fgross(like APT) via photosynthesis at specific hour calculate growth per day at three times per day, morning, middday, and evenning and this amount is weighted based on how much light is hitting hte plant based on the latitude of your study site
          dtga[hr] <- fgross[hr]*wgaus[hr]          #weights fgross for specific time of day
        }
        
        
        dtgaCollapsed <- sum(dtga)*twlvg            #calculates total biomass gained across plant (twlvg is amount of leaves/green matter) --> you feed the model total biomass and then from that we determine how much leaf mass there is and so then basically an average of how much that average leaf will produce multiplied by the number of leaves, this is assuming that all leaves are mature                            
        assimilatedCH20 <- dtgaCollapsed*Light$daylength[day[j]]  #Total biomass for day length
        gphot <- assimilatedCH20*(30/44)            #converts carbohydrates to glucose where photosysnthesis unit is glucose and then we later convert that glucose to biomass in another section 
        
        dailyphoto[j] <- gphot
        
      } 
      
      # if then statement ends glucose generation at end of growing season
      if (day[j]>= endGrow){
        gphot <- 0
      }
      
      #############################################################################################################################
      #This section calculates death and mortality of the plant
      #############################################################################################################################
      
      #Mortality function 
      #how much growth (lv, st, rt) are going to die each day based on time of year and temperature, with assumption that there are more leaves dying as it gets hotter 
      #Death rate based on julian day and temp - Best and Boyd A24-25
      if (day[j]>leafDieOff){
        tempi <- c(-100, 19, 19.01, 30, 30.01, 40, 40.01, 100)
        dRi <- c(0.021, 0.021, 0.042, 0.042, 0.084, 0.084, 1, 1)
        death <- approx(tempi,dRi,Light$meantemp[day[j]]) #approx returns a list with components x and y, containing n coordinates which interpolate the given 
        #data points according to the method (and rule) desired
        deathFrac <- as.numeric(death[2])               #as.numeric converts simply the index part of factor into numeric.
      } else {
        
        deathFrac <- 0.0
      }
      
      
      #amount of leaves, stems and roots that will die each day where day is morning, midday, and evenning combined 
      #plant death including temperature related death and probablistic death from binomial risk equation (don't have risk equation), 
      #ideath = 1, risk = 1, dt_balance = 1 (time step),#twlvdNEW = 1 * twlvg * deathFrac(as above)^1 * 1, so this statement is actually total weight  X deathFrac.
      twlvdNEW <- twlvg*deathFrac^risk*dt_balance #total weight of dead matter produced per timestep
      twstdNEW <- twstg*deathFrac^risk*dt_balance #total weight of dead matter produced per timestep
      twrtdNEW <- twrtg*deathFrac^risk*dt_balance #total weight of dead matter produced per timestep
      
      
      #plant growth
      growthWeight <- (gphot-respMaint)/glocseReg     #total weight of dry matter produced per day
      
      
      if (twlvg<minSize){ #if plant is less than minimum size, it won't grow
        growthWeight <- 0
      }
      
      #distribute the life and death into the growing part with a floor of 0
      #FracDM variables are constants of the amount in percentage of how much growth goes towards the three plant parts 
     
      twstg <- twstg + FracDM_STG*growthWeight - twstdNEW
      twlvg <- twlvg + FracDM_LVG*growthWeight - twlvdNEW
      twrtg <- twrtg + FracDM_RTG*growthWeight - twrtdNEW
      
      #if plant is below minimum size, it won't grow and will loose green weight 
      if (twstg<0.01){
        twstg <- 0
      }
      
      if (twlvg<0.01){
        twlvg <- 0
      }
      
      if (twrtg<0.01){
        twrtg <- 0
      }
      
      #growth when plant is too small for photosyntehsis - 10% of root grows to stems and leaves
      reDist <- 0.1*twrtg
      
      if ((twlvg>=minSize) | (twlvg==0.1)){
        reDist <- 0
      }
      
      twstg <- twstg + (reDist*(FracDM_STG/(FracDM_STG+FracDM_LVG))) #redistributes biomass to stems
      twlvg <- twlvg + (reDist*(FracDM_LVG/(FracDM_STG+FracDM_LVG))) #redistributes biomass to leaves
      twrtg <- twrtg - reDist                                        #bookkeeping to clear redistributed biomass from the roots
      
      #constrain the plant biomass
      aboveBiomass <- twlvg + twstg
      extra <- aboveBiomass-maxBiomass
      
      if ((aboveBiomass<maxBiomass) | (twrtg<0)){
        extra <- 0
      }
      
      extraLv <-extra*(FracDM_LVG/(FracDM_STG + FracDM_LVG))
      extraSt <- extra*(FracDM_STG/(FracDM_STG + FracDM_LVG))
      twlvg <- twlvg-extraLv
      twstg <- twstg-extraSt
      
      #if root is below minimum weight then entire plant dies and clear out dead material to reset any potential plant age 
      #problems (i.e., including dead material from previous plant)
      
      # distribute the new death into the dead part
      twstd <- twstd + twstdNEW
      twlvd <- twlvd + twlvdNEW
      twrtd <- twrtd + twrtdNEW
      
      # need to update biomass variable  
      biomass <- twlvg + twstg
      fracbiomass <- biomass/maxBiomass
      
  
      if (twstg<(minSize*FracDM_STG)){
        twstg <- 0
      }
      
      if (twlvg<(minSize*FracDM_LVG)){
        twlvg <- 0
      }
      
      if (twrtg<(minSize*FracDM_RTG)){
        twrtg <- 0
      }
      
  
      # need to update biomass var
      biomass <- twlvg+twstg
      fracbiomass <- biomass/maxBiomass
      
      
      #Distribute final biomass values into plants, stems, roots
      plant_dens <- twstg/tillerDensity*maxTillerWeight
      
      if (plant_dens>maxDensity){
        plant_dens <- maxDensity
      }
      
      till_dens <- twstg/(plant_dens*tillerDensity)
      till_dens[is.na(till_dens)] <- 0
      
      if (till_dens<0){
        till_dens <- 0
      }
      
      till_ht <- till_dens*rWeightTillerHeight/100
      
      if (till_ht>maxTillerHeight){
        till_ht <- maxTillerHeight/100
      }
      
      root_lg <- RSratio*till_ht
      
      if (root_lg>maxRootLength){
        root_lg <- maxRootLength
      }
      
      st_dia <- 0
      
      if (till_dens>0){
        st_dia <- 2*(TCA/pi)^0.5/100
      }
      
      
      Density<-till_dens  #Density value 
      Height <- till_ht   #Height value 
      Roots <- root_lg    #Roots value 
      Dia <- st_dia       #Stem Diameter value 
      AGBiomass <- twlvg + twlvd + twstg + twstd
      TOTBiomass <- twlvg + twlvd + twstg + twstd + twrtd + twrtg
    }
    #Daily Results outputted to vectors
    
    dailyleafg[j] <- twlvg
    dailyrootg[j] <- twrtg
    dailystemg[j] <- twstg
    dailyleafd[j] <- twlvd
    dailyrootd[j] <- twrtd
    dailystemd[j] <- twstd
    dailyAGbiomass[j] <- twlvg + twlvd + twstg + twstd
    dailyTOTbiomass [j] <- twlvg + twlvd + twstg + twstd + twrtd + twrtg 
    dailyRTbiomass [j] <- twrtd +twrtg
    # dailyAGbiomass1 [j] <- twlvg + twlvd + twstg + twstd
    
  }
  
  
  #Creates output file
  VegOutput <- cbind(day, Light$meantemp, dailyleafg, dailyrootg, dailystemg, dailyleafd, dailyrootd, dailystemd, dailyAGbiomass, dailyTOTbiomass)
  colnames(VegOutput) <- c("day", "Mean_Temp", "lvsg","rootsg", "stemg", "lvsd", "rootd", "stemd", "AGB", "TOTB"  )
  VegOutput <- data.frame(VegOutput)
  write.csv(VegOutput, file="Veg_output.csv", row.names=FALSE)
  
  
