#What does this script do?
#Mechanistic model of box experiment
#7th iteration


# Set up environment ------------------------------------------------------
rm(list = ls())
Sys.setenv(lang = "en_US")

local({
  Sys.setenv(LANGUAGE="en")
})

## Libraries ---------------------------------------------------------------

#Data processing
library(readxl)
library(readr)
library(dplyr)
library(zoo)
#Fitting likelihood and parallelisation
library(optimParallel)
library(parallel)
library(numDeriv)
library(readr)


## Data --------------------------------------------------------------------
setwd("C:/Users/admin/Desktop/MSThesis/1-Data")

dataExperimentBox_df <- read_delim("./AttemptSides_21-05-24.csv", 
                                   escape_double = FALSE, trim_ws = TRUE)
colnames(dataExperimentBox_df) <- c("id", "phase", "opening", "trial", "innovationL", "innovationP", "attemptL", "attemptP")

#For Canteloup version:
# dataExperimentBox_df <- read_excel("./CanteloupData/CanteloupMechModInput_23-07-24.xlsx")

## Functions ---------------------------------------------------------------

estimateProbability <- function(
    #Now accounts for innovation
  #Differentiates between switching when no lock was available and when one was
  #Added whether the first pulling is in itself the innovation or not (and then the first switch is the innovation)
  data,
  SS,
  SC,
  IS,
  IC,
  E,
  countFirstSwitchAsInnovation = FALSE
){
  if(countFirstSwitchAsInnovation){
    data <- data %>%
      mutate(
        innovationSimple = innovationSimple2
      )
  }
  
  #Help for writing formula  
  # SS
  # IS
  # IC
  # SC
  # E
  # attempt
  # switching
  # innovationSimple
  # innovationComplex
  # previousPull
  # previousLift
  # 
  # ((**())*((1 - )**()))**()
  
  data <- data %>% 
    rowwise() %>% 
    mutate(
      probability = NA,
      probability = ifelse(
        opening != "a" & phase == 1,
        (((1 - IS)**(1 - switching))*(IS**(switching)))**(1 - innovationSimple)*
          (((1 - SS)**(1 - switching))*((SS)**(switching)))**(innovationSimple)*
          (((1 - E)**(attempt))*((E)**(1 - attempt)))**(1 - (1 - innovationSimple)*switching),
        probability
      ),
      probability = ifelse(
        opening == "a" & phase == 1,
        (((1 - IS)**(1 - 1))*(IS**(1)))**(1 - innovationSimple)*
          (((1 - SS)**(1 - 1))*((SS)**(1)))**(innovationSimple)*
          (((1 - E)**(attempt))*((E)**(1 - attempt)))**(1 - (1 - innovationSimple)*1) +
          
          
          (((1 - IS)**(1 - 0))*(IS**(0)))**(1 - innovationSimple)*
          (((1 - SS)**(1 - 0))*((SS)**(0)))**(innovationSimple)*
          (((1 - E)**(attempt))*((E)**(1 - attempt)))**(1 - (1 - innovationSimple)*0),
        probability
      ),
      
      
      
      #Phase 2
      probability = ifelse(
        opening != "a" & phase == 2,
        ((SC**(switching))*((1 - SC)**(1 - switching)))**(innovationSimple)*
          ((IS**(switching))*((1 - IS)**(1 - switching)))**(1 - innovationSimple)*
          ((0**(attempt))*(1**(1 - attempt)))**((1-previousPull)*(1 - innovationSimple)*switching)*
          (
            ((E**(1 - attempt))*((1 - E)**(attempt)))**(1 - (1 - innovationL)*max(previousPull*switching, (1 - previousPull)*(1 - switching)))*
              ((IC**(1 - attempt))*((1 - IC)**(attempt)))**((1 - innovationL)*max(previousPull*switching, (1 - previousPull)*(1 - switching)))
          )**(1 - ((1-previousPull)*(1 - innovationSimple)*switching)),
        probability
      ),
      probability = ifelse(
        opening == "a" & phase == 2,
        
        ((SC**(0))*((1 - SC)**(1 - 0)))**(innovationSimple)*
          ((IS**(0))*((1 - IS)**(1 - 0)))**(1 - innovationSimple)*
          ((0**(attempt))*(1**(1 - attempt)))**((1-previousPull)*(1 - innovationSimple)*0)*
          (
            ((E**(1 - attempt))*((1 - E)**(attempt)))**(1 - (1 - innovationL)*max(previousPull*0, (1 - previousPull)*(1 - 0)))*
              ((IC**(1 - attempt))*((1 - IC)**(attempt)))**((1 - innovationL)*max(previousPull*0, (1 - previousPull)*(1 - 0)))
          )**(1 - ((1-previousPull)*(1 - innovationSimple)*0)) +
          
          
          ((SC**(1))*((1 - SC)**(1 - 1)))**(innovationSimple)*
          ((IS**(1))*((1 - IS)**(1 - 1)))**(1 - innovationSimple)*
          ((0**(attempt))*(1**(1 - attempt)))**((1-previousPull)*(1 - innovationSimple)*1)*
          (
            ((E**(1 - attempt))*((1 - E)**(attempt)))**(1 - (1 - innovationL)*max(previousPull*1, (1 - previousPull)*(1 - 1)))*
              ((IC**(1 - attempt))*((1 - IC)**(attempt)))**((1 - innovationL)*max(previousPull*1, (1 - previousPull)*(1 - 1)))
          )**(1 - ((1-previousPull)*(1 - innovationSimple)*1)),
        probability
      ),
      #Phase 3
      probability = ifelse(
        opening != "a" & phase == 3,
        ((SC**(switching))*((1 - SC)**(1 - switching)))**(innovationSimple)*
          ((IS**(switching))*((1 - IS)**(1 - switching)))**(1 - innovationSimple)*
          ((E**(1 - attempt))*((1 - E)**(attempt)))**(max(c(innovationL*switching*(1 - previousLift),innovationL*(1 - switching)*previousLift,
                                                            innovationP * switching * previousLift, innovationP * (1 - switching) * previousPull)))*
          ((IC**(1 - attempt))*((1 - IC)**(attempt)))**(1 - max(c(innovationL*switching*(1 - previousLift),innovationL*(1 - switching)*previousLift,
                                                                  innovationP * switching * previousLift, innovationP * (1 - switching) * previousPull))),
        probability
      ),
      probability = ifelse(
        opening == "a" & phase == 3,
        ((SC**(0))*((1 - SC)**(1 - 0)))**(innovationSimple)*
          ((IS**(0))*((1 - IS)**(1 - 0)))**(1 - innovationSimple)*
          ((E**(1 - attempt))*((1 - E)**(attempt)))**(max(c(innovationL*0*(1 - previousLift),innovationL*(1 - 0)*previousLift,
                                                            innovationP * 0 * previousLift, innovationP * (1 - 0) * previousPull)))*
          ((IC**(1 - attempt))*((1 - IC)**(attempt)))**(1 - max(c(innovationL*0*(1 - previousLift),innovationL*(1 - 0)*previousLift,
                                                                  innovationP * 0 * previousLift, innovationP * (1 - 0) * previousPull))) +
          
          ((SC**(1))*((1 - SC)**(1 - 1)))**(innovationSimple)*
          ((IS**(1))*((1 - IS)**(1 - 1)))**(1 - innovationSimple)*
          ((E**(1 - attempt))*((1 - E)**(attempt)))**(max(c(innovationL*1*(1 - previousLift),innovationL*(1 - 1)*previousLift,
                                                            innovationP * 1 * previousLift, innovationP * (1 - 1) * previousPull)))*
          ((IC**(1 - attempt))*((1 - IC)**(attempt)))**(1 - max(c(innovationL*1*(1 - previousLift),innovationL*(1 - 1)*previousLift,
                                                                  innovationP * 1 * previousLift, innovationP * (1 - 1) * previousPull))),
        probability
      )
    )
  
  return(data$probability)
}


#Added ES vs EC compared to v1
estimateProbability2 <- function(
    #Now accounts for innovation
  #Differentiates between switching when no lock was available and when one was
  #Added whether the first pulling is in itself the innovation or not (and then the first switch is the innovation)
  data,
  SS,
  SC,
  IS,
  IC,
  ES,
  EC,
  countFirstSwitchAsInnovation = FALSE
){
  if(countFirstSwitchAsInnovation){
    data <- data %>%
      mutate(
        innovationSimple = innovationSimple2
      )
  }
  
  #Help for writing formula  
  # SS
  # IS
  # IC
  # SC
  # E
  # attempt
  # switching
  # innovationSimple
  # innovationComplex
  # previousPull
  # previousLift
  # 
  # ((**())*((1 - )**()))**()
  
  data <- data %>% 
    rowwise() %>% 
    mutate(
      probability = NA,
      probability = ifelse(
        opening != "a" & phase == 1,
        (((1 - IS)**(1 - switching))*(IS**(switching)))**(1 - innovationSimple)*
          (((1 - SS)**(1 - switching))*((SS)**(switching)))**(innovationSimple)*
          (((1 - ES)**(attempt))*((ES)**(1 - attempt)))**(1 - (1 - innovationSimple)*switching),
        probability
      ),
      probability = ifelse(
        opening == "a" & phase == 1,
        (((1 - IS)**(1 - 1))*(IS**(1)))**(1 - innovationSimple)*
          (((1 - SS)**(1 - 1))*((SS)**(1)))**(innovationSimple)*
          (((1 - SE)**(attempt))*((ES)**(1 - attempt)))**(1 - (1 - innovationSimple)*1) +
          (((1 - IS)**(1 - 0))*(IS**(0)))**(1 - innovationSimple)*
          (((1 - SS)**(1 - 0))*((SS)**(0)))**(innovationSimple)*
          (((1 - ES)**(attempt))*((ES)**(1 - attempt)))**(1 - (1 - innovationSimple)*0),
        probability
      ),
      #Phase 2
      probability = ifelse(
        opening != "a" & phase == 2,
        ((SC**(switching))*((1 - SC)**(1 - switching)))**(innovationSimple)*
          ((IS**(switching))*((1 - IS)**(1 - switching)))**(1 - innovationSimple)*
          ((0**(attempt))*(1**(1 - attempt)))**((1-previousPull)*(1 - innovationSimple)*switching)*
          (
            ((ifelse(previousPull*(1-switching) == 1 | previousLift*switching == 1, ES, EC)**(1 - attempt))*((1 - ifelse(previousPull*(1-switching) == 1 | previousLift*switching == 1, ES, EC))**(attempt)))**(1 - (1 - innovationL)*max(previousPull*switching, (1 - previousPull)*(1 - switching)))*
              ((IC**(1 - attempt))*((1 - IC)**(attempt)))**((1 - innovationL)*max(previousPull*switching, (1 - previousPull)*(1 - switching)))
          )**(1 - ((1-previousPull)*(1 - innovationSimple)*switching)),
        probability
      ),
      probability = ifelse(
        opening == "a" & phase == 2,
        ((SC**(0))*((1 - SC)**(1 - 0)))**(innovationSimple)*
          ((IS**(0))*((1 - IS)**(1 - 0)))**(1 - innovationSimple)*
          ((0**(attempt))*(1**(1 - attempt)))**((1-previousPull)*(1 - innovationSimple)*0)*
          (
            ((ifelse(previousPull*(1-0) == 1 | previousLift*0 == 1, ES, EC)**(1 - attempt))*((1 - ifelse(previousPull*(1-0) == 1 | previousLift*0 == 1, ES, EC))**(attempt)))**(1 - (1 - innovationL)*max(previousPull*0, (1 - previousPull)*(1 - 0)))*
              ((IC**(1 - attempt))*((1 - IC)**(attempt)))**((1 - innovationL)*max(previousPull*0, (1 - previousPull)*(1 - 0)))
          )**(1 - ((1-previousPull)*(1 - innovationSimple)*0)) +
          ((SC**(1))*((1 - SC)**(1 - 1)))**(innovationSimple)*
          ((IS**(1))*((1 - IS)**(1 - 1)))**(1 - innovationSimple)*
          ((0**(attempt))*(1**(1 - attempt)))**((1-previousPull)*(1 - innovationSimple)*1)*
          (
            ((ifelse(previousPull*(1-1) == 1 | previousLift*1 == 1, ES, EC)**(1 - attempt))*((1 - ifelse(previousPull*(1-1) == 1 | previousLift*1 == 1, ES, EC))**(attempt)))**(1 - (1 - innovationL)*max(previousPull*1, (1 - previousPull)*(1 - 1)))*
              ((IC**(1 - attempt))*((1 - IC)**(attempt)))**((1 - innovationL)*max(previousPull*1, (1 - previousPull)*(1 - 1)))
          )**(1 - ((1-previousPull)*(1 - innovationSimple)*1)),
        probability
      ),
      
      #Phase 3
      probability = ifelse(
        opening != "a" & phase == 3,
        ((SC**(switching))*((1 - SC)**(1 - switching)))**(innovationSimple)*
          ((IS**(switching))*((1 - IS)**(1 - switching)))**(1 - innovationSimple)*
          ((EC**(1 - attempt))*((1 - EC)**(attempt)))**(max(c(innovationL*switching*(1 - previousLift),innovationL*(1 - switching)*previousLift,
                                                              innovationP * switching * previousLift, innovationP * (1 - switching) * previousPull)))*
          ((IC**(1 - attempt))*((1 - IC)**(attempt)))**(1 - max(c(innovationL*switching*(1 - previousLift),innovationL*(1 - switching)*previousLift,
                                                                  innovationP * switching * previousLift, innovationP * (1 - switching) * previousPull))),
        probability
      ),
      probability = ifelse(
        opening == "a" & phase == 3,
        ((SC**(0))*((1 - SC)**(1 - 0)))**(innovationSimple)*
          ((IS**(0))*((1 - IS)**(1 - 0)))**(1 - innovationSimple)*
          ((EC**(1 - attempt))*((1 - EC)**(attempt)))**(max(c(innovationL*0*(1 - previousLift),innovationL*(1 - 0)*previousLift,
                                                              innovationP * 0 * previousLift, innovationP * (1 - 0) * previousPull)))*
          ((IC**(1 - attempt))*((1 - IC)**(attempt)))**(1 - max(c(innovationL*0*(1 - previousLift),innovationL*(1 - 0)*previousLift,
                                                                  innovationP * 0 * previousLift, innovationP * (1 - 0) * previousPull))) +
          ((SC**(1))*((1 - SC)**(1 - 1)))**(innovationSimple)*
          ((IS**(1))*((1 - IS)**(1 - 1)))**(1 - innovationSimple)*
          ((EC**(1 - attempt))*((1 - EC)**(attempt)))**(max(c(innovationL*1*(1 - previousLift),innovationL*(1 - 1)*previousLift,
                                                              innovationP * 1 * previousLift, innovationP * (1 - 1) * previousPull)))*
          ((IC**(1 - attempt))*((1 - IC)**(attempt)))**(1 - max(c(innovationL*1*(1 - previousLift),innovationL*(1 - 1)*previousLift,
                                                                  innovationP * 1 * previousLift, innovationP * (1 - 1) * previousPull))),
        probability
      )
    )
  
  return(data$probability)
}

#Added to v2 the SC should only be when one innovation complex occurred
estimateProbability3 <- function(
    #Now accounts for innovation
  #Differentiates between switching when no lock was available and when one was
  #Added whether the first pulling is in itself the innovation or not (and then the first switch is the innovation)
  data,
  SS,
  SC,
  IS,
  IC,
  ES,
  EC,
  countFirstSwitchAsInnovation = FALSE
){
  if(countFirstSwitchAsInnovation){
    data <- data %>%
      mutate(
        innovationSimple = innovationSimple2
      )
  }
  
  #Help for writing formula  
  # SS
  # IS
  # IC
  # SC
  # E
  # attempt
  # switching
  # innovationSimple
  # innovationComplex
  # previousPull
  # previousLift
  # 
  # ((**())*((1 - )**()))**()
  
  data <- data %>% 
    rowwise() %>% 
    mutate(
      probability = NA,
      probability = ifelse(
        opening != "a" & phase == 1,
        (((1 - IS)**(1 - switching))*(IS**(switching)))**(1 - innovationSimple)*
          (((1 - SS)**(1 - switching))*((SS)**(switching)))**(innovationSimple)*
          (((1 - ES)**(attempt))*((ES)**(1 - attempt)))**(1 - (1 - innovationSimple)*switching),
        probability
      ),
      probability = ifelse(
        opening == "a" & phase == 1,
        (((1 - IS)**(1 - 1))*(IS**(1)))**(1 - innovationSimple)*
          (((1 - SS)**(1 - 1))*((SS)**(1)))**(innovationSimple)*
          (((1 - SE)**(attempt))*((ES)**(1 - attempt)))**(1 - (1 - innovationSimple)*1) +
          (((1 - IS)**(1 - 0))*(IS**(0)))**(1 - innovationSimple)*
          (((1 - SS)**(1 - 0))*((SS)**(0)))**(innovationSimple)*
          (((1 - ES)**(attempt))*((ES)**(1 - attempt)))**(1 - (1 - innovationSimple)*0),
        probability
      ),
      #Phase 2
      probability = ifelse(
        opening != "a" & phase == 2,
        ((ifelse(innovationL == 1 | innovationP == 1, SC, SS)**(switching))*((1 - ifelse(innovationL == 1 | innovationP == 1, SC, SS))**(1 - switching)))**(innovationSimple)*
          ((IS**(switching))*((1 - IS)**(1 - switching)))**(1 - innovationSimple)*
          ((0**(attempt))*(1**(1 - attempt)))**((1-previousPull)*(1 - innovationSimple)*switching)*
          (
            ((ifelse(previousPull*(1-switching) == 1 | previousLift*switching == 1, ES, EC)**(1 - attempt))*((1 - ifelse(previousPull*(1-switching) == 1 | previousLift*switching == 1, ES, EC))**(attempt)))**(1 - (1 - innovationL)*max(previousPull*switching, (1 - previousPull)*(1 - switching)))*
              ((IC**(1 - attempt))*((1 - IC)**(attempt)))**((1 - innovationL)*max(previousPull*switching, (1 - previousPull)*(1 - switching)))
          )**(1 - ((1-previousPull)*(1 - innovationSimple)*switching)),
        probability
      ),
      probability = ifelse(
        opening == "a" & phase == 2,
        ((ifelse(innovationL == 1 | innovationP == 1, SC, SS)**(0))*((1 - ifelse(innovationL == 1 | innovationP == 1, SC, SS))**(1 - 0)))**(innovationSimple)*
          ((IS**(0))*((1 - IS)**(1 - 0)))**(1 - innovationSimple)*
          ((0**(attempt))*(1**(1 - attempt)))**((1-previousPull)*(1 - innovationSimple)*0)*
          (
            ((ifelse(previousPull*(1-0) == 1 | previousLift*0 == 1, ES, EC)**(1 - attempt))*((1 - ifelse(previousPull*(1-0) == 1 | previousLift*0 == 1, ES, EC))**(attempt)))**(1 - (1 - innovationL)*max(previousPull*0, (1 - previousPull)*(1 - 0)))*
              ((IC**(1 - attempt))*((1 - IC)**(attempt)))**((1 - innovationL)*max(previousPull*0, (1 - previousPull)*(1 - 0)))
          )**(1 - ((1-previousPull)*(1 - innovationSimple)*0)) +
          ((ifelse(innovationL == 1 | innovationP == 1, SC, SS)**(1))*((1 - ifelse(innovationL == 1 | innovationP == 1, SC, SS))**(1 - 1)))**(innovationSimple)*
          ((IS**(1))*((1 - IS)**(1 - 1)))**(1 - innovationSimple)*
          ((0**(attempt))*(1**(1 - attempt)))**((1-previousPull)*(1 - innovationSimple)*1)*
          (
            ((ifelse(previousPull*(1-1) == 1 | previousLift*1 == 1, ES, EC)**(1 - attempt))*((1 - ifelse(previousPull*(1-1) == 1 | previousLift*1 == 1, ES, EC))**(attempt)))**(1 - (1 - innovationL)*max(previousPull*1, (1 - previousPull)*(1 - 1)))*
              ((IC**(1 - attempt))*((1 - IC)**(attempt)))**((1 - innovationL)*max(previousPull*1, (1 - previousPull)*(1 - 1)))
          )**(1 - ((1-previousPull)*(1 - innovationSimple)*1)),
        probability
      ),
      
      #Phase 3
      probability = ifelse(
        opening != "a" & phase == 3,
        ((ifelse(innovationL == 1 | innovationP == 1, SC, SS)**(switching))*((1 - ifelse(innovationL == 1 | innovationP == 1, SC, SS))**(1 - switching)))**(innovationSimple)*
          ((IS**(switching))*((1 - IS)**(1 - switching)))**(1 - innovationSimple)*
          ((EC**(1 - attempt))*((1 - EC)**(attempt)))**(max(c(innovationL*switching*(1 - previousLift),innovationL*(1 - switching)*previousLift,
                                                              innovationP * switching * previousLift, innovationP * (1 - switching) * previousPull)))*
          ((IC**(1 - attempt))*((1 - IC)**(attempt)))**(1 - max(c(innovationL*switching*(1 - previousLift),innovationL*(1 - switching)*previousLift,
                                                                  innovationP * switching * previousLift, innovationP * (1 - switching) * previousPull))),
        probability
      ),
      probability = ifelse(
        opening == "a" & phase == 3,
        ((ifelse(innovationL == 1 | innovationP == 1, SC, SS)**(0))*((1 - ifelse(innovationL == 1 | innovationP == 1, SC, SS))**(1 - 0)))**(innovationSimple)*
          ((IS**(0))*((1 - IS)**(1 - 0)))**(1 - innovationSimple)*
          ((EC**(1 - attempt))*((1 - EC)**(attempt)))**(max(c(innovationL*0*(1 - previousLift),innovationL*(1 - 0)*previousLift,
                                                              innovationP * 0 * previousLift, innovationP * (1 - 0) * previousPull)))*
          ((IC**(1 - attempt))*((1 - IC)**(attempt)))**(1 - max(c(innovationL*0*(1 - previousLift),innovationL*(1 - 0)*previousLift,
                                                                  innovationP * 0 * previousLift, innovationP * (1 - 0) * previousPull))) +
          ((ifelse(innovationL == 1 | innovationP == 1, SC, SS)**(1))*((1 - ifelse(innovationL == 1 | innovationP == 1, SC, SS))**(1 - 1)))**(innovationSimple)*
          ((IS**(1))*((1 - IS)**(1 - 1)))**(1 - innovationSimple)*
          ((EC**(1 - attempt))*((1 - EC)**(attempt)))**(max(c(innovationL*1*(1 - previousLift),innovationL*(1 - 1)*previousLift,
                                                              innovationP * 1 * previousLift, innovationP * (1 - 1) * previousPull)))*
          ((IC**(1 - attempt))*((1 - IC)**(attempt)))**(1 - max(c(innovationL*1*(1 - previousLift),innovationL*(1 - 1)*previousLift,
                                                                  innovationP * 1 * previousLift, innovationP * (1 - 1) * previousPull))),
        probability
      )
    )
  
  return(data$probability)
}

#Added to v1 the information, in attempt, of the attempted solution (lift and/or pull)
estimateProbability4 <- function(
    #Now accounts for innovation
  #Differentiates between switching when no lock was available and when one was
  #Added whether the first pulling is in itself the innovation or not (and then the first switch is the innovation)
  data,
  SS,
  SC,
  IS,
  IC,
  E,
  countFirstSwitchAsInnovation = FALSE
){
  if(countFirstSwitchAsInnovation){
    data <- data %>%
      mutate(
        innovationSimple = innovationSimple2
      )
  }
  
  #Help for writing formula  
  # SS
  # IS
  # IC
  # SC
  # E
  # attempt
  # switching
  # innovationSimple
  # innovationComplex
  # previousPull
  # previousLift
  # 
  # ((**())*((1 - )**()))**()
  
  data <- data %>% 
    rowwise() %>% 
    mutate(
      probability = NA,
      probability = ifelse(
        opening != "a" & phase == 1,
        (((1 - IS)**(1 - switching))*(IS**(switching)))**(1 - innovationSimple)*
          (((1 - SS)**(1 - switching))*((SS)**(switching)))**(innovationSimple)*
          (((1 - E)**(attempt))*((E)**(1 - attempt)))**(1 - (1 - innovationSimple)*switching),
        probability
      ),
      #Case when attempted, knowing it attempted only on one side (or no knowledge of which side)
      probability = ifelse(
        opening == "a" & phase == 1 & attemptBoth != 1,
        
        ((((1 - IS)**(1 - 1))*(IS**(1)))**(1 - innovationSimple)*
           (((1 - SS)**(1 - 1))*((SS)**(1)))**(innovationSimple)*
           (((1 - E)**(attempt))*((E)**(1 - attempt)))**(1 - (1 - innovationSimple)*1))**max(previousPull*attemptL, previousLift*attemptP) +
          
          ((((1 - IS)**(1 - 0))*(IS**(0)))**(1 - innovationSimple)*
             (((1 - SS)**(1 - 0))*((SS)**(0)))**(innovationSimple)*
             (((1 - E)**(attempt))*((E)**(1 - attempt)))**(1 - (1 - innovationSimple)*0))**max(previousLift*attemptL, previousPull*attemptP),
        probability
      ),
      #Case when attempted, knowing it attempted on the two sides
      probability = ifelse(
        opening == "a" & phase == 1 & attemptBoth == 1,
        
        ((((1 - IS)**(1 - 1))*(IS**(1)))**(1 - innovationSimple)*
           (((1 - SS)**(1 - 1))*((SS)**(1)))**(innovationSimple)*
           (((1 - E)**(attempt))*((E)**(1 - attempt)))**(1 - (1 - innovationSimple)*1))**max(previousPull*attemptL, previousLift*attemptP) *#Here it is now multiplied
          
          ((((1 - IS)**(1 - 0))*(IS**(0)))**(1 - innovationSimple)*
             (((1 - SS)**(1 - 0))*((SS)**(0)))**(innovationSimple)*
             (((1 - E)**(attempt))*((E)**(1 - attempt)))**(1 - (1 - innovationSimple)*0))**max(previousLift*attemptL, previousPull*attemptP),
        probability
      ),
      
      #Phase 2
      probability = ifelse(
        opening != "a" & phase == 2,
        ((SC**(switching))*((1 - SC)**(1 - switching)))**(innovationSimple)*
          ((IS**(switching))*((1 - IS)**(1 - switching)))**(1 - innovationSimple)*
          ((0**(attempt))*(1**(1 - attempt)))**((1-previousPull)*(1 - innovationSimple)*switching)*
          (
            ((E**(1 - attempt))*((1 - E)**(attempt)))**(1 - (1 - innovationL)*max(previousPull*switching, (1 - previousPull)*(1 - switching)))*
              ((IC**(1 - attempt))*((1 - IC)**(attempt)))**((1 - innovationL)*max(previousPull*switching, (1 - previousPull)*(1 - switching)))
          )**(1 - ((1-previousPull)*(1 - innovationSimple)*switching)),
        probability
      ),
      #Case when attempted, knowing it attempted only on one side (or no knowledge of which side)
      probability = ifelse(
        opening == "a" & phase == 2 & attemptBoth != 1,
        
        (((SC**(0))*((1 - SC)**(1 - 0)))**(innovationSimple)*
           ((IS**(0))*((1 - IS)**(1 - 0)))**(1 - innovationSimple)*
           ((0**(attempt))*(1**(1 - attempt)))**((1-previousPull)*(1 - innovationSimple)*0)*
           (
             ((E**(1 - attempt))*((1 - E)**(attempt)))**(1 - (1 - innovationL)*max(previousPull*0, (1 - previousPull)*(1 - 0)))*
               ((IC**(1 - attempt))*((1 - IC)**(attempt)))**((1 - innovationL)*max(previousPull*0, (1 - previousPull)*(1 - 0)))
           )**(1 - ((1-previousPull)*(1 - innovationSimple)*0)))**max(previousLift*attemptL, previousPull*attemptP) +
          
          (((SC**(1))*((1 - SC)**(1 - 1)))**(innovationSimple)*
             ((IS**(1))*((1 - IS)**(1 - 1)))**(1 - innovationSimple)*
             ((0**(attempt))*(1**(1 - attempt)))**((1-previousPull)*(1 - innovationSimple)*1)*
             (
               ((E**(1 - attempt))*((1 - E)**(attempt)))**(1 - (1 - innovationL)*max(previousPull*1, (1 - previousPull)*(1 - 1)))*
                 ((IC**(1 - attempt))*((1 - IC)**(attempt)))**((1 - innovationL)*max(previousPull*1, (1 - previousPull)*(1 - 1)))
             )**(1 - ((1-previousPull)*(1 - innovationSimple)*1)))**max(previousPull*attemptL, previousLift*attemptP),
        probability
      ),
      #Case when attempted, knowing it attempted on the two sides
      probability = ifelse(
        opening == "a" & phase == 2 & attemptBoth == 1,
        
        (((SC**(0))*((1 - SC)**(1 - 0)))**(innovationSimple)*
           ((IS**(0))*((1 - IS)**(1 - 0)))**(1 - innovationSimple)*
           ((0**(attempt))*(1**(1 - attempt)))**((1-previousPull)*(1 - innovationSimple)*0)*
           (
             ((E**(1 - attempt))*((1 - E)**(attempt)))**(1 - (1 - innovationL)*max(previousPull*0, (1 - previousPull)*(1 - 0)))*
               ((IC**(1 - attempt))*((1 - IC)**(attempt)))**((1 - innovationL)*max(previousPull*0, (1 - previousPull)*(1 - 0)))
           )**(1 - ((1-previousPull)*(1 - innovationSimple)*0)))**max(previousLift*attemptL, previousPull*attemptP) *#Here it is now multiplied
          
          (((SC**(1))*((1 - SC)**(1 - 1)))**(innovationSimple)*
             ((IS**(1))*((1 - IS)**(1 - 1)))**(1 - innovationSimple)*
             ((0**(attempt))*(1**(1 - attempt)))**((1-previousPull)*(1 - innovationSimple)*1)*
             (
               ((E**(1 - attempt))*((1 - E)**(attempt)))**(1 - (1 - innovationL)*max(previousPull*1, (1 - previousPull)*(1 - 1)))*
                 ((IC**(1 - attempt))*((1 - IC)**(attempt)))**((1 - innovationL)*max(previousPull*1, (1 - previousPull)*(1 - 1)))
             )**(1 - ((1-previousPull)*(1 - innovationSimple)*1)))**max(previousPull*attemptL, previousLift*attemptP),
        probability
      ),
      
      #Phase 3
      probability = ifelse(
        opening != "a" & phase == 3,
        ((SC**(switching))*((1 - SC)**(1 - switching)))**(innovationSimple)*
          ((IS**(switching))*((1 - IS)**(1 - switching)))**(1 - innovationSimple)*
          ((E**(1 - attempt))*((1 - E)**(attempt)))**(max(c(innovationL*switching*(1 - previousLift),innovationL*(1 - switching)*previousLift,
                                                            innovationP * switching * previousLift, innovationP * (1 - switching) * previousPull)))*
          ((IC**(1 - attempt))*((1 - IC)**(attempt)))**(1 - max(c(innovationL*switching*(1 - previousLift),innovationL*(1 - switching)*previousLift,
                                                                  innovationP * switching * previousLift, innovationP * (1 - switching) * previousPull))),
        probability
      ),
      #Case when attempted, knowing it attempted only on one side (or no knowledge of which side)
      probability = ifelse(
        opening == "a" & phase == 3 & attemptBoth != 1,
        
        (((SC**(0))*((1 - SC)**(1 - 0)))**(innovationSimple)*
           ((IS**(0))*((1 - IS)**(1 - 0)))**(1 - innovationSimple)*
           ((E**(1 - attempt))*((1 - E)**(attempt)))**(max(c(innovationL*0*(1 - previousLift),innovationL*(1 - 0)*previousLift,
                                                             innovationP * 0 * previousLift, innovationP * (1 - 0) * previousPull)))*
           ((IC**(1 - attempt))*((1 - IC)**(attempt)))**(1 - max(c(innovationL*0*(1 - previousLift),innovationL*(1 - 0)*previousLift,
                                                                   innovationP * 0 * previousLift, innovationP * (1 - 0) * previousPull))))**max(previousLift*attemptL, previousPull*attemptP) +
          
          (((SC**(1))*((1 - SC)**(1 - 1)))**(innovationSimple)*
             ((IS**(1))*((1 - IS)**(1 - 1)))**(1 - innovationSimple)*
             ((E**(1 - attempt))*((1 - E)**(attempt)))**(max(c(innovationL*1*(1 - previousLift),innovationL*(1 - 1)*previousLift,
                                                               innovationP * 1 * previousLift, innovationP * (1 - 1) * previousPull)))*
             ((IC**(1 - attempt))*((1 - IC)**(attempt)))**(1 - max(c(innovationL*1*(1 - previousLift),innovationL*(1 - 1)*previousLift,
                                                                     innovationP * 1 * previousLift, innovationP * (1 - 1) * previousPull))))**max(previousPull*attemptL, previousLift*attemptP),
        probability
      ),
      #Case when attempted, knowing it attempted on the two sides
      probability = ifelse(
        opening == "a" & phase == 3 & attemptBoth == 1,
        
        (((SC**(0))*((1 - SC)**(1 - 0)))**(innovationSimple)*
           ((IS**(0))*((1 - IS)**(1 - 0)))**(1 - innovationSimple)*
           ((E**(1 - attempt))*((1 - E)**(attempt)))**(max(c(innovationL*0*(1 - previousLift),innovationL*(1 - 0)*previousLift,
                                                             innovationP * 0 * previousLift, innovationP * (1 - 0) * previousPull)))*
           ((IC**(1 - attempt))*((1 - IC)**(attempt)))**(1 - max(c(innovationL*0*(1 - previousLift),innovationL*(1 - 0)*previousLift,
                                                                   innovationP * 0 * previousLift, innovationP * (1 - 0) * previousPull))))**max(previousLift*attemptL, previousPull*attemptP) *#Here it is now multiplied
          
          (((SC**(1))*((1 - SC)**(1 - 1)))**(innovationSimple)*
             ((IS**(1))*((1 - IS)**(1 - 1)))**(1 - innovationSimple)*
             ((E**(1 - attempt))*((1 - E)**(attempt)))**(max(c(innovationL*1*(1 - previousLift),innovationL*(1 - 1)*previousLift,
                                                               innovationP * 1 * previousLift, innovationP * (1 - 1) * previousPull)))*
             ((IC**(1 - attempt))*((1 - IC)**(attempt)))**(1 - max(c(innovationL*1*(1 - previousLift),innovationL*(1 - 1)*previousLift,
                                                                     innovationP * 1 * previousLift, innovationP * (1 - 1) * previousPull))))**max(previousPull*attemptL, previousLift*attemptP),
        probability
      ),
    )
  
  return(data$probability)
}


#Similar to v4, however, correcting for the error in attempts: in phase 1 it was always as a function of e/1-e, which should not be the case, and then it could include probability 0 for a fail attempt.
#Now it consider Innovation for the first switch simple, and innovation if this is then successful.
estimateProbability5 <- function(
    #Now accounts for innovation
  #Differentiates between switching when no lock was available and when one was
  #Added whether the first pulling is in itself the innovation or not (and then the first switch is the innovation)
  data,
  SS,
  SC,
  IS,
  IC,
  E,
  countFirstSwitchAsInnovation = FALSE
){
  if(countFirstSwitchAsInnovation){
    data <- data %>%
      mutate(
        innovationSimple = innovationSimple2
      )
  }
  
  #Help for writing formula  
  # SS
  # IS
  # IC
  # SC
  # E
  # attempt
  # switching
  # innovationSimple
  # innovationComplex
  # previousPull
  # previousLift
  # 
  # ((**())*((1 - )**()))**()
  
  data <- data %>% 
    rowwise() %>% 
    mutate(
      probability = NA,
      probability = ifelse(
        opening != "a" & phase == 1,
        (((1 - ifelse(hasattemptBothAlready == 1, SS, IS))**(1 - switching))*(ifelse(hasattemptBothAlready == 1, SS, IS)**(switching)))**(1 - innovationSimple)*
          (((1 - SS)**(1 - switching))*((SS)**(switching)))**(innovationSimple)*
          (((1 - E)**(attempt))*((E)**(1 - attempt)))**(1 - (1 - innovationSimple)*switching)*
          (((1 - IS)**(attempt))*((IS)**(1 - attempt)))**((1 - innovationSimple)*switching),
        probability
      ),
      #Case when attempted, knowing it attempted only on one side (or no knowledge of which side)
      probability = ifelse(
        opening == "a" & phase == 1 & attemptBoth != 1,
        ((((1 - ifelse(hasattemptBothAlready == 1, SS, IS))**(1 - 1))*(ifelse(hasattemptBothAlready == 1, SS, IS)**(1)))**(1 - innovationSimple)*
           (((1 - SS)**(1 - 1))*((SS)**(1)))**(innovationSimple)*
           (((1 - E)**(attempt))*((E)**(1 - attempt)))**(1 - (1 - innovationSimple)*1)*
           (((1 - IS)**(attempt))*((IS)**(1 - attempt)))**((1 - innovationSimple)*1))**max(previousPull*attemptL, previousLift*attemptP) +
          
          ((((1 - ifelse(hasattemptBothAlready == 1, SS, IS))**(1 - 0))*(ifelse(hasattemptBothAlready == 1, SS, IS)**(0)))**(1 - innovationSimple)*
             (((1 - SS)**(1 - 0))*((SS)**(0)))**(innovationSimple)*
             (((1 - E)**(attempt))*((E)**(1 - attempt)))**(1 - (1 - innovationSimple)*0)*
             (((1 - IS)**(attempt))*((IS)**(1 - attempt)))**((1 - innovationSimple)*0))**max(previousLift*attemptL, previousPull*attemptP),
        probability
      ),
      #Case when attempted, knowing it attempted on the two sides
      probability = ifelse(
        opening == "a" & phase == 1 & attemptBoth == 1,
        
        ((((1 - ifelse(hasattemptBothAlready == 1, SS, IS))**(1 - 1))*(ifelse(hasattemptBothAlready == 1, SS, IS)**(1)))**(1 - innovationSimple)*
           (((1 - SS)**(1 - 1))*((SS)**(1)))**(innovationSimple)*
           (((1 - E)**(attempt))*((E)**(1 - attempt)))**(1 - (1 - innovationSimple)*1)*
           (((1 - IS)**(attempt))*((IS)**(1 - attempt)))**((1 - innovationSimple)*1))**max(previousPull*attemptL, previousLift*attemptP) *#Here it is now multiplied
          
          ((((1 - ifelse(hasattemptBothAlready == 1, SS, IS))**(1 - 0))*(ifelse(hasattemptBothAlready == 1, SS, IS)**(0)))**(1 - innovationSimple)*
             (((1 - SS)**(1 - 0))*((SS)**(0)))**(innovationSimple)*
             (((1 - E)**(attempt))*((E)**(1 - attempt)))**(1 - (1 - innovationSimple)*0)*
             (((1 - IS)**(attempt))*((IS)**(1 - attempt)))**((1 - innovationSimple)*0))**max(previousLift*attemptL, previousPull*attemptP),
        probability
      ),
      
      #Phase 2
      probability = ifelse(
        opening != "a" & phase == 2,
        ((SC**(switching))*((1 - SC)**(1 - switching)))**(innovationSimple)*
          ((ifelse(hasattemptBothAlready == 1, SS, IS)**(switching))*((1 - ifelse(hasattemptBothAlready == 1, SS, IS))**(1 - switching)))**(1 - innovationSimple)*
          (((1 - IS)**(attempt))*(IS**(1 - attempt)))**((1-previousPull)*(1 - innovationSimple)*switching)*
          (
            ((E**(1 - attempt))*((1 - E)**(attempt)))**(1 - (1 - innovationL)*max(previousPull*switching, (1 - previousPull)*(1 - switching)))*
              ((IC**(1 - attempt))*((1 - IC)**(attempt)))**((1 - innovationL)*max(previousPull*switching, (1 - previousPull)*(1 - switching)))
          )**(1 - ((1-previousPull)*(1 - innovationSimple)*switching)),
        probability
      ),
      #Case when attempted, knowing it attempted only on one side (or no knowledge of which side)
      probability = ifelse(
        opening == "a" & phase == 2 & attemptBoth != 1,
        
        (((SC**(0))*((1 - SC)**(1 - 0)))**(innovationSimple)*
           ((ifelse(hasattemptBothAlready == 1, SS, IS)**(0))*((1 - ifelse(hasattemptBothAlready == 1, SS, IS))**(1 - 0)))**(1 - innovationSimple)*
           (((1 - IS)**(attempt))*(IS**(1 - attempt)))**((1-previousPull)*(1 - innovationSimple)*0)*
           (
             ((E**(1 - attempt))*((1 - E)**(attempt)))**(1 - (1 - innovationL)*max(previousPull*0, (1 - previousPull)*(1 - 0)))*
               ((IC**(1 - attempt))*((1 - IC)**(attempt)))**((1 - innovationL)*max(previousPull*0, (1 - previousPull)*(1 - 0)))
           )**(1 - ((1-previousPull)*(1 - innovationSimple)*0)))**max(previousLift*attemptL, previousPull*attemptP) +
          
          (((SC**(1))*((1 - SC)**(1 - 1)))**(innovationSimple)*
             ((ifelse(hasattemptBothAlready == 1, SS, IS)**(1))*((1 - ifelse(hasattemptBothAlready == 1, SS, IS))**(1 - 1)))**(1 - innovationSimple)*
             (((1 - IS)**(attempt))*(IS**(1 - attempt)))**((1-previousPull)*(1 - innovationSimple)*1)*
             (
               ((E**(1 - attempt))*((1 - E)**(attempt)))**(1 - (1 - innovationL)*max(previousPull*1, (1 - previousPull)*(1 - 1)))*
                 ((IC**(1 - attempt))*((1 - IC)**(attempt)))**((1 - innovationL)*max(previousPull*1, (1 - previousPull)*(1 - 1)))
             )**(1 - ((1-previousPull)*(1 - innovationSimple)*1)))**max(previousPull*attemptL, previousLift*attemptP),
        probability
      ),
      #Case when attempted, knowing it attempted on the two sides
      probability = ifelse(
        opening == "a" & phase == 2 & attemptBoth == 1,
        
        (((SC**(0))*((1 - SC)**(1 - 0)))**(innovationSimple)*
           ((ifelse(hasattemptBothAlready == 1, SS, IS)**(0))*((1 - ifelse(hasattemptBothAlready == 1, SS, IS))**(1 - 0)))**(1 - innovationSimple)*
           (((1 - IS)**(attempt))*(IS**(1 - attempt)))**((1-previousPull)*(1 - innovationSimple)*0)*
           (
             ((E**(1 - attempt))*((1 - E)**(attempt)))**(1 - (1 - innovationL)*max(previousPull*0, (1 - previousPull)*(1 - 0)))*
               ((IC**(1 - attempt))*((1 - IC)**(attempt)))**((1 - innovationL)*max(previousPull*0, (1 - previousPull)*(1 - 0)))
           )**(1 - ((1-previousPull)*(1 - innovationSimple)*0)))**max(previousLift*attemptL, previousPull*attemptP)*#Here it is now multiplied
          
          (((SC**(1))*((1 - SC)**(1 - 1)))**(innovationSimple)*
             ((ifelse(hasattemptBothAlready == 1, SS, IS)**(1))*((1 - ifelse(hasattemptBothAlready == 1, SS, IS))**(1 - 1)))**(1 - innovationSimple)*
             (((1 - IS)**(attempt))*(IS**(1 - attempt)))**((1-previousPull)*(1 - innovationSimple)*1)*
             (
               ((E**(1 - attempt))*((1 - E)**(attempt)))**(1 - (1 - innovationL)*max(previousPull*1, (1 - previousPull)*(1 - 1)))*
                 ((IC**(1 - attempt))*((1 - IC)**(attempt)))**((1 - innovationL)*max(previousPull*1, (1 - previousPull)*(1 - 1)))
             )**(1 - ((1-previousPull)*(1 - innovationSimple)*1)))**max(previousPull*attemptL, previousLift*attemptP),
        probability
      ),
      
      #Phase 3
      probability = ifelse(
        opening != "a" & phase == 3,
        ((SC**(switching))*((1 - SC)**(1 - switching)))**(innovationSimple)*
          ((ifelse(hasattemptBothAlready == 1, SS, IS)**(switching))*((1 - ifelse(hasattemptBothAlready == 1, SS, IS))**(1 - switching)))**(1 - innovationSimple)*
          ((E**(1 - attempt))*((1 - E)**(attempt)))**(max(c(innovationL*switching*(1 - previousLift),innovationL*(1 - switching)*previousLift,
                                                            innovationP * switching * previousLift, innovationP * (1 - switching) * previousPull)))*
          ((IC**(1 - attempt))*((1 - IC)**(attempt)))**(1 - max(c(innovationL*switching*(1 - previousLift),innovationL*(1 - switching)*previousLift,
                                                                  innovationP * switching * previousLift, innovationP * (1 - switching) * previousPull))),
        probability
      ),
      #Case when attempted, knowing it attempted only on one side (or no knowledge of which side)
      probability = ifelse(
        opening == "a" & phase == 3 & attemptBoth != 1,
        
        (((SC**(0))*((1 - SC)**(1 - 0)))**(innovationSimple)*
           ((ifelse(hasattemptBothAlready == 1, SS, IS)**(0))*((1 - ifelse(hasattemptBothAlready == 1, SS, IS))**(1 - 0)))**(1 - innovationSimple)*
           ((E**(1 - attempt))*((1 - E)**(attempt)))**(max(c(innovationL*0*(1 - previousLift),innovationL*(1 - 0)*previousLift,
                                                             innovationP * 0 * previousLift, innovationP * (1 - 0) * previousPull)))*
           ((IC**(1 - attempt))*((1 - IC)**(attempt)))**(1 - max(c(innovationL*0*(1 - previousLift),innovationL*(1 - 0)*previousLift,
                                                                   innovationP * 0 * previousLift, innovationP * (1 - 0) * previousPull))))**max(previousLift*attemptL, previousPull*attemptP) +
          
          (((SC**(1))*((1 - SC)**(1 - 1)))**(innovationSimple)*
             ((ifelse(hasattemptBothAlready == 1, SS, IS)**(1))*((1 - ifelse(hasattemptBothAlready == 1, SS, IS))**(1 - 1)))**(1 - innovationSimple)*
             ((E**(1 - attempt))*((1 - E)**(attempt)))**(max(c(innovationL*1*(1 - previousLift),innovationL*(1 - 1)*previousLift,
                                                               innovationP * 1 * previousLift, innovationP * (1 - 1) * previousPull)))*
             ((IC**(1 - attempt))*((1 - IC)**(attempt)))**(1 - max(c(innovationL*1*(1 - previousLift),innovationL*(1 - 1)*previousLift,
                                                                     innovationP * 1 * previousLift, innovationP * (1 - 1) * previousPull))))**max(previousPull*attemptL, previousLift*attemptP),
        probability
      ),
      #Case when attempted, knowing it attempted on the two sides
      probability = ifelse(
        opening == "a" & phase == 3 & attemptBoth == 1,
        
        (((SC**(0))*((1 - SC)**(1 - 0)))**(innovationSimple)*
           ((ifelse(hasattemptBothAlready == 1, SS, IS)**(0))*((1 - ifelse(hasattemptBothAlready == 1, SS, IS))**(1 - 0)))**(1 - innovationSimple)*
           ((E**(1 - attempt))*((1 - E)**(attempt)))**(max(c(innovationL*0*(1 - previousLift),innovationL*(1 - 0)*previousLift,
                                                             innovationP * 0 * previousLift, innovationP * (1 - 0) * previousPull)))*
           ((IC**(1 - attempt))*((1 - IC)**(attempt)))**(1 - max(c(innovationL*0*(1 - previousLift),innovationL*(1 - 0)*previousLift,
                                                                   innovationP * 0 * previousLift, innovationP * (1 - 0) * previousPull))))**max(previousLift*attemptL, previousPull*attemptP) *#Here it is now multiplied
          
          (((SC**(1))*((1 - SC)**(1 - 1)))**(innovationSimple)*
             ((ifelse(hasattemptBothAlready == 1, SS, IS)**(1))*((1 - ifelse(hasattemptBothAlready == 1, SS, IS))**(1 - 1)))**(1 - innovationSimple)*
             ((E**(1 - attempt))*((1 - E)**(attempt)))**(max(c(innovationL*1*(1 - previousLift),innovationL*(1 - 1)*previousLift,
                                                               innovationP * 1 * previousLift, innovationP * (1 - 1) * previousPull)))*
             ((IC**(1 - attempt))*((1 - IC)**(attempt)))**(1 - max(c(innovationL*1*(1 - previousLift),innovationL*(1 - 1)*previousLift,
                                                                     innovationP * 1 * previousLift, innovationP * (1 - 1) * previousPull))))**max(previousPull*attemptL, previousLift*attemptP),
        probability
      ),
    )
  
  return(data$probability)
}


#Similar to v5, however, correcting for the switch simple when both attempts was done in the previous trial.
#Now, consider no switching for the trial following a both attempt trial
#Remove IS by making innovation simple always equal to 1.
#Remove any observations with unknown attempts
#Fixed the computation that in one side attempt, it does not add 1
estimateProbability6 <- function(
    data,
    SS,
    SC,
    IS,
    IC,
    E,
    countFirstSwitchAsInnovation = FALSE
){
  if(countFirstSwitchAsInnovation){
    data <- data %>%
      mutate(
        innovationSimple = innovationSimple2
      )
  }
  
  #Help for writing formula  
  # SS
  # IS
  # IC
  # SC
  # E
  # attempt
  # switching
  # innovationSimple
  # innovationComplex
  # previousPull
  # previousLift
  # 
  # ((**())*((1 - )**()))**()
  
  data <- data %>% 
    rowwise() %>% 
    mutate(
      probability = NA,
      probability = ifelse(
        opening != "a" & phase == 1 & !is.na(switching),
        (((1 - ifelse(hasattemptBothAlready == 1, SS, IS))**(1 - switching))*(ifelse(hasattemptBothAlready == 1, SS, IS)**(switching)))**(1 - innovationSimple)*
          (((1 - SS)**(1 - switching))*((SS)**(switching)))**(innovationSimple)*
          (((1 - E)**(attempt))*((E)**(1 - attempt)))**(1 - (1 - innovationSimple)*switching)*
          (((1 - IS)**(attempt))*((IS)**(1 - attempt)))**((1 - innovationSimple)*switching),
        probability
      ),
      #Case when attempted, knowing it attempted only on one side (or no knowledge of which side)
      probability = ifelse(
        opening == "a" & phase == 1 & attemptBoth != 1 & !is.na(switching),
        ((((1 - ifelse(hasattemptBothAlready == 1, SS, IS))**(1 - 1))*(ifelse(hasattemptBothAlready == 1, SS, IS)**(1)))**(1 - innovationSimple)*
           (((1 - SS)**(1 - 1))*((SS)**(1)))**(innovationSimple)*
           (((1 - E)**(attempt))*((E)**(1 - attempt)))**(1 - (1 - innovationSimple)*1)*
           (((1 - IS)**(attempt))*((IS)**(1 - attempt)))**((1 - innovationSimple)*1))*max(previousPull*attemptL, previousLift*attemptP) +
          
          ((((1 - ifelse(hasattemptBothAlready == 1, SS, IS))**(1 - 0))*(ifelse(hasattemptBothAlready == 1, SS, IS)**(0)))**(1 - innovationSimple)*
             (((1 - SS)**(1 - 0))*((SS)**(0)))**(innovationSimple)*
             (((1 - E)**(attempt))*((E)**(1 - attempt)))**(1 - (1 - innovationSimple)*0)*
             (((1 - IS)**(attempt))*((IS)**(1 - attempt)))**((1 - innovationSimple)*0))*max(previousLift*attemptL, previousPull*attemptP),
        probability
      ),
      #Case when attempted, knowing it attempted on the two sides
      probability = ifelse(
        opening == "a" & phase == 1 & attemptBoth == 1 & !is.na(switching), 
        
        ((((1 - ifelse(hasattemptBothAlready == 1, SS, IS))**(1 - 1))*(ifelse(hasattemptBothAlready == 1, SS, IS)**(1)))**(1 - innovationSimple)*
           (((1 - SS)**(1 - 1))*((SS)**(1)))**(innovationSimple)*
           (((1 - E)**(attempt))*((E)**(1 - attempt)))**(1 - (1 - innovationSimple)*1)*
           (((1 - IS)**(attempt))*((IS)**(1 - attempt)))**((1 - innovationSimple)*1))**max(previousPull*attemptL, previousLift*attemptP) *#Here it is now multiplied
          
          ((((1 - ifelse(hasattemptBothAlready == 1, SS, IS))**(1 - 0))*(ifelse(hasattemptBothAlready == 1, SS, IS)**(0)))**(1 - innovationSimple)*
             (((1 - SS)**(1 - 0))*((SS)**(0)))**(innovationSimple)*
             (((1 - E)**(attempt))*((E)**(1 - attempt)))**(1 - (1 - innovationSimple)*0)*
             (((1 - IS)**(attempt))*((IS)**(1 - attempt)))**((1 - innovationSimple)*0))**max(previousLift*attemptL, previousPull*attemptP),
        probability
      ),
      #Case when attempted, previous attempt is both sides
      probability = ifelse(
        phase == 1 & is.na(switching), 
        (((1 - E)**(attempt))*((E)**(1 - attempt)))**(2**attemptBoth),
        probability
      ),
      
      #Phase 2
      probability = ifelse(
        opening != "a" & phase == 2 & !is.na(switching),
        ((SC**(switching))*((1 - SC)**(1 - switching)))**(innovationSimple)*
          ((ifelse(hasattemptBothAlready == 1, SS, IS)**(switching))*((1 - ifelse(hasattemptBothAlready == 1, SS, IS))**(1 - switching)))**(1 - innovationSimple)*
          (((1 - IS)**(attempt))*(IS**(1 - attempt)))**((1-previousPull)*(1 - innovationSimple)*switching)*
          (
            ((E**(1 - attempt))*((1 - E)**(attempt)))**(1 - (1 - innovationL)*max(previousPull*switching, (1 - previousPull)*(1 - switching)))*
              ((IC**(1 - attempt))*((1 - IC)**(attempt)))**((1 - innovationL)*max(previousPull*switching, (1 - previousPull)*(1 - switching)))
          )**(1 - ((1-previousPull)*(1 - innovationSimple)*switching)),
        probability
      ),
      #Case when attempted, knowing it attempted only on one side (or no knowledge of which side)
      probability = ifelse(
        opening == "a" & phase == 2 & attemptBoth != 1 & !is.na(switching),
        
        (((SC**(0))*((1 - SC)**(1 - 0)))**(innovationSimple)*
           ((ifelse(hasattemptBothAlready == 1, SS, IS)**(0))*((1 - ifelse(hasattemptBothAlready == 1, SS, IS))**(1 - 0)))**(1 - innovationSimple)*
           (((1 - IS)**(attempt))*(IS**(1 - attempt)))**((1-previousPull)*(1 - innovationSimple)*0)*
           (
             ((E**(1 - attempt))*((1 - E)**(attempt)))**(1 - (1 - innovationL)*max(previousPull*0, (1 - previousPull)*(1 - 0)))*
               ((IC**(1 - attempt))*((1 - IC)**(attempt)))**((1 - innovationL)*max(previousPull*0, (1 - previousPull)*(1 - 0)))
           )**(1 - ((1-previousPull)*(1 - innovationSimple)*0)))*max(previousLift*attemptL, previousPull*attemptP) +
          
          (((SC**(1))*((1 - SC)**(1 - 1)))**(innovationSimple)*
             ((ifelse(hasattemptBothAlready == 1, SS, IS)**(1))*((1 - ifelse(hasattemptBothAlready == 1, SS, IS))**(1 - 1)))**(1 - innovationSimple)*
             (((1 - IS)**(attempt))*(IS**(1 - attempt)))**((1-previousPull)*(1 - innovationSimple)*1)*
             (
               ((E**(1 - attempt))*((1 - E)**(attempt)))**(1 - (1 - innovationL)*max(previousPull*1, (1 - previousPull)*(1 - 1)))*
                 ((IC**(1 - attempt))*((1 - IC)**(attempt)))**((1 - innovationL)*max(previousPull*1, (1 - previousPull)*(1 - 1)))
             )**(1 - ((1-previousPull)*(1 - innovationSimple)*1)))*max(previousPull*attemptL, previousLift*attemptP),
        probability
      ),
      #Case when attempted, knowing it attempted on the two sides
      probability = ifelse(
        opening == "a" & phase == 2 & attemptBoth == 1 & !is.na(switching),
        
        (((SC**(0))*((1 - SC)**(1 - 0)))**(innovationSimple)*
           ((ifelse(hasattemptBothAlready == 1, SS, IS)**(0))*((1 - ifelse(hasattemptBothAlready == 1, SS, IS))**(1 - 0)))**(1 - innovationSimple)*
           (((1 - IS)**(attempt))*(IS**(1 - attempt)))**((1-previousPull)*(1 - innovationSimple)*0)*
           (
             ((E**(1 - attempt))*((1 - E)**(attempt)))**(1 - (1 - innovationL)*max(previousPull*0, (1 - previousPull)*(1 - 0)))*
               ((IC**(1 - attempt))*((1 - IC)**(attempt)))**((1 - innovationL)*max(previousPull*0, (1 - previousPull)*(1 - 0)))
           )**(1 - ((1-previousPull)*(1 - innovationSimple)*0)))**max(previousLift*attemptL, previousPull*attemptP)*#Here it is now multiplied
          
          (((SC**(1))*((1 - SC)**(1 - 1)))**(innovationSimple)*
             ((ifelse(hasattemptBothAlready == 1, SS, IS)**(1))*((1 - ifelse(hasattemptBothAlready == 1, SS, IS))**(1 - 1)))**(1 - innovationSimple)*
             (((1 - IS)**(attempt))*(IS**(1 - attempt)))**((1-previousPull)*(1 - innovationSimple)*1)*
             (
               ((E**(1 - attempt))*((1 - E)**(attempt)))**(1 - (1 - innovationL)*max(previousPull*1, (1 - previousPull)*(1 - 1)))*
                 ((IC**(1 - attempt))*((1 - IC)**(attempt)))**((1 - innovationL)*max(previousPull*1, (1 - previousPull)*(1 - 1)))
             )**(1 - ((1-previousPull)*(1 - innovationSimple)*1)))**max(previousPull*attemptL, previousLift*attemptP),
        probability
      ),
      #Case when  previous trial is attempt on both sides
      probability = ifelse(
        phase == 2 & is.na(switching),
        (
          ((E**(1 - attempt))*((1 - E)**(attempt)))**(sum(c(attemptP, ifelse(opening == "p", 1, 0),
                                                            innovationL*ifelse(opening == "l", 1, 0), innovationL*attemptL), na.rm = TRUE))*
            ((IC**(1 - attempt))*((1 - IC)**(attempt))) ** ((1 - innovationL) * max(attemptL, ifelse(opening == "l", 1, 0), na.rm = TRUE))
        ),
        probability
      ),
      
      #Phase 3
      probability = ifelse(
        opening != "a" & phase == 3 & !is.na(switching),
        ((SC**(switching))*((1 - SC)**(1 - switching)))**(innovationSimple)*
          ((ifelse(hasattemptBothAlready == 1, SS, IS)**(switching))*((1 - ifelse(hasattemptBothAlready == 1, SS, IS))**(1 - switching)))**(1 - innovationSimple)*
          ((E**(1 - attempt))*((1 - E)**(attempt)))**(max(c(innovationL*switching*(1 - previousLift),innovationL*(1 - switching)*previousLift,
                                                            innovationP * switching * previousLift, innovationP * (1 - switching) * previousPull)))*
          ((IC**(1 - attempt))*((1 - IC)**(attempt)))**(1 - max(c(innovationL*switching*(1 - previousLift),innovationL*(1 - switching)*previousLift,
                                                                  innovationP * switching * previousLift, innovationP * (1 - switching) * previousPull))),
        probability
      ),
      #Case when attempted, knowing it attempted only on one side (or no knowledge of which side)
      probability = ifelse(
        opening == "a" & phase == 3 & attemptBoth != 1 & !is.na(switching),
        
        (((SC**(0))*((1 - SC)**(1 - 0)))**(innovationSimple)*
           ((ifelse(hasattemptBothAlready == 1, SS, IS)**(0))*((1 - ifelse(hasattemptBothAlready == 1, SS, IS))**(1 - 0)))**(1 - innovationSimple)*
           ((E**(1 - attempt))*((1 - E)**(attempt)))**(max(c(innovationL*0*(1 - previousLift),innovationL*(1 - 0)*previousLift,
                                                             innovationP * 0 * previousLift, innovationP * (1 - 0) * previousPull)))*
           ((IC**(1 - attempt))*((1 - IC)**(attempt)))**(1 - max(c(innovationL*0*(1 - previousLift),innovationL*(1 - 0)*previousLift,
                                                                   innovationP * 0 * previousLift, innovationP * (1 - 0) * previousPull))))*max(previousLift*attemptL, previousPull*attemptP) +
          
          (((SC**(1))*((1 - SC)**(1 - 1)))**(innovationSimple)*
             ((ifelse(hasattemptBothAlready == 1, SS, IS)**(1))*((1 - ifelse(hasattemptBothAlready == 1, SS, IS))**(1 - 1)))**(1 - innovationSimple)*
             ((E**(1 - attempt))*((1 - E)**(attempt)))**(max(c(innovationL*1*(1 - previousLift),innovationL*(1 - 1)*previousLift,
                                                               innovationP * 1 * previousLift, innovationP * (1 - 1) * previousPull)))*
             ((IC**(1 - attempt))*((1 - IC)**(attempt)))**(1 - max(c(innovationL*1*(1 - previousLift),innovationL*(1 - 1)*previousLift,
                                                                     innovationP * 1 * previousLift, innovationP * (1 - 1) * previousPull))))*max(previousPull*attemptL, previousLift*attemptP),
        probability
      ),
      #Case when attempted, knowing it attempted on the two sides
      probability = ifelse(
        opening == "a" & phase == 3 & attemptBoth == 1 & !is.na(switching),
        
        (((SC**(0))*((1 - SC)**(1 - 0)))**(innovationSimple)*
           ((ifelse(hasattemptBothAlready == 1, SS, IS)**(0))*((1 - ifelse(hasattemptBothAlready == 1, SS, IS))**(1 - 0)))**(1 - innovationSimple)*
           ((E**(1 - attempt))*((1 - E)**(attempt)))**(max(c(innovationL*0*(1 - previousLift),innovationL*(1 - 0)*previousLift,
                                                             innovationP * 0 * previousLift, innovationP * (1 - 0) * previousPull)))*
           ((IC**(1 - attempt))*((1 - IC)**(attempt)))**(1 - max(c(innovationL*0*(1 - previousLift),innovationL*(1 - 0)*previousLift,
                                                                   innovationP * 0 * previousLift, innovationP * (1 - 0) * previousPull))))**max(previousLift*attemptL, previousPull*attemptP) *#Here it is now multiplied
          
          (((SC**(1))*((1 - SC)**(1 - 1)))**(innovationSimple)*
             ((ifelse(hasattemptBothAlready == 1, SS, IS)**(1))*((1 - ifelse(hasattemptBothAlready == 1, SS, IS))**(1 - 1)))**(1 - innovationSimple)*
             ((E**(1 - attempt))*((1 - E)**(attempt)))**(max(c(innovationL*1*(1 - previousLift),innovationL*(1 - 1)*previousLift,
                                                               innovationP * 1 * previousLift, innovationP * (1 - 1) * previousPull)))*
             ((IC**(1 - attempt))*((1 - IC)**(attempt)))**(1 - max(c(innovationL*1*(1 - previousLift),innovationL*(1 - 1)*previousLift,
                                                                     innovationP * 1 * previousLift, innovationP * (1 - 1) * previousPull))))**max(previousPull*attemptL, previousLift*attemptP),
        probability
      ),
      #Case when phase 3, previous trial had both attempts
      probability = ifelse(
        phase == 3 & is.na(switching),
        ((E**(1 - attempt))*((1 - E)**(attempt)))**(sum(c(innovationL * ifelse(opening == "l", 1, 0), innovationL * ifelse(is.na(attemptL), 0, attemptL),
                                                          innovationP * ifelse(opening == "p", 1, 0), innovationP * ifelse(is.na(attemptP), 0, attemptP)), na.rm = TRUE))*
          ((IC**(1 - attempt))*((1 - IC)**(attempt)))**(sum(c((1 - innovationL) * ifelse(opening == "l", 1, 0), (1 - innovationL) * ifelse(is.na(attemptL), 0, attemptL),
                                                              (1 - innovationP) * ifelse(opening == "p", 1, 0), (1 - innovationP) * ifelse(is.na(attemptP), 0, attemptP)), na.rm = TRUE)),
        probability
      )
    )
  
  return(data$probability)
}


#Similar to v6, however, adds IS back in
#Allow first trial back in

estimateProbability7 <- function(
    data,
    SS,
    SC,
    IS,
    IC,
    E,
    countFirstSwitchAsInnovation = FALSE
){
  if(countFirstSwitchAsInnovation){
    data <- data %>%
      mutate(
        innovationSimple = innovationSimple2
      )
  }
  
  #Help for writing formula  
  # SS
  # IS
  # IC
  # SC
  # E
  # attempt
  # switching
  # innovationSimple
  # innovationComplex
  # previousPull
  # previousLift
  # 
  # ((**())*((1 - )**()))**()
  
  data <- data %>% 
    rowwise() %>% 
    mutate(
      probability = NA,
      probability = ifelse(
        opening != "a" & phase == 1 & !is.na(switching),
        (((1 - SS)**(1 - switching))*((SS)**(switching)))*
          (((1 - E)**(attempt))*((E)**(1 - attempt)))**(max(innovationLS*ifelse(opening == "l", 1, 0), innovationPS*ifelse(opening == "p", 1, 0), na.rm = TRUE))*
          (((1 - IS)**(attempt))*((IS)**(1 - attempt)))**(1 - (max(innovationLS*ifelse(opening == "l", 1, 0), innovationPS*ifelse(opening == "p", 1, 0), na.rm = TRUE))),
        probability
      ),
      #Case when attempted, knowing it attempted only on one side
      probability = ifelse(
        opening == "a" & phase == 1 & attemptBoth != 1 & !is.na(switching),
        ((((1 - SS)**(max(previousLift*attemptL, previousPull*attemptP)))*((SS)**(max(previousPull*attemptL, previousLift*attemptP))))*
           (((1 - E)**(attempt))*((E)**(1 - attempt)))**(max(innovationLS*attemptL, innovationPS*attemptP, na.rm = TRUE))*
           (((1 - IS)**(attempt))*((IS)**(1 - attempt)))**(max((1 - innovationLS)*attemptL, (1 - innovationPS)*attemptP, na.rm = TRUE))),
        probability
      ),
      #Case when attempted, knowing it attempted on the two sides
      probability = ifelse(
        opening == "a" & phase == 1 & attemptBoth == 1 & !is.na(switching), 
        (((1 - SS)*(SS))*
           (((1 - E)**(attempt))*((E)**(1 - attempt)))**(sum(innovationLS*attemptL, innovationPS*attemptP, na.rm = TRUE))*
           (((1 - IS)**(attempt))*((IS)**(1 - attempt)))**(sum((1 - innovationLS)*attemptL, (1 - innovationPS)*attemptP, na.rm = TRUE))),
        probability
      ),
      #Case when attempted, previous attempt is both sides
      probability = ifelse(
        phase == 1 & is.na(switching), 
        (((1 - E)**(attempt))*((E)**(1 - attempt)))**(sum(innovationLS*attemptL, innovationLS*ifelse(opening == "l", 1, 0), 
                                                          innovationPS*attemptP, innovationPS*ifelse(opening == "p", 1, 0), na.rm = TRUE))*
          (((1 - IS)**(attempt))*((IS)**(1 - attempt)))**(sum((1 - innovationLS)*attemptL, (1 - innovationLS)*ifelse(opening == "l", 1, 0), 
                                                              (1 - innovationPS)*attemptP, (1 - innovationPS)*ifelse(opening == "p", 1, 0), na.rm = TRUE)),
        probability
      ),
      
      #Phase 2
      probability = ifelse(
        opening != "a" & phase == 2 & !is.na(switching),
        ((SC**(switching))*((1 - SC)**(1 - switching)))*
          (((1 - IS)**(attempt))*(IS**(1 - attempt)))**(max((1 - innovationLS)*ifelse(opening == "l", 1, 0), (1 - innovationPS)*ifelse(opening == "p", 1, 0), na.rm = TRUE))*
          ((E**(1 - attempt))*((1 - E)**(attempt)))**(max(innovationLS*ifelse(opening == "l", 1, 0), innovationPS*ifelse(opening == "p", 1, 0), na.rm = TRUE))*
          ((IC**(1 - attempt))*((1 - IC)**(attempt)))**((1-innovationL)*ifelse(opening == "l", 1, 0)),
        probability
      ),
      #Case when attempted, knowing it attempted only on one side
      probability = ifelse(
        opening == "a" & phase == 2 & attemptBoth != 1 & !is.na(switching),
        ((SC**(switching))*((1 - SC)**(1 - switching)))*
          (((1 - IS)**(attempt))*(IS**(1 - attempt)))**(max((1 - innovationLS)*attemptL, (1 - innovationPS)*attemptP, na.rm = TRUE))*
          ((E**(1 - attempt))*((1 - E)**(attempt)))**(max(innovationLS*attemptL, innovationPS*attemptP, na.rm = TRUE))*
          ((IC**(1 - attempt))*((1 - IC)**(attempt)))**((1-innovationL)*attemptL),
        probability
      ),
      #Case when attempted, knowing it attempted on the two sides
      probability = ifelse(
        opening == "a" & phase == 2 & attemptBoth == 1 & !is.na(switching),
        ((SC**(switching))*((1 - SC)**(1 - switching)))*
          (((1 - IS)**(attempt))*(IS**(1 - attempt)))**(sum((1 - innovationLS)*attemptL, (1 - innovationPS)*attemptP, na.rm = TRUE))*
          ((E**(1 - attempt))*((1 - E)**(attempt)))**(sum(innovationLS*attemptL, innovationPS*attemptP, na.rm = TRUE))*
          ((IC**(1 - attempt))*((1 - IC)**(attempt)))**((1-innovationL)*attemptL),
        probability
      ),
      #Case when  previous trial is attempt on both sides
      probability = ifelse(
        phase == 2 & is.na(switching),
        (
          ((E**(1 - attempt))*((1 - E)**(attempt)))**(sum(c(innovationPS*attemptP, innovationPS*ifelse(opening == "p", 1, 0),
                                                            innovationLS*innovationL*ifelse(opening == "l", 1, 0), innovationLS*innovationL*attemptL), na.rm = TRUE))*
            ((IS**(1 - attempt))*((1 - IS)**(attempt)))**(sum(c((1 - innovationPS)*attemptP, 
                                                                (1 - innovationPS)*ifelse(opening == "p", 1, 0),
                                                                (1 - innovationLS)*(1 - innovationL)*ifelse(opening == "l", 1, 0), 
                                                                (1 - innovationLS)*(1 - innovationL)*attemptL), na.rm = TRUE))*
            ((IC**(1 - attempt))*((1 - IC)**(attempt))) ** ((1 - innovationL) * max(attemptL, ifelse(opening == "l", 1, 0), na.rm = TRUE))
        ),
        probability
      ),
      
      #Phase 3
      probability = ifelse(
        opening != "a" & phase == 3 & !is.na(switching),
        ((SC**(switching))*((1 - SC)**(1 - switching)))*
          (((1 - IS)**(attempt))*(IS**(1 - attempt)))**(max((1 - innovationLS)*ifelse(opening == "l", 1, 0), (1 - innovationPS)*ifelse(opening == "p", 1, 0), na.rm = TRUE))*
          ((E**(1 - attempt))*((1 - E)**(attempt)))**(max(innovationL*innovationLS*ifelse(opening == "l", 1, 0), innovationP*innovationPS*ifelse(opening == "p", 1, 0), na.rm = TRUE))*
          ((IC**(1 - attempt))*((1 - IC)**(attempt)))**(max((1-innovationL)*ifelse(opening == "l", 1, 0), (1-innovationP)*ifelse(opening == "p", 1, 0))),
        probability
      ),
      #Case when attempted, knowing it attempted only on one side (or no knowledge of which side)
      probability = ifelse(
        opening == "a" & phase == 3 & attemptBoth != 1 & !is.na(switching),
        ((SC**(switching))*((1 - SC)**(1 - switching)))*
          (((1 - IS)**(attempt))*(IS**(1 - attempt)))**(max((1 - innovationLS)*attemptL, (1 - innovationPS)*attemptP, na.rm = TRUE))*
          ((E**(1 - attempt))*((1 - E)**(attempt)))**(max(innovationL*innovationLS*attemptL, innovationP*innovationPS*attemptP, na.rm = TRUE))*
          ((IC**(1 - attempt))*((1 - IC)**(attempt)))**(max((1-innovationL)*attemptL, (1-innovationP)*attemptP)),
        probability
      ),
      #Case when attempted, knowing it attempted on the two sides
      probability = ifelse(
        opening == "a" & phase == 3 & attemptBoth == 1 & !is.na(switching),
        ((SC**(switching))*((1 - SC)**(1 - switching)))*
          (((1 - IS)**(attempt))*(IS**(1 - attempt)))**(sum((1 - innovationLS)*attemptL, (1 - innovationPS)*attemptP, na.rm = TRUE))*
          ((E**(1 - attempt))*((1 - E)**(attempt)))**(sum(innovationL*innovationLS*attemptL, innovationP*innovationPS*attemptP, na.rm = TRUE))*
          ((IC**(1 - attempt))*((1 - IC)**(attempt)))**(sum((1-innovationL)*attemptL, (1-innovationP)*attemptP)),
        probability
      ),
      #Case when phase 3, previous trial had both attempts
      probability = ifelse(
        phase == 3 & is.na(switching),
        ((E**(1 - attempt))*((1 - E)**(attempt)))**(sum(c(innovationL * ifelse(opening == "l", 1, 0), innovationL * ifelse(is.na(attemptL), 0, attemptL),
                                                          innovationP * ifelse(opening == "p", 1, 0), innovationP * ifelse(is.na(attemptP), 0, attemptP)), na.rm = TRUE))*
          ((IS**(1 - attempt))*((1 - IS)**(attempt)))**(sum(c((1 - innovationPS)*attemptP, 
                                                              (1 - innovationPS)*ifelse(opening == "p", 1, 0),
                                                              (1 - innovationLS)*(1 - innovationL)*ifelse(opening == "l", 1, 0), 
                                                              (1 - innovationLS)*(1 - innovationL)*attemptL), na.rm = TRUE))*
          ((IC**(1 - attempt))*((1 - IC)**(attempt)))**(sum(c((1 - innovationL) * ifelse(opening == "l", 1, 0), (1 - innovationL) * ifelse(is.na(attemptL), 0, attemptL),
                                                              (1 - innovationP) * ifelse(opening == "p", 1, 0), (1 - innovationP) * ifelse(is.na(attemptP), 0, attemptP)), na.rm = TRUE)),
        probability
      )
    )
  
  return(data$probability)
}

estimateLL <- function(parameters_v,
                       data){
  # print(parameters_v)
  
  parameters_v = ifelse(parameters_v < 0.00001, 0.00001, parameters_v)
  parameters_v = ifelse(parameters_v > 0.99999, 0.99999, parameters_v)
  
  probability_v <- estimateProbability(
    data = data,
    SS = parameters_v[1],
    SC = parameters_v[2],
    IS = parameters_v[3],
    IC = parameters_v[4],
    E = parameters_v[5])
  # print(probability_v[8:9])
  likelihood <- sum(log(probability_v))
  # print(likelihood)
  return(-likelihood)
}

estimateLL2 <- function(parameters_v,
                        data){
  # print("test")
  # print(parameters_v)
  
  parameters_v = ifelse(parameters_v < 0.00001, 0.00001, parameters_v)
  parameters_v = ifelse(parameters_v > 0.99999, 0.99999, parameters_v)
  
  probability_v <- estimateProbability2(
    data = data,
    SS = parameters_v[1],
    SC = parameters_v[2],
    IS = parameters_v[3],
    IC = parameters_v[4],
    ES = parameters_v[5],
    EC = parameters_v[6])
  # print(probability_v[8:9])
  # print(probability_v)
  likelihood <- sum(log(probability_v))
  # print(likelihood)
  return(-likelihood)
}

estimateLL3 <- function(parameters_v,
                        data){
  # print(parameters_v)
  
  parameters_v = ifelse(parameters_v < 0.00001, 0.00001, parameters_v)
  parameters_v = ifelse(parameters_v > 0.99999, 0.99999, parameters_v)
  
  probability_v <- estimateProbability3(
    data = data,
    SS = parameters_v[1],
    SC = parameters_v[2],
    IS = parameters_v[3],
    IC = parameters_v[4],
    ES = parameters_v[5],
    EC = parameters_v[6])
  # print(probability_v[8:9])
  likelihood <- sum(log(probability_v))
  # print(likelihood)
  return(-likelihood)
}

estimateLL4 <- function(parameters_v,
                        data){
  # print(parameters_v)
  
  parameters_v = ifelse(parameters_v < 0.00001, 0.00001, parameters_v)
  parameters_v = ifelse(parameters_v > 0.99999, 0.99999, parameters_v)
  
  probability_v <- estimateProbability4(
    data = data,
    SS = parameters_v[1],
    SC = parameters_v[2],
    IS = parameters_v[3],
    IC = parameters_v[4],
    ES = parameters_v[5],
    EC = parameters_v[6])
  # print(probability_v[8:9])
  likelihood <- sum(log(probability_v))
  # print(likelihood)
  return(-likelihood)
}

estimateLL5 <- function(parameters_v,
                        data){
  # print(parameters_v)
  
  parameters_v = ifelse(parameters_v < 0.00001, 0.00001, parameters_v)
  parameters_v = ifelse(parameters_v > 0.99999, 0.99999, parameters_v)
  
  probability_v <- estimateProbability5(
    data = data,
    SS = parameters_v[1],
    SC = parameters_v[2],
    IS = parameters_v[3],
    IC = parameters_v[4],
    E = parameters_v[5])
  # print(probability_v[8:9])
  likelihood <- sum(log(probability_v))
  # print(likelihood)
  return(-likelihood)
}

estimateLL6 <- function(parameters_v,
                        data){
  # print(parameters_v)
  
  parameters_v = ifelse(parameters_v < 0.00001, 0.00001, parameters_v)
  parameters_v = ifelse(parameters_v > 0.99999, 0.99999, parameters_v)
  
  probability_v <- estimateProbability6(
    data = data,
    SS = parameters_v[1],
    SC = parameters_v[2],
    IS = parameters_v[3],
    IC = parameters_v[4],
    E = parameters_v[5])
  # print(probability_v[8:9])
  likelihood <- sum(log(probability_v))
  # print(likelihood)
  return(-likelihood)
}

estimateLL7 <- function(parameters_v,
                        data){
  # print(parameters_v)
  
  parameters_v = ifelse(parameters_v < 0.00001, 0.00001, parameters_v)
  parameters_v = ifelse(parameters_v > 0.99999, 0.99999, parameters_v)
  
  probability_v <- estimateProbability7(
    data = data,
    SS = parameters_v[1],
    SC = parameters_v[2],
    IS = parameters_v[3],
    IC = parameters_v[4],
    E = parameters_v[5])
  # print(probability_v[8:9])
  likelihood <- sum(log(probability_v))
  # print(likelihood)
  return(-likelihood)
}

#Estimate 95%CI of parameters
est95FromHessian <- function(hess, par, whichNotFitted){
  
  if(length(whichNotFitted) > 0){
    prop_sigma <- sqrt(diag(solve(hess[-whichNotFitted, -whichNotFitted])))#I minimise the log(likelihood, so take the hess directly)
    # 95% CIs are calculated as: parameter estimates +/- 1.96 * prop_sigm
    est95 <- cbind(
      par[-whichNotFitted],
      par[-whichNotFitted] - 1.96*prop_sigma,
      par[-whichNotFitted] + 1.96*prop_sigma
    ) %>% as.data.frame()
  }else{
    prop_sigma <- sqrt(diag(solve(hess)))#I minimise the log(likelihood, so take the hess directly)
    # 95% CIs are calculated as: parameter estimates +/- 1.96 * prop_sigm
    est95 <- cbind(
      par,
      par - 1.96*prop_sigma,
      par + 1.96*prop_sigma
    ) %>% as.data.frame()
  }
  
  colnames(est95) <- c("est", "estLwr", "estUpr")
  return(est95)
}


# Analysis ----------------------------------------------------------------

## Preparing the model -----------------------------------------------------

### Extract complementary descriptive covariates ------------------------------------

complementaryCovariates_df <- dataExperimentBox_df %>% 
  arrange(id, phase, trial) %>% 
  group_by(id) %>% 
  mutate(
    previousOpening = lag(opening),
    firstInnovationSimple = ifelse((opening == "l" & previousOpening == "p") | (opening == "l" & previousOpening == "l"), TRUE, FALSE),
    trialBeforeFirstInnovationSimple = which(firstInnovationSimple)[1],
    firstInnovationL = ifelse(innovationL == 0 & opening == "ul", TRUE, FALSE), 
    trialBeforeFirstInnovationL = which(firstInnovationL & phase > 1)[1],
    firstInnovationP = ifelse(innovationP == 0 & opening == "up", TRUE, FALSE),
    trialBeforeFirstInnovationP = which(firstInnovationP & phase > 2)[1]
  ) %>%  
  group_by(id, phase) %>% 
  summarise(
    totTrials = length(opening),
    successRate = length(opening[opening %in% c("l", "ul", "p", "up")])/totTrials,
    maxFrequency = ifelse(length(opening[opening %in% c("p", "up")]) >  length(opening[opening %in% c("l", "ul")]), 
                          length(opening[opening %in% c("p", "up") ])/(length(opening[opening %in% c("p", "up")]) + length(opening[opening %in% c("l", "ul")])), 
                          length(opening[opening %in% c("l", "ul")])/(length(opening[opening %in% c("p", "up")]) + length(opening[opening %in% c("l", "ul")]))),
    whichMaxFrequency = ifelse(length(opening[opening %in% c("p", "up")]) > length(opening[opening %in% c("l", "ul")]), "p", "l"),
    whichMaxFrequency = ifelse(successRate == 0, NA, whichMaxFrequency),
    
    #Simple innovation from lifting to pulling or pulling to lifting, without the lock
    innovationSimpleOccurred = ifelse(phase == 1 & "p" %in% opening & "l" %in% opening, TRUE, FALSE),
    innovationSimpleOccurred = ifelse(phase == 2 & FALSE %in% innovationSimpleOccurred[phase == 1] &
                                        !("l" %in% opening[phase == 1]) & "p" %in% opening[phase == 2], TRUE, innovationSimpleOccurred),
    
    #Accounting for the lock innovation (lifting lock)
    innovationComplexLOccurred = ifelse(TRUE %in% firstInnovationL, TRUE, FALSE),
    
    #Accounting for the lock innovation (pulling lock)
    innovationComplexPOccurred = ifelse(TRUE %in% firstInnovationP, TRUE, FALSE),
    
    howManyTrialsBeforeInnovationSimple = ifelse(TRUE %in% firstInnovationSimple, trialBeforeFirstInnovationSimple, NA),
    howManyTrialsBeforeInnovationInnovationL = ifelse(TRUE %in% firstInnovationL, trialBeforeFirstInnovationL, NA),
    howManyTrialsBeforeInnovationInnovationP = ifelse(TRUE %in% firstInnovationP, trialBeforeFirstInnovationP, NA)
    
  ) %>% 
  unique()

### Preparing the data to be modelled --------------------------------------
dataExperimentBox_df <- dataExperimentBox_df %>% 
  arrange(id, phase, trial) %>% 
  group_by(id) %>% 
  mutate(
    rowNumber = 1:length(opening),
    
    #Estimate first switch "normal"
    openingBinary = ifelse(opening == "p" | opening == "l", opening, NA), 
    openingBinary = ifelse(openingBinary == ifelse(!is.null(openingBinary[!is.na(openingBinary)]), openingBinary[!is.na(openingBinary)][1], "impossible"), 0, 1),
    #Solve the problem for those only having attempt otherwise does not compute
    openingBinary = ifelse(trial == 1 & is.na(openingBinary), 0, openingBinary),
    openingBinary = na.locf(openingBinary),
    openingBinary = cumsum(openingBinary),
    openingBinary = ifelse(openingBinary == 0, 0, 1),

    innovationLS = ifelse(+row_number() > min(which(opening == "l")) | phase > min(phase[opening == "l"]), 1, 0), 
    innovationPS = ifelse(+row_number() > min(which(opening == "p")) | phase > min(phase[opening == "p"]), 1, 0),

    #Extract innovation simple (first switch)
    opening = ifelse(opening == "ul", "l", opening),
    opening = ifelse(opening == "up", "p", opening),
    previousOpening = lag(opening, default = "l"),
    previousOpening = ifelse(previousOpening != "p" & previousOpening != "l", NA, previousOpening),
    previousOpening = na.locf(previousOpening),
    attempt = ifelse(opening == "a", 1, 0),
    previousLift = ifelse(previousOpening == "l", 1, 0),
    previousPull = ifelse(previousOpening == "p", 1, 0),
    
    #Process the attempts
    attemptBoth = ifelse(attemptL == 1 & attemptP == 1, 1, 0),
    #For those with no knowledge (attemptL and attemptP are 0), replace to "1" as if tried both
    attemptUnknown = ifelse(attemptL == 0 & attemptP == 0, 1, 0),
    attemptL = ifelse(attemptUnknown == 1, 1, attemptL),
    attemptP = ifelse(attemptUnknown == 1, 1, attemptP),
    previousOpening = ifelse(mapply(identical, lag(attemptBoth), 1), NA, 
                             ifelse(mapply(identical, lag(attemptL), 1) & mapply(identical, lag(attemptP), 0), "l",
                                    ifelse( mapply(identical, lag(attemptL), 0) & mapply(identical, lag(attemptP), 1), "p", previousOpening))),
    switching = ifelse(opening == previousOpening, 0, 1),
    switching = ifelse(mapply(identical, lag(attemptBoth), 1), NA, switching),
    
    previousAttemptL = lag(attemptL),
    previousAttemptP = lag(attemptP),

    attemptBoth = ifelse(is.na(attemptBoth), 0, attemptBoth),

  ) %>% 
  # filter(phase != 1 | trial != 1) %>%
  filter(!mapply(identical, attemptUnknown, 1)) %>%
  #only include those with at least 10 phase 1 trials
  filter(!id %in% (dataExperimentBox_df %>%
                     group_by(id, phase) %>%
                     filter(phase == 1) %>%
                     summarise(phasesum = n()) %>%
                     filter(phasesum < 10))$id) %>%
  dplyr::select(-rowNumber, -openingBinary)

#allow first phase to be included without switching
dataExperimentBox_df$switching[dataExperimentBox_df$trial == 1 & dataExperimentBox_df$phase == 1] <- NA


#both sides should be counted in every attempt
dataExperimentBox_df$attemptBoth[dataExperimentBox_df$opening == "a"] <- 1
dataExperimentBox_df$attemptL[dataExperimentBox_df$opening == "a"] <- 1
dataExperimentBox_df$attemptP[dataExperimentBox_df$opening == "a"] <- 1


## Fitting the model -------------------------------------------------------

uniqueId_v <- unique(dataExperimentBox_df$id)

# estimateLL(
#   parameters_v = c(0.25,0.25, 0.5, 0.5, 0.5),
#   data = dataOfIndividualInterest_df
# )
# 
# estimateLL(
#   parameters_v = c(0.5,0.5, 0.5, 0.5, 0.5),
#   data = dataOfIndividualInterest_df
# )

#Do it by individual

library(pbapply)

fittedModels_l <- pblapply(1:length(uniqueId_v), function(whichIndiv){
  #print(whichIndiv)
  dataOfIndividualInterest_df <- dataExperimentBox_df %>% filter(id == uniqueId_v[whichIndiv])

  
  optimisedLL <- optim(
    par =  rep(0.25, times = 5),
    fn = estimateLL7,
    data = dataOfIndividualInterest_df, 
    hessian = TRUE,
    lower = rep(0.00001, times = 5),
    upper = rep(0.99999, times = 5),
    method = "L-BFGS-B"#,
    #control = control
  )
  
  #Estimate hessian matrix and 95% CI
  hessianMat <- hessian(estimateLL7,
                        x = optimisedLL$par,
                        method = "Richardson",
                        data = dataExperimentBox_df %>% filter(id == uniqueId_v[whichIndiv])
  )
  
  #Prepare CI95%
  ci95 <- data.frame(
    est = optimisedLL$par,
    estLwr = NA,
    estUpr = NA
  )
  
  #Check which parameters were not fitted (e.g., the individual never switched simple because never innovated simple early enough)
  whichNotFitted <- which(optimisedLL$par == 0.25)
  
  tryCatch(
    {
      
      if(length(whichNotFitted) > 0){
        ci95[-whichNotFitted,] <- est95FromHessian(hess = hessianMat, par = optimisedLL$par, whichNotFitted = whichNotFitted) %>% 
          mutate(
            estLwr = ifelse(estLwr < 0, 0, estLwr),
            estUpr = ifelse(estUpr > 1, 1, estUpr)
          )
      }else{
        ci95 <- est95FromHessian(hess = hessianMat, par = optimisedLL$par, whichNotFitted = whichNotFitted) %>% 
          mutate(
            estLwr = ifelse(estLwr < 0, 0, estLwr),
            estUpr = ifelse(estUpr > 1, 1, estUpr)
          )
      }
      
      
      ci95[-whichNotFitted,] <- est95FromHessian(hess = hessianMat, par = optimisedLL$par, whichNotFitted = whichNotFitted) %>% 
        mutate(
          estLwr = ifelse(estLwr < 0, 0, estLwr),
          estUpr = ifelse(estUpr > 1, 1, estUpr)
        )
      hasPhase2 <- ifelse(length(unique(dataOfIndividualInterest_df$phase)) >= 2, TRUE, FALSE)
      
      if(!hasPhase2){
        ci95[c(2,4),] <- NA
      }
      
    },
    error = function(e){
      hasPhase2 <- ifelse(length(unique(dataOfIndividualInterest_df$phase)) >= 2, TRUE, FALSE)
      
      if(!hasPhase2){
        ci95[c(2,4),] <<- NA
      }
      
    }
  )
  
  return(list(optimisedLL, hessianMat, ci95))
})


# Summarize results fits in one table -------------------------------------

summaryFit_df <- lapply(1:length(uniqueId_v), function(whichIndiv){
  estimates_df <- fittedModels_l[[whichIndiv]][3][[1]] %>% as.data.frame()
  estimates_df$what <- c("SS", "SC", "IS", "IC", "E")
  estimates_df$id <- uniqueId_v[whichIndiv]
  return(estimates_df) 
})
summaryFit_df <- do.call("rbind", summaryFit_df) %>% 
  as.data.frame() %>% 
  #Remove estimates equal to initialisation
  mutate(
    est = ifelse(est == 0.25, NA, est)
  )
