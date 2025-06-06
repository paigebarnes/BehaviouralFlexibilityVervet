#What does this script do?
#Analysis of PB masters thesis - check accuracy of our mechanistic model
# Posterior predictive check
# compare observed and simulated parameters (100 repetitions)
# parameter error based on different phase lengths


# Set up environment ------------------------------------------------------
rm(list = ls())
Sys.setenv(lang = "en_US")

local({
  Sys.setenv(LANGUAGE="en")
})

## Libraries ---------------------------------------------------------------

#Data processing
library(readr)
library(glmmTMB)
library(dplyr)
library(zoo)
library(pbapply)
library(ggplot2)
library(gridExtra)
library(grid)
#Fitting likelihood and parallelisation
library(optimParallel)
library(parallel)
library(numDeriv)
library(readr)
library(pbapply)


## Data --------------------------------------------------------------------
setwd("C:/Users/admin/Desktop/MSThesis")

dataExperimentBox_df <- read_delim("./1-Data/AttemptSides_21-05-24.csv",
                                   escape_double = FALSE, trim_ws = TRUE)
# dataExperimentBox_df <- read_delim("C:/Users/paige/Downloads/AttemptSides_21-05-24.csv", 
                                   # escape_double = FALSE, trim_ws = TRUE)
colnames(dataExperimentBox_df) <- c("id", "phase", "opening", "trial", "innovationL", "innovationP", "attemptL", "attemptP")

source("./2-Scripts/0_parameters.R")
# load("C:/Users/paige/Downloads/ParameterFit_bothsideattempts_24-06-24.RData")

## Functions ---------------------------------------------------------------

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


# Posterior Predictive Check
posteriorPredictiveCheck <- function(data, params, ph1count = NA, ph2count = NA, ph3count = NA, ph1successcount = NA) {
  
  #take real first choice from each monkey in phase 1
  first_trials <- data %>%
    group_by(id) %>%
    arrange(phase, trial) %>%
    slice(1)
  
  simulated_data <- data.frame(matrix(NA, nrow = 0, ncol = 20))
  names(simulated_data) <- c(names(data), "probability")
  
  counts <- data %>%
    group_by(id, phase) %>%
    summarise(count = n())
  
  for (i in 1:nrow(first_trials)){
    
    #set up counts of phase 2 and 3 based on real data
    if(is.na(ph2count)){
      ph2count <- counts$count[counts$phase == 2 & counts$id == first_trials$id[i]]
      ph2count<- ifelse(length(ph2count) == 0, 0, ph2count)
    }
    if(is.na(ph3count)){
      ph3count <- counts$count[counts$phase == 3 & counts$id == first_trials$id[i]]
      ph3count<- ifelse(length(ph3count) == 0, 0, ph3count)
    }
    
    # ph2count <- counts$count[counts$phase == 2 & counts$id == first_trials$id[i]]
    # ph3count <- counts$count[counts$phase == 3 & counts$id == first_trials$id[i]]
    # ph2count<- ifelse(length(ph2count) == 0, 0, ph2count)
    # ph3count<- ifelse(length(ph3count) == 0, 0, ph3count)
    
    #take only the first line
    data <- first_trials[i,]
    
    #extract parameters from real summary
    IS <- params$est[params$id == data$id[1] & params$what == "IS"]
    IC <- params$est[params$id == data$id[1] & params$what == "IC"]
    IC <- ifelse(!is.na(IC), IC, 0)
    SS <- params$est[params$id == data$id[1] & params$what == "SS"]
    SC <- params$est[params$id == data$id[1] & params$what == "SC"]
    SC <- ifelse(!is.na(SC), SC, 0)
    E <- params$est[params$id == data$id[1] & params$what == "E"]
    
    #set up possible trials (l, p, a)
    possible_trials <- data.frame(matrix(NA, nrow = 3, ncol = 20))
    names(possible_trials) <- c(names(data), "probability")
    possible_trials$phase <- 1
    
    while(possible_trials$phase[1] != 0){
      possible_trials$id <- data$id[nrow(data)]
      
      if(is.na(ph1count) & is.na(ph1successcount)){
        possible_trials$phase <- ifelse(length(data$opening[data$opening != "a" & data$phase == 1]) < 10 & length(data$opening[data$phase == 1]) < 22, 1,
                                        ifelse(length(data$opening[data$phase == 2]) < ph2count, 2,
                                               ifelse(length(data$opening[data$phase == 3]) < ph3count, 3, 0)))
      } else if ( !is.na(ph1count) & is.na(ph1successcount)) {
        possible_trials$phase <- ifelse(length(data$opening[data$phase == 1]) < ph1count, 1,
                                        ifelse(length(data$opening[data$phase == 2]) < ph2count, 2,
                                               ifelse(length(data$opening[data$phase == 3]) < ph3count, 3, 0)))
      } else {
        possible_trials$phase <- ifelse(length(data$opening[data$opening != "a" & data$phase == 1]) < ph1successcount & length(data$opening[data$phase == 1]) < 22, 1,
                                        ifelse(length(data$opening[data$phase == 2]) < ph2count, 2,
                                               ifelse(length(data$opening[data$phase == 3]) < ph3count, 3, 0)))
      }
      
      possible_trials$trial <- data$trial[nrow(data)] + 1
      
      possible_trials$innovationLS <- ifelse(data$innovationLS[nrow(data)] == 1 |
                                               (data$phase[nrow(data)] == 1 & data$opening[nrow(data)] == "l" ), 1, 0)
      possible_trials$innovationPS <- ifelse(data$innovationPS[nrow(data)] == 1 |
                                               (data$phase[nrow(data)] != 3 & data$opening[nrow(data)] == "p" ), 1, 0)
      possible_trials$innovationL <- ifelse(data$innovationL[nrow(data)] == 1 |
                                              (data$phase[nrow(data)] != 1 & data$opening[nrow(data)] == "l" ), 1, 0)
      possible_trials$innovationP <- ifelse(data$innovationP[nrow(data)] == 1 |
                                              (data$phase[nrow(data)] == 3 & data$opening[nrow(data)] == "p" ), 1, 0)
      
      possible_trials$previousOpening <- ifelse(data$opening[nrow(data)] == "l", "l", 
                                                ifelse(data$opening[nrow(data)] == "p", "p",
                                                       data$previousOpening[nrow(data)]))
      possible_trials$previousLift <- ifelse(data$opening[nrow(data)] == "l", 1, 0)
      possible_trials$previousPull <- ifelse(data$opening[nrow(data)] == "p", 1, 0)
      possible_trials$previousAttemptL <- ifelse(data$attemptL[nrow(data)] == "l", 1, 0)
      possible_trials$previousAttemptP <- ifelse(data$attemptP[nrow(data)] == "p", 1, 0)
      
      #potential trials:
      # success l
      # success p
      # attempt
      possible_trials$opening[1] <- "l"
      possible_trials$attempt[1] <- 0
      possible_trials$attemptBoth[1] <- 0
      possible_trials$switching[1] <- ifelse(possible_trials$previousPull[1] == 1, 1, 0)
      
      possible_trials$opening[2] <- "p"
      possible_trials$attempt[2] <- 0
      possible_trials$attemptBoth[2] <- 0
      possible_trials$switching[2] <- ifelse(possible_trials$previousLift[2] == 1, 1, 0)
      
      possible_trials$opening[3] <- "a"
      possible_trials$attempt[3] <- 1
      possible_trials$attemptBoth[3] <- 1
      possible_trials$attemptL[3] <- 1
      possible_trials$attemptP[3] <- 1
      
      possible_trials <- possible_trials %>% 
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
      #random probability
      p_trial <- runif(1,0,1)
      #check which possible trial the random probability corresponds to
      data <- if(possible_trials$phase[1] != 0){
        rbind(data, possible_trials[ifelse(p_trial > 0 & p_trial < possible_trials$probability[1], 1,
                                           ifelse(p_trial > possible_trials$probability[1] & p_trial < possible_trials$probability[1] + possible_trials$probability[2], 2,
                                                  3)), ])
      } else{data}
    }
    #update list of simulated trials
    simulated_data <- rbind(simulated_data, data)
    
  }
  #return simulated data
  return(simulated_data)
}


simulate_parameters <- function(data, params, repeats, ph1count = NA, ph2count = NA, ph3count = NA, ph1successcount =  NA) {
  
  #set up dataframe
  simulated_params <- data.frame(matrix(NA, nrow = 0, ncol = 5))
  names(simulated_params) <- names(params)
  
  #while loop to repeat process specified amounts of time
  j <- 1
  while(j <= repeats){
    #simulate new data
    simulated_data <- posteriorPredictiveCheck(data, params, ph1count, ph2count, ph3count, ph1successcount)
    
    
    ### Extract complementary descriptive covariates ------------------------------------
    
    complementaryCovariates_df <- simulated_data %>% 
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
    
    
    ### Fitting the model -------------------------------------------------------
    
    uniqueId_v <- unique(simulated_data$id)
    
    fittedModels_l <- pblapply(1:length(uniqueId_v), function(whichIndiv){
      #print(whichIndiv)
      dataOfIndividualInterest_df <- simulated_data %>% filter(id == uniqueId_v[whichIndiv])
      
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
                            data = simulated_data %>% filter(id == uniqueId_v[whichIndiv])
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
    #update the list of simulated parameters
    simulated_params <- rbind(simulated_params, summaryFit_df)
    j <- j+1 #for while loop count
  }
  return(simulated_params)
  
}


## Preparing the data to be modeled --------------------------------------

set.seed(42)

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
    
    #Reset standard lebelling of opening independent of innovation
    opening = ifelse(opening == "ul", "l", opening),
    opening = ifelse(opening == "up", "p", opening),
    #Record previous move (whether it is p or l, discard attempts a) after relabelling
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
# dataExperimentBox_df$innovationSimple <- 1
# dataExperimentBox_df$innovationSimple2 <- 1

#both sides should be counted in every attempt
dataExperimentBox_df$attemptBoth[dataExperimentBox_df$opening == "a"] <- 1
dataExperimentBox_df$attemptL[dataExperimentBox_df$opening == "a"] <- 1
dataExperimentBox_df$attemptP[dataExperimentBox_df$opening == "a"] <- 1




## Run the PPC simulated data through the mechanistic model #########################

# run simulations
sim_param_test <- simulate_parameters(data_test_trials, summaryFit_trail_testing, 100)
# save(sim_param_test, file = "C:/Users/admin/Desktop/MSThesis/2-Scripts/MechMod_trialsTest_DATE.RData")

## Run trial count simulations -----------------------------------------------------
# for trail testing
data_test_trials <- dataExperimentBox_df %>%
  filter(id %in% c("Ange", "Dra", "Gin")) %>%
  filter(phase == 1 & trial == 1)
summaryFit_trial_testing <- summaryFit_df %>%
  filter(id %in% c("Ange", "Dra", "Gin")) %>%
  mutate(est = 0.5) %>%
  select(est, what, id)

# running simulations for each type of start trial (l, p, a) with all params 0.5
# repeat with different trial counts for phases 1 and 2
sim_param_trialTestp1p2 <- data.frame(matrix(NA, nrow = 0, ncol = 8))
names(sim_param_trialTestp1p2) <- c(names(summaryFit_df), "ph1count", "ph2count", "ph3count") 

for (i in 1:15){
  for (j in 1:10) {
    sim_param_current <- simulate_parameters(data_test_trials, summaryFit_trial_testing, 50, i, j, 0)
    # assign(paste("sim_param_test_", i, "_", j, "_", 0), sim_param_test)
    sim_param_current$ph1count <- i
    sim_param_current$ph2count <- j
    sim_param_current$ph3count <- 0
    sim_param_trialTestp1p2 <- rbind(sim_param_trialTestp1p2, sim_param_current)
  }
}
# save(sim_param_trialTestp1p2, file = "C:/Users/paige/Desktop/MechModTrialTesting/sim_param_trialTestp1p2_13-09-24.RData")

# Success count variations in phase 1 and trial count variations phase 2
sim_param_p1successp2trial <- data.frame(matrix(NA, nrow = 0, ncol = 8))
names(sim_param_p1successp2trial) <- c(names(summaryFit_df), "ph1success", "ph2count", "ph3count") 
for (i in 1:15){
  for (j in 1:10) {
    sim_param_current <- simulate_parameters(data_test_trials, summaryFit_trial_testing, 50, NA, j, 0, i)
    # assign(paste("sim_param_test_", i, "_", j, "_", 0), sim_param_test)
    sim_param_current$ph1success <- i
    sim_param_current$ph2count <- j
    sim_param_current$ph3count <- 0
    sim_param_p1successp2trial <- rbind(sim_param_p1successp2trial, sim_param_current)
  }
}
# save(sim_param_p1successp2trial, file = "C:/Users/paige/Desktop/MechModTrialTesting/sim_param_p1successp2trial_13-09-24.RData")


# running simulations for each type of start trial (l, p, a) with all params 0.5
# repeat with different trial counts for phases 2 and 3
sim_param_trialTestp2p3 <- data.frame(matrix(NA, nrow = 0, ncol = 8))
names(sim_param_trialTestp2p3) <- c(names(summaryFit_df), "ph1count", "ph2count", "ph3count") 
for (j in 1:10){
  for (k in 1:10) {
    sim_param_current <- simulate_parameters(data_test_trials, summaryFit_trial_testing, 50, NA, j, k)
    sim_param_current$ph1count <- "original rule"
    sim_param_current$ph2count <- j
    sim_param_current$ph3count <- k
    sim_param_trialTestp2p3 <- rbind(sim_param_trialTestp2p3, sim_param_current)
  }
}
# save(sim_param_trialTestp2p3, file = "C:/Users/paige/Desktop/MechModTrialTesting/sim_param_trialTestp2p3_13-09-24.RData")



### sample 15 with init param set to observed mean for each one #######################################################


# for trail testing
data_test_trials <- dataExperimentBox_df %>%
  filter(id %in% c("Ange", "Dra", "Gin")) %>%
  filter(phase == 1 & trial == 1)
summaryFitMean_trial_testing <- summaryFit_df %>%
  group_by(what) %>%
  mutate(est = mean(est, na.rm = TRUE)) %>%
  filter(id %in% c("Ange", "Dra", "Gin")) %>%
  select(est, what, id)

# running simulations for each type of start trial (l, p, a) with all params their mean
# repeat with different trial counts for phases 1 and 2
sim_paramMean_trialTestp1p2 <- data.frame(matrix(NA, nrow = 0, ncol = 8))
names(sim_paramMean_trialTestp1p2) <- c(names(summaryFit_df), "ph1count", "ph2count", "ph3count") 

for (i in 1:15){
  for (j in 1:10) {
    sim_param_current <- simulate_parameters(data_test_trials, summaryFitMean_trial_testing, 5, i, j, 0)
    # assign(paste("sim_param_test_", i, "_", j, "_", 0), sim_param_test)
    sim_param_current$ph1count <- i
    sim_param_current$ph2count <- j
    sim_param_current$ph3count <- 0
    sim_paramMean_trialTestp1p2 <- rbind(sim_paramMean_trialTestp1p2, sim_param_current)
  }
}
# save(sim_paramMean_trialTestp1p2, file = "C:/Users/paige/Desktop/MechModTrialTesting/sim_paramMean_trialTestp1p2_24-09-24.RData")

# Success count variations in phase 1 and trial count variations phase 2
sim_paramMean_p1successp2trial <- data.frame(matrix(NA, nrow = 0, ncol = 8))
names(sim_paramMean_p1successp2trial) <- c(names(summaryFit_df), "ph1success", "ph2count", "ph3count") 
for (i in 1:15){
  for (j in 1:10) {
    sim_param_current <- simulate_parameters(data_test_trials, summaryFitMean_trial_testing, 5, NA, j, 0, i)
    # assign(paste("sim_param_test_", i, "_", j, "_", 0), sim_param_test)
    sim_param_current$ph1success <- i
    sim_param_current$ph2count <- j
    sim_param_current$ph3count <- 0
    sim_paramMean_p1successp2trial <- rbind(sim_paramMean_p1successp2trial, sim_param_current)
  }
}
# save(sim_paramMean_p1successp2trial, file = "C:/Users/paige/Desktop/MechModTrialTesting/sim_paramMean_p1successp2trial_24-09-24.RData")


# running simulations for each type of start trial (l, p, a) with all params 0.5
# repeat with different trial counts for phases 2 and 3
sim_paramMean_trialTestp2p3 <- data.frame(matrix(NA, nrow = 0, ncol = 8))
names(sim_paramMean_trialTestp2p3) <- c(names(summaryFit_df), "ph1count", "ph2count", "ph3count") 
for (j in 1:10){
  for (k in 1:10) {
    sim_param_current <- simulate_parameters(data_test_trials, summaryFitMean_trial_testing, 5, NA, j, k)
    sim_param_current$ph1count <- "original rule"
    sim_param_current$ph2count <- j
    sim_param_current$ph3count <- k
    sim_paramMean_trialTestp2p3 <- rbind(sim_paramMean_trialTestp2p3, sim_param_current)
  }
}
# save(sim_paramMean_trialTestp2p3, file = "C:/Users/paige/Desktop/MechModTrialTesting/sim_paramMean_trialTestp2p3_24-09-24.RData")


## model simulated vs. obs difference  -------------------------------------------

#Load in data if you start here
load("./MechModTest100_06-08-24.RData")
load("./ParameterFit_bothsideattempts_24-06-24.RData")

load("./MechModTest100_06-08-24.RData")
load("./ParameterFit_bothsideattempts_24-06-24.RData")


#compare the mean of each simulated parameter for each monkey to obs
compareParams <- sim_param_test %>%
  group_by(id, what) %>%
  summarise(est = mean(est, na.rm = TRUE),
            estLwr = mean(estLwr, na.rm = TRUE),
            estUpr = mean(estUpr, na.rm = TRUE)) %>%
  mutate(source = "sim") %>%
  rbind(summaryFit_df %>% mutate(source = "obs")) %>%
  merge(raiding_df %>% select("ID", "Group", "Age"), by.x = "id", by.y = "ID", all.x = TRUE) %>%
  mutate(id = ifelse(Age == "a", toupper(id), tolower(id)),
         Age = ifelse(Age == "a", "a", "j"),
         id = ifelse(id == "curiousboy", "cur",
                     ifelse(id == "cutepinky", "cute",
                            ifelse(id == "dentedchin", "den",
                                   ifelse(id == "goateegirl", "goat",
                                          ifelse(id == "spottedchin", "spot",
                                                 ifelse(id == "whitedot", "wdo",
                                                        ifelse(id == "whitepatches", "wpa", id)))))))) %>%
  arrange(Group, Age)

# check if there is a difference for each parameter between obs and sim
#IS
model_is <- glmmTMB(estc ~ source + (1|id), data = compareParams %>% 
                      filter(what == "IS") %>%
                      filter(!is.na(est)) %>%
                      group_by(id) %>%
                      mutate(n = n()) %>%
                      filter(n > 1) %>%
                      ungroup() %>%
                      mutate(estc = (est*(nrow(.) - 1) + 0.5)/nrow(.)),
                    family = beta_family(link = "logit"))
summary(model_is) # sig
#SS
model_ss <- glmmTMB(estc ~ source + (1|id), data = compareParams %>% 
                      filter(what == "SS") %>%
                      filter(!is.na(est)) %>%
                      group_by(id) %>%
                      mutate(n = n()) %>%
                      filter(n > 1) %>%
                      ungroup() %>%
                      mutate(estc = (est*(nrow(.) - 1) + 0.5)/nrow(.)),
                    family = beta_family(link = "logit"))
summary(model_ss) # sig
#IC
model_ic <- glmmTMB(estc ~ source + (1|id), data = compareParams %>% 
                      filter(what == "IC") %>%
                      filter(!is.na(est)) %>%
                      group_by(id) %>%
                      mutate(n = n()) %>%
                      filter(n > 1) %>%
                      ungroup() %>%
                      mutate(estc = (est*(nrow(.) - 1) + 0.5)/nrow(.)),
                    family = beta_family(link = "logit"))
summary(model_ic) # not sig
#SC
model_sc <- glmmTMB(estc ~ source + (1|id), data = compareParams %>% 
                      filter(what == "SC") %>%
                      filter(!is.na(est)) %>%
                      group_by(id) %>%
                      mutate(n = n()) %>%
                      filter(n > 1) %>%
                      ungroup() %>%
                      mutate(estc = (est*(nrow(.) - 1) + 0.5)/nrow(.)),
                    family = beta_family(link = "logit"))
summary(model_sc) # sig
#E
model_e <- glmmTMB(estc ~ source + (1|id), data = compareParams %>% 
                     filter(what == "E") %>%
                     filter(!is.na(est)) %>%
                     group_by(id) %>%
                     mutate(n = n()) %>%
                     filter(n > 1) %>%
                     ungroup() %>%
                     mutate(estc = (est*(nrow(.) - 1) + 0.5)/nrow(.)),
                   family = beta_family(link = "logit"))
summary(model_e) # sig


## plots -------------------------------------------------------------------------

#Load in data if you start here
load("./MechModTest100_06-08-24.RData")
load("./ParameterFit_bothsideattempts_24-06-24.RData")


### density plots ----------------------------------------------------------------

#example density graph
ggplot(sim_param_test[sim_param_test$id == "Ange" & sim_param_test$what == "E",], aes(x=est)) + 
  geom_density() +
  geom_vline(aes(xintercept=summaryFit_df$est[summaryFit_df$what == "E" & summaryFit_df$id == "Ange"]),
             color="blue", linetype="dashed", size=1) +
  theme_bw()

ecdf(sim_param_test$est[sim_param_test$id == "Ange" & sim_param_test$what == "E"])(summaryFit_df$est[summaryFit_df$what == "E" & summaryFit_df$id == "Ange"])

# make a dataframe comparing all of the density plots
densities <- summaryFit_df %>%
  filter(!is.na(est)) %>%
  rowwise() %>%
  mutate(sim_est_mean = mean(sim_param_test$est[sim_param_test$what == what & sim_param_test$id == id], na.rm = TRUE),
         sim_low_mean = mean(sim_param_test$estLwr[sim_param_test$what == what & sim_param_test$id == id], na.rm = TRUE),
         sim_high_mean = mean(sim_param_test$estUpr[sim_param_test$what == what & sim_param_test$id == id], na.rm = TRUE),
         density = ifelse(any(!is.na(sim_param_test$est[sim_param_test$id == id & sim_param_test$what == what])),
                          ecdf(sim_param_test$est[sim_param_test$id == id & sim_param_test$what == what & !is.na(sim_param_test$est)])(est),
                          NA))

compareParamsMeans <- compareParams %>%
  mutate(what = ifelse(what == "E", "R",
                       ifelse(what == "IS", "IS",
                              ifelse(what == "IC", "IT",
                                     ifelse(what == "SS", "SS",
                                            "ST"))))) %>%
  group_by(what, source) %>%
  dplyr::summarise(est = mean(est, na.rm = TRUE),
                   estLwr = mean(estLwr, na.rm = TRUE),
                   estUpr = mean(estUpr, na.rm = TRUE))


ggplot(compareParamsMeans, aes(x=source, y=est, color=source)) + 
  geom_point(aes(col = source)) +  
  geom_errorbar(aes(ymin = estLwr, ymax = estUpr), width = 0) +
  facet_wrap(. ~ what, ncol = 5) +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill="grey97")) +
  labs(y= "Estimate value") + 
  ggtitle("Parameters") +
  scale_color_manual(labels = c("empirical mean and CI95%", "simulated mean and CI95%"),
                     values = c("#88CCEE", "#882255"))

### plot real vs. sim param for each individual ---------------------------------

#IS
ggplot(compareParams[compareParams$what == "IS",], aes(x=source, y=est, color=source)) + 
  geom_point(aes(col = source)) +  
  geom_errorbar(aes(ymin = estLwr, ymax = estUpr), width = 0) +
  facet_wrap(. ~ fct_inorder(id), ncol=9) +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill="grey97"),
        legend.position = c(1, 0),
        legend.justification = c(1, 0))+
  labs(y= "Estimate value", x = NULL) + 
  ggtitle("IS per individual") +
  scale_color_manual(labels = c("empirical mean and CI95%", "simulated mean and CI95%"),
                     values = c("#88CCEE", "#882255"))

#SS
ggplot(compareParams[compareParams$what == "SS",], aes(x=source, y=est, color=source)) + 
  geom_point(aes(col = source)) +  
  geom_errorbar(aes(ymin = estLwr, ymax = estUpr), width = 0) +
  facet_wrap(. ~ fct_inorder(id), ncol=9) +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill="grey97"),
        legend.position = c(1, 0),
        legend.justification = c(1, 0))+
  labs(y= "Estimate value", x = NULL) + 
  ggtitle("SS per individual") +
  scale_color_manual(labels = c("empirical mean and CI95%", "simulated mean and CI95%"),
                     values = c("#88CCEE", "#882255"))

#IC
ggplot(compareParams[compareParams$what == "IC",], aes(x=source, y=est, color=source)) + 
  geom_point(aes(col = source)) +  
  geom_errorbar(aes(ymin = estLwr, ymax = estUpr), width = 0) +
  facet_wrap(. ~ fct_inorder(id), ncol=9) +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill="grey97"),
        legend.position = c(1, 0),
        legend.justification = c(1, 0))+
  labs(y= "Estimate value", x = NULL) + 
  ggtitle("IT per individual") +
  scale_color_manual(labels = c("empirical mean and CI95%", "simulated mean and CI95%"),
                     values = c("#88CCEE", "#882255"))

#SC
ggplot(compareParams[compareParams$what == "SC",], aes(x=source, y=est, color=source)) + 
  geom_point(aes(col = source)) +  
  geom_errorbar(aes(ymin = estLwr, ymax = estUpr), width = 0) +
  facet_wrap(. ~ fct_inorder(id), ncol=9) +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill="grey97"),
        legend.position = c(1, 0),
        legend.justification = c(1, 0))+
  labs(y= "Estimate value", x = NULL) + 
  ggtitle("ST per individual") +
  scale_color_manual(labels = c("empirical mean and CI95%", "simulated mean and CI95%"),
                     values = c("#88CCEE", "#882255"))

#E
ggplot(compareParams[compareParams$what == "E",], aes(x=source, y=est, color=source)) + 
  geom_point(aes(col = source)) +  
  geom_errorbar(aes(ymin = estLwr, ymax = estUpr), width = 0) +
  facet_wrap(. ~ fct_inorder(id), ncol=9) +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill="grey97"),
        legend.position = c(1, 0),
        legend.justification = c(1, 0))+
  labs(y= "Estimate value", x = NULL) + 
  ggtitle("R per individual") +
  scale_color_manual(labels = c("empirical mean and CI95%", "simulated mean and CI95%"),
                     values = c("#88CCEE", "#882255"))


### Heatmap of trial lengths -----------------------------------------------------

load("C:/Users/paige/Desktop/MechModTrialTesting/sim_param_trialTestp1p2_13-09-24.RData")
load("C:/Users/paige/Desktop/MechModTrialTesting/sim_param_trialTestp2p3_13-09-24.RData")
load("C:/Users/paige/Desktop/MechModTrialTesting/sim_param_p1successp2trial_13-09-24.RData")

#### Phases 1 and 2 compared ------------------------------------------------------

heat_df <- sim_param_trialTestp1p2 %>%
  group_by(id, what, ph1count, ph2count) %>%
  slice (1:33) %>%
  ungroup() %>%
  group_by(what, ph1count, ph2count) %>%
  mutate(est = ifelse(is.na(est), 0, est)) %>%
  summarise(err = mean(abs(0.5 - est), na.rm = TRUE),
            errNoAbs = mean(0.5-est, na.rm = TRUE)) %>%
  mutate(paramTitle = factor(ifelse(what == "IS", "IS",
                                    ifelse(what == "SS", "SS",
                                           ifelse(what == "IC", "IT",
                                                  ifelse(what == "SC", "ST",
                                                         ifelse(what == "E", "R",)))))))

#IS 
is12 <- ggplot(heat_df[heat_df$what == "IS",], aes(x = ph1count, y = ph2count, fill = err)) +
  geom_tile() +
  scale_fill_gradient(high = "#D81B60", low = "#1E88E5") +
  labs(x = "Phase 1 trials", y = "Phase 2 trials", title = "IS") +
  geom_hline(yintercept = 6.5) +
  geom_vline(xintercept = 10.5) +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  coord_fixed()+
  theme_bw()

#SS 
ss12 <- ggplot(heat_df[heat_df$what == "SS",], aes(x = ph1count, y = ph2count, fill = err)) +
  geom_tile() +
  scale_fill_gradient(high = "#D81B60", low = "#1E88E5") +
  labs(x = "Phase 1 trials", y = "Phase 2 trials", title = "SS") +
  geom_hline(yintercept = 6.5) +
  geom_vline(xintercept = 10.5) +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  coord_fixed()+
  theme_bw()

#IC
ic12 <- ggplot(heat_df[heat_df$what == "IC",], aes(x = ph1count, y = ph2count, fill = err)) +
  geom_tile() +
  scale_fill_gradient(high = "#D81B60", low = "#1E88E5") +
  labs(x = "Phase 1 trials", y = "Phase 2 trials", title = "IT") +
  geom_hline(yintercept = 6.5) +
  geom_vline(xintercept = 10.5) +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  coord_fixed()+
  theme_bw()

#SC 
sc12 <- ggplot(heat_df[heat_df$what == "SC",], aes(x = ph1count, y = ph2count, fill = err)) +
  geom_tile() +
  scale_fill_gradient(high = "#D81B60", low = "#1E88E5") +
  labs(x = "Phase 1 trials", y = "Phase 2 trials", title = "ST") +
  geom_hline(yintercept = 6.5) +
  geom_vline(xintercept = 10.5) +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  coord_fixed()+
  theme_bw()

#E
e12 <- ggplot(heat_df[heat_df$what == "E",], aes(x = ph1count, y = ph2count, fill = err)) +
  geom_tile() +
  scale_fill_gradient(high = "#D81B60", low = "#1E88E5") +
  labs(x = "Phase 1 trials", y = "Phase 2 trials", title = "R") +
  geom_hline(yintercept = 6.5) +
  geom_vline(xintercept = 10.5) +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  coord_fixed()+
  theme_bw()

# #all parameters
ggplot(heat_df, aes(x = ph1count, y = ph2count, fill = err)) +
  geom_tile() +
  scale_fill_gradient(high = "#D81B60", low = "#1E88E5") +
  facet_wrap(. ~ factor(paramTitle, c("SS", "IS", "ST", "IT", "R"))) +
  labs(x = "Phase 1 trial count",
       y = "Phase 2 trial count",
       title = "Phases 1 and 2 trial count heatmap",
       fill = "Difference from\ninitial value\n") +
  geom_hline(yintercept = 6.5) +
  geom_vline(xintercept = 10.5) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  coord_fixed()+
  theme_bw() +
  theme(strip.background = element_rect(fill="grey97"))

# #overall params
# any_param_heat <- sim_param_trialTestp1p2 %>%
#   group_by(id, what, ph1count, ph2count) %>%
#   slice (1:33) %>%
#   ungroup() %>%
#   group_by(ph1count, ph2count) %>%
#   mutate(est = ifelse(is.na(est), 0, est)) %>%
#   summarise(est_mean = mean(est, na.rm = TRUE),
#             err = abs(0.5 - est_mean))
# # 
# all12 <-ggplot(any_param_heat, aes(x = ph1count, y = ph2count, fill = err)) +
#   geom_tile() +
#   scale_fill_gradient(high = "#D81B60", low = "#1E88E5") +
#   labs(x = "Phase 1 trials", y = "Phase 2 trials", title = "Mean of all parameters")  +
#   geom_hline(yintercept = 6.5) +
#   geom_vline(xintercept = 10.5) +
#   coord_fixed() +
#   theme_bw()

grid.arrange(is12, ss12, ic12, sc12, e12, #all12,
             top = textGrob("Phases 1 and 2",gp=gpar(fontsize=20)))

#For mean init parameters (not 0.5)
heatMean_df <- sim_paramMean_trialTestp1p2 %>%
  group_by(what, ph1count, ph2count) %>%
  merge((summaryFitMean_trial_testing %>% mutate(origEst = est) %>% dplyr::select(origEst, what)), by = "what", all.x = TRUE) %>%
  mutate(est = ifelse(is.na(est), 0, est)) %>%
  group_by(what, ph1count, ph2count) %>%
  summarise(err = mean(abs(0.5 - est), na.rm = TRUE),
            errNoAbs = mean(0.5-est, na.rm = TRUE)) %>%
  mutate(paramTitle = factor(ifelse(what == "IS", "IS",
                                    ifelse(what == "SS", "SS",
                                           ifelse(what == "IC", "IT",
                                                  ifelse(what == "SC", "ST",
                                                         ifelse(what == "E", "R",)))))))

ggplot(heatMean_df, aes(x = ph1count, y = ph2count, fill = err)) +
  geom_tile() +
  scale_fill_gradient(high = "#D81B60", low = "#1E88E5") +
  facet_wrap(. ~ factor(paramTitle, c("SS", "IS", "ST", "IT", "R"))) +
  labs(x = "Phase 1 trials", y = "Phase 2 trials", title = "Phases 1 and 2 heatmap") +
  geom_hline(yintercept = 6.5) +
  geom_vline(xintercept = 10.5) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  coord_fixed()+
  theme_bw()



#### Phase 1 successes and phase 2 trials compared -------------------------------

heat_1s2_df <- sim_param_p1successp2trial %>%
  group_by(id, what, ph1success, ph2count) %>%
  ungroup() %>%
  group_by(what, ph1success, ph2count) %>%
  mutate(est = ifelse(is.na(est), 0, est)) %>%
  summarise(err = mean(abs(0.5 - est), na.rm = TRUE),
            errNoAbs = mean(0.5-est, na.rm = TRUE)) %>%
  mutate(paramTitle = factor(ifelse(what == "IS", "IS",
                                    ifelse(what == "SS", "SS",
                                           ifelse(what == "IC", "IT",
                                                  ifelse(what == "SC", "ST",
                                                         ifelse(what == "E", "R",)))))))

#IS 
is1s2 <- ggplot(heat_1s2_df[heat_1s2_df$what == "IS",], aes(x = ph1success, y = ph2count, fill = err)) +
  geom_tile() +
  scale_fill_gradient(high = "#D81B60", low = "#1E88E5") +
  labs(x = "Phase 1 successes", y = "Phase 2 trials", title = "IS") +
  geom_hline(yintercept = 6.5) +
  geom_vline(xintercept = 10.5) +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  coord_fixed()+
  theme_bw()

#SS 
ss1s2 <- ggplot(heat_1s2_df[heat_1s2_df$what == "SS",], aes(x = ph1success, y = ph2count, fill = err)) +
  geom_tile() +
  scale_fill_gradient(high = "#D81B60", low = "#1E88E5") +
  labs(x = "Phase 1 successes", y = "Phase 2 trials", title = "SS") +
  geom_hline(yintercept = 6.5) +
  geom_vline(xintercept = 10.5) +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  coord_fixed()+
  theme_bw()

#IC
ic1s2 <- ggplot(heat_1s2_df[heat_1s2_df$what == "IC",], aes(x = ph1success, y = ph2count, fill = err)) +
  geom_tile() +
  scale_fill_gradient(high = "#D81B60", low = "#1E88E5") +
  labs(x = "Phase 1 successes", y = "Phase 2 trials", title = "IT") +
  geom_hline(yintercept = 6.5) +
  geom_vline(xintercept = 10.5) +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  coord_fixed()+
  theme_bw()

#SC 
sc1s2 <- ggplot(heat_1s2_df[heat_1s2_df$what == "SC",], aes(x = ph1success, y = ph2count, fill = err)) +
  geom_tile() +
  scale_fill_gradient(high = "#D81B60", low = "#1E88E5") +
  labs(x = "Phase 1 successes", y = "Phase 2 trials", title = "ST") +
  geom_hline(yintercept = 6.5) +
  geom_vline(xintercept = 10.5) +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  coord_fixed()+
  theme_bw()

#E
e1s2 <- ggplot(heat_1s2_df[heat_1s2_df$what == "E",], aes(x = ph1success, y = ph2count, fill = err)) +
  geom_tile() +
  scale_fill_gradient(high = "#D81B60", low = "#1E88E5") +
  labs(x = "Phase 1 successes", y = "Phase 2 trials", title = "R") +
  geom_hline(yintercept = 6.5) +
  geom_vline(xintercept = 10.5) +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  coord_fixed()+
  theme_bw()

# #all parameters
ggplot(heat_1s2_df, aes(x = ph1success, y = ph2count, fill = err)) +
  geom_tile() +
  scale_fill_gradient(high = "#D81B60", low = "#1E88E5") +
  facet_wrap(. ~ factor(paramTitle, c("SS", "IS", "ST", "IT", "R"))) +
  labs(x = "Phase 1 success requirement",
       y = "Phase 2 trial count",
       title = "Phases 1 success and 2 trial count heatmap",
       fill = "Difference from\ninitial value\n") +
  geom_hline(yintercept = 6.5) +
  geom_vline(xintercept = 10.5) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  coord_fixed()+
  theme_bw() +
  theme(strip.background = element_rect(fill="grey97"))

# # overall params
# any_param_heat_1s2 <- sim_param_p1successp2trial %>%
#   group_by(id, what, ph1success, ph2count) %>%
#   slice (1:33) %>%
#   ungroup() %>%
#   group_by(ph1success, ph2count) %>%
#   mutate(est = ifelse(is.na(est), 0, est)) %>%
#   summarise(est_mean = mean(est, na.rm = TRUE),
#             err = abs(0.5 - est_mean))
# 
# all1s2 <-ggplot(any_param_heat_1s2, aes(x = ph1success, y = ph2count, fill = err)) +
#   geom_tile() +
#   scale_fill_gradient(high = "#D81B60", low = "#1E88E5") +
#   labs(x = "Phase 1 successes", y = "Phase 2 trials", title = "Mean of all parameters")  +
#   geom_hline(yintercept = 6.5) +
#   geom_vline(xintercept = 10.5) +
#   coord_fixed() +
#   theme_bw()

grid.arrange(is1s2, ss1s2, ic1s2, sc1s2, e1s2, #all1s2,
             top = textGrob("Phases 1 successes and 2 trials",gp=gpar(fontsize=20)))

#For mean init parameters
heatMean_1s2_df <- sim_paramMean_p1successp2trial %>%
  group_by(id, what, ph1success, ph2count) %>%
  ungroup() %>%
  group_by(what, ph1success, ph2count) %>%
  mutate(est = ifelse(is.na(est), 0, est)) %>%
  summarise(err = mean(abs(0.5 - est), na.rm = TRUE),
            errNoAbs = mean(0.5-est, na.rm = TRUE)) %>%
  mutate(paramTitle = factor(ifelse(what == "IS", "IS",
                                    ifelse(what == "SS", "SS",
                                           ifelse(what == "IC", "IT",
                                                  ifelse(what == "SC", "ST",
                                                         ifelse(what == "E", "R",)))))))


ggplot(heatMean_1s2_df, aes(x = ph1success, y = ph2count, fill = err)) +
  geom_tile() +
  scale_fill_gradient(high = "#D81B60", low = "#1E88E5") +
  facet_wrap(. ~ factor(paramTitle, c("SS", "IS", "ST", "IT", "R"))) +
  labs(x = "Phase 1 trials", y = "Phase 2 trials", title = "Phases 1 and 2 heatmap") +
  geom_hline(yintercept = 6.5) +
  geom_vline(xintercept = 10.5) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  coord_fixed()+
  theme_bw()




#### Phases 2 and 3 compared -----------------------------------------------------

heat2_3_df <- sim_param_trialTestp2p3 %>%
  group_by(what, ph2count, ph3count) %>%
  mutate(est = ifelse(is.na(est), 0, est)) %>%
  summarise(err = mean(abs(0.5 - est), na.rm = TRUE),
            errNoAbs = mean(0.5-est, na.rm = TRUE)) %>%
  mutate(paramTitle = factor(ifelse(what == "IS", "IS",
                                    ifelse(what == "SS", "SS",
                                           ifelse(what == "IC", "IT",
                                                  ifelse(what == "SC", "ST",
                                                         ifelse(what == "E", "R",)))))))

#IS 
is23 <- ggplot(heat2_3_df[heat2_3_df$what == "IS",], aes(x = ph2count, y = ph3count, fill = err)) +
  geom_tile() +
  scale_fill_gradient(high = "#D81B60", low = "#1E88E5") +
  labs(x = "Phase 2 trials", y = "Phase 3 trials", title = "IS") +
  geom_hline(yintercept = 6.5) +
  geom_vline(xintercept = 6.5) +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  coord_fixed()+
  theme_bw()

#SS 
ss23 <- ggplot(heat2_3_df[heat2_3_df$what == "SS",], aes(x = ph2count, y = ph3count, fill = err)) +
  geom_tile() +
  scale_fill_gradient(high = "#D81B60", low = "#1E88E5") +
  labs(x = "Phase 2 trials", y = "Phase 3 trials", title = "SS") +
  geom_hline(yintercept = 6.5) +
  geom_vline(xintercept = 6.5) +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  coord_fixed()+
  theme_bw()

#IC
ic23 <- ggplot(heat2_3_df[heat2_3_df$what == "IC",], aes(x = ph2count, y = ph3count, fill = err)) +
  geom_tile() +
  scale_fill_gradient(high = "#D81B60", low = "#1E88E5") +
  labs(x = "Phase 2 trials", y = "Phase 3 trials", title = "IT") +
  geom_hline(yintercept = 6.5) +
  geom_vline(xintercept = 6.5) +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  coord_fixed()+
  theme_bw()

#SC 
sc23 <- ggplot(heat2_3_df[heat2_3_df$what == "SC",], aes(x = ph2count, y = ph3count, fill = err)) +
  geom_tile() +
  scale_fill_gradient(high = "#D81B60", low = "#1E88E5") +
  labs(x = "Phase 2 trials", y = "Phase 3 trials", title = "ST") +
  geom_hline(yintercept = 6.5) +
  geom_vline(xintercept = 6.5) +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  coord_fixed()+
  theme_bw()

#E
e23 <- ggplot(heat2_3_df[heat2_3_df$what == "E",], aes(x = ph2count, y = ph3count, fill = err)) +
  geom_tile() +
  scale_fill_gradient(high = "#D81B60", low = "#1E88E5") +
  labs(x = "Phase 2 trials", y = "Phase 3 trials", title = "R") +
  geom_hline(yintercept = 6.5) +
  geom_vline(xintercept = 6.5) +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  coord_fixed()+
  theme_bw()

# #all parameters
ggplot(heat2_3_df, aes(x = ph2count, y = ph3count, fill = err)) +
  geom_tile() +
  scale_fill_gradient(high = "#D81B60", low = "#1E88E5") +
  facet_wrap(. ~ factor(paramTitle, c("SS", "IS", "ST", "IT", "R"))) +
  labs(x = "Phase 2 trial count",
       y = "Phase 3 trial count",
       title = "Phases 2 and 3 trial count heatmap",
       fill = "Difference from\ninitial value\n") +
  geom_hline(yintercept = 6.5) +
  geom_vline(xintercept = 6.5) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  coord_fixed()+
  theme_bw() +
  theme(strip.background = element_rect(fill="grey97"))

# overall params
# any_param_heat2_3 <- sim_param_trialTestp2p3 %>%
#   group_by(id, what, ph2count, ph3count) %>%
#   slice (1:33) %>%
#   ungroup() %>%
#   group_by(ph2count, ph3count) %>%
#   mutate(est = ifelse(is.na(est), 0, est)) %>%
#   summarise(est_mean = mean(est, na.rm = TRUE),
#             err = abs(0.5 - est_mean))
# all23 <- ggplot(any_param_heat2_3, aes(x = ph2count, y = ph3count, fill = err)) +
#   geom_tile() +
#   scale_fill_gradient(high = "#D81B60", low = "#1E88E5") +
#   labs(x = "Phase 2 trials", y = "Phase 3 trials", title = "Mean of all parameters") +
#   geom_hline(yintercept = 6.5) +
#   geom_vline(xintercept = 6.5) +
#   coord_fixed()+
#   theme_bw()

grid.arrange(is23, ss23, ic23, sc23, e23, #all23,
             top = textGrob("Phases 2 and 3",gp=gpar(fontsize=20)))


#For mean init parameters

heatMean2_3_df <- sim_paramMean_trialTestp2p3 %>%
  group_by(what, ph2count, ph3count) %>%
  mutate(est = ifelse(is.na(est), 0, est)) %>%
  summarise(err = mean(abs(0.5 - est), na.rm = TRUE),
            errNoAbs = mean(0.5-est, na.rm = TRUE)) %>%
  mutate(paramTitle = factor(ifelse(what == "IS", "IS",
                                    ifelse(what == "SS", "SS",
                                           ifelse(what == "IC", "IT",
                                                  ifelse(what == "SC", "ST",
                                                         ifelse(what == "E", "R",)))))))

ggplot(heatMean2_3_df, aes(x = ph2count, y = ph3count, fill = err)) +
  geom_tile() +
  scale_fill_gradient(high = "#D81B60", low = "#1E88E5") +
  facet_wrap(. ~ factor(paramTitle, c("SS", "IS", "ST", "IT", "R"))) +
  labs(x = "Phase 2 trials", y = "Phase 3 trials", title = "Phases 2 and 3 heatmap") +
  geom_hline(yintercept = 6.5) +
  geom_vline(xintercept = 6.5) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  coord_fixed()+
  theme_bw()
