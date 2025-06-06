#What does this script do?
#Analysis of PB et al. manuscript - box success vs. parameters and human food metric

rm(list = ls())
Sys.setenv(lang = "en_US")

# Set up environment -------------------------------------------------------


## Libraries ---------------------------------------------------------------

#Statistical modelling
library(AICcmodavg)
library(glmmTMB)
library(effects)
library(ggplot2)
library(forcats)
library(DHARMa)
library("extrafont")
library(plot3D)
library(performance)
library(sjPlot)
library(librarian)
library(car)
library("ggpubr")

#Data processing
library(dplyr)
library("readxl")
library(tidyr)
library(readr)


## Functions -----------------------------------------------------------------

diagnostics.plot.dharma <-
  function(mod.res,
           col = grey(level = 0.25, alpha = 0.5),
           breaks.histo = 20,
           quantreg = TRUE) {
    old.par = par(no.readonly = TRUE)
    par(mfrow = c(2, 2))
    par(mar = c(3, 3, 3, 0.5))
    hist(
      residuals(mod.res),
      probability = T,
      xlab = "",
      ylab = "",
      main = "",
      breaks = breaks.histo
    )
    mtext(text = "Histogram of residuals",
          side = 3,
          line = 0)
    x = seq(min(residuals(mod.res)), max(residuals(mod.res)), length.out =
              100)
    lines(x, dnorm(x, mean = 0, sd = sd(residuals(mod.res))))
    
    library(DHARMa)
    simulationOutput <-
      simulateResiduals(fittedModel = mod.res, plot = FALSE)
    
    plotQQunif(simulationOutput) # left plot in plot.DHARMa()
    plotResiduals(simulationOutput, quantreg = quantreg)
    #Old way without dharma, from Roger Mundry
    # qqnorm(residuals(mod.res), main="", pch=19)
    # qqline(residuals(mod.res))
    # mtext(text="qq-plot of residuals", side=3, line=0)
    # plot(fitted(mod.res), residuals(mod.res), pch=19, col=col)
    # abline(h=0, lty=2)
    # mtext(text="residuals against fitted values", side=3, line=0)
    par(old.par)
  }

#Function to check linear model validity
check_model <- function(#Still in construction not tested for all model types
  model,
  hasInteraction = FALSE,
  modelWithoutInteractionForVif = NA,
  computedfBetas = FALSE,
  fittedWith = c("lm", "lmer", "glm", "glmer", "glmmTMB")
){
  
  shelf(performance)
  shelf(lme4)
  # shelf(influence.ME)
  
  #Convergence
  singularityTest <- check_singularity(model)
  convergenceTest <- check_convergence(model)
  
  if(fittedWith %in% c("lm", "lmer")){
    family = "gaussian"
  }else{
    family = family(model)$family
  }
  #Overdispersion
  if(family[1] == "poisson"){
    overdispersionTest <- check_overdispersion(model)
  }else{
    overdispersionTest <- NA
  }
  
  #Autocorr
  autocorrTest <- check_autocorrelation(model)
  
  #Check correlation var (VIF)
  if(hasInteraction){
    collinearityTest <- check_collinearity(modelWithoutInteractionForVif)
  }else{
    collinearityTest <- check_collinearity(model)
  }
  
  #Check dfbetas
  if(computedfBetas){
    if(fittedWith == "lm" | fittedWith == "glm"){
      dfbetas <- 
        round(cbind(coefficients(model), 
                    coefficients(model)+t(apply(X=dfbeta(model), 
                                                MARGIN=2, FUN=range))), 5)
      #Check outliers
      outlierTest <- check_outliers(model)
    }else if(fittedWith == "lmer" | fittedWith == "glmer"){
      dfbetas <- round(dfbetas(influence.ME::influence(model, obs=TRUE)), 2) 
      dfbetas <- as.data.frame(cbind(
        rownames(summary(model)$coefficients),
        summary(model)$coefficients[,1],
        summary(model)$coefficients[,1] + apply(dfbetas, 2, min),
        summary(modelDistance)$coefficients[,1] + apply(dfbetas, 2, max)
      )
      )
      colnames(dfbetas) <- c("Estimate", "Min", "Max")
      print("Dfbetas are computed at the singular value level. You may want to run it at the random level too with the influence() function of the lme4 package")
      #Check outliers
      outlierTest <- check_outliers(model)
    }else{
      randomEffects <- ifelse(length(ranef(model)$cond) > 0, TRUE, FALSE)
      
      # testData <- data.frame(
      #   output = rbinom(100, 1, 0.3),
      #   pred1 = runif(100, 0, 1),
      #   random = rep(c("A", "B"), times = 50)
      # )
      # library(glmmTMB)
      # model <- glmmTMB(output ~ pred1 + (1|random), data = testData, family = binomial)
      # 
      if(randomEffects){
        print("Dfbetas are computed at the singular value level. You may want to run it at the random level too with the influence() function of the lme4 package")
        dfToTest <- model$frame
        if(family == "gaussian"){
          modelnonglmmTMB <- lmer(formula = formula(model), data = dfToTest)
        }else{
          modelnonglmmTMB <- glmer(formula = formula(model), data = dfToTest, family = family(model))
        }
        #Check outliers
        outlierTest <- check_outliers(modelnonglmmTMB)
        dfbetas <- round(dfbetas(influence.ME::influence(modelnonglmmTMB, data = dfToTest, obs=TRUE)), 2) 
        dfbetas <- as.data.frame(cbind(
          rownames(summary(modelnonglmmTMB)$coefficients),
          summary(modelnonglmmTMB)$coefficients[,1],
          summary(modelnonglmmTMB)$coefficients[,1] + apply(dfbetas, 2, min),
          summary(modelnonglmmTMB)$coefficients[,1] + apply(dfbetas, 2, max)
        )
        )
      }else{
        if(family == "gaussian"){
          modelnonglmmTMB <- lm(formula = formula(model), data = model$frame)
        }else{
          modelnonglmmTMB <- glm(formula = formula(model), data = model$frame, family = family(model))
        }
        #Check outliers
        outlierTest <- check_outliers(modelnonglmmTMB)
        dfbetas <- 
          round(cbind(coefficients(modelnonglmmTMB), 
                      coefficients(modelnonglmmTMB)+t(apply(X=dfbeta(modelnonglmmTMB), 
                                                            MARGIN=2, FUN=range))), 2)
      }
    }
    colnames(dfbetas) <- c("Estimate", "Min", "Max")
  }else{
    dfbetas <- NA
  }
  
  
  #Check assumptions
  dharmaOutput <- diagnostics.plot.dharma(model)
  if(family == "gaussian"){
    heteroskTest <- check_heteroscedasticity(model)
  }else{
    heteroskTest <- NA
  }
  
  output <- list(
    singularityTest,
    convergenceTest,
    overdispersionTest,
    autocorrTest,
    collinearityTest,
    outlierTest,
    dfbetas,
    heteroskTest,
    dharmaOutput
  )
  names(output) <- 
    c(
      "singularityTest",
      "convergenceTest",
      "overdispersionTest",
      "autocorrTest",
      "collinearityTest",
      "outlierTest",
      "dfbetas",
      "heteroskTest",
      "dharmaOutput"
    )
  return(
    output
  )
}


## Data --------------------------------------------------------------------

#success rates will come from experiment data
setwd("C:/Users/admin/Desktop/MSThesis/1-Data")

dataExperimentBox_df <- read_delim("./AttemptSides_21-05-24.csv", 
                                   escape_double = FALSE, trim_ws = TRUE)
colnames(dataExperimentBox_df) <- c("id", "phase", "opening", "trial", "innovationL", "innovationP", "attemptL", "attemptP")

setwd("C:/Users/admin/Desktop/MSThesis/2-Scripts")
source("./0_parameters.R")
# load("./HFadlibrates_28-06-24.RData")
#Human food data comes from march 5th - june 20th
load("./HFadlibrates_03-12-24.RData")
hf_rate <- hfrate_df
#human food data comes from July 11 2023 to Sep 30 2024

load("./ParameterFit_bothsideattempts_24-06-24.RData")



# Analysis ----------------------------------------------------------------

## Data processing ---------------------------------------------------------

#if not an attempt, it's a success
dataExperimentBox_df$success <- 1
dataExperimentBox_df$success[dataExperimentBox_df$opening == "a"] <- 0

#success rates for each phase
popsuc1 <- sum(dataExperimentBox_df$success[dataExperimentBox_df$phase == 1]) /
  length(dataExperimentBox_df$success[dataExperimentBox_df$phase == 1])
popsuc2 <- sum(dataExperimentBox_df$success[dataExperimentBox_df$phase == 2]) /
  length(dataExperimentBox_df$success[dataExperimentBox_df$phase == 2])
popsuc3 <- sum(dataExperimentBox_df$success[dataExperimentBox_df$phase == 3]) /
  length(dataExperimentBox_df$success[dataExperimentBox_df$phase == 3])

#weighting different phases based on group success rates
w1 <- 1
w2 <- popsuc1 / popsuc2
w3 <- popsuc1 / popsuc3

#scale so sum of all sums to 1
w1scale <- w1 / sum(w1, w2, w3)
w2scale <- w2 / sum(w1, w2, w3)
w3scale <- w3 / sum(w1, w2, w3)

#scale so sum of all sums to 1, excluding phase3
w1scale_p1p2 <- w1 / sum(w1, w2)
w2scale_p1p2 <- w2 / sum(w1, w2)

#change phase name to make summaries easier
dataExperimentBox_df$phase[dataExperimentBox_df$phase == 1] <- "phase1"
dataExperimentBox_df$phase[dataExperimentBox_df$phase == 2] <- "phase2"
dataExperimentBox_df$phase[dataExperimentBox_df$phase == 3] <- "phase3"


#success rates for each individual
indivsuc <- dataExperimentBox_df %>%
  group_by(id, phase) %>%
  summarise(sucrateph = sum(success) / length(success)) %>%
  ungroup() %>%
  spread(key = phase, value = sucrateph) %>%
  group_by(id) %>%
  summarise(boxsucrate = (w1scale * ifelse(is.na(phase1), 0, phase1)) +
              (w2scale * ifelse(is.na(phase2), 0, phase2)) +
              (w3scale * ifelse(is.na(phase3), 0, phase3)),
            boxsucrate_p1p2 = (w1scale_p1p2 * ifelse(is.na(phase1), 0, phase1)) +
              (w2scale_p1p2 * ifelse(is.na(phase2), 0, phase2)))

# exposuresimple_df <- exposure_df %>%
#   group_by(id) %>%
#   filter(trial == min(trial, na.rm = TRUE))

sucCompiled_df <- summaryFit_df %>%
  dplyr::select(id, what, est) %>%
  spread(key = what, value = est) %>%
  merge(raiding_df, by.x = "id", by.y = "ID", all.x = TRUE) %>%
  merge(exposure_df, by = "id", all.x = TRUE) %>%
  merge(rank_df, by = "id", all.x = TRUE) %>%
  merge(indivsuc, by = "id", all.x = TRUE) %>%
  #merge(humanFoodrates, by = "id", all.x = TRUE) %>%
  merge(hf_rate, by.x = "id", by.y = "ID", all.x = TRUE)



### Fit model ---------------------------------------------------------------

#add 0s where the individual never got human food during focal
# sucCompiled_df$HF_dietpercent[is.na(sucCompiled_df$HF_dietpercent)] <- 0
# sucCompiled_df$HF_durrate[is.na(sucCompiled_df$HF_durrate)] <- 0
# sucCompiled_df$HF_typeeventrate[is.na(sucCompiled_df$HF_typeeventrate)] <- 0
# sucCompiled_df$HF_dayeventrate[is.na(sucCompiled_df$HF_dayeventrate)] <- 0

# sucCompiled_df$hf_rate <- datawizard::standardise(sucCompiled_df$hf_rate)
# sucCompiled_df$IC <- datawizard::standardise(sucCompiled_df$IC)
# sucCompiled_df$SS <- datawizard::standardise(sucCompiled_df$SS)
# sucCompiled_df$E <- datawizard::standardise(sucCompiled_df$E)

# model <- glmmTMB(boxsucratec ~ IC + SS + hf_rate, data = sucCompiled_df %>% 
#                    filter(!is.na(IC)) %>% 
#                    mutate(SSc = (SS*(nrow(.) - 1) + 0.5)/nrow(.), 
#                           ICc = (IC*(nrow(.) - 1) + 0.5)/nrow(.),
#                           boxsucratec = (boxsucrate*(nrow(.) - 1) + 0.5)/nrow(.),
#                           hf_ratec = (hf_rate*(nrow(.) - 1) + 0.5)/nrow(.)), 
#                  family = beta_family(link = "logit"))
# summary(model)
# 
# 
# model2 <- glmmTMB(boxsucrate_p1p2c ~ IC + SS + hf_rate, data = sucCompiled_df %>% 
#                    filter(!is.na(IC)) %>% 
#                    mutate(SSc = (SS*(nrow(.) - 1) + 0.5)/nrow(.), 
#                           ICc = (IC*(nrow(.) - 1) + 0.5)/nrow(.),
#                           boxsucrate_p1p2c = (boxsucrate_p1p2*(nrow(.) - 1) + 0.5)/nrow(.),
#                           hf_ratec = (hf_rate*(nrow(.) - 1) + 0.5)/nrow(.)), 
#                  family = beta_family(link = "logit"))
# summary(model2)
# 
# 
# model3 <- glmmTMB(boxsucrate_p1p2c ~ 1, data = sucCompiled_df %>% 
#                     filter(!is.na(IC)) %>% 
#                     mutate(SSc = (SS*(nrow(.) - 1) + 0.5)/nrow(.), 
#                            ICc = (IC*(nrow(.) - 1) + 0.5)/nrow(.),
#                            boxsucrate_p1p2c = (boxsucrate_p1p2*(nrow(.) - 1) + 0.5)/nrow(.),
#                            hf_ratec = (hf_rate*(nrow(.) - 1) + 0.5)/nrow(.)), 
#                   family = beta_family(link = "logit"))
# summary(model3)
# 
# anova(model3, model2)

### New analyses ###

#should these just compare individuals that we have a IC score for?

#pearson correlation between human food and P1P2 box success
cor.test(sucCompiled_df$hf_rate, sucCompiled_df$boxsucrate_p1p2, method ="pearson", exact = FALSE)

#rank version of correlation
sucCompiled_df$hf_rank <- order(sucCompiled_df$hf_rate)
sucCompiled_df$boxsucrate_p1p2_rank <- order(sucCompiled_df$boxsucrate_p1p2)
cor.test(sucCompiled_df$hf_rank, sucCompiled_df$boxsucrate_p1p2_rank, method ="spearman")

#box success model
boxSucmodel <- glmmTMB(boxsucratec ~ IC + poly(SS, 2), data = sucCompiled_df %>% 
                   filter(!is.na(IC)) %>% 
                   mutate(boxsucratec = (boxsucrate*(nrow(.) - 1) + 0.5)/nrow(.)), 
                 family = beta_family(link = "logit"))
summary(boxSucmodel)

nullmodBS <- glmmTMB(boxsucratec ~ 1, data = sucCompiled_df %>% 
                       filter(!is.na(IC)) %>% 
                       mutate(boxsucratec = (boxsucrate*(nrow(.) - 1) + 0.5)/nrow(.)), 
                     family = beta_family(link = "logit"))
anova(boxSucmodel, nullmodBS)



#human food success model
HFmodel <- glmmTMB(hf_ratec ~ IC + poly(SS, 2), data = sucCompiled_df %>% 
                         filter(!is.na(IC)) %>% 
                         mutate(hf_ratec = (hf_rate*(nrow(.) - 1) + 0.5)/nrow(.)), 
                       family = beta_family(link = "logit"))
summary(HFmodel)

nullmodHF <- glmmTMB(hf_ratec ~ 1, data = sucCompiled_df %>% 
                       filter(!is.na(IC)) %>% 
                       mutate(hf_ratec = (hf_rate*(nrow(.) - 1) + 0.5)/nrow(.)), 
                     family = beta_family(link = "logit"))
anova(HFmodel, nullmodHF)



## Plot and check models ####

### Box success ####
variables <- c("IC", "poly(SS, 2)1", "poly(SS, 2)2")

var_weightBS <- data.frame(var = variables,
                         weight = rep(NA, 3),
                         lowerCI = rep(NA, 3),
                         upperCI = rep(NA, 3))

Nobs <- nrow(sucCompiled_df)
summodbS <- summary(boxSucmodel)
#do this with every variable
for (i in (1:length(variables))){
  var_weightBS$weight[var_weightBS$var == variables[i]] <- as.numeric(summodbS$coefficients$cond[i+1])
  var_weightBS$lowerCI[var_weightBS$var == variables[i]] <- as.numeric(summodbS$coefficients$cond[i+1] + qt(0.05/2, Nobs - 1)*summodbS$coefficients$cond[i + 6])
  var_weightBS$upperCI[var_weightBS$var == variables[i]] <- as.numeric(summodbS$coefficients$cond[i+1] - qt(0.05/2, Nobs - 1)*summodbS$coefficients$cond[i + 6])
  
}


var_weightBS$category <- factor(c("IT", 
                         "SS", "SS2"), ordered = TRUE)


#plot the weights with confidence intervals
pm_colors <- c("#6e3323", "lightskyblue3")

pm <- plot_models(boxSucmodel, HFmodel,
                  transform = NULL,
                  show.p = TRUE,
                  p.shape = FALSE,
                  show.values = TRUE,
                  m.labels=c("Problem-solving success", "Anthropogenic food consumption"),
                  legend.title = "Models",
                  legend.pval.title = NULL,
                  # title = "Parameters on success models",
                  colors = pm_colors) +
  theme_bw() +
  theme(text = element_text(size = 14)) +
  geom_hline(yintercept= 0, colour = "black")
#change names
pm$data$term <- factor(pm$data$term)
pm <- pm + scale_x_discrete(labels = c(expression(SS^2), "SS", "IT"))
pm

#sjPlot plot_model

tab_model(boxSucmodel, HFmodel, show.p = FALSE)

ggplot(var_weightBS, aes(x = fct_inorder(category), y = weight, color = category)) +
  geom_point(aes(size = 3)) +  
  geom_errorbar(aes(ymin = lowerCI, ymax = upperCI), width = 0) +
  geom_hline(yintercept = 0) + 
  theme_set(theme_bw(base_family = 'Calibri')) + 
  theme(axis.text = element_text(size=20),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 20),
        legend.position="none") +
  scale_size(guide = 'none') +
  ylab("Variable estimate") +
  scale_color_manual(values = c("IT" = "#6c719c",
                                "SS" = "#b5739d",
                                "SS2" = "#6e3323",
                                "R" = "#81B29A")) +
  coord_flip() +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 20))


library(sjPlot)
plot_model(boxSucmodel)


checksBS <- check_model(boxSucmodel,
                      hasInteraction = FALSE,
                      modelWithoutInteractionForVif = NA,
                      computedfBetas = FALSE,
                      fittedWith = "glmmTMB"
)

vifmodBS <- data.frame(check_collinearity(boxSucmodel))
# vifmodBS$Term <- factor("IT", "poly(SS, 2)", "R", ordered = TRUE)
tab_df(vifmodBS)
ggtexttable((vifmodBS[,1:2]), rows = NULL, theme = ttheme("light", base_size = 10)) %>%
  table_cell_font(row = 2:3, column = 1, face = "italic", size = 10)



### Human food ####
var_weightHF <- data.frame(var = variables,
                           weight = rep(NA, 3),
                           lowerCI = rep(NA, 3),
                           upperCI = rep(NA, 3))

summodHF <- summary(HFmodel)
#do this with every variable
for (i in (1:length(variables))){
  var_weightHF$weight[var_weightHF$var == variables[i]] <- as.numeric(summodHF$coefficients$cond[i+1])
  var_weightHF$lowerCI[var_weightHF$var == variables[i]] <- as.numeric(summodHF$coefficients$cond[i+1] + qt(0.05/2, Nobs - 1)*summodHF$coefficients$cond[i + 5])
  var_weightHF$upperCI[var_weightHF$var == variables[i]] <- as.numeric(summodHF$coefficients$cond[i+1] - qt(0.05/2, Nobs - 1)*summodHF$coefficients$cond[i + 5])
  
}


var_weightHF$category <- factor(c("IT", 
                                  "SS", "SS2",
                                  "R"), ordered = TRUE)


#plot the weights with confidence intervals

ggplot(var_weightHF, aes(x = fct_inorder(category), y = weight, color = category)) +
  geom_point(aes(size = 3)) +  
  geom_errorbar(aes(ymin = lowerCI, ymax = upperCI), width = 0) +
  geom_hline(yintercept = 0) + 
  theme_set(theme_bw(base_family = 'Calibri')) + 
  theme(axis.text = element_text(size=20),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 20),
        legend.position="none") +
  scale_size(guide = 'none') +
  ylab("Variable estimate") +
  scale_color_manual(values = c("IT" = "#6c719c",
                                "SS" = "#b5739d",
                                "SS2" = "#6e3323",
                                "R" = "#81B29A")) +
  coord_flip() +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 20))


library(sjPlot)
plot_model(HFmodel)


checksHF <- check_model(HFmodel,
                        hasInteraction = FALSE,
                        modelWithoutInteractionForVif = NA,
                        computedfBetas = FALSE,
                        fittedWith = "glmmTMB"
)

vifmodHF <- data.frame(check_collinearity(HFmodel))
# vifmodHF$Term <- factor("IT", "SS", "SS2", "R", ordered = TRUE)
ggtexttable((vifmodHF[,1:2]), rows = NULL, theme = ttheme("light", base_size = 10)) %>%
  table_cell_font(row = 2:3, column = 1, face = "italic", size = 10)

#both together
tab_dfs(list(vifmodBS[,1:2], vifmodHF[,1:2]))







#--- End of Script ---#

# sucmod1.1 <- lm(data = sucCompiled_df %>% filter(!is.na(IC)), boxsucrate ~ poly(IC, 2) + HF_dietpercent)
# summary(sucmod1.1)
# drop1(sucmod1.1, test = "Chisq")
# 
# diagnostics.plot.dharma(sucmod1.1)
# 
# sucmod2.1 <- lm(data = sucCompiled_df %>% filter(!is.na(IC)), boxsucrate ~ IC + SS)
# summary(sucmod2.1)
# drop1(sucmod2.1, test = "Chisq")
# #just E is significant
# 
# diagnostics.plot.dharma(sucmod2.1)
# 
# 
# #Phases 1, 2, and 3 weighted success
# sucmod1.1 <- lm(data = sucCompiled_df, boxsucrate ~ IC * SS + HF_dietpercent)
# summary(sucmod1.1)
# drop1(sucmod1.1, test = "Chisq")
# #no sig
# 
# sucmod1.2 <- lm(data = sucCompiled_df, boxsucrate ~ IC * SS + HF_durrate)
# summary(sucmod1.2)
# drop1(sucmod1.2, test = "Chisq")
# #HF is sig
# 
# sucmod1.3 <- lm(data = sucCompiled_df, boxsucrate ~ IC * SS + HF_typeeventrate)
# summary(sucmod1.3)
# drop1(sucmod1.3, test = "Chisq")
# #no sig
# 
# sucmod1.4 <- lm(data = sucCompiled_df, boxsucrate ~ IC * SS + HF_dayeventrate)
# summary(sucmod1.4)
# drop1(sucmod1.4, test = "Chisq")
# #no sig
# 
# sucmod1.5 <- lm(data = sucCompiled_df, boxsucrate ~ IC * SS)
# summary(sucmod1.5)
# drop1(sucmod1.5, test = "Chisq")
# #no sig
# 
# 
# sucmod2.1 <- lm(data = sucCompiled_df, boxsucrate ~ E + SS + HF_dietpercent)
# summary(sucmod2.1)
# drop1(sucmod2.1, test = "Chisq")
# #E sig
# 
# sucmod2.2 <- lm(data = sucCompiled_df, boxsucrate ~ E + SS + HF_durrate)
# summary(sucmod2.2)
# drop1(sucmod2.2, test = "Chisq")
# #E sig
# 
# sucmod2.3 <- lm(data = sucCompiled_df, boxsucrate ~ E + SS + HF_typeeventrate)
# summary(sucmod2.3)
# drop1(sucmod2.3, test = "Chisq")
# #E sig
# 
# sucmod2.4 <- lm(data = sucCompiled_df, boxsucrate ~ E + SS + HF_dayeventrate)
# summary(sucmod2.4)
# drop1(sucmod2.4, test = "Chisq")
# #E sig
# 
# sucmod2.5 <- lm(data = sucCompiled_df, boxsucrate ~ E + SS)
# summary(sucmod2.5)
# drop1(sucmod2.5, test = "Chisq")
# #E sig
# 
# 
# sucmod3.1 <- lm(data = sucCompiled_df, boxsucrate ~ IC + SS + HF_dietpercent)
# summary(sucmod3.1)
# drop1(sucmod3.1, test = "Chisq")
# #IC is sig
# 
# sucmod3.2 <- lm(data = sucCompiled_df, boxsucrate ~ IC + SS + HF_durrate)
# summary(sucmod3.2)
# drop1(sucmod3.2, test = "Chisq")
# #IC andHF is sig
# 
# sucmod3.3 <- lm(data = sucCompiled_df, boxsucrate ~ IC + SS + HF_typeeventrate)
# summary(sucmod3.3)
# drop1(sucmod3.3, test = "Chisq")
# #IC is sig
# 
# sucmod3.4 <- lm(data = sucCompiled_df, boxsucrate ~ IC + SS + HF_dayeventrate)
# summary(sucmod3.4)
# drop1(sucmod3.4, test = "Chisq")
# #IC is sig
# 
# sucmod3.5 <- lm(data = sucCompiled_df, boxsucrate ~ IC + SS)
# summary(sucmod3.5)
# drop1(sucmod3.5, test = "Chisq")
# #IC is sig
# 
# 
# sucmodnull1 <- lm(data = sucCompiled_df[!is.na(sucCompiled_df$IC),], boxsucrate ~ 1)
# summary(sucmodnull1)
# sucmodnull2 <- lm(data = sucCompiled_df[!is.na(sucCompiled_df$E),], boxsucrate ~ 1)
# summary(sucmodnull2)
# 
# 
# anova(sucmod1.3, sucmodnull1)
# # not significantly different
# anova(sucmod2.3, sucmodnull2)
# # not sig
# 
# 
# #Phases 1, and 2 weighted success
# sucmod_p1p2_1.1 <- lm(data = sucCompiled_df, boxsucrate_p1p2 ~ IC * SS + HF_dietpercent)
# summary(sucmod_p1p2_1.1)
# drop1(sucmod_p1p2_1.1, test = "Chisq")
# #no sig
# 
# sucmod_p1p2_1.2 <- lm(data = sucCompiled_df, boxsucrate_p1p2 ~ IC * SS + HF_durrate)
# summary(sucmod_p1p2_1.2)
# drop1(sucmod_p1p2_1.2, test = "Chisq")
# #no sig
# 
# sucmod_p1p2_1.3 <- lm(data = sucCompiled_df, boxsucrate_p1p2 ~ IC * SS + HF_typeeventrate)
# summary(sucmod_p1p2_1.3)
# drop1(sucmod_p1p2_1.3, test = "Chisq")
# #no sig
# 
# sucmod_p1p2_1.4 <- lm(data = sucCompiled_df, boxsucrate_p1p2 ~ IC * SS + HF_dayeventrate)
# summary(sucmod_p1p2_1.4)
# drop1(sucmod_p1p2_1.4, test = "Chisq")
# #no sig
# 
# sucmod_p1p2_1.5 <- lm(data = sucCompiled_df, boxsucrate_p1p2 ~ IC * SS)
# summary(sucmod_p1p2_1.5)
# drop1(sucmod_p1p2_1.5, test = "Chisq")
# #no sig
# 
# 
# sucmod_p1p2_2.1 <- lm(data = sucCompiled_df, boxsucrate_p1p2 ~ E + SS + HF_dietpercent)
# summary(sucmod_p1p2_2.1)
# drop1(sucmod_p1p2_2.1, test = "Chisq")
# #no sig
# 
# sucmod_p1p2_2.2 <- lm(data = sucCompiled_df, boxsucrate_p1p2 ~ E + SS + HF_durrate)
# summary(sucmod_p1p2_2.2)
# drop1(sucmod_p1p2_2.2, test = "Chisq")
# #no sig
# 
# sucmod_p1p2_2.3 <- lm(data = sucCompiled_df, boxsucrate_p1p2 ~ E + SS + HF_typeeventrate)
# summary(sucmod_p1p2_2.3)
# drop1(sucmod_p1p2_2.3, test = "Chisq")
# #no sig
# 
# sucmod_p1p2_2.4 <- lm(data = sucCompiled_df, boxsucrate_p1p2 ~ E + SS + HF_dayeventrate)
# summary(sucmod_p1p2_2.4)
# drop1(sucmod_p1p2_2.4, test = "Chisq")
# #no sig
# 
# sucmod_p1p2_2.5 <- lm(data = sucCompiled_df, boxsucrate_p1p2 ~ E + SS )
# summary(sucmod_p1p2_2.5)
# drop1(sucmod_p1p2_2.5, test = "Chisq")
# #no sig
# 
# 
# sucmod_p1p2_3.1 <- lm(data = sucCompiled_df, boxsucrate_p1p2 ~ IC + SS + HF_dietpercent)
# summary(sucmod_p1p2_3.1)
# drop1(sucmod_p1p2_3.1, test = "Chisq")
# #no sig
# 
# sucmod_p1p2_3.2 <- lm(data = sucCompiled_df, boxsucrate_p1p2 ~ IC + SS + HF_durrate)
# summary(sucmod_p1p2_3.2)
# drop1(sucmod_p1p2_3.2, test = "Chisq")
# #no sig
# 
# sucmod_p1p2_3.3 <- lm(data = sucCompiled_df, boxsucrate_p1p2 ~ IC + SS + HF_typeeventrate)
# summary(sucmod_p1p2_3.3)
# drop1(sucmod_p1p2_3.3, test = "Chisq")
# #no sig
# 
# sucmod_p1p2_3.4 <- lm(data = sucCompiled_df, boxsucrate_p1p2 ~ IC + SS + HF_dayeventrate)
# summary(sucmod_p1p2_3.4)
# drop1(sucmod_p1p2_3.4, test = "Chisq")
# #no sig
# 
# sucmod_p1p2_3.5 <- lm(data = sucCompiled_df, boxsucrate_p1p2 ~ IC + SS)
# summary(sucmod_p1p2_3.5)
# drop1(sucmod_p1p2_3.5, test = "Chisq")
# #no sig
# 
# 
# 
# sucmodnull_p1p2_1 <- lm(data = sucCompiled_df[!is.na(sucCompiled_df$IC),], boxsucrate_p1p2 ~ 1)
# summary(sucmodnull_p1p2_1)
# sucmodnull_p1p2_2 <- lm(data = sucCompiled_df[!is.na(sucCompiled_df$E),], boxsucrate_p1p2 ~ 1)
# summary(sucmodnull_p1p2_2)
# 
# 
# anova(sucmod_p1p2_1.3, sucmodnull_p1p2_1)
# # not significantly different
# anova(sucmod_p1p2_2.3, sucmodnull_p1p2_2)
# # not sig



### Check model -------------------------------------------------------------





### Plot model --------------------------------------------------------------




##--##--orientation = 
## END OF SCRIPT

