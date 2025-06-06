#What does this script do?
#Analysis of PB masters thesis - innovation complex


#To do:
#Latency up to the first IC to the order that they were tested?
#or IC score vs order of tested.

#Try to separate house and bin data -> also we have more data now (June 20, 2024)
#Diversity of food types
#Repeat IC analysis with SS

#if include rank, scale for rank in different groups (then divide by the standard dev of both the groups together)


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
library(plotly)
library(stringr)
library(sjPlot)

library(FactoMineR)
library(vcd)
library(factoextra)
library(librarian)
library(car)
library(patchwork)
library(gridGraphics)
library(grid)
library(ggpubr)

#Data processing
library(dplyr)
library("readxl")
library(tidyr)


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

setwd("C:/Users/admin/Desktop/MSThesis/2-Scripts")
source("./0_parameters.R")

# Analysis ----------------------------------------------------------------

## Data processing ---------------------------------------------------------
compiled_df <- summaryFit_df %>%
  mutate(IC_cat = ifelse(what == "IC" & estUpr > IC_cutoff, 1, NA)) %>%
  group_by(id) %>%
  mutate(IC_cat = ifelse(sum(IC_cat, na.rm = TRUE) >=1, 1,0)) %>%
  dplyr::select(id, what, est, IC_cat) %>%
  spread(key = what, value = est) %>%
  merge(raiding_df, by.x = "id", by.y = "ID", all.x = TRUE) %>%
  merge(exposure_df, by.x = "id", all.x = TRUE) %>%
  merge(rank_df, by = "id", all.x = TRUE)

compiled_df$Age <- factor(compiled_df$Age, 
                          levels = c("1y", "2y", "3y", "4y", "a"),
                          ordered = TRUE)
compiled_df$Sex <- factor(compiled_df$Sex)
compiled_df$Group <- factor(compiled_df$Group)

compiled_df$rank <- as.numeric(compiled_df$rank)


## Fitting correlation model ----------------------------------------------

model <- glmmTMB(ICc ~ E, data = compiled_df %>% 
                   filter(!is.na(IC)) %>% 
                   mutate(Ec = (E*(nrow(.) - 1) + 0.5)/nrow(.),
                          ICc = (IC*(nrow(.) - 1) + 0.5)/nrow(.)), family = beta_family(link = "logit"))
summary(model)
model <- glmmTMB(Ec ~ IC, data = compiled_df %>%
                   filter(!is.na(IC)) %>%
                   mutate(ICc = (IC*(nrow(.) - 1) + 0.5)/nrow(.),
                          Ec = (E*(nrow(.) - 1) + 0.5)/nrow(.)), family = beta_family(link = "logit"))
summary(model)
confint(model, 'IC', level=0.95)

model <- glmmTMB(ICc ~ SS, data = compiled_df %>% 
                   filter(!is.na(IC)) %>% 
                   mutate(ICc = (IC*(nrow(.) - 1) + 0.5)/nrow(.)), family = beta_family(link = "logit"))
summary(model)
confint(model, 'SS', level=0.95)

model <- glmmTMB(SSc ~ IC, data = compiled_df %>% 
                   filter(!is.na(IC)) %>% 
                   mutate(SSc = (SS*(nrow(.) - 1) + 0.5)/nrow(.)), family = beta_family(link = "logit"))
summary(model)
confint(model, 'IC', level=0.95)

model <- glmmTMB(Ec ~ SS, data = compiled_df %>% 
                   filter(!is.na(IC)) %>% 
                   mutate(Ec = (E*(nrow(.) - 1) + 0.5)/nrow(.)), family = beta_family(link = "logit"))
summary(model)
confint(model, 'SS', level=0.95)

model <- glmmTMB(SSc ~ E, data = compiled_df %>% 
                   filter(!is.na(IC)) %>% 
                   mutate(SSc = (SS*(nrow(.) - 1) + 0.5)/nrow(.)), family = beta_family(link = "logit"))
summary(model)
confint(model, 'E', level=0.95)





library(performance)
check_outliers(model)
summary(model)


### Fit model ---------------------------------------------------------------

#do this before scaling
compiled_df$first_is_exposure_sqrt <- sqrt(compiled_df$first_is_exposure)
compiled_df$first_is_expos2week_sqrt <- sqrt(compiled_df$first_is_expos2week)
compiled_df$raidrate_sqrt <- sqrt(compiled_df$raidrate)
compiled_df$ph2inn_sqrt <- sqrt(compiled_df$ph2inn)
compiled_df$ph2exossucrate_sqrt <- sqrt(compiled_df$ph2exossucrate)
compiled_df$ph2attend5m_sqrt <- sqrt(compiled_df$ph2attend5m)
compiledTInit_df <- compiled_df 
compiled_df <- compiledTInit_df 
compiled_df <- compiled_df %>%
  mutate(first_is_expos2week_sqrt = datawizard::standardise(first_is_expos2week_sqrt),
         raidrate_sqrt = datawizard::standardise(raidrate_sqrt),
         Age = as.numeric(as.factor(ifelse(Age != "a", "j", "a"))),
         Sex = as.numeric(as.factor(Sex)),
         first_is_exposure_sqrt = datawizard::standardise(first_is_exposure_sqrt),
         ph2inn_sqrt = datawizard::standardise(ph2inn_sqrt),
         ph2exossucrate_sqrt = datawizard::standardise(ph2exossucrate_sqrt),
         ph2attend5m_sqrt = datawizard::standardise(ph2attend5m_sqrt)
  )

compiled_df$rank[compiled_df$Group == "Acacia"] <- datawizard::standardise(compiled_df$rank[compiled_df$Group == "Acacia"])
compiled_df$rank[compiled_df$Group == "Savanna"] <- datawizard::standardise(compiled_df$rank[compiled_df$Group == "Savanna"])
compiled_df$rank <- compiled_df$rank / sd(compiled_df$rank)


#Try a PCA for the social exposure component
#I used a FAMD, which is a generalization of PCS that can include categorical variables

#check
library(corrplot)
library("Hmisc")


#Social metric
# pca_df <- (na.omit(compiled_df %>% select(expos_sidesolution.l, expos_sidesolution.p, ph2trialcount_sqrt, ph2inn)))
soc_pca_df <- na.omit(compiled_df %>% filter(!is.na(IC)) %>% select(first_is_expos2week_sqrt, ph2inn_sqrt, ph2attend5m_sqrt)) #, ph2exossucrate_sqrt)), ph2trialcount_sqrt))
socrescorr <- rcorr(as.matrix(soc_pca_df))
# corrplot(socrescorr$r, method = "number")
colnames(socrescorr$r) <- c("Attended same \nsolution trial \ncount Phase 1", "Attended innovation \nsuccess count \nPhase 2", "Attended within \n5m trial count \nPhase 2")
corrplot(socrescorr$r, order = 'AOE', addCoef.col = 'black', tl.pos = 'd', tl.col = 'black',
         cl.pos = 'n', col = COL1('Blues'))
names(soc_pca_df) <- c("Attended same \nsolution trial \ncount Phase 1", "Attended innovation \nsuccess count \nPhase 2", "Attended within \n5m trial count \nPhase 2")
soc_pca <- prcomp(soc_pca_df, scale = TRUE)
summary(soc_pca)
# fviz_pca_var(soc_pca, col.ind = "cos2",
#               gradient.cols = c("blue", "orange", "red"),
#               repel = TRUE)
fviz_pca_var(soc_pca, col.var = "contrib", repel=T) +
  scale_color_gradient2(low = "white", high = "darkblue") +
  theme_bw() +
  ggtitle("Social dimension variables - PCA")


#Demographic metric
# pca_df <- (na.omit(compiled_df %>% select(expos_sidesolution.l, expos_sidesolution.p, ph2trialcount_sqrt, ph2inn)))
innate_pca_df <- na.omit(compiled_df %>% filter(!is.na(IC)) %>% select(rank, Sex, Age))
names(innate_pca_df) <- c("Rank", "Sex", "Age")
innaterescorr <- rcorr(as.matrix(innate_pca_df))
# corrplot(innaterescorr$r)
corrplot(innaterescorr$r, order = 'AOE', addCoef.col = 'black', tl.pos = 'd', tl.col = 'black',
         cl.pos = 'n', col = COL1('Blues'))
innate_pca <- prcomp(innate_pca_df, scale = TRUE)
summary(innate_pca)
# fviz_pca_var(innate_pca, col.ind = "cos2",
#              gradient.cols = c("blue", "orange", "red"),
#              repel = TRUE)
fviz_pca_var(innate_pca, col.var = "contrib", repel=T) +
  scale_color_gradient2(low = "white", high = "darkblue") +
  theme_bw() +
  ggtitle("Demographic dimension variables - PCA")

#enviro metric
# pca_df <- (na.omit(compiled_df %>% select(expos_sidesolution.l, expos_sidesolution.p, ph2trialcount_sqrt, ph2inn)))
# enviro_pca_df <- na.omit(compiled_df %>% filter(!is.na(IC)) %>% select(garbageraidrate_sqrt, houseraidrate_sqrt))
# envirorescorr <- rcorr(as.matrix(enviro_pca_df))
# corrplot(envirorescorr$r, method = "number")
# #highly correlated
# enviro_pca <- prcomp(enviro_pca_df, scale = TRUE)
# summary(enviro_pca)
# fviz_pca_var(enviro_pca, col.ind = "cos2",
#              gradient.cols = c("blue", "orange", "red"),
#              repel = TRUE)



compiled_df$socpca[!is.na(compiled_df$IC)] <- datawizard::standardise(soc_pca$x[1:(length(soc_pca$x)/3)])
compiled_df$innatepca[!is.na(compiled_df$IC)] <- datawizard::standardise(innate_pca$x[1:(length(innate_pca$x)/3)])

table(compiled_df$Age, compiled_df$Sex)
#Vector of formulas
vectorFormulas_v <- c(
  as.formula(IC_cat ~ 1),
  as.formula(IC_cat ~ innatepca + raidrate_sqrt),
  as.formula(IC_cat ~ innatepca + socpca),
  as.formula(IC_cat ~ socpca + raidrate_sqrt),
  as.formula(IC_cat ~ innatepca),
  as.formula(IC_cat ~ socpca),
  as.formula(IC_cat ~ raidrate_sqrt),
  as.formula(IC_cat ~ innatepca + socpca + raidrate_sqrt)
)

#Run models
modelsToCompare_l <- lapply(1:length(vectorFormulas_v), function(whichModel){
  #Fit the seasonal model centered on the week of interest
  model <- glmmTMB(
    vectorFormulas_v[whichModel][[1]],
    data = compiled_df %>% mutate(IC = (IC*(nrow(.) - 1) + 0.5)/nrow(.)),
    family =  "binomial" #beta_family(link = "logit")
  )
  diagnostics.plot.dharma(model)
  check_collinearity(model)
  
  #Return the model
  return(model)
})
names(modelsToCompare_l) <- c("null",
                             # "Demographic, Social, Experience",
                              "Demographic, Experience",
                              "Demographic, Social",
                              "Social, Experience",
                              "Demographic",
                              "Social",
                              "Experience",
                              "Demographic, Social, Experience"
                             )

#Extract model ranking
selectionTable <- aictab(modelsToCompare_l)
selectionTable


selectionTablePrint <- selectionTable %>%
         mutate('Model names' = Modnames,
                AICc = round(AICc, 3),
                'ΔAICc' = round(Delta_AICc, 3),
                Weight = round(AICcWt, 3),
                'Log-likelihood' = round(LL, 3),
                'Cumulative weight' = round(Cum.Wt, 3)) %>%
         select('Model names', K, AICc, 'ΔAICc', 'Log-likelihood', Weight, 'Cumulative weight')

tab_df(selectionTablePrint)


variables <- c("socpca", "innatepca", "raidrate_sqrt")
               # "first_is_exposure_sqrt", "ph2inn1",
               # #"ph2trialcount_sqrt", "ph2inn",
               # "Sexm", "Agej", "raidrate_sqrt")
var_weight <- data.frame(var = variables,
                         weight = rep(NA, 3),
                         lowerCI = rep(NA, 3),
                         upperCI = rep(NA, 3))

#do this with every variable
for (i in (1:length(variables))){
  var_weight$weight[var_weight$var == variables[i]] <- as.numeric(modavg(modelsToCompare_l,
                                                                         parm = variables[i])[3][[1]])
  var_weight$lowerCI[var_weight$var == variables[i]] <- as.numeric(modavg(modelsToCompare_l,
                                                                          parm = variables[i])[6][[1]])
  var_weight$upperCI[var_weight$var == variables[i]] <- as.numeric(modavg(modelsToCompare_l,
                                                                          parm = variables[i])[7][[1]])
  
}


for (i in (1:length(variables))){
  var_weight$weight[var_weight$var == variables[i]] <- as.numeric(modavgShrink(modelsToCompare_l,
                                                                         parm = variables[i])[3][[1]])
  var_weight$lowerCI[var_weight$var == variables[i]] <- as.numeric(modavgShrink(modelsToCompare_l,
                                                                          parm = variables[i])[6][[1]])
  var_weight$upperCI[var_weight$var == variables[i]] <- as.numeric(modavgShrink(modelsToCompare_l,
                                                                          parm = variables[i])[7][[1]])
  
}



var_weight$category <- c("Social exposure to experiment", #"Social", 
                         "Demographic traits", #"Demographic",
                         "Anthropogenic raiding tendency")



### Check model -------------------------------------------------------------

#IC
fullmodelIC <- model <- glmmTMB(
  vectorFormulas_v[5][[1]],
  data = compiled_df[!is.na(compiled_df$innatepca),],
  family = "binomial"
)
summary(fullmodelIC)

plot(allEffects(fullmodelIC))

# ICinnate_df <- compiled_df[!is.na(compiled_df$innatepca),] %>%
#   select(IC_cat, innatepca, Sex, Age, rank)
ICinnate_df <- compiled_df[!is.na(compiled_df$IC),] %>%
  select(IC_cat, Sex, Age, rank)
ICinnate_df$IC_cat <- as.factor(ICinnate_df$IC_cat)
ICinnate_df$Sex <- as.factor(ICinnate_df$Sex)
ICinnate_df$age <- "a"
ICinnate_df$age[ICinnate_df$Age != "a"] <- "j"
ICinnate_df$Age <- as.factor(ICinnate_df$Age)
# ggplot(ICinnate_df, aes(x = IC_cat, y = innatepca, group = IC_cat)) +
#   geom_boxplot()

ICinnSex <- ICinnate_df %>%
  group_by(Sex) %>%
  summarise(IC_catmean = mean(as.numeric(IC_cat)))

ICinnAge <- ICinnate_df %>%
  group_by(Age) %>%
  summarise(IC_catmean = mean(as.numeric(IC_cat)))

ggplot(ICinnate_df, aes(x = Sex, y = IC_cat, group = Sex)) +
  geom_violin()
ggplot(ICinnate_df, aes(x = age, y = IC_cat, group = age)) +
  geom_violin()
ggplot(ICinnate_df, aes(x = rank, y = IC_cat)) +
  geom_point() +
  theme_set(theme_bw(base_family = 'Calibri')) + 
  theme(axis.text = element_text(size=20),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 20),
        legend.position="none") +
  stat_summary(fun.data=mean_cl_normal) + 
  geom_smooth(method=lm, na.rm = TRUE, fullrange= TRUE,
              aes(group=1),colour="black")+
  scale_size(guide = 'none') +
  ylab("Variable estimate") 



### Plot model --------------------------------------------------------------



#plot the weights with confidence intervals

MS_forest <- ggplot(var_weight, aes(x = fct_inorder(category), y = weight, color = category)) +
  geom_point(aes(size = 3)) +  
  geom_errorbar(aes(ymin = lowerCI, ymax = upperCI), width = 0) +
  geom_hline(yintercept = 0) + 
  theme_set(theme_bw(base_family = 'Calibri')) + 
  theme(axis.text = element_text(size=15),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 15),
        legend.position="none") +
  scale_size(guide = 'none') +
  ylab("Variable estimate") +
  scale_color_manual(values = c("Social exposure to experiment" = "#B48168",
                                "Anthropogenic raiding tendency" = "#81B29A",
                                "Demographic traits" = "#F2CC8F")) +
  coord_flip() +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 15))

  

#x axis: variables
#y axis 1: probability of being selected (sum of the normalized 
#weights of the model including the variable)
#y axis 2: absolute average estimate (mean weighted estimate of
#all the estimates of each model in the 95% confidence 
#set of best models)

group_var <- c("Social", "Experience", "Demographic")
plot_df <- data.frame(group_var = group_var,
                      p_select_ic = rep(NA, length(group_var)),
                      avg_est_ic = rep(NA, length(group_var)))

for(i in (1:length(group_var))){
  plot_df$p_select_ic[i] <- sum(as.data.frame(selectionTable)$AICcWt[grep(plot_df$group_var[i], 
                                                                         as.data.frame(selectionTable)$Modnames)])
  plot_df$avg_est_ic[i] <- mean(abs(as.data.frame(selectionTable)$Delta_AICc[grep(plot_df$group_var[i], 
                                                                                  as.data.frame(selectionTable)$Modnames)]))
}

plot_df$group_var <- c("Social exposure to experiment", "Anthropogenic raiding tendency", "Demographic traits")

MS_bar <- ggplot() +
  geom_bar(data=plot_df, aes(x=group_var, y=p_select_ic, fill = group_var), stat="identity", color = "black") + 
  theme_set(theme_bw(base_family = 'Calibri')) + 
  scale_fill_manual(values = c("Social exposure to experiment" = "#B48168",
                               "Anthropogenic raiding tendency" = "#96916A",
                               "Demographic traits" = "#F2CC8F")) +
  ylab("Model selection weight") +
  theme(axis.text = element_text(size=15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        legend.position="none") +
  ylim(c(0,0.75)) +
  scale_size(guide = 'none') +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10))

MS_bar + MS_forest



## Repeat with SS ####

vectorFormulas_SS <- c(
  as.formula(SS ~ 1),
  as.formula(SS ~ innatepca + raidrate_sqrt),
  as.formula(SS ~ innatepca + socpca),
  as.formula(SS ~ socpca + raidrate_sqrt),
  as.formula(SS ~ innatepca),
  as.formula(SS ~ socpca),
  as.formula(SS ~ raidrate_sqrt),
  as.formula(SS ~ innatepca + socpca + raidrate_sqrt)
)

#Run models
modelsToCompare_SS <- lapply(1:length(vectorFormulas_SS), function(whichModel){
  #Fit the seasonal model centered on the week of interest
  model <- glmmTMB(
    vectorFormulas_SS[whichModel][[1]],
    data = compiled_df %>% mutate(SS = (SS*(nrow(.) - 1) + 0.5)/nrow(.)),
    family =  beta_family(link = "logit")
  )
  diagnostics.plot.dharma(model)
  check_collinearity(model)
  
  #Return the model
  return(model)
})
names(modelsToCompare_SS) <- c("null",
                              "Demographic, Experience",
                              "Demographic, Social",
                              "Social, Experience",
                              "Demographic",
                              "Social",
                              "Experience",
                              "Demographic, Social, Experience"
)

#Extract model ranking
selectionTableSS <- aictab(modelsToCompare_SS)
selectionTableSS


selectionTablePrintSS <- selectionTableSS %>%
  mutate('Model names' = Modnames,
         AICc = round(AICc, 3),
         'ΔAICc' = round(Delta_AICc, 3),
         Weight = round(AICcWt, 3),
         'Log-likelihood' = round(LL, 3),
         'Cumulative weight' = round(Cum.Wt, 3)) %>%
  select('Model names', K, AICc, 'ΔAICc', 'Log-likelihood', Weight, 'Cumulative weight')

tab_df(selectionTablePrintSS)


variables <- c("socpca", "innatepca", "raidrate_sqrt")
# "first_is_exposure_sqrt", "ph2inn1",
# #"ph2trialcount_sqrt", "ph2inn",
# "Sexm", "Agej", "raidrate_sqrt")
var_weightSS <- data.frame(var = variables,
                         weight = rep(NA, 3),
                         lowerCI = rep(NA, 3),
                         upperCI = rep(NA, 3))

#do this with every variable
for (i in (1:length(variables))){
  var_weightSS$weight[var_weightSS$var == variables[i]] <- as.numeric(modavg(modelsToCompare_SS,
                                                                         parm = variables[i])[3][[1]])
  var_weightSS$lowerCI[var_weightSS$var == variables[i]] <- as.numeric(modavg(modelsToCompare_SS,
                                                                          parm = variables[i])[6][[1]])
  var_weightSS$upperCI[var_weightSS$var == variables[i]] <- as.numeric(modavg(modelsToCompare_SS,
                                                                          parm = variables[i])[7][[1]])
  
}


for (i in (1:length(variables))){
  var_weightSS$weight[var_weightSS$var == variables[i]] <- as.numeric(modavgShrink(modelsToCompare_SS,
                                                                               parm = variables[i])[3][[1]])
  var_weightSS$lowerCI[var_weightSS$var == variables[i]] <- as.numeric(modavgShrink(modelsToCompare_SS,
                                                                                parm = variables[i])[6][[1]])
  var_weight$upperCI[var_weight$var == variables[i]] <- as.numeric(modavgShrink(modelsToCompare_l,
                                                                                parm = variables[i])[7][[1]])
  
}



var_weightSS$category <- c("Social exposure to experiment", #"Social", 
                         "Demographic traits", #"Demographic",
                         "Anthropogenic raiding tendency")


### Check model -------------------------------------------------------------

#IC
fullmodelSS <- model <- glmmTMB(
  vectorFormulas_SS[8][[1]],
  data = compiled_df[!is.na(compiled_df$innatepca),],
  family = beta_family(link = "logit")
)
summary(fullmodelSS)

plot(allEffects(fullmodelSS))

### Plot model --------------------------------------------------------------

#plot the weights with confidence intervals

MS_forestSS <- ggplot(var_weightSS, aes(x = fct_inorder(category), y = weight, color = category)) +
  geom_point(aes(size = 3)) +  
  geom_errorbar(aes(ymin = lowerCI, ymax = upperCI), width = 0) +
  geom_hline(yintercept = 0) + 
  theme_set(theme_bw(base_family = 'Calibri')) + 
  theme(axis.text = element_text(size=15),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 15),
        legend.position="none") +
  scale_size(guide = 'none') +
  ylab("Variable estimate") +
  scale_color_manual(values = c("Social exposure to experiment" = "#B48168",
                                "Anthropogenic raiding tendency" = "#81B29A",
                                "Demographic traits" = "#F2CC8F")) +
  coord_flip() +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 15))



#x axis: variables
#y axis 1: probability of being selected (sum of the normalized 
#weights of the model including the variable)
#y axis 2: absolute average estimate (mean weighted estimate of
#all the estimates of each model in the 95% confidence 
#set of best models)

group_var <- c("Social", "Experience", "Demographic")
plot_dfSS <- data.frame(group_var = group_var,
                      p_select_ic = rep(NA, length(group_var)),
                      avg_est_ic = rep(NA, length(group_var)))

for(i in (1:length(group_var))){
  plot_dfSS$p_select_ic[i] <- sum(as.data.frame(selectionTableSS)$AICcWt[grep(plot_dfSS$group_var[i], 
                                                                          as.data.frame(selectionTableSS)$Modnames)])
  plot_dfSS$avg_est_ic[i] <- mean(abs(as.data.frame(selectionTableSS)$Delta_AICc[grep(plot_dfSS$group_var[i], 
                                                                                  as.data.frame(selectionTableSS)$Modnames)]))
}

plot_dfSS$group_var <- c("Social exposure to experiment", "Anthropogenic raiding tendency", "Demographic traits")

MS_barSS <- ggplot() +
  geom_bar(data=plot_dfSS, aes(x=group_var, y=p_select_ic, fill = group_var), stat="identity", color = "black") + 
  theme_set(theme_bw(base_family = 'Calibri')) + 
  scale_fill_manual(values = c("Social exposure to experiment" = "#B48168",
                               "Anthropogenic raiding tendency" = "#96916A",
                               "Demographic traits" = "#F2CC8F")) +
  ylab("Model selection weight") +
  theme(axis.text = element_text(size=15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        legend.position="none") +
  ylim(c(0,0.75)) +
  scale_size(guide = 'none') +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10))

MS_barSS + MS_forestSS

MS_bar + MS_forest + MS_barSS + MS_forestSS

plot_dfSS$para <- "Simple switch\n tendency"
plot_df$para <- "Technical innovation"

plot_df_combo <- rbind(plot_df, plot_dfSS)
ggplot(data=plot_df_combo, aes(col = para, x=group_var, y=p_select_ic, group = para, fill = group_var)) +
  geom_bar(stat="identity",  position = "dodge", size = 2) + 
  theme_bw() + 
  scale_fill_manual(values = c("Social exposure to experiment" = "#B48168",
                               "Anthropogenic raiding tendency" =  "#81B29A", #"#96916A",
                               "Demographic traits" = "#F2CC8F")) +
  scale_color_manual(values = c("#B5739D", "#6C719C")) +
  ylab("Model selection weight") +
  theme(axis.text = element_text(size=15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        legend.position = c(.25, .8),
        legend.text = element_text(size = 10)) +
  guides(fill = "none", col = guide_legend(title = "Parameter")) +
  ylim(c(0,0.75)) +
  scale_size(guide = 'none') +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10))





### model the mechanistic model parameter summaries per individual ########
summaryFit_df <- summaryFit_df %>%
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
  arrange(desc(Age), desc(id))


#SS
ggplot(summaryFit_df %>% filter(!is.na(est) & what == "SS"), aes(x = fct_inorder(id), y = est, color = Group)) +
  geom_point(aes(color = Group), size = 2) +  
  geom_errorbar(aes(ymin = estLwr, ymax = estUpr, color = Group), width = 0) +
  theme_set(theme_bw(base_family = 'Calibri')) + 
  theme(axis.text = element_text(size=12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        title = element_text(size = 12),
        legend.position="none") +
  ggtitle("Switch Tendency Simple (SS)") +
  scale_size(guide = 'none') +
  ylab("Parameter estimate") +
  xlab("ID") +
  coord_flip() +
  ylim(0,1) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 15)) +
  scale_color_manual(values = c("Savanna" = "darkorange2",
                                "Acacia" = "darkolivegreen4"))


#IS
ggplot(summaryFit_df %>% filter(!is.na(est) & what == "IS"), aes(x = fct_inorder(id), y = est), color = Group) +
  geom_point(aes(color = Group), size = 2) +  
  geom_errorbar(aes(ymin = estLwr, ymax = estUpr, color = Group), width = 0) +
  theme_set(theme_bw(base_family = 'Calibri')) + 
  theme(axis.text = element_text(size=12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        title = element_text(size = 12),
        legend.position="none") +
  ggtitle("Innovation Simple (IS)") +
  scale_size(guide = 'none') +
  ylab("Parameter estimate") +
  xlab("ID") +
  coord_flip() +
  ylim(0,1) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 15)) +
  scale_color_manual(values = c("Savanna" = "darkorange2",
                                "Acacia" = "darkolivegreen4"))

#SC
ggplot(summaryFit_df %>% filter(!is.na(est) & what == "SC"), aes(x = fct_inorder(id), y = est, color = Group)) +
  geom_point(aes(color = Group), size = 2) +  
  geom_errorbar(aes(ymin = estLwr, ymax = estUpr, color = Group), width = 0) +
  theme_set(theme_bw(base_family = 'Calibri')) + 
  theme(axis.text = element_text(size=12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        title = element_text(size = 12),
        legend.position="none") +
  ggtitle("Technical Switch Tendency (ST)") +
  scale_size(guide = 'none') +
  ylab("Parameter estimate") +
  xlab("ID") +
  coord_flip() +
  ylim(0,1) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 15)) +
  scale_color_manual(values = c("Savanna" = "darkorange2",
                                "Acacia" = "darkolivegreen4"))


#IC
ggplot(summaryFit_df %>% filter(!is.na(est) & what == "IC"), aes(x = fct_inorder(id), y = est, color = Group)) +
  geom_point(aes(color = Group), size = 2) +  
  geom_errorbar(aes(ymin = estLwr, ymax = estUpr, color = Group), width = 0) +
  theme_set(theme_bw(base_family = 'Calibri')) + 
  theme(axis.text = element_text(size=12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        title = element_text(size = 12),
        legend.position="none") +
  ggtitle("Technical Innovation (IT)") +
  scale_size(guide = 'none') +
  ylab("Parameter estimate") +
  xlab("ID") +
  coord_flip() +
  ylim(0,1) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 15)) +
  scale_color_manual(values = c("Savanna" = "darkorange2",
                                "Acacia" = "darkolivegreen4"))


#E
ggplot(summaryFit_df %>% filter(!is.na(est) & what == "E"), aes(x = fct_inorder(id), y = est, color = Group)) +
  geom_point(aes(color = Group), size = 2) +  
  geom_errorbar(aes(ymin = estLwr, ymax = estUpr, color = Group), width = 0) +
  theme_set(theme_bw(base_family = 'Calibri')) + 
  theme(axis.text = element_text(size=12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        title = element_text(size = 12),
        legend.position="none") +
  ggtitle("Learning Sensitivity (V)") +
  scale_size(guide = 'none') +
  ylab("Parameter estimate") +
  xlab("ID") +
  coord_flip() +
  ylim(0,1) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 15)) +
  scale_color_manual(values = c("Savanna" = "darkorange2",
                                "Acacia" = "darkolivegreen4"))


### Check most complex model convergence ##########
  

#### IT model checking ####
#pull out most complex model 
model8 <- glmmTMB(IC_cat ~ innatepca + socpca + raidrate_sqrt,
                  data = compiled_df %>% mutate(IC = (IC*(nrow(.) - 1) + 0.5)/nrow(.)),
                  family =  "binomial" #beta_family(link = "logit")
)

check_convergence(model8) #TRUE, this is fine
AICcmodavg::checkConv(model8)
sapply(modelsToCompare_l, AICcmodavg::checkConv) #all converged
check_singularity(model8) #FALSE, good
check_outliers(model8) #not able

modelForCheck <- glm(formula = formula(model8), data = model8$frame, family = "binomial")
check_outliers(modelForCheck) #fine
extractCN(modelForCheck) #not able -> model should be refitted to match the unmarkFit class; check the unmarkfit package and the vignette were extractCN is used
CN(model8) #look into this?
extractCN(model8) #not able
vifmod8 <- data.frame(check_collinearity(model8))
vifmod8$Term <- c("Demographic traits", "Social exposure to experiment", "Anthropogenic raiding tendency")
vifmod8 <- vifmod8 %>%
  arrange(Term)
ggtexttable((vifmod8[,1:2]), rows = NULL, theme = ttheme("light", base_size = 10))


checks <- check_model(model8,
  hasInteraction = FALSE,
  modelWithoutInteractionForVif = NA,
  computedfBetas = FALSE,
  fittedWith = "glmmTMB"
)


#### SS model checking ####
#pull out most complex model 
modelSS <- glmmTMB(SS ~ innatepca + socpca + raidrate_sqrt,
                  data = compiled_df %>% mutate(IC = (IC*(nrow(.) - 1) + 0.5)/nrow(.)),
                  family =  beta_family(link = "logit")
)

check_convergence(modelSS) #TRUE, this is fine
AICcmodavg::checkConv(modelSS)
sapply(modelsToCompare_SS, AICcmodavg::checkConv) #all converged
check_singularity(modelSS) #FALSE, good
check_outliers(modelSS) #not able

modelForCheckSS <- glm(formula = formula(modelSS), data = modelSS$frame, family = beta_family(link = "logit"))
check_outliers(modelForCheckSS) #fine
extractCN(modelForCheckSS) #not able -> model should be refitted to match the unmarkFit class; check the unmarkfit package and the vignette were extractCN is used
CN(modelSS) #look into this?
extractCN(modelSS) #not able
vifmodSS <- data.frame(check_collinearity(modelSS))
vifmodSS$Term <- c("Demographic traits", "Social exposure to experiment", "Anthropogenic raiding tendency")
vifmodSS <- vifmodSS %>%
  arrange(Term)
ggtexttable((vifmodSS[,1:2]), rows = NULL, theme = ttheme("light", base_size = 10))


checks <- check_model(model8,
                      hasInteraction = FALSE,
                      modelWithoutInteractionForVif = NA,
                      computedfBetas = FALSE,
                      fittedWith = "glmmTMB"
)

##--##--orientation = 
## END OF SCRIPT


library(readr)
dataExperimentBox_df <- read_delim("C:/Users/admin/Desktop/MSThesis/AttemptSides_21-05-24.csv", 
                                   escape_double = FALSE, trim_ws = TRUE)
colnames(dataExperimentBox_df) <- c("id", "phase", "opening", "trial", "innovationL", "innovationP", "attemptL", "attemptP")


toMatch_df <- dataExperimentBox_df %>% 
  filter(opening != "a") %>% 
  group_by(id) %>% 
  summarise(
    phaseMax = max(phase, na.rm = TRUE)
  )

toPlot_df<- compiled_df %>% 
  left_join(toMatch_df, by = c("id"))

library(ggplot2)
ggplot(toPlot_df, mapping = aes(x = phaseMax, y = IC)) +
  geom_point(pch = 19)

ggplot(toPlot_df, mapping = aes(x = phaseMax, y = SC)) +
  geom_point(pch = 19)

plot(compiled_df$IC, compiled_df$SC)


#Social metric
# pca_df <- (na.omit(compiled_df %>% select(expos_sidesolution.l, expos_sidesolution.p, ph2trialcount_sqrt, ph2inn)))
soc_pca_df <- na.omit(compiled_df %>% filter(!is.na(SS)) %>% select(first_is_expos2week_sqrt, ph2trialcount_sqrt, ph2inn))
socrescorr <- rcorr(as.matrix(soc_pca_df))
corrplot(socrescorr$r)
soc_pca <- prcomp(soc_pca_df, scale = TRUE)
summary(soc_pca)
fviz_pca_var(soc_pca, col.ind = "cos2",
             gradient.cols = c("blue", "orange", "red"),
             repel = TRUE)

#Demographic metric
# pca_df <- (na.omit(compiled_df %>% select(expos_sidesolution.l, expos_sidesolution.p, ph2trialcount_sqrt, ph2inn)))
innate_pca_df <- na.omit(compiled_df %>% filter(!is.na(SS)) %>% select(rank, Sex, Age))
innaterescorr <- rcorr(as.matrix(innate_pca_df))
corrplot(innaterescorr$r)
innate_pca <- prcomp(innate_pca_df, scale = TRUE)
summary(innate_pca)
fviz_pca_var(innate_pca, col.ind = "cos2",
             gradient.cols = c("blue", "orange", "red"),
             repel = TRUE)


compiled_df$SS_cat <- ifelse(compiled_df$SS > 0.5, 1,0)


compiled_df$socpca[!is.na(compiled_df$IC)] <- datawizard::standardise(soc_pca$x[1:(length(soc_pca$x)/3)])
compiled_df$innatepca[!is.na(compiled_df$IC)] <- datawizard::standardise(innate_pca$x[1:(length(innate_pca$x)/3)])

table(compiled_df$Age, compiled_df$Sex)
#Vector of formulas
vectorFormulas_v <- c(
  as.formula(SS_cat ~ 1),
  as.formula(SS_cat ~ innatepca + raidrate_sqrt),
  as.formula(SS_cat ~ innatepca + socpca),
  as.formula(SS_cat ~ socpca + raidrate_sqrt),
  as.formula(SS_cat ~ innatepca),
  as.formula(SS_cat ~ socpca),
  as.formula(SS_cat ~ raidrate_sqrt)
)

#Run models
modelsToCompare_l <- lapply(1:length(vectorFormulas_v), function(whichModel){
  #Fit the seasonal model centered on the week of interest
  model <- glmmTMB(
    vectorFormulas_v[whichModel][[1]],
    data = compiled_df %>% mutate(IC = (IC*(nrow(.) - 1) + 0.5)/nrow(.)),
    family =  "binomial" #beta_family(link = "logit")
  )
  diagnostics.plot.dharma(model)
  check_collinearity(model)
  
  #Return the model
  return(model)
})
names(modelsToCompare_l) <- c("null",
                              # "Demographic, Social, Experience",
                              "Demographic, Experience",
                              "Demographic, Social",
                              "Social, Experience",
                              "Demographic",
                              "Social",
                              "Experience")

#Extract model ranking
selectionTable <- aictab(modelsToCompare_l)
selectionTable
