#sensitivity testing

# Set up environment ------------------------------------------------------
rm(list = ls())
Sys.setenv(lang = "en_US")

local({
  Sys.setenv(LANGUAGE="en")
})

## Libraries ---------------------------------------------------------------
#library
library(ggplot2)
library(ggpubr)
library(AICcmodavg)
library(stringr)
library(glmmTMB)
library(patchwork)
library(forcats)
library(dplyr)
library(extrafont)
library(sjPlot)
library(Hmisc)
library(kableExtra)


#load in data
setwd("C:/Users/admin/OneDrive - Universität Zürich UZH/sensitivityTesting")
# setwd("C:/Users/paige/OneDrive - Universität Zürich UZH/sensitivityTesting")

#model selection table for social exposure (SE) days variations for IT
load("./modSel7days.RData")
load("./modSel10days.RData")
load("./modSel12days.RData")
load("./modSel14days.RData")
load("./modSel16days.RData")
load("./modSel18days.RData")
load("./modSel21days.RData")

#model selection table for social exposure (SE) days variations for SS
load("./modSelSS7days.RData")
load("./modSelSS10days.RData")
load("./modSelSS12days.RData")
load("./modSelSS14days.RData")
load("./modSelSS16days.RData")
load("./modSelSS18days.RData")
load("./modSelSS21days.RData")

#mod selection table for IC cutoff variations
load("./modSelIT3.RData")
load("./modSelIT4.RData")
load("./modSelIT5.RData")
load("./modSelIT6.RData")
load("./modSelIT7.RData")



# Social Exposure day sensitivity ######################################################


#sjPlot to print tables

selectToPrint <- function(modelsToCompare){
  selectionTable <- aictab(modelsToCompare)
  selectionTablePrint <- selectionTable %>%
    mutate('days' = "",
           'Model names' = Modnames,
           'ΔAICc' = Delta_AICc,
           Weight = AICcWt,
           'Log-likelihood' = LL,
           'Cumulative weight' = Cum.Wt) %>%
    dplyr::select('days', 'Model names', K, AICc, 'ΔAICc', 'Log-likelihood', Weight, 'Cumulative weight')
  return(selectionTablePrint)
}

#IT social exposure
tab_dfs(list(selectToPrint(modelsToCompare_l07),
             selectToPrint(modelsToCompare_l10),
             selectToPrint(modelsToCompare_l12),
             selectToPrint(modelsToCompare_l14),
             selectToPrint(modelsToCompare_l16),
             selectToPrint(modelsToCompare_l18),
             selectToPrint(modelsToCompare_l21)))

#SS social exposure
tab_dfs(list(selectToPrint(modelsToCompare_lSS07),
             selectToPrint(modelsToCompare_lSS10),
             selectToPrint(modelsToCompare_lSS12),
             selectToPrint(modelsToCompare_lSS14),
             selectToPrint(modelsToCompare_lSS16),
             selectToPrint(modelsToCompare_lSS18),
             selectToPrint(modelsToCompare_lSS21)))

#IT cutoff
tab_dfs(list(selectToPrint(modelsToCompare_lIT3),
             selectToPrint(modelsToCompare_lIT4),
             selectToPrint(modelsToCompare_lIT5),
             selectToPrint(modelsToCompare_lIT6),
             selectToPrint(modelsToCompare_lIT7)))


## prepare data for visualizations  ----------------------------


#list models to compare with different number of days for social exposure
daysExpoMS <- list(modelsToCompare_l07, modelsToCompare_l10, modelsToCompare_l12,
                   modelsToCompare_l14, modelsToCompare_l16, modelsToCompare_l18,
                   modelsToCompare_l21)
names(daysExpoMS) <- c("modelsToCompare_l07", "modelsToCompare_l10", "modelsToCompare_l12",
                       "modelsToCompare_l14", "modelsToCompare_l16", "modelsToCompare_l18",
                       "modelsToCompare_l21")

#set up forest plot data frame
sensitivitySocExp <- data.frame(days = NULL,
                                var = NULL,
                                weight = NULL,
                                lowerCI = NULL,
                                upperCI = NULL,
                                category =NULL)

#extract values for each model selection
for (j in 1:7){
  variables <- c("socpca", "innatepca", "raidrate_sqrt")
  var_weight <- data.frame(days = rep(str_sub(deparse((names(daysExpoMS)[j])), -3, -2), 3),
                           var = variables,
                           weight = rep(NA, 3),
                           lowerCI = rep(NA, 3),
                           upperCI = rep(NA, 3))
  
  #do this with every variable
  for (i in (1:length(variables))){
    var_weight$weight[var_weight$var == variables[i]] <- as.numeric(modavg(daysExpoMS[[j]],
                                                                           parm = variables[i])[3][[1]])
    var_weight$lowerCI[var_weight$var == variables[i]] <- as.numeric(modavg(daysExpoMS[[j]],
                                                                            parm = variables[i])[6][[1]])
    var_weight$upperCI[var_weight$var == variables[i]] <- as.numeric(modavg(daysExpoMS[[j]],
                                                                            parm = variables[i])[7][[1]])
    
  }
  
  for (i in (1:length(variables))){
    var_weight$weight[var_weight$var == variables[i]] <- as.numeric(modavgShrink(daysExpoMS[[j]],
                                                                                 parm = variables[i])[3][[1]])
    var_weight$lowerCI[var_weight$var == variables[i]] <- as.numeric(modavgShrink(daysExpoMS[[j]],
                                                                                  parm = variables[i])[6][[1]])
    var_weight$upperCI[var_weight$var == variables[i]] <- as.numeric(modavgShrink(daysExpoMS[[j]],
                                                                                  parm = variables[i])[7][[1]])
    
  }
  var_weight$category <- c("Social exposure\n to experiment", #"Social", 
                           "Demographic traits", #"Demographic",
                           "Anthropogenic raiding\n experience")
  sensitivitySocExp <- rbind(sensitivitySocExp, var_weight)
  
  #prepare sensitivity tables to be printed
  # nam <- paste("selectionTable", str_sub(deparse((names(daysExpoMS)[j])), -3, -2), "days", sep = "")
  # assign(nam, as.data.frame(aictab((daysExpoMS)[[j]])) %>%
  #          mutate('Model names' = Modnames,
  #                 AICc = round(AICc, 3),
  #                 'ΔAICc' = round(Delta_AICc, 3),
  #                 Weight = round(AICcWt, 3),
  #                 'Log-likelihood' = round(LL, 3),
  #                 'Cumulative weight' = round(Cum.Wt, 3)) %>%
  #          select('Model names', K, AICc, 'ΔAICc', 'Log-likelihood', Weight, 'Cumulative weight'))
  
}

## plot tables ----------------------------
# ggarrange(plotlist = list(print(ggtexttable((selectionTable07days), rows = NULL, theme = ttheme("light", base_size = 10)) %>% table_cell_font(row = 2:7, column = 1:7, face = "bold", size = 10)),
#                           print(ggtexttable(selectionTable10days, rows = NULL, theme = ttheme("light", base_size = 10)) %>% table_cell_font(row = 2:6, column = 1:7, face = "bold", size = 10)),
#                           print(ggtexttable(selectionTable12days, rows = NULL, theme = ttheme("light", base_size = 10)) %>% table_cell_font(row = 2:6, column = 1:7, face = "bold", size = 10)),
#                           print(ggtexttable(selectionTable14days, rows = NULL, theme = ttheme("light", tbody.style = tbody_style(color = "black", fill = "#e8f3de", size = 10))) %>%
#                                   table_cell_font(row = 2:6, column = 1:6, face = "bold", size = 10)) ),
#           labels = c("7 days", "10 days", "12 days", "14 days"),
#           font.label = list(size = 11),
#           nrow = 4
# )
# 
# ggarrange(plotlist = list(print(ggtexttable(selectionTable16days, rows = NULL, theme = ttheme("light", base_size = 10)) %>% table_cell_font(row = 2:6, column = 1:7, face = "bold", size = 10)),
#                           print(ggtexttable(selectionTable18days, rows = NULL, theme = ttheme("light", base_size = 10)) %>% table_cell_font(row = 2:6, column = 1:7, face = "bold", size = 10)),
#                           print(ggtexttable(selectionTable21days, rows = NULL, theme = ttheme("light", base_size = 10)) %>% table_cell_font(row = 2:6, column = 1:7, face = "bold", size = 10))),
#           labels = c("16 days", "18 days", "21 days"),
#           font.label = list(size = 11),
#           nrow = 3
# )


## forest plot  ----------------------------

sensitivitySocExp <- sensitivitySocExp %>%
  mutate(days = as.numeric(days),
         weight = as.numeric(weight)) %>%
  arrange(category)

#forest plot
SE_forest <- ggplot(data = sensitivitySocExp, aes(x = days, y = weight)) +
  geom_point(aes(x = days, y = weight, colour = factor(category)),size = 1) +  
  geom_errorbar(aes(ymin = lowerCI, ymax = upperCI, colour = factor(category)), width = 0) +
  geom_hline(yintercept = 0) + 
  theme_bw() + 
  scale_colour_manual(values = c("#81B29A", "#F2CC8F", "#B48168")) +
  theme(strip.text = element_text(size = 10),
        axis.text = element_text(size=10),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        legend.position="none",
        strip.background = element_rect(fill = "grey94")) +
  facet_wrap(. ~ fct_inorder(category), ncol = 1) +
  scale_size(guide = 'none') +
  ylab("Variable estimate") +
  xlab("Days of social exposure") +
  coord_flip() +
  scale_x_continuous(breaks = seq(6, 22, by = 2))


## bar plot -------------------------------------------


#list models to compare with different number of days for social exposure
daysExpoST <- list(aictab(modelsToCompare_l07), aictab(modelsToCompare_l10), aictab(modelsToCompare_l12),
                   aictab(modelsToCompare_l14), aictab(modelsToCompare_l16), aictab(modelsToCompare_l18),
                   aictab(modelsToCompare_l21))
names(daysExpoST) <- c("selectionTable07days", "selectionTable10days", "selectionTable12days",
                       "selectionTable14days", "selectionTable16days", "selectionTable18days",
                       "selectionTable21days")

sensitivitySocExpST <- data.frame(days = NULL,
                                  group_var = NULL,
                                  p_select_ic = NULL,
                                  avg_est_ic = NULL,
                                  group_var = NULL)


for (j in 1:7){
  group_var <- c("Social", "Experience", "Demographic")
  plot_df <- data.frame(days = rep(str_sub(deparse((names(daysExpoST)[j])), 16, -6), 3),
                        group_var = group_var,
                        p_select_ic = rep(NA, length(group_var)),
                        avg_est_ic = rep(NA, length(group_var)))

  for(i in (1:length(group_var))){
    plot_df$p_select_ic[i] <- sum(as.data.frame(daysExpoST[[j]])$AICcWt[grep(plot_df$group_var[i],
                                                                            as.data.frame(daysExpoST[[j]])$'Modnames')])
    plot_df$avg_est_ic[i] <- mean(abs(as.data.frame(daysExpoST[[j]])$'Delta_AICc'[grep(plot_df$group_var[i],
                                                                                    as.data.frame(daysExpoST[[j]])$'Modnames')]))
  }

  plot_df$group_var <- c("Social exposure\n to experiment", "Anthropogenic raiding\n experience", "Demographic traits")

  sensitivitySocExpST <- rbind(sensitivitySocExpST, plot_df)
}

sensitivitySocExpST <- sensitivitySocExpST %>%
  mutate(days = as.numeric(days),
         p_select_ic = as.numeric(p_select_ic),
         group_var = factor(group_var, levels = c("Anthropogenic raiding\n experience",  "Demographic traits", "Social exposure\n to experiment"),
                            ordered = TRUE)) %>%
  arrange(group_var)

SE_bar <- ggplot(sensitivitySocExpST, aes(x=fct_inorder(group_var), y=p_select_ic, fill = group_var)) +
  geom_bar(stat = "identity", color = "black") + 
  facet_wrap(. ~ days, ncol = 3,
             labeller = labeller(days = 
                                   c("7" = "7 days",
                                     "10" = "10 days",
                                     "12" = "12 days",
                                     "14" = "14 days",
                                     "16" = "16 days",
                                     "18" = "18 days",
                                     "21" = "21 days"))) +
  theme_bw() + 
  ylab("Model selection weight") +
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=10),
        strip.background = element_rect(fill = "grey94"),
        legend.position = c(1, 0),
        legend.justification = c(1, 0),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("Social exposure\n to experiment" = "#B48168",
                                "Anthropogenic raiding\n experience" = "#81B29A",
                                "Demographic traits" = "#F2CC8F")) +
  guides(fill = guide_legend(title = "Dimension")) +
  xlab("Models") +
  scale_size(guide = 'none') +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10))
  # ggtitle("IT model: sensitivity test for days of Social exposure\n to experiment")

## version with models instead of variables
# sensitivitySocExpST <- data.frame(days = NULL,
#                                   model_name = NULL,
#                                   p_select_ic = NULL,
#                                   avg_est_ic = NULL)
# #for each model selection
# for (j in 1:7){
#   model_name <- data.frame(daysExpoST[[j]])$Modnames
#   plot_df <- data.frame(days = rep(str_sub(deparse((names(daysExpoST)[j])), 16, -6), length((daysExpoST[[j]])[,1])),
#                         model_name = model_name,
#                         p_select_ic = rep(NA, length((daysExpoST[[j]])[,1])),
#                         avg_est_ic = rep(NA, length((daysExpoST[[j]])[,1])))
#   
#   for(i in (1:length(model_name))){
#     plot_df$p_select_ic[i] <- data.frame(daysExpoST[[j]])$AICcWt[i]
#     plot_df$avg_est_ic[i] <- data.frame(daysExpoST[[j]])$Delta_AICc[i]
#   }
#   
#   sensitivitySocExpST <- rbind(sensitivitySocExpST, plot_df)
# }
# 
# sensitivitySocExpST <- sensitivitySocExpST %>%
#   mutate(days = as.numeric(days),
#          p_select_ic = as.numeric(p_select_ic),
#          model_name = factor(model_name, levels = c("Experience", "Demographic, Experience", 
#                              "Demographic", "Demographic, Social", "Social", "Social, Experience", 
#                              "Demographic, Social, Experience", "null"), ordered = TRUE)) %>%
#   arrange(model_name)
# 
# SE_bar <- ggplot(sensitivitySocExpST, aes(x=fct_inorder(model_name), y=p_select_ic, fill = model_name)) +
#   geom_bar(stat = "identity", color = "black") + 
#   facet_wrap(. ~ days, ncol = 3,
#              labeller = labeller(days = 
#                                    c("7" = "7 days",
#                                      "10" = "10 days",
#                                      "12" = "12 days",
#                                      "14" = "14 days",
#                                      "16" = "16 days",
#                                      "18" = "18 days",
#                                      "21" = "21 days"))) +
#   theme_set(theme_bw(base_family = 'Calibri')) + 
#   scale_fill_brewer("Models", palette="Dark2") +
#   ylab("Model selection weight") +
#   theme(axis.text.y = element_text(size=10),
#         axis.text.x = element_blank(),
#         axis.title.x = element_blank(),
#         axis.title.y = element_text(size=10),
#         strip.background = element_rect(fill = "grey94"),
#         legend.position = c(1, 0),
#         legend.justification = c(1, 0),
#         axis.ticks.x = element_blank()) +
#   xlab("Models") +
#   scale_size(guide = 'none') +
#   scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
#   ggtitle("Sensitivity test for days of Social exposure\n to experiment")

SE_bar + SE_forest + plot_layout(widths = c(3, 1))

# "#B48168"
# "#9B9A81"
# "#81B29A"
# "#BABF95"
# "#F2CC8F"
# "#D3A77C"
# "#BFAA84"



# SS model days exposure ######################################################

## prepare data for visualizations  ----------------------------


#list models to compare with different number of days for social exposure
daysExpoMS_SS <- list(modelsToCompare_lSS07, modelsToCompare_lSS10, modelsToCompare_lSS12,
                   modelsToCompare_lSS14, modelsToCompare_lSS16, modelsToCompare_lSS18,
                   modelsToCompare_lSS21)
names(daysExpoMS_SS) <- c("modelsToCompare_lSS07", "modelsToCompare_lSS10", "modelsToCompare_lSS12",
                       "modelsToCompare_lSS14", "modelsToCompare_lSS16", "modelsToCompare_lSS18",
                       "modelsToCompare_lSS21")

#set up forest plot data frame
sensitivitySocExpSS <- data.frame(days = NULL,
                                var = NULL,
                                weight = NULL,
                                lowerCI = NULL,
                                upperCI = NULL,
                                category =NULL)

#extract values for each model selection
for (j in 1:7){
  variables <- c("socpca", "innatepca", "raidrate_sqrt")
  var_weight <- data.frame(days = rep(str_sub(deparse((names(daysExpoMS_SS)[j])), -3, -2), 3),
                           var = variables,
                           weight = rep(NA, 3),
                           lowerCI = rep(NA, 3),
                           upperCI = rep(NA, 3))
  
  #do this with every variable
  for (i in (1:length(variables))){
    var_weight$weight[var_weight$var == variables[i]] <- as.numeric(modavg(daysExpoMS_SS[[j]],
                                                                           parm = variables[i])[3][[1]])
    var_weight$lowerCI[var_weight$var == variables[i]] <- as.numeric(modavg(daysExpoMS_SS[[j]],
                                                                            parm = variables[i])[6][[1]])
    var_weight$upperCI[var_weight$var == variables[i]] <- as.numeric(modavg(daysExpoMS_SS[[j]],
                                                                            parm = variables[i])[7][[1]])
    
  }
  
  for (i in (1:length(variables))){
    var_weight$weight[var_weight$var == variables[i]] <- as.numeric(modavgShrink(daysExpoMS_SS[[j]],
                                                                                 parm = variables[i])[3][[1]])
    var_weight$lowerCI[var_weight$var == variables[i]] <- as.numeric(modavgShrink(daysExpoMS_SS[[j]],
                                                                                  parm = variables[i])[6][[1]])
    var_weight$upperCI[var_weight$var == variables[i]] <- as.numeric(modavgShrink(daysExpoMS_SS[[j]],
                                                                                  parm = variables[i])[7][[1]])
    
  }
  var_weight$category <- c("Social exposure\n to experiment", #"Social", 
                           "Demographic traits", #"Demographic",
                           "Anthropogenic raiding\n experience")
  sensitivitySocExpSS <- rbind(sensitivitySocExpSS, var_weight)
  
}


## forest plot  ----------------------------

sensitivitySocExpSS <- sensitivitySocExpSS %>%
  mutate(days = as.numeric(days),
         weight = as.numeric(weight)) %>%
  arrange(category)

#forest plot
SE_forestSS <- ggplot(data = sensitivitySocExpSS, aes(x = days, y = weight)) +
  geom_point(aes(x = days, y = weight, colour = factor(category)),size = 1) +  
  geom_errorbar(aes(ymin = lowerCI, ymax = upperCI, colour = factor(category)), width = 0) +
  geom_hline(yintercept = 0) + 
  theme_bw() + 
  scale_colour_manual(values = c("#81B29A", "#F2CC8F", "#B48168")) +
  theme(strip.text = element_text(size = 10),
        axis.text = element_text(size=10),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        legend.position="none",
        strip.background = element_rect(fill = "grey94")) +
  facet_wrap(. ~ fct_inorder(category), ncol = 1) +
  scale_size(guide = 'none') +
  ylab("Variable estimate") +
  xlab("Days of social exposure") +
  coord_flip() +
  scale_x_continuous(breaks = seq(6, 22, by = 2))


## bar plot -------------------------------------------


#list models to compare with different number of days for social exposure
daysExpoST_SS <- list(aictab(modelsToCompare_lSS07), aictab(modelsToCompare_lSS10), aictab(modelsToCompare_lSS12),
                   aictab(modelsToCompare_lSS14), aictab(modelsToCompare_lSS16), aictab(modelsToCompare_lSS18),
                   aictab(modelsToCompare_lSS21))
names(daysExpoST_SS) <- c("selectionTableSS07days", "selectionTableSS10days", "selectionTableSS12days",
                       "selectionTable14SSdays", "selectionTableSS16days", "selectionTableSS18days",
                       "selectionTable21SSdays")

sensitivitySocExpSTSS <- data.frame(days = NULL,
                                  group_var = NULL,
                                  p_select_ic = NULL,
                                  avg_est_ic = NULL,
                                  group_var = NULL)


for (j in 1:7){
  group_var <- c("Social", "Experience", "Demographic")
  plot_df <- data.frame(days = rep(str_sub(deparse((names(daysExpoST)[j])), 16, -6), 3),
                        group_var = group_var,
                        p_select_ic = rep(NA, length(group_var)),
                        avg_est_ic = rep(NA, length(group_var)))
  
  for(i in (1:length(group_var))){
    plot_df$p_select_ic[i] <- sum(as.data.frame(daysExpoST_SS[[j]])$AICcWt[grep(plot_df$group_var[i],
                                                                             as.data.frame(daysExpoST_SS[[j]])$'Modnames')])
    plot_df$avg_est_ic[i] <- mean(abs(as.data.frame(daysExpoST_SS[[j]])$'Delta_AICc'[grep(plot_df$group_var[i],
                                                                                       as.data.frame(daysExpoST_SS[[j]])$'Modnames')]))
  }
  
  plot_df$group_var <- c("Social exposure\n to experiment", "Anthropogenic raiding\n experience", "Demographic traits")
  
  sensitivitySocExpSTSS <- rbind(sensitivitySocExpSTSS, plot_df)
}

sensitivitySocExpSTSS <- sensitivitySocExpSTSS %>%
  mutate(days = as.numeric(days),
         p_select_ic = as.numeric(p_select_ic),
         group_var = factor(group_var, levels = c("Anthropogenic raiding\n experience",  "Demographic traits", "Social exposure\n to experiment"),
                            ordered = TRUE)) %>%
  arrange(group_var)

SE_barSS <- ggplot(sensitivitySocExpSTSS, aes(x=fct_inorder(group_var), y=p_select_ic, fill = group_var)) +
  geom_bar(stat = "identity", color = "black") + 
  facet_wrap(. ~ days, ncol = 3,
             labeller = labeller(days = 
                                   c("7" = "7 days",
                                     "10" = "10 days",
                                     "12" = "12 days",
                                     "14" = "14 days",
                                     "16" = "16 days",
                                     "18" = "18 days",
                                     "21" = "21 days"))) +
  theme_bw() + 
  ylab("Model selection weight") +
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=10),
        strip.background = element_rect(fill = "grey94"),
        legend.position = c(1, 0),
        legend.justification = c(1, 0),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("Social exposure\n to experiment" = "#B48168",
                               "Anthropogenic raiding\n experience" = "#81B29A",
                               "Demographic traits" = "#F2CC8F")) +
  guides(fill = guide_legend(title = "Dimension")) +
  xlab("Models") +
  scale_size(guide = 'none') +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10))
  # ggtitle("SS model: sensitivity test for days of Social exposure\n to experiment")


SE_barSS + SE_forestSS + plot_layout(widths = c(3, 1))






# IC cutoff sensitivity ######################################################

#list models to compare with different number of days for social exposure
icCutMS <- list(modelsToCompare_lIT3, modelsToCompare_lIT4, modelsToCompare_lIT5,
                modelsToCompare_lIT6, modelsToCompare_lIT7)
names(icCutMS) <- c("modelsToCompare_l3", "modelsToCompare_l4", "modelsToCompare_l5",
                    "modelsToCompare_l6", "modelsToCompare_l7")

#set up forest plot data frame
sensitivityIC <- data.frame(cutoff = NULL,
                            var = NULL,
                            weight = NULL,
                            lowerCI = NULL,
                            upperCI = NULL,
                            category =NULL)

#extract values for each model selection
for (j in 1:5){
  variables <- c("socpca", "innatepca", "raidrate_sqrt")
  var_weight <- data.frame(days = rep(str_sub(deparse((names(icCutMS)[j])), -3, -2), 3),
                           var = variables,
                           weight = rep(NA, 3),
                           lowerCI = rep(NA, 3),
                           upperCI = rep(NA, 3))
  
  #do this with every variable
  for (i in (1:length(variables))){
    var_weight$weight[var_weight$var == variables[i]] <- as.numeric(modavg(icCutMS[[j]],
                                                                           parm = variables[i])[3][[1]])
    var_weight$lowerCI[var_weight$var == variables[i]] <- as.numeric(modavg(icCutMS[[j]],
                                                                            parm = variables[i])[6][[1]])
    var_weight$upperCI[var_weight$var == variables[i]] <- as.numeric(modavg(icCutMS[[j]],
                                                                            parm = variables[i])[7][[1]])
    
  }
  
  for (i in (1:length(variables))){
    var_weight$weight[var_weight$var == variables[i]] <- as.numeric(modavgShrink(icCutMS[[j]],
                                                                                 parm = variables[i])[3][[1]])
    var_weight$lowerCI[var_weight$var == variables[i]] <- as.numeric(modavgShrink(icCutMS[[j]],
                                                                                  parm = variables[i])[6][[1]])
    var_weight$upperCI[var_weight$var == variables[i]] <- as.numeric(modavgShrink(icCutMS[[j]],
                                                                                  parm = variables[i])[7][[1]])
    
  }
  var_weight$category <- c("Social exposure\n to experiment", #"Social", 
                           "Demographic traits", #"Demographic",
                           "Anthropogenic raiding\n experience")
  sensitivitySocExp <- rbind(sensitivitySocExp, var_weight)
  
  #prepare sensitivity tables to be printed
  # namIC <- paste("selectionTable0", str_sub(deparse((names(icCutMS)[j])), -2, -2), sep = "")
  # assign(namIC, as.data.frame(aictab((icCutMS)[[j]])) %>%
  #          mutate('Model names' = Modnames,
  #                 AICc = round(AICc, 3),
  #                 'ΔAICc' = round(Delta_AICc, 3),
  #                 Weight = round(AICcWt, 3),
  #                 'Log-likelihood' = round(LL, 3),
  #                 'Cumulative weight' = round(Cum.Wt, 3)) %>%
  #          select('Model names', K, AICc, 'ΔAICc', 'Log-likelihood', Weight, 'Cumulative weight'))
  
}

## plot tables ----------------------------
# ggarrange(plotlist = list(print(ggtexttable(selectionTable03, rows = NULL, theme = ttheme("light", base_size = 10)) %>% table_cell_font(row = 2:7, column = 1:7, face = "bold", size = 10)),
#                           print(ggtexttable(selectionTable04, rows = NULL, theme = ttheme("light", base_size = 10)) %>% table_cell_font(row = 2:7, column = 1:7, face = "bold", size = 10)),
#                           print(ggtexttable(selectionTable05, rows = NULL, theme = ttheme("light", tbody.style = tbody_style(color = "black", fill = "#e8f3de", size = 10))) %>%
#                                   table_cell_font(row = 2:6, column = 1:7, face = "bold", size = 10)) ),
#           labels = c("IC > 0.3", "IC > 0.4", "IC > 0.5"),
#           font.label = list(size = 11),
#           nrow = 3
# )
# 
# ggarrange(plotlist = list(print(ggtexttable(selectionTable06, rows = NULL, theme = ttheme("light", base_size = 10)) %>% table_cell_font(row = 2:7, column = 1:7, face = "bold", size = 10)),
#                           print(ggtexttable(selectionTable07, rows = NULL, theme = ttheme("light", base_size = 10)) %>% table_cell_font(row = 2:8, column = 1:7, face = "bold", size = 10))),
#           labels = c("IC > 0.6", "IC > 0.7"),
#           font.label = list(size = 11),
#           nrow = 2
# )
# 
# ggarrange(plotlist = list(print(ggtexttable(selectionTable03, rows = NULL, theme = ttheme("light", base_size = 10)) %>% table_cell_font(row = 2:7, column = 1:7, face = "bold", size = 10)),
#                           print(ggtexttable(selectionTable04, rows = NULL, theme = ttheme("light", base_size = 10)) %>% table_cell_font(row = 2:7, column = 1:7, face = "bold", size = 10)),
#                           print(ggtexttable(selectionTable05, rows = NULL, theme = ttheme("light", tbody.style = tbody_style(color = "black", fill = "#e8f3de", size = 10))) %>%
#                                   table_cell_font(row = 2:6, column = 1:7, face = "bold", size = 10)),
#                           print(ggtexttable(selectionTable06, rows = NULL, theme = ttheme("light", base_size = 10)) %>% table_cell_font(row = 2:7, column = 1:7, face = "bold", size = 10)),
#                           print(ggtexttable(selectionTable07, rows = NULL, theme = ttheme("light", base_size = 10)) %>% table_cell_font(row = 2:8, column = 1:7, face = "bold", size = 10))),
#           labels = c("IC > 0.3", "IC > 0.4", "IC > 0.5", "IC > 0.6", "IC > 0.7"),
#           font.label = list(size = 11),
#           nrow = 5
# )


## make forest plot and bar graph  ----------------------------

#extract values for each model selection
for (j in 1:5){
  variables <- c("socpca", "innatepca", "raidrate_sqrt")
  var_weight <- data.frame(cutoff = rep(str_sub(deparse((names(icCutMS)[j])), -2, -2), 3),
                           var = variables,
                           weight = rep(NA, 3),
                           lowerCI = rep(NA, 3),
                           upperCI = rep(NA, 3))
  
  #do this with every variable
  for (i in (1:length(variables))){
    var_weight$weight[var_weight$var == variables[i]] <- as.numeric(modavg(icCutMS[[j]],
                                                                           parm = variables[i])[3][[1]])
    var_weight$lowerCI[var_weight$var == variables[i]] <- as.numeric(modavg(icCutMS[[j]],
                                                                            parm = variables[i])[6][[1]])
    var_weight$upperCI[var_weight$var == variables[i]] <- as.numeric(modavg(icCutMS[[j]],
                                                                            parm = variables[i])[7][[1]])
    
  }
  
  for (i in (1:length(variables))){
    var_weight$weight[var_weight$var == variables[i]] <- as.numeric(modavgShrink(icCutMS[[j]],
                                                                                 parm = variables[i])[3][[1]])
    var_weight$lowerCI[var_weight$var == variables[i]] <- as.numeric(modavgShrink(icCutMS[[j]],
                                                                                  parm = variables[i])[6][[1]])
    var_weight$upperCI[var_weight$var == variables[i]] <- as.numeric(modavgShrink(icCutMS[[j]],
                                                                                  parm = variables[i])[7][[1]])
    
  }
  var_weight$category <- c("Social exposure\n to experiment", #"Social", 
                           "Demographic traits", #"Demographic",
                           "Anthropogenic raiding\n experience")
  sensitivityIC <- rbind(sensitivityIC, var_weight)
}

sensitivityIC <- sensitivityIC %>%
  mutate(cutoff = as.numeric(cutoff),
         cutoff = ifelse(cutoff == 3, 0.3,
                         ifelse(cutoff == 4, 0.4,
                                ifelse(cutoff == 5, 0.5,
                                      ifelse(cutoff == 6, 0.6, 0.7)))),
         weight = as.numeric(weight)) %>%
  arrange(category)


#forest plot
IC_forest <- ggplot(data = sensitivityIC, aes(x = cutoff, y = weight)) +
  geom_point(aes(x = cutoff, y = weight, colour = factor(category)),size = 1) +  
  geom_errorbar(aes(ymin = lowerCI, ymax = upperCI, colour = factor(category)), width = 0) +
  geom_hline(yintercept = 0) + 
  theme_bw() + 
  scale_colour_manual(values = c("#81B29A", "#F2CC8F", "#B48168")) +
  theme(strip.text = element_text(size = 10),
        axis.text = element_text(size=10),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        legend.position="none",
        strip.background = element_rect(fill = "grey94")) +
  facet_wrap(. ~ fct_inorder(category), ncol = 1) +
  scale_size(guide = 'none') +
  ylab("Variable estimate") +
  xlab("IT confidence interval cutoff") +
  coord_flip()
  # scale_x_continuous(breaks = seq(6, 22, by = 2))


## bar plot -------------------------------------------

#list models to compare with different number of days for social exposure
icCutST <- list(aictab(modelsToCompare_lIT3), aictab(modelsToCompare_lIT4), aictab(modelsToCompare_lIT5),
                aictab(modelsToCompare_lIT6), aictab(modelsToCompare_lIT7))
names(icCutST) <- c("selectionTable03", "selectionTable04", "selectionTable05",
                    "selectionTable06", "selectionTable07")

sensitivityIC_ST <- data.frame(cutoff = NULL,
                               model_name = NULL,
                               p_select_ic = NULL,
                               avg_est_ic = NULL)

for (j in 1:5){
  group_var <- c("Social", "Experience", "Demographic")
  plot_df <- data.frame(cutoff = rep(str_sub(deparse((names(icCutST)[j])), 16, -2), 3),
                        group_var = group_var,
                        p_select_ic = rep(NA, length(group_var)),
                        avg_est_ic = rep(NA, length(group_var)))
  
  for(i in (1:length(group_var))){
    plot_df$p_select_ic[i] <- sum(as.data.frame(icCutST[[j]])$AICcWt[grep(plot_df$group_var[i],
                                                                             as.data.frame(icCutST[[j]])$'Modnames')])
    plot_df$avg_est_ic[i] <- mean(abs(as.data.frame(icCutST[[j]])$'Delta_AICc'[grep(plot_df$group_var[i],
                                                                                       as.data.frame(icCutST[[j]])$'Modnames')]))
  }
  
  plot_df$group_var <- c("Social exposure\n to experiment", "Anthropogenic raiding\n experience", "Demographic traits")
  
  sensitivityIC_ST <- rbind(sensitivityIC_ST, plot_df)
}

sensitivityIC_ST <- sensitivityIC_ST %>%
  mutate(p_select_ic = as.numeric(p_select_ic),
         group_var = factor(group_var, levels = c("Anthropogenic raiding\n experience",  "Demographic traits", "Social exposure\n to experiment"),
                            ordered = TRUE)) %>%
  arrange(group_var)




IC_bar <- ggplot(sensitivityIC_ST, aes(x=fct_inorder(group_var), y=p_select_ic, fill = group_var)) +
  geom_bar(stat = "identity", color = "black") + 
  facet_wrap(. ~ cutoff, ncol = 3,
             labeller = labeller(cutoff = c("03" = "IT > 0.3",
                                            "04" = "IT > 0.4",
                                            "05" = "IT > 0.5",
                                            "06" = "IT > 0.6",
                                            "07" = "IT > 0.7"))) +
  theme_bw() + 
  ylab("Model selection weight") +
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=10),
        strip.background = element_rect(fill = "grey94"),
        legend.position = c(1, 0),
        legend.justification = c(1, 0),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("Social exposure\n to experiment" = "#B48168",
                               "Anthropogenic raiding\n experience" = "#81B29A",
                               "Demographic traits" = "#F2CC8F")) +
  guides(fill = guide_legend(title = "Dimension")) +
  xlab("Models") +
  scale_size(guide = 'none') +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) 
  # ggtitle("Sensitivity test for IT confidence interval cutoff")

IC_bar + IC_forest + plot_layout(widths = c(3, 1))

# #for each model selection
# for (j in 1:5){
#   model_name <- data.frame(icCutST[[j]])$Modnames
#   plot_df <- data.frame(cutoff = rep(str_sub(deparse((names(icCutST)[j])), 16, -2), length((icCutST[[j]])[,1])),
#                         model_name = model_name,
#                         p_select_ic = rep(NA, length((icCutST[[j]])[,1])),
#                         avg_est_ic = rep(NA, length((icCutST[[j]])[,1])))
#   
#   for(i in (1:length(model_name))){
#     plot_df$p_select_ic[i] <- data.frame(icCutST[[j]])$AICcWt[i]
#     plot_df$avg_est_ic[i] <- data.frame(icCutST[[j]])$Delta_AICc[i]
#   }
#   
#   sensitivityIC_ST <- rbind(sensitivityIC_ST, plot_df)
# }
# sensitivityIC_ST <- sensitivityIC_ST %>%
#   mutate(p_select_ic = as.numeric(p_select_ic), 
#          model_name = factor(model_name, levels = c("Experience", "Demographic, Experience", 
#                                                                                            "Demographic", "Demographic, Social", "Social", "Social, Experience", 
#                                                                                            "Demographic, Social, Experience", "null"), ordered = TRUE)) %>%
#   arrange(model_name)
# IC_bar <- ggplot(sensitivityIC_ST, aes(x=fct_inorder(model_name), y=p_select_ic, fill = model_name)) +
#   geom_bar(stat = "identity", color = "black") + 
#   facet_wrap(. ~ cutoff, ncol = 3,
#              labeller = labeller(cutoff = 
#                                    c("03" = "IC > 0.3",
#                                      "04" = "IC > 0.4",
#                                      "05" = "IC > 0.5",
#                                      "06" = "IC > 0.6",
#                                      "07" = "IC > 0.7"))) +
#   theme_bw() + 
#   scale_fill_brewer("Models", palette="Dark2") +
#   ylab("Model selection weight") +
#   theme(axis.text.y = element_text(size=10),
#         axis.text.x = element_blank(),
#         axis.title.x = element_blank(),
#         axis.title.y = element_text(size=10),
#         strip.background = element_rect(fill = "grey94"),
#         legend.position = c(1, 0),
#         legend.justification = c(1, 0),
#         axis.ticks.x = element_blank()) +
#   xlab("Models") +
#   scale_size(guide = 'none') +
#   scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
#   ggtitle("Sensitivity test for IC confidence interval cutoff")


