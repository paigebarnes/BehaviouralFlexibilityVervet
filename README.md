# BehaviouralFlexibilityVervet

This repository contains codes and data for reproducing simulations of the manuscript: *"Behavioural flexibility and its drivers in semi-urban vervet monkeys"* by Paige Barnes, Benjamin Robira, Stephanie Mercier, and Sofia Forss.

The associated manuscript has been submitted and is currently under review. Hence, the code has gone through only a limited cleaning so as to allow readability. It might still include commented lines for explored - but not included in the manuscript - parts, in case it is useful for the reviewing process. Hopefully, once the manuscript is accepted, it will go through a substantial "layout" editing before being added to a repository.

The associated data can be available upon request, but present here are the scripts to show how the measures of interest were computed from the data, as well as the subsequent analyses. These are as follows:

+ __probabalisticModel.R__: Uses experimental trial data to determine the scores for each individual of: simple/technical switching tendency, simple/technical innovativeness, and learning sensitivity. This is done using a probabilistic model (version 7 in the script, this took some fine tuning).
+ __fitModels.R__: Fits the models for questions i and ii.
  + (i) This models innovativeness and switch tendency as a function of raiding rates, a pca to summarize demographic traits, and a pca to summarize the social exposure to the box experiment.
  + (ii) Additionally, models innovativeness, switch tendency, and learning sensitivity in beta regression models to check for the presence of a behavioral flexibility syndrome.
+ __successModels.R__: Models the task success and human food consumption each as a function of innovativeness and switch tendency.
+ __simulationProbMod.R__: For model validation (in the supplementary material, not the main text). Tests the probabilistic model by simulating trial sequences based on parameter output, then running these through the model for empirical vs. simulated parameter comparisons.
+ __sensitivityTesting.R__: For model validation (in the supplementary material, not the main text). Compares outcomes when "arbitrary decisions" varied in our models (i.e. whether we counted 14 days as a salient number of days for social exposure to the box or increased/decreased this number).

For questions or comments, I can be reached at paigebarnes1115@gmail.com.
