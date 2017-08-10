# TUSCAN workflow

**data** contains the data files used to construct the model

- Indel.txt: feature matrix for the Indel dataset (Chari et al 2015)
- FlowCytometry.txt: feature matrix for the FC dataset (Doench et al 2014)
- Resistance.txt: feautre matrix for the RES dataset (Doench et al 2016)

**code** contains the scripts used to construct the models 

- classifier_workflow.R: workflow for construction of the classification model
- regression_workflow.R: workflow for construction of the regression model
- random_state_seed.RData: loading this file sets the Random Seed, ensuring that the results are always fully reproducible

