# Hepatitis C infectious disease model
Model of hepatitis C transmission among people who inject drugs and immigrants in Norway used in Whittaker, Midtbø and Kløvstad (submitted)

## How to run
The model is run by standing in the model/ directory and sourcing `run_hcv_model_odin.R`.
You will need to look at both that file and also `plots.R` to see which libraries are required in order to run. For the `odin/dust/odin.dust/mcstate` package suite, we recommend installing the latest main-branch versions straight from their respective github repositories (see readme in those repos).

## Data availability
All data required to run the model is included, with one notable exception: The country-specific estimates for prevalence of chronic hepatitis C are from a non-open-access publication[^1]. Therefore, we inlcude here only a file with mock numbers to demonstrate the format, and it is up to the user to obtain the correct numbers and include them in the model.

[^1]: The Polaris Observatory HCV Collaborators, Global change in hepatitis C virus prevalence and cascade of care between 2015 and 2020: a modelling study. The Lancet Gastroenterology & Hepatology (2022), https://doi.org/10.1016/S2468-1253(21)00472-6