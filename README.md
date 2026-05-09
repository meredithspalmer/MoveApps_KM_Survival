# KM Survival App 

Github repository: 

https://github.com/meredithspalmer/MoveApps_Survival

## Description

Perform basic Kaplan-Meier survival analyses and optional group comparisons via the log-rank test.

These analyses can be performed across an entire dataset, within defined time periods, and/or across data subsets. 

## Documentation

This app implements fundamental Kaplan-Meier (KM) survival estimation functions. It produces life tables, survival curves, cumulative hazard plots, and, if applicable, statistical comparisons of per-group survival estimation for different attributes.  

**Kaplan-Meier Survival Estimation:** The KM estimator is a non-parametric method used to estimate the survival function, that is, the probability that an individual survives past time t, from lifetime data. This analysis allows for:
- *Right-censoring*, where the exact time of death is unknown for some individuals because they are still alive at the end of the study period, lost to follow-up (e.g., collar failure), or exit the study period alive for other reasons. 
- *Staggered entry* (also called left truncation or delayed entry), where individuals enter the study at different times rather than all starting at the same baseline. 

**KM Survival Curves:** These plots depict the estimated probability of survival (e.g., remaining mortality-free) over time. They include visual representations of when events (deaths) occur. 

**Cumulative Hazard Plots:** These plots depict the the total accumulated risk (e.g., the expected number of mortalities) the population has experienced up to a specific time.

**Log-rank Test:** Log-rank tests (also known as Mantel–Cox tests) compare survival distributions across groups. The log-rank test assesses whether there are statistically significant differences in survival between two or more independent groups. The grouping variable is user-defined: any column name present in the data can be entered (e.g., sex, life stage, attachment type, tag model, survival year).


Data subsetting: 
- Data can be subset up to two times based on specific variables (e.g., only females or only adult females) - this then allows the users to perform comparisons within this smaller subset (e.g., comparing survival of adult females with different collar types). 
- Users can enter additional information to define survival years and animal life stages for comparison across these variables. For life stages, data must contain the column `animal_birth_hatch_year` and a data frame mapping age to life stage must be uploaded (see template in this repro). 
- Users can also define a study period, censor data to exclude post-capture mortality events, and specify how to handle missing time-stamp information. 

Data pre-processing includes:  
- Removing or updating empty or invalid data (according to user specification)
- Flagging and handling marked outliers and test data  
- Checking for logical errors (e.g., start date after end date)  
- Subsetting data to defined study period 
- Subsetting data based on user-selected attributes

*Optional* The app then summarizes the data into:
- A figure depicting per-individual duration within the study period with a survival event indicator 
- A population-level figure depicting monthly mortality across the study period. 

Finally, the app performs Kaplan-Meier survival estimation on the cleaned dataset, generating: 
- Life table (summary statistics) 
- A KM survival curve (estimated survival probability over time, with optional 95% confidence intervals)
- A cumulative hazard plot (total accumulated risk over time, with optional 95% confidence intervals)

When group comparison is enabled, the app additionally produces:  
- A table of log-rank test results  
- A comparative KM survival curve plot contrasting the groups (with optional 95% CIs)
- A comparative cumulative hazard plot contrasting the groups (with optional 95% CIs)

Users have the option to zoom in along the y-axis for KM survival curves and cumulative hazard plots to improve figure readability. 

### Application scope

#### Generality of App usability

This app has currently been tested on mammals and birds with N_individuals > 50, N_deaths > 1, and study duration >1 year, but should be applicable to any dataset containing mortality data of sufficient sample size (see below). 

This app allows for staggered entry during the defined study period and for censored data (individuals that are "lost" from the study, e.g., due to equipment failure and individuals that survive the study period). 


#### Required data properties

**Events**: This app can only be used if mortality information (indication of event along with an associated end date) is captured in one of the following columns: death_comments, deployment_end_comments, deployment_end_type, mortality_location_filled. 

The app will produce a warning and terminate if none of the individuals in the study experienced a mortality event during the study period or data subset. 

**Sample size**: The required sample size for a KM survival analysis depends primarily on the total number of deaths, rather than the total number of individuals. In general, the lower the event (death) rate, the higher the number of required individuals. 

This app does *not* perform a power analysis prior to performing the survival analyzes. However, if fewer than 10 mortality events are detected, the app will generate a warning that the model may have low statistical power, potentially resulting in unreliable estimates and poor predictive power. 

Note that a larger sample size is required for comparison across groups. 

**Subsetting data**: Users can select what subset of individuals to include in the study (only females, only individuals alive within a specific survival year, individuals wearing certain collar types, etc.). Data can be subset based on up to two attributes. Available subsetting variables are: sex, life stage, survival year, attachment type, and tag model.

**Comparison types**: Users can compare survival across any categorical grouping variable present in their data by entering the exact column name. The grouping variable is not restricted to a predefined list. Note that group labels are not automatically standardized, so spelling, capitalization, and formatting must be consistent within the data (e.g., "male" and "Male" will be treated as separate groups).

**Survival years**: Users can define a "survival year" for analyses, different from a calendar year in that this period extends from when an animal is typically born to the end of their first year.  

**Life stages**: If the user wants to calculate life stage, the input data must contain the column "animal_birth_hatch_year". Users must also upload auxiliary information (see template: https://github.com/meredithspalmer/MoveApps_Survival/blob/master/animal_birth_hatch_year_table.csv) that links individual animal ages to species-specific life stages. 


### Input type

`move2::move2_loc`

### Output type

`move2::move2_loc`


### Artefacts

**Tracking history**: Figure (`tracking_history.png`) detailing the start and end dates of each individual during the tracking period (with gaps in collaring noted), along with an indicator of how each individual was terminated (i.e., death, censored, survived). (*Optional*) 

**Mortality plot:** Diagnostic plot (`monthly_mortality.png`) depicting mortality rate per month across the study perio (*Optional*). 

**Life table:** Output of KM survival analysis; table (`life_table.csv`) with the time, number of individuals at risk, number of events, survival, standard error, and 95% confidence intervals.

**Annual survival by year**: When survival years are defined, a table (`annual_survival_by_year.csv`) reporting the end-of-year survival probability, standard error, and 95% confidence intervals estimated from a separate KM fit for each survival year, along with the number of individuals at risk.

**KM survival curve:** Plot (`km_survival_curve.png`), depicting survival probability over time, with median survival time indicated and optional 95% confidence intervals.

**Cumulative hazard plot:** Plot (`cumulative_hazard_plot.png`), depicting total accumulated risk (expected number of accumulated deaths) over time, with optional 95% confidence intervals.

**Log-rank test:** Output of comparing survival curves between groups; table (`logrank_table_statistics.csv`) with group-level sample sizes, observed and expected events, O/E ratio, chi-square statistic, degrees of freedom, and p-value.

**Comparison KM curves:** Plot (`km_comparison_curves.png`), depicting survival curves by selected group, with optional 95% confidence intervals.

**Comparison hazard plots:** Plot (`cumulative_hazard_comparison_plot.png`), depicting cumulative hazard by selected group, with optional 95% confidence intervals.

### Settings 

`Start date` and `End date`: Temporal limits of the study period. If left null (default), the analysis will encompass the entire study period. Useful for defining a study year. Unit: `date`. 

`Fix empty start times` and `Fix empty end times`: Defines how the app handles NA dates (`deploy_on_timestamp` and`deploy_off_timestamp`) at the beginning and end of the study. Options include: 
- Replacing with the first/last recorded timestamp (default)
- Replacing with the current date (for end times only)
- Removing missing data 

`Censor capture-related mortality`: If capture-related mortality is a concern, this setting allows users to define a number of days post-capture to exclude from the overall analysis. Default is no censoring. Unit: `days`. 

`Perform analysis on subset of data` and `Define subset condition`: If user wants to perform analyses on a subset of the data, user can define the grouping parameter they want to split the data on (subset condition; e.g., "sex") and the group of interest they wish to retain in the study (subset definition; e.g., "f"). Available subsetting variables are: sex, life stage, survival year, attachment type, and tag model. Default is to include all data. Up to two sequential subsets can be defined. 

`Groups for comparison`: To compare survival across groups, enter the exact column name of the categorical grouping variable of interest. The variable must be present in the data. For survival year or life stage comparisons, survival year start dates and/or auxiliary life stage data must also be provided. Default is no group comparisons.

`Survival year start date`: If comparing across survival years (see above), the user can define the day and month that each 'survival year' begins. The code assumes a year runs 365(6) days. Default is null. Unit: `date`. 

`Animal birth/hatch year definitions`: Optional auxiliary file a user can upload if they are running analyses by survival year and wish to classify individual life stage within a specific year. This file maps animal age to user-defined life stages. A template can be found here: https://github.com/meredithspalmer/MoveApps_Survival/blob/master/animal_birth_hatch_year_table.csv. App expects files in `.csv` format. 

`Life table length`: How often to generate statistical output in the life table (e.g., every # days). Default value is 30 days. Unit: `days`. 

`Plot monthly mortality`: Optionally generate a heatmap showing the number of confirmed mortality events per month across the study period. Default is no plot.

`Plot tracking history`: Optionally generate a figure showing the collaring history of each individual, including deployment start and end dates and gaps in tracking. Default is no plot.

`Add confidence intervals`: Add 95% confidence bands to KM survival and cumulative hazard plots (both overall and comparison plots). Default is no confidence intervals.

`Zoom to plot`: Zoom the y-axis of KM survival and cumulative hazard plots to the range of the data, rather than displaying the full scale. Useful for visualizing subtle differences between groups. Default is no zoom.


### Most common errors

Please document and send errors to mspalmer.zool@gmail.com. 


### Null or error handling

See null handling outlined in *Settings*. 

The app will log warnings quantifying the reductions in sample size due to data cleaning and censoring, as well as indicating how many null timestamps were updated or removed.  

The app will log warnings when data includes fewer than 10 mortality events and will terminate if there are no mortality events. 

Errors may arise due to how mortality events are recorded (what metadata flag mortality and whether mortality is associated with a specific end date). If mortality is logged in a column not mentioned in this documentation, it will not be recognized. 
