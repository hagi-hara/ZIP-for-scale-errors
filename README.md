ZIP-for-scale-errors
====


## Overview
This repository includes the analysis codes, graphing codes, and dataset for the following paper:

[Blinded for peer review]
<!-- Hagihara, H., Ishibashi, M., Moriguchi Y., & Shinya Y. (2022). ZLarge-scale data demystify children’s scale errors: A meta-analytic approach using the Zero-Inflated Poisson models. [Manuscript in progress] -->

## Dataset and R scripts
- **01_DescriptiveStats.R**
  - R scripts to analyze the sample characteristics and experimental manipulations.
- **02_SchematicFigures.R**
  - R scripts to create Figure 1 and Supplementary Figure s1.
- **03_Analysis1_InLab.R**
  - R scripts to perform Analysis 1 (In-Lab data).
- **04_Analysis1_Classroom.R**
  - R scripts to perform Analysis 1 (Classroom data).
- **05_Analysis2_TotalVcb.R**
  - R scripts to perform Analysis 2 (Total vocabulary size).
- **06_Analysis2_Noun.R**
  - R scripts to perform Analysis 2 (Noun vocabulary size).
- **07_Analysis2_Verb.R**
  - R scripts to perform Analysis 2 (Verb vocabulary size).
- **08_Analysis2_Adj.R**
  - R scripts to perform Analysis 2 (Adjective vocabulary size).
- **09_SupAnalysis1_SpWords.R**
  - R scripts to perform Supplementary Analysis 1.
- **10_SupAnalysis2_DurationGender.R**
  - R scripts to perform Supplementary Analysis 2.
- **11_SupAnalysis3_DurationInhibitoryControl.R**
  - R scripts to perform Supplementary Analysis 3.
- **12_SupAnalysis4_BSplines.R**
  - R scripts to perform Supplementary Analysis 4.
- **sedata.csv**
  - Dataset used in the study.
- **CDI_data.xlsx**
  - Raw data for vocabulary measures (CDI) used in the study.
- **Models (folder)**
  - Stan codes to perform Bayesian statistics.


## Data Structure
- sedata.csv
- The original studies are as follows (see also the explanation of "groupid"):
  - **A2020**: Arterberry, M. E., Hespos, S. J., Walsh, C. A., & Daniels, C. I. (2020). Integration of thought and action continued: Scale errors and categorization in toddlers. *Infancy*, *25*(6), 851–870. https://doi.org/10.1111/infa.12364
  - **G2019**: Grzyb, B. J., Cangelosi, A., Cattani, A., & Floccia, C. (2019). Children’s scale errors: A by‐product of lexical development? *Developmental Science*, *22*(2), e12741. https://doi.org/10.1111/desc.12741
  - **H2022**: Hagihara, H., Ishibashi, M., Moriguchi, Y., & Shinya, Y. (2022b). Object labeling activates young children’s scale errors at an early stage of verb vocabulary growth. *Journal of Experimental Child Psychology*, *222*, 105471. https://doi.org/10.1016/j.jecp.2022.105471
  - **I2021JP** / **I2021UK**: Ishibashi, M., Twomey, K. E., Westermann, G., & Uehara, I. (2021). Children’s scale errors and object processing: Early evidence for cross-cultural differences. *Infant Behavior and Development*, *65*, 101631. https://doi.org/10.1016/j.infbeh.2021.101631 [JP = Japan sample, UK = UK sample]
  - **IM2017**: Ishibashi, M., & Moriguchi, Y. (2017). Understanding why children commit scale errors: Scale error and its relation to action planning and inhibitory control, and the concept of size. *Frontiers in Psychology*, *8*, 826. https://doi.org/10.3389/fpsyg.2017.00826
  - **IU2020**: Ishibashi, M., & Uehara, I. (2020). The relationship between children’s scale error production and play patterns including pretend play. *Frontiers in Psychology*, *11*, 1776. https://doi.org/10.3389/fpsyg.2020.01776
  - **R2009**: Rosengren, K. S., Carmichael, C., Schein, S. S., Anderson, K. N., & Gutiérrez, I. T. (2009). A method for eliciting scale errors in preschool classrooms. *Infant Behavior and Development*, *32*(3), 286–290. https://doi.org/10.1016/j.infbeh.2009.03.001
  - **R2010**: Rosengren, K. S., Schein, S. S., & Gutiérrez, I. T. (2010). Individual differences in children’s production of scale errors. *Infant Behavior and Development*, *33*(3), 309–313. https://doi.org/10.1016/j.infbeh.2010.03.011
  - **R2010U**: Rosengren, K. S. (n.d.). Variability in young children’s interactions with scale replicas: Exploratory play, general play, pretense, and scale errors. *Unpublished data*.

| Column Name | Variable                | Example   | Explanation                                                                                       |
| ----        | ----                    | ----      | ----                                                                                              |
| paridraw    |qualitative              | NA        | Original participant ID (all values are converted to NA)                                          |    
| parid       |qualitative              | H2022_01  | Participant ID                                                                                    |
| groupid     |quntitative              | H2022     | Dataset ID                                                                                        | 
| gender      |qualitative              | f         | Participants' gender ("f" = female; "m" = male)                                                   |
| age         |quantitative (integer)   | 21        | Participants' age in months                                                                       | 
| se_occ      |quantitative (binary)    | 0         | Whether scale errors were observed or not ("1" = a participant produced at least one scale error) |
| se_sum      |quantitative (integer)   | 2         | The number of scale errors observed                                                               |
| se_dur      |quantitative (continuous)| 55.7      | The duration of scale errors ("inlab" = seconds unite; "classroom" = minutes unite)               |
| country     |quantitative             | JP        | Country where the task was performed ("JP" or "UK")                                               |
| setting     |quantitative             | inlab     | task settings ("inlab" or "classroom")                                                            |
| session     |quantitative (integer)   | 1         | Time when data were collected                                                                     |
| n_obj       |quantitative (integer)   | 3         | Nuber of target objects                                                                           |
| dur         |quantitative (integer)   | 5         | Duration of the task period (min)                                                                 |
| date        |quantitative (date)      | 2009/3/9  | Date when observations were performed                                                             |
| ec          |quantitative (continuous)| 1         | Effortful control score                                                                           |
| vcball      |quantitative (integer)   | 259       | Participants' total vocabulary size                                                               |
| noun        |quantitative (integer)   | 201       | Participants' noun vocabulary size                                                                |
| verb        |quantitative (integer)   | 18        | Participants' verb vocabulary size                                                                |
| adj         |quantitative (integer)   | 2         | Participants' adjective vocabulary size                                                           |
| noun_sp     |quantitative (continuous)| 0.8       | Participants' noun lexical score (see Supplementary Analysis 1)                                   |
| verb_sp     |quantitative (continuous)| 0.75      | Participants' verb lexical score (see Supplementary Analysis 1)                                   |
| adj_sp      |quantitative (continuous)| 0.5       | Participants' adjective lexical score (see Supplementary Analysis 1)                              |
| sounds1     |quantitative (integer)   | 2         | Participants' sub-vocabulary size (same for the columns between animals and others)               |


## Software & Package Versions
- R: 4.2.2
- cmdstanr: 0.5.3
- ggmcmc: 1.5.1.1
- here: 1.0.1
- loo: 2.5.1
- rstan: 2.26.13
- tidyverse: 2.0.0
- dplyr: 1.1.1
- splines2: 0.5.0
- 
  
## Authors of This Repository
[Blinded for peer review]

<!-- If you have any questions, please email at **hiromichi.h(AT)gmail.com** (please replace **(AT)** with **@**).	- [Hiromichi Hagihara](https://github.com/hagi-hara) -->
