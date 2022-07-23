ZIP-for-scale-errors
====


## Overview
This repository includes the analysis codes, graphing codes, and dataset for the following paper:

[Blinded for peer review]
<!-- Hagihara, H., Ishibashi, M., Moriguchi Y., & Shinya Y. (2022). Zero-inflated Poisson models for developmental change in young children’s scale errors. [Manuscript in progress] -->

## Dataset and R scripts
- **1_DescriptiveStats.R**
  - R scripts to analyze the sample characteristics and experimental manipulations.
- **2_SchematicFigures.R**
  - R scripts to create Figures 1 and 2.
- **3_ZipPolynomial.R**
  - R scripts to perform Study 1 (Models 1-4).
- **4_ZipBspline.R**
  - R scripts to perform Study 1 (Model 5).
- **5_Poisson.R**
  - R scripts to perform Study 1 (Models 6-8).
- **6_ZipPolynomial_vcb.R**
  - R scripts to perform Study 2 (Models 1-4 and 6-9).
- **7_ZipBspline_vcb.R**
  - R scripts to perform Study 2 (Models 5 and 10).
- **sedata.csv**
  - Dataset used in the study.
- **Models (folder)**
  - Stan codes to perform Bayesian statistics.


## Data Structure
- sedata.csv
- The original studies are as follows (see also the explanation of "groupid"):
  - **H2022F** / **H2022S**: Hagihara, H., Ishibashi, M., Moriguchi, Y., & Shinya, Y. (2022b). Object labeling activates young children’s scale errors at an early stage of verb vocabulary growth. *Journal of Experimental Child Psychology*, *222*, 105471. https://doi.org/10.1016/j.jecp.2022.105471 [F = 1st session, S = 2nd session]
  - **I2021JP** / **I2021UK**: Ishibashi, M., Twomey, K. E., Westermann, G., & Uehara, I. (2021). Children’s scale errors and object processing: Early evidence for cross-cultural differences. *Infant Behavior and Development*, *65*, 101631. https://doi.org/10.1016/j.infbeh.2021.101631 [JP = Japan sample, UK = UK sample]
  - **IM2017**: Ishibashi, M., & Moriguchi, Y. (2017). Understanding why children commit scale errors: Scale error and its relation to action planning and inhibitory control, and the concept of size. *Frontiers in Psychology*, *8*, 826. https://doi.org/10.3389/fpsyg.2017.00826
  - **IU2020**: Ishibashi, M., & Uehara, I. (2020). The relationship between children’s scale error production and play patterns including pretend play. *Frontiers in Psychology*, *11*, 1776. https://doi.org/10.3389/fpsyg.2020.01776

| Column Name | Variable             | Example   | Explanation                                                                                       |
| ----        | ----                 | ----      | ----                                                                                              |
| paridraw    |qualitative           | NA        | Original participant ID (all values are converted to NA)                                          |    
| parid       |qualitative           | H2022F_01 | Participant ID                                                                                    |
| groupid     |quntitative           | H2022F    | Dataset ID                                                                                        | 
| gender      |qualitative           | f         | Participants' gender ("f" = female; "m" = male)                                                   |
| age         |quntitative (integer) | 21        | Participants' age in months                                                                       | 
| se_occ      |quantitative (binary) | 0         | Whether scale errors were observed or not ("1" = a participant produced at least one scale error) |
| se_sum      |quantitative (integer)| 2         | The number of scale errors observed                                                               |
| country     |quantitative          | JP        | Country where the task was performed ("JP" or "UK")                                               |
| session     |qualitative           | first     | Time when data were collected (i.e., 1st or 2nd session)                                          |
| vcb         |quntitative (integer) | 259       | Participants' total vocabulary size                                                               |


## Software & Package Versions
- RStudio: 1.4.1106
- R: 4.04
- cmdstanr: 0.3.0
- ggmcmc: 1.5.1.1
- here: 1.0.1
- loo: 2.4.1
- rstan: 2.21.2
- tidyverse: 1.3.0

  
## Authors of This Repository
[Blinded for peer review]

<!-- If you have any questions, please email at **hiromichi.h(AT)gmail.com** (please replace **(AT)** with **@**).	- [Hiromichi Hagihara](https://github.com/hagi-hara) -->
