# Bayesian IV Models with Variable Selection in R

![R](https://img.shields.io/badge/R-%25276696C3.svg?style=for-the-badge&logo=r&logoColor=white)
![Bayesian](https://img.shields.io/badge/Bayesian-Research-009687.svg?style=for-the-badge)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=for-the-badge)

## 🌱 Preface

Until I organized the code and published them in my github, it's been a long time. if there are anything wrong please let me know! Although maybe in a long time i could not fix them on time, i will try my best! it is my first try not only writting so many code but also truly realizing how is a research project going, it's a gift for my master graduation and maybe my last time being a student. i will appreciate my first project that can be genuinely regarded as a "work" and never forget this experience! Thanks! 

---

This repository contains the complete code and resources for my Master's research project:  
**"Bayesian Analysis of Linear Models with Endogeneity using Instrumental Variables"**

The project implements, compares, and applies several Bayesian models in **R** to tackle two classic econometric challenges:  
- **Endogeneity bias**  
- **Variable selection in sparse settings**

The proposed **SimplifiedBayesianIV** and **BayesianLassoIV** methods demonstrate superior performance in both simulation studies and a real-world application to nutritional epidemiology (NHANES data).

---

## 🔍 Project Overview

The project implements and rigorously compares four main approaches:

1. **Naive**: Standard Bayesian linear regression that ignores endogeneity (baseline).  
2. **FullBayesianIV**: Full Bayesian IV model that jointly samples all parameters.  
3. **SimplifiedBayesianIV**: A simplified, more computationally efficient Bayesian IV model (**key contribution**).  
4. **BayesianLassoIV**: Extension of the simplified model with Bayesian Lasso priors for automatic variable selection.  

---

## 📦 Repository Structure

```
Bayesian-IV-Variable-Selection-R/
├── README.md                 # Project documentation
├── LICENSE                   # MIT License
├── .gitignore                # Standard R .gitignore
├── Bayesian-IV-Variable-Selection.Rproj     # RStudio Project File
├── R/                                       # Core R function libraries
│   ├── funs_models.R                        # 4 Models fitting functions (Gibbs samplers)
│   ├── funs_sim_endogeneity.R               # Functions to simulate data for endogeneity study
│   ├── funs_sim_variable_selection.R        # Functions to simulate data for variable selection study
│   └── funs_appli_variable_selection.R      # Functions to application for variable selection study (inclduing functions of plots and tables)
├── scripts/                                 # Analysis scripts (can call functions in R/)
│   ├── 01_sim_endogeneity.R                 # Run simulation study for endogeneity
│   ├── 02_sim_variable_selection.R          # Run simulation study for variable selection (based on best model)
│   └── 03_appli_variable_selection.R        # Run application study for variable selection
└── results_report/
    ├── Thesis_Lin_Dongyuan.pdf  # Full thesis manuscript
    └── slide.pdf                # Slides for final presentation in USYD
```

---


## ℹ️ Usage Notes

### Data Availability
No primary data is included in this repository.  

- **Simulated data**: Can be generated using functions in `R/funs_sim_*.R`.  
- **NHANES data**: Must be obtained separately due to restrictions.  
  - Download via: `nhanes` R package  
  - Or from: CDC NHANES Website  

### Code Purpose
The provided code demonstrates:

- Implementation of Bayesian IV models with Gibbs sampling  
- Simulation studies for method validation  
- Basic application framework for real data  

> Note: The application to endogenous settings is straightforward and not elaborated here. For complete results and interpretation, please refer to the thesis document.

---

## 🚀 Getting Started

### Prerequisites
- R (>= 4.1.0)  
- Required packages: `mvtnorm`, `MASS`, `MCMCpack`, `brms`, `doParallel`, `foreach`, `coda`, `ggplot2`, `dplyr`  

### Basic Usage
1. Clone the repository  
2. Open the RStudio project file  
3. Install required packages  
4. Run simulation scripts to generate data and test models  
5. For NHANES application: download data separately and adapt scripts as needed  

Example workflow:

```r
# Fit models using Gibbs sampler
source("R/funs_models.R")
result <- bayesian_iv(data$Y, data$X, data$Z,
                      method = "M3_Lasso", niter = 2000, burn = 500)
```
---

## 👨‍💻 Author
**Dongyuan Lin**  
- Email: *lin_dongyuan@foxmail.com*  
- GitHub: [@dongyuanlin](https://github.com/dongyuanlin)
- LinkedIn: [@dongyuan lin](optional)*(https://www.linkedin.com/in/dongyuan-lin-094237297?utm_source=share&utm_campaign=share_via&utm_content=profile&utm_medium=ios_app) 

---

## 🙏 Acknowledgments
I am deeply grateful to my supervisors, **Associate Professor Clara Grazian** and **Dr. Linh Nghiem**, for their invaluable guidance and support. I had some assistance from smart tools to refine my code. 

This work was completed at **The University of Sydney**.  

---
## 📜 License
This project is licensed under the [MIT License](LICENSE).  
*Note: This license does **not** cover the NHANES data itself.*
