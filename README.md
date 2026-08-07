# Tweedie GLMM Practical: MAIHDA Approach in R

[![R](https://img.shields.io/badge/R-%E2%89%A5%_4.1-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)
[![glmmTMB](https://img.shields.io/badge/Package-glmmTMB-blue)](#)
[![synthpop](https://img.shields.io/badge/Privacy-synthpop-success)](#)

A practical masterclass script for modeling zero-inflated, right-skewed continuous data using **Tweedie Generalized Linear Mixed Models (GLMMs)** within a **MAIHDA** (Multilevel Analysis of Individual Heterogeneity and Discriminatory Accuracy) framework.

## 🎯 The Core Problem Solved
Standard OLS and log(y+1) transformations fail when analyzing data with a massive spike at zero combined with a long right tail (e.g., weekly alcohol units consumed). 

While traditional two-part (Hurdle) models require two sets of coefficients and complex covariance terms that frequently fail to converge in a MAIHDA multilevel context, the **Tweedie distribution** handles this seamlessly. It provides ONE modeling part, ONE set of random effects, and highly interpretable Variance Partition Coefficients (VPC).

## 📊 About the Data
To comply with safe data requirements and preserve privacy, this repository uses **synthetic data** generated via the `synthpop` package. It learns the joint distribution of variables from the real Scottish Health Survey (SHeS 2017) and generates a new dataset (`SHeS_Masterclass_Synthetic.csv`) that preserves the statistical structure without containing any real individual records.

## 🚀 Usage (Live Coding Script)

The practical is divided into two parts:
* **PART 1 (Synthetic Data Generation):** Shows how the synthetic dataset was generated using CART methods to perfectly mirror the real SHeS data.
* **PART 2 (Student Live Coding):** The main tutorial script walking through:
  * Data exploration and identifying zero-inflation
  * Building the MAIHDA intersectional strata
  * Fitting the Null MAIHDA model and Full MAIHDA model using `glmmTMB`
  * Extracting Variance Partition Coefficients (VPC) and Proportional Change in Variance (PCV)
  * Checking GLM assumptions using `DHARMa`

### Dependencies
You will need the following R packages installed before running the scripts:
```R
install.packages(c("dplyr", "glmmTMB", "easystats", "synthpop", "haven", "DHARMa"))
