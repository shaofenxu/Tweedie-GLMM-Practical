# ==============================================================================
# Tweedie GLMM (using MAIHDA) Masterclass - R Scripts
# Author: Shaofen Xu
# Description: Two-part script for PGT Methods Masterclass
#   PART 1 (for Shaofen only): Generate synthetic SHeS data
#   PART 2 (for students): Live coding script
# ==============================================================================




# ==============================================================================
# PART 1: SYNTHETIC DATA GENERATION
# Prepare data for students
# Students should only use Part 2 for practice
# ==============================================================================

# --- 1.1 Install / load required packages ------------------------------------

library(synthpop)   # Privacy-preserving synthetic data generation
library(dplyr)
library(haven)
library(glmmTMB)   
library(easystats)

# --- 1.2 Load and prepare real data -------------------------------------

shes_17i_archive_v1 <- read_dta("shes_17i_archive_v1.dta")

age_breaks <- c(15, 24, 34, 44, 54, 64, 74, Inf)
age_labels <- c("16-24", "25-34", "35-44", "45-54", "55-64", "65-74", "75+")

real_data <- shes_17i_archive_v1 %>%
  filter(drating >= 0, age >= 16) %>%
  mutate(
    drating    = as.numeric(drating),
    age_group  = cut(age, breaks = age_breaks, labels = age_labels, right = TRUE),
    Sex        = haven::as_factor(Sex),
    SIMD16_RPa = factor(as.numeric(SIMD16_RPa),
                        levels = 1:5,
                        labels = c("1 (least deprived)", "2", "3", "4", "5 (most deprived)"))
  ) %>%
  select(drating, Sex, age_group, SIMD16_RPa) %>%
  filter(complete.cases(.))

# --- 1.3 Fit a quick model to firstly check real p (psi or ψ) and phi (φ) -------------------------

m_check <- glmmTMB(
  drating ~ Sex + age_group + SIMD16_RPa + (1 | Sex:age_group:SIMD16_RPa), # MAIHDA approach: level 1 -- Individual/observation, level 2 -- intersectional strata
  data   = real_data,
  family = tweedie()  # Tweedie log-link function
)

# Extract variance power (p or ψ), back-transforming from internal psi parameter
psi_est      <- m_check$fit$par["psi"]
variance_power <- plogis(psi_est) + 1   # Should be between 1 and 2
cat("Real data variance power (p):", round(variance_power, 3), "\n")

# Extract dispersion parameter (φ) - stored in log scale
phi_est <- sigma(m_check)
cat("Real data dispersion (phi):  ", round(phi_est, 3), "\n")

# Note the zero proportion in real data for later comparison
zero_pct_real <- mean(real_data$drating == 0)
cat("Zero proportion in real data:", round(zero_pct_real * 100, 1), "%\n")

# --- 1.4 Generate synthetic data using synthpop ------------------------------
# synthpop learns the joint distribution of all variables from the real data
# and generates a new dataset that preserves that structure without
# containing any real individual records. This satisfies safe data requirements.

set.seed(20260501)  # Reproducible; change if preferring a different seed

syn_result <- syn(
  real_data %>% select(Sex, age_group, SIMD16_RPa, drating),
  method = c("sample", "cart", "cart", "cart"),  # CART for all four variables
  m      = 1,       # Generate one synthetic dataset
  k      = 4000     # Target sample size for teaching (manageable for students)
)

synthetic_data <- syn_result$syn

# --- 1.5 Validate: compare key distributional features -----------------------

cat("\n--- Validation: Zero proportion ---\n")
cat("Real data:      ", round(mean(real_data$drating == 0)          * 100, 1), "%\n")
cat("Synthetic data: ", round(mean(synthetic_data$drating == 0)     * 100, 1), "%\n")

cat("\n--- Validation: Mean among drinkers ---\n")
cat("Real data:      ", round(mean(real_data$drating[real_data$drating > 0]), 2), "units\n")
cat("Synthetic data: ", round(mean(synthetic_data$drating[synthetic_data$drating > 0]), 2), "units\n")

cat("\n--- Validation: Variable distributions ---\n")
cat("Sex distribution (real):\n");      print(prop.table(table(real_data$Sex)))
cat("Sex distribution (synthetic):\n"); print(prop.table(table(synthetic_data$Sex)))

# --- 1.6 Final tidy and export -----------------------------------------------

teaching_data <- synthetic_data %>%
  mutate(
    drating    = round(as.numeric(drating), 1),
    Sex        = factor(Sex),
    age_group  = factor(age_group,  levels = age_labels),
    SIMD16_RPa = factor(SIMD16_RPa,
                        levels = c("1 (least deprived)", "2", "3",
                                   "4", "5 (most deprived)"))
  )

write.csv(teaching_data, "SHeS_Masterclass_Synthetic.csv", row.names = FALSE)

cat("\nSynthetic dataset saved: SHeS_Masterclass_Synthetic.csv\n")
cat("Rows:", nrow(teaching_data), " | Strata structure: 2 x 7 x 5 = 70\n")







# ==============================================================================
# PART 2: LIVE CODING SCRIPT
# ==============================================================================

# ==============================================================================
# Tweedie MAIHDA Masterclass - Student Live Coding Script
#
# Data: "SHeS_Masterclass_Synthetic.csv" (Distributed to students by Shaofen)
#       (Synthetic dataset based on Scottish Health Survey 2017,
#        generated to preserve distributional features.)
#
# Packages needed - please install before the session:
#   install.packages(c("dplyr", "glmmTMB", "easystats"))
# ==============================================================================

# --- Load libraries -----------------------------------------------------------
library(dplyr)
library(glmmTMB)     # for constructing generalised linear models
library(easystats)   # Loads parameters, performance, insight, etc.

# --- 0. Load data -------------------------------------------------------------
shes_data <- read.csv(file.choose())     #Select "SHeS_Masterclass_Synthetic.csv" where you stored it
head(shes_data)
str(shes_data)

# --- 1. Explore the data: the core problem ------------------------------------

# Look at the distribution of weekly alcohol units
# You will see a large spike at zero AND a long right tail.
# This is why standard OLS and log(y+1) transformations fail.

hist(shes_data$drating,
     breaks = 50,
     main   = "Weekly Alcohol Unit Distribution (Synthetic SHeS 2017)",
     xlab   = "Units per week",
     col    = "#3182bd",
     border = "white")

# What proportion are exactly zero?
cat("Zero proportion:", round(mean(shes_data$drating == 0) * 100, 1), "%\n")

# This is continuous (not integer counts), so ZIP / ZINB fails
cat("Non-integer values among drinkers:",
    sum(shes_data$drating[shes_data$drating > 0] !=
        floor(shes_data$drating[shes_data$drating > 0])), "\n")

# --- 2. Data preparation ------------------------------------------------------
shes_data <- shes_data %>%
  mutate(
    drating    = as.numeric(drating),
    Sex        = factor(Sex),
    age_group  = factor(age_group,
                        levels = c("16-24","25-34","35-44",
                                   "45-54","55-64","65-74","75+")),
    SIMD16_RPa = factor(SIMD16_RPa,
                        levels = c("1 (least deprived)", "2", "3",
                                   "4", "5 (most deprived)"))
  )

# Build the MAIHDA intersectional strata variable (2 x 7 x 5 = 70 strata)
shes_data$strata <- factor(paste0(shes_data$Sex, ", ",
                                   shes_data$age_group, ", ",
                                   shes_data$SIMD16_RPa))

cat("Number of strata:", nlevels(shes_data$strata), "\n")  # Should be 70

# --- 3. Model 1: Null MAIHDA model -------------------------------------------
# No fixed effects controlled.
#The random intercept variance captures total inequalities between strata
#
# Key: family = tweedie()
#   - Handles the spike at zero (via Poisson component)
#   - Handles the right-skewed continuous tail (via Gamma component)
#   - Log-link ensures predicted means are always positive
#   - Automatically estimates variance power p (psi or ψ) (should be 1 < p < 2)

cat("Fitting Null Model...\n")

m_null <- glmmTMB(
  drating ~ 1 + (1 | strata),
  data   = shes_data,
  family = tweedie()
)

# View fixed and random effect summary
model_parameters(m_null)

# The Variance Partition Coefficient (VPC / ICC):
# What proportion of total variance (inequalities) in alcohol units sits between strata?
# A VPC > 5% is evidence that intersectional inequality (discriminatory accuracy) matters.
# We will compare this VPC with Full Model's VPC later on.
icc(m_null)

# --- 4. Model 2: Full MAIHDA model -----------------------------------
# Control Sex, age_group, SIMD as additive fixed effects.
#The random intercept variance now captures interaction-associated inequalities
#(among sex, age and deprivation) only between strata

cat("Fitting Full Model...\n")

m_main <- glmmTMB(
  drating ~ Sex + age_group + SIMD16_RPa + (1 | strata),
  data   = shes_data,
  family = tweedie()
)

# How much variance remains between strata after accounting for additive effects for Null Model?
# Therefore, Full Model's VPC is evidence of interaction effects (also, > 5% means it matters).
cat("\n--- VPC comparison ---\n")
cat("Null Model VPC:       ", round(icc(m_null)$ICC_adjusted, 3), "\n")
cat("Full Model VPC:     ", round(icc(m_main)$ICC_adjusted, 3), "\n")

# Proportional Change in Variance (PCV): how much did the additive main effects explain?
v_null <- get_variance(m_null)
v_main <- get_variance(m_main)
pcv    <- (v_null$var.random - v_main$var.random) / v_null$var.random
cat("PCV (Null Model --> Full Model):   ", round(pcv * 100, 1), "%\n")

# --- 5. Extract and inspect the variance power (p or ψ or psi) ---------------------------
# p is estimated on logit-scale automatically by glmmTMB using "qlogis(psi - 1)".
# Backtransform using "plogis(psi) + 1" to confirm it sits between 1 and 2.
# Note: Different packages estimate p differently, check package manuals for backtransforming.

psi            <- m_null$fit$par["psi"]
variance_power <- plogis(psi) + 1
cat("Null Model's variance power p:", round(variance_power, 5), "\n")

psi            <- m_main$fit$par["psi"]
variance_power <- plogis(psi) + 1
cat("Full Model's variance power p:", round(variance_power, 5), "\n")
# Interpretation: 1 < p < 2 confirms valid Tweedie (compound Poisson-Gamma) structure.

# --- 6. Interpret coefficients on the original (response) scale --------------
# Because we used log-link, log-scale coefficients or intercepts must be exponentiated to Rate Ratios or Means.
# Rate Ratio: how much the expected alcohol units are multiplied
#             for a one-unit / one-level change in that predictor.

model_parameters(m_main, exponentiate = TRUE)
compare_parameters(m_null, m_main, exponentiate = TRUE)

# Example interpretation based on m_main (Full Model) output:
# 
# Intercept is the 8.29 alcohol units, which is baseline mean
#
# Rate Ratio > 1: Look at 'Sex [Male]' = 2.20
# Males are expected to consume 2.20 times the alcohol units per week 
# (i.e., a 120% increase) compared to females, holding age and SIMD constant.
#
# Rate Ratio < 1: Look at 'age group [75+]' = 0.49
# Individuals aged 75 and over are expected to consume 51% fewer alcohol units
# (1 - 0.49 = 0.51) compared to the baseline age group (16-24), holding others constant.


# --- 7. Why NOT a two-part (Hurdle) model here? -----------------
# In a standard Hurdle model you would need:
#   Step 1: Logistic regression on whether drating is 0  -->  one set of coefficients
#   Step 2: Gamma/Log-Normal regression when drating > 0  -->  another set
#
# In MAIHDA you would then need TWO sets of coefficients for both Null Model and Full Model,
# plus a covariance term between them. Models of this complexity frequently
# fail to converge, and computing interpretable VPC and PCV becomes intractable.
#
# Tweedie solves all of this with ONE modelling part, ONE set of random effects,
# and ONE interpretable set of VPC and PCV. This is the core benefit.


#BONUS: You can always use "DHARMa" package for checking GLM assumptions: in our case, they are uniform distribution, linearity, and homoscedasticity.
#       Even that standard linear regression assumes normality, this can be checked via uniform distribution, because standard linear regression belongs to GLMs.
#       However, for multilevel models (GLMMs) like MAIHDA, other assumptions (e.g., random intercept normality) need to be separately checked.

library(DHARMa)
simulateResiduals(m_null, plot = TRUE)
simulateResiduals(m_main, plot = TRUE)

#Note for students: With a large sample size, "DHARMa" tests are highly sensitive. Focus on overall visual checks rather than red warnings.

# ==============================================================================
# END OF STUDENT SCRIPT
# For questions after the session, contact: x.xu.2@research.gla.ac.uk
# ==============================================================================
