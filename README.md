# LAMBDA

(Data & Codes for:)<br>
**Evolution of virulence in emerging epidemics: from theory to experimental evolution and back**

Wakinyan Benhamou, François Blanquart, Marc Choisy, Thomas W. Berngruber, Rémi Choquet and Sylvain Gandon

## Data (source and short descriptions)

This is previous experimental time series from [[Berngruber *et al.*, 2013](https://doi.org/10.1371/journal.ppat.1003209)] or simulated data from this study. All data used in the scripts are in the 'Data' folder and, within, simulated datasets are in the 'Simulated_data' subfolder.

The evolution experiment designed in [[Berngruber *et al.*, 2013](https://doi.org/10.1371/journal.ppat.1003209)] monitored the competition between two strains of phage $\lambda$ with distinct life-history strategies in continuous cultures of *Escherichia coli*. The first strain is the wildtype which is known to have a relatively large lysogenization rate and low reactivation rate. The second strain is the $\lambda$ cI857 variant, which is known to be more virulent and transmitted mostly horizontally through lytic cycles. There are 8 chemostats: 4 replicates x 2 treatments (epidemic, with an initial prevalence around 1% *vs.* endemic, with an initial prevalence around 99%). Viruses are first introduced as prophages with initial ratio 1:1.

### data_FACS.csv

&#11169;&emsp; Experimental times series of the (logit-)prevalence and viral strain (logit-)frequencies among infected cells.

### data_qPCR.csv

&#11169;&emsp; Experimental times series of the viral strain (logit-)frequencies in the culture medium.

### Experimental_data_mean_and_figures.RData

&#11169;&emsp; This R Data file (generated by **Experimental_data_mean_and_figures.R**, in the root directory of this repositery) stores:
- The R dataframes corresponding to **data_FACS.csv** and **data_qPCR.csv**, respectively;
- The R dataframes containing mean values accross chemostats for each treatment;
- Plots of experimental times series (logit-prevalence, logit-frequency infected by virulent phage and logit-frequency free virulent phage over time).

### Simulated_data/

&#11169;&emsp; Simulated datasets, as generated by **Simulate_Data.R** (root directory of this repositery), for 8 chemostats: 4 replicates x 2 treatments (epidemic *vs.* endemic). The file **parms_simul.csv** stores the parameter values used for the simulations. Files whose name contains '*Simul_data_FACS*' simulate times series of the (logit-)prevalence and viral strain (logit-)frequencies among infected cells; files whose name contains '*Simul_data_qPCR*' simulate times series of the viral strain (logit-)frequencies in the free virus stage. Eventually, the meaning of the code number (between 1 and 4) just before the extension '*.csv*' is as follows:
|                                        | Samplings every 0.1 h | Samplings every 1 h |
| :---                                   |         :---:         |         :---:       |
| Almost no measurement errors (SD=0.01) |           1           |          3          |
| Greater measurement errors (SD=0.5)    |           2           |          4          |



## R codes

We developed an inference approach to estimate the parameters of the model. We implemented a two-step MLE (Maximum Likelihood Estimation) procedure: (i) we first obtained point estimates of the rates of prophage reactivation $\alpha_w$ and $\alpha_m$ of both viral strains, (ii) then we fixed $\alpha_w$ and $\alpha_m$ to their point estimates and we ran non-linear optimizations to infer the remaining parameters of the model (except the burst size $B$ which had to be fixed because of an identifiability issue).

### Inference_rates_prophage_reactivation/

Scripts:
- **STEP1_functions.R** (script with the custom functions used in the first step)
- **STEP1_Script_UK_data.Rmd** (analyses using the UK data)
- **STEP1_Identifiability_profiles.Rmd** (build and plot identifiability profiles of the model)

### Non_linear_optimizations/

Scripts:
- **STEP2_functions.R** (script with the custom functions used in the second step)
- **STEP2_Script_England_data.Rmd** (analyses - including the mixed-effects models - using the data from England)

### LogLikelihood_landscape_b_B/

---

In addition, in the root directory:

- **Selection_coefficient_vs_Stringency_Index.R** explores the correlation between the selection coefficient of the Alpha variant and the Stringency Index in the UK; the script generates Fig. S2
- **Fig_2step_analysis_with_real_data.R** generates Fig. 1

---

## Outputs

The 'Outputs_csv' folder hosts a part of what was generated by the scripts (csv and rds files). Most of them were generated after too long running times not to be stored - e.g. non-linear optimizations - and/or are reused elsewhere - e.g. a file generated in the first step and imported in the second step.

Non-linear optimizations

- Initial values: 
**Initial_values_v4_gamma01_kappa02_pS09.csv**

- Parameter estimates: 
**Estim_tab_v4_gamma01_kappa02_pS09_R025.csv**

- Best estimates from the previous file: 
**Best_estimates_phase1_v4_gamma01_kappa02_pS09_R025.csv**

Wild bootstraps

- Using Mammen's 2-points distribution:

     (estimates) **bootstrap_wild_v4_Mammen2_gamma01_kappa02_pS09_R025.csv**

     (simulations) **simul_obs_bootstrap_wild_v4_Mammen2_gamma01_kappa02_pS09_R025.csv**

- Using Rademacher distribution:

     (estimates) **bootstrap_wild_v4_Rademacher_gamma01_kappa02_pS09_R025.csv**

     (simulations) **simul_obs_bootstrap_wild_v4_Rademacher_gamma01_kappa02_pS09_R025.csv**

Identifiability profiles

- Parameter E(t0)/N: 
**Identifiability_profile_pE_pS09_gamma01_kappa02_R025_norm_noise_sd5000_500_5.csv**
- Parameter omega: 
**Identifiability_profile_omega_pS09_gamma01_kappa02_R025_norm_noise_sd5000_500_5.csv**
- Parameter p: 
**Identifiability_profile_p_pS09_gamma01_kappa02_R025_norm_noise_sd5000_500_5.csv**
- Parameter alpha: 
**Identifiability_profile_alpha_pS09_gamma01_kappa02_R025_norm_noise_sd5000_500_5.csv**
- Parameter k: 
**Identifiability_profile_k_pS09_gamma01_kappa02_R025_norm_noise_sd5000_500_5.csv**
- Parameter a: 
**Identifiability_profile_a_pS09_gamma01_kappa02_R025_norm_noise_sd5000_500_5.csv**
- Parameter eta: 
**Identifiability_profile_eta_pS09_gamma01_kappa02_R025_norm_noise_sd5000_500_5.csv**
- Parameter mu: 
**Identifiability_profile_mu_pS09_gamma01_kappa02_R025_norm_noise_sd5000_500_5.csv**

Effects of small variations in the fixed parameters (phase 1)

- Varying gamma: 
**Optim_v4_vary_gamma_nstarts500.csv**
- Varying kappa: 
**Optim_v4_vary_kappa_nstarts500.csv**
- Varying R0 (basic reproduction number): 
**Optim_v4_vary_R0_nstarts500.csv**
- Varying S(t0)/N: 
**Optim_v4_vary_SN0_nstarts500.csv**
- Overview: 
**Estimates_with_pertubed_parameters_phase1.csv**
- Figure (R file): 
**Fig_optim_var_step1.rds**
