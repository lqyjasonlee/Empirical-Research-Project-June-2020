# Documentation for ERP Data and Code Files
Empirical Research Project  
Qiyuan Li 2019 SMU Ph.D. student
***

This repository contains MATLAB and Stata code and dataset that replicates the results from the simulation and empirical results. Please ensure that the code and the called files are in the same directory.

## Files for simulation illustrations
### **ERP_simulation.m**
the MATLAB code used for generating results in FIGURE 1 and FIGURE 2.
- In line 17:
    - set `design = 1` for exact sparsity case in FIGURE 1;
    - set `design = 22` for approximate sparsity case in FIGURE 2.
- The results are stored in **simulation_result.mat** and **simulation_result_sparse.mat**.
- This function calls the following files:
    - **Heteroskedastic_se.m**
    - **LassoShooting.m**
    - **MC_TE_Design_New.m**
    - **MC_TE_FixedDesign_Heteroskedastic_Lasso_RedForm.m**
    - **MC_TE_GetCoef_RedForm.m**
    - **MC_TE_GetSupport.m**
    - **MC_TE_LassoHeteroskedastic_SplitSample.m**
    - **MC_TE_LassoHeteroskedastic_unpenalized.m**
    - **MC_TE_PostEstimator.m**
    - **MC_TE_SimulateLambda.m**
    - **process_options.m**

## Files for empirical results
### **ERP_data_qiyuan.dta** 
the dataset used for generating results in TABLE 1, use `summarize(varname)`.
### **ERP_code.do**
the Stata code for generating results in TABLE 2.
- The results are stored in **erp_result_new.txt** in the same directory.
- This file calls the following files:
    - **ERP_data_qiyuan.dta**
    - **lassoShooting.ado**
