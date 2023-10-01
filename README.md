# hidden-controller-framework
Codes and data for the paper "Data-Driven Framework for Uncovering Hidden Control Strategies in Evolutionary Analysis"

## Notes: 
In the paper, two sets of observed COVID-19 data were utilized:
  ```
  (i) Data from five prefectures in Japan, provided by JX PRESS Corporation, and
  (ii) Data from nine countries, sourced from the "worldwide epidemiological database for COVID-19" (Guidotti, 2022) and the "COVID-19 data hub" (Guidotti & Ardia, 2020).
  ```

According to the contract with JX PRESS Corporation, the first dataset (i) cannot be made public. 
Therefore, this GitHub repository includes only the second dataset (ii), along with the two codes, HiddenControllerFramework_for_data_ii.m and Heatmap_for_data_ii.r, designated for this dataset.

## **HiddenControllerFramework_for_data_ii.m** (Matlab):
  Estimate the parameters of A, B, Q, and R in our hidden controller framework.
  ```
  INPUT: data from nine countries: 
    AUS.csv		CHL.csv		CZE.csv		JPN.csv		ZAF.csv
    BRA.csv		COL.csv		DEU.csv		LTU.csv   
  OUTPUT: 
    parameters of A, B, Q, and R
  ```

## **Heatmap_for_data_ii.r** (R):
  Generate a heatmap illustrating the z-scores derived from Generalized Additive Models ("mgcv" package) for predictor covariates. 
  These covariates should be associated with the number of individuals who are infected, deceased, and recovered across the nine countries under study.
  ```
  INPUT: 
    output_interior-point_original_i1_d1_{AUS, BRA, CHL, COL, CZE, DEU, JPN, LTU, ZAF}/*.csv  
  OUTPUT:
    heatmap (Figure 8) and dendrogram (Figure 9)
  ```
    
  
