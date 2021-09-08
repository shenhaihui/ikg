# Source Codes of Integrated Knowledge Gradient

This repository contains the source codes used in the following paper:
* Ding, Liang, L. Jeff Hong, Haihui Shen<sup>*</sup>, Xiaowei Zhang (2021). Technical Note--Knowledge gradient for selection with covariates: Consistency and computation. Accepted in *Naval Research Logistics*.

### Disclaimer
LICENSE: Redistribution and use in source and binary forms, with or without modification, are permitted, under the terms of the [BSD license](./BSD_License.txt).

WARNING: These codes are written only for the purpose of demonstration and verification. While the correctness has been carefully checked, the quality such as standardability,
clarity, generality, and efficiency has not been well considered.

## 1 Introduction
The `\IKG2021` folder contains the MATLAB implementations of the numerical experiments in Ding et al. (2021).
* `\IKG2021\Main`, contains codes for numerical experiments in Ding et al. (2021, §5.1).
The codes for the other compared polices are included in the folder `\IKG2021\Main\OtherPolicy`.
* The function `\IKG2021\Main\trandn.m` is a third-party function from Botev (2020).
* `\IKG2021\Addtional`, contains codes for numerical experiments in Ding et al. (2021, Appendix I).

## 2 Installation
The codes were written and run in MATLAB R2018b, on Windows 7 Enterprise 64-bit OS,
with Intel i9-9900K CPU @ 3.60 GHz, 16 GB RAM.

To install the MATLAB codes, just copy the entire folder `\R&S-C2020` into your MATLAB directory, or change the path of MATLAB to the folder `\R&S-C2020`

## 3 Details on Numerical Experiments and Case Study
### 3.1 Numerical Experiments in Shen et al. (2020, §6, EC.6)
Get into folder `\NumericalExperiments`.
* To run Procedure TS and TS+ for all the 13 problems with pre-specified target
<img src="https://latex.codecogs.com/svg.latex?{\text{PCS}_{\text{E}}\geq{1-\alpha}}">, run script `PCSE.m`, which is self-explanatory.
 Scripts `find_h_hom_PCSE.m` and `find_h_hom_PCSE_SA.m` are functions to calculates the constant <img src="https://latex.codecogs.com/svg.latex?{h}"> of Procedure TS for homoscedastic errors, where the the former uses numerical integration (`integral` or `trapz`) and MATLAB built-in root finding function `fzero`, while the latter uses stochastic approximation method.
 Similarly, scripts `find_h_het_PCSE.m` and `find_h_het_PCSE_SA.m` are functions to calculates the constant <img src="https://latex.codecogs.com/svg.latex?{h_{\text{Het}}}"> of Procedure TS+ for heteroscedastic errors.
 
* To run Procedure TS and TS+ for all the 13 problems with pre-specified target
<img src="https://latex.codecogs.com/svg.latex?{\text{PCS}_{\text{min}}\geq{1-\alpha}}">, run script `PCSmin.m`, which is self-explanatory.
 Scripts `find_h_hom_PCSmin.m` and `find_h_het_PCSmin.m` are functions to calculates the constant <img src="https://latex.codecogs.com/svg.latex?{h}"> of Procedure TS and constant <img src="https://latex.codecogs.com/svg.latex?{h_{\text{Het}}}"> of Procedure TS+ , respectively.

### 3.2 Robustness Test in Shen et al. (2020, §7.2)
Get into folder `\RobustnessTest`, and run script `RobustnessTest.m`, which is self-explanatory.

### 3.3 Case Study in Shen et al. (2020, §8)
Get into folder `\EsophagealCancer`.
* The script `EsophagealCancerSim.m` is the implementation of a discrete-time Markov chain model that simulates the progression of Barrett’s esophagus (BE) to esophageal
adenocarcinoma (EAC), which was developed by Hur et al. (2004) and Choi et al. (2014). See more details [here](https://simopt.github.io/ECSim).

* The script `BruteForceSim.m` runs the simulation model for 10^6 replications on a grid of points to approximate the true response surfaces, whose results are saved in `QALY_true.mat`.

* To select personalized treatment regimen with Procedure TS+ and compare the performance with traditional approach, run script `PersonalizedTreatment.m`, which
is self-explanatory.

## References
* Choi, Sung Eun, Katherine E. Perzan, Angela C. Tramontano, Chung Yin Kong, and Chin Hur (2014). Statins and aspirin for chemoprevention in Barrett’s esophagus: Results
of a cost-effectiveness analysis. *Cancer Prevention Research* **7**(3), 341–350.

* Hur, Chin, Norman S. Nishioka, and G. Scott Gazelle (2004). Cost-effectiveness of aspirin chemoprevention for Barrett’s esophagus. *Journal of the National Cancer Institute* **96**(4), 316–325.
