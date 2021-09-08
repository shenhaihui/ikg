# Source Codes of Integrated Knowledge Gradient

This repository contains the source codes used in the following paper:
* Ding, Liang, L. Jeff Hong, Haihui Shen<sup>*</sup>, Xiaowei Zhang (2021). Technical Note--Knowledge gradient for selection with covariates: Consistency and computation. Accepted in *Naval Research Logistics*.

### Disclaimer
LICENSE: Redistribution and use in source and binary forms, with or without modification, are permitted, under the terms of the [BSD license](./BSD_License.txt).

WARNING: These codes are written only for the purpose of demonstration and verification. While the correctness has been carefully checked, the quality such as standardability,
clarity, generality, and efficiency has not been well considered.

## 1 Introduction
The `\IKG2021` folder contains the MATLAB implementations of the numerical experiments in Ding et al. (2021).
* `\IKG2021\Main`, contains codes for numerical experiments in Ding et al. (2021, §5).
The codes for the other compared polices are included in the folder `\IKG2021\Main\OtherPolicy`.
* The function `\IKG2021\Main\trandn.m` is a third-party function from Botev (2020).
* `\IKG2021\Addtional`, contains codes for numerical experiments in Ding et al. (2021, Appendix I).

## 2 Installation
The codes were written and run in MATLAB R2018b, on Windows 7 Enterprise 64-bit OS,
with Intel i9-9900K CPU @ 3.60 GHz, 16 GB RAM.

To install the MATLAB codes, just copy the entire folder `\IKG2021` into your MATLAB directory, or change the path of MATLAB to the folder `\IKG2021`

## 3 Details on Numerical Experiments
### 3.1 Numerical Experiments in Ding et al. (2021, §5.1)
Get into folder `\Main`. Run the first and second sections of script `IKG_Experiments.m`. 
* Set parameter `density_type=1` for P1, and set parameter `density_type=2` for P2.
* The compared policies, *IKGwRC*, *BSE*, *PRS* are implemented in folder `\OtherPolicy`.
* The BSE policy needs to tune a parameter
<img src="https://latex.codecogs.com/svg.latex?{M}">,
which is conducted in script `TuneMforBSE.m`,
and intermediate results are saved in `BSE_regret_reduction.mat`. 

### 3.2 Numerical Experiments in Ding et al. (2021, §5.2)
Run the third section of script `IKG_Experiments.m` for P3.

### 3.3 Numerical Experiments in Ding et al. (2021, Appendix I)
Get into folder `\Addtional`.
* For the **Computational cost comparison** part, get into folder `\Addtional\CompuationCost` and open script `CompareWithSAA.m`, which is self-explanatory. This part involves manually and roughly tuning the sample size in step (i) of method 2, for each
<img src="https://latex.codecogs.com/svg.latex?{d=1,\ldots,7}">.
To do that, first fix a value of
<img src="https://latex.codecogs.com/svg.latex?{d}">,
then run the first and second sections with different value of `sample_size_ratio`, as explained there.
To directly check the computation time with our tuned values of `sample_size_ratio`, run the third section of `CompareWithSAA.m` to record the time and run the fourth section to plot.
Note that different computer configuration may require different time, but the comparison should be similar.
* 
* 
* The script `EsophagealCancerSim.m` is the implementation of a discrete-time Markov chain model that simulates the progression of Barrett’s esophagus (BE) to esophageal
adenocarcinoma (EAC), which was developed by Hur et al. (2004) and Choi et al. (2014). See more details [here](https://simopt.github.io/ECSim).

* The script `BruteForceSim.m` runs the simulation model for 10^6 replications on a grid of points to approximate the true response surfaces, whose results are saved in `QALY_true.mat`.

* To select personalized treatment regimen with Procedure TS+ and compare the performance with traditional approach, run script `PersonalizedTreatment.m`, which
is self-explanatory.

## References
* Choi, Sung Eun, Katherine E. Perzan, Angela C. Tramontano, Chung Yin Kong, and Chin Hur (2014). Statins and aspirin for chemoprevention in Barrett’s esophagus: Results
of a cost-effectiveness analysis. *Cancer Prevention Research* **7**(3), 341–350.

* Hur, Chin, Norman S. Nishioka, and G. Scott Gazelle (2004). Cost-effectiveness of aspirin chemoprevention for Barrett’s esophagus. *Journal of the National Cancer Institute* **96**(4), 316–325.
