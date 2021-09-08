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

* For the **Estimated sampling variance** part, get into folder `\Addtional\EstimateVar` and open scripts `IKG_EstimateVar.m` and `noise_var.m`.
To consider the sampling variance (1) that is constant and sampling variance (2) that is varying, respectively,
just comment and uncomment the corresponding parts in `noise_var.m`, then run `IKG_EstimateVar.m` to see the comparison for
<img src="https://latex.codecogs.com/svg.latex?{d=1}">
and
<img src="https://latex.codecogs.com/svg.latex?{d=3}">.

## References
* Botev, Zdravko (2020). Truncated Normal Generator. *MATLAB Central File Exchange*.
Retrieved June 29, 2020. https://www.mathworks.com/matlabcentral/fileexchange/53180-truncated-normal-generator
