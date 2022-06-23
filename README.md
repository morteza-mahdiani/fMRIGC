# Connectivity
_Investigating the connectivity of human brain's recorded signals_

The first version of the code was implemented by [Jon Walbrin](https://orcid.org/0000-0001-9740-4471) and [Morteza Mahdiani](https://morteza-mahdiani.github.io/) at [PROACTION Laboratory](https://proactionlab.fpce.uc.pt/).

## Description
This repository provides an object oriented code for investigating the directed connectivity between signals obtained from human brain fMRI data.

### Dependencies
* [CoSMoMVPA](https://cosmomvpa.org/)
* [The Multivariate Granger Causality (MVGC) Toolbox
](https://www.mathworks.com/matlabcentral/fileexchange/78727-the-multivariate-granger-causality-mvgc-toolbox)
* [Statistics and Machine Learning Toolbox
](https://www.mathworks.com/products/statistics.html)


### Quick Guide
First, you should create an instance of the GCParameters class to set the parameters that are required for granger cuasality analysis obtained by MVGC toolbox as bellow:

```bash
gc_instance = GCParameters(ntrials,nobs, regmode, icregmode, morder, momax, tstat, alpha, mhtc, seed);
```

The parameters are:

- **ntrials**<br>
	-number of trials
- **nobs**
	-number of observations per trial
- **regmode**----------VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
- **icregmode**----------information criteria regression mode ('OLS', 'LWR' or empty for default)
- **morder**-------------model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
- **momax**--------------maximum model order for model order estimation
- **tstat**--------------statistical test for MVGC:  'chi2' for Geweke's chi2 test (default) or'F' for Granger's F-test
- **alpha**--------------significance level for significance test
- **mhtc**---------------multiple hypothesis test correction (see routine 'significance')
- **seed**---------------random seed (0 for unseeded)

For example we can set them like this:

```bash
gcInst = GCParameters(1,1100,'OLS', 'LWR', 'AIC', 20, 'F', 0.05, 'FDR', 0)
```

Secondly we need to create an instance for Connectivity class by calling its constructor as below:

```bash
inst = Connectivity(pathOfData, pathOfMasks, outputPath, fID, lID, ROIs, path_to_MVGC)
```

The parameters are:

- **pathOfData**		% path to the nifti files of fMRI data
- **pathOfMasks**		% path to the masks that we want to apply on fMRI data
- **outputPath**		% path to a directory for storing the outputs
- **fID**     		% the index of the first subject in the nifti files
- **lID** 			% the index of the last subject in the nifti files
- **ROIs** 			% the resions that we want to use for masking input data
- **path_to_MVGC** 	% path to the MVGC toolbox root in our local device

For example we can set them like this:

```bash
inst = Connectivity('/Documents/ResidualTimeCourse_THBFP_FIR/','/Documents/SubjReg_SearchSpaces_GM_ASMasked/','/Documents/out/', 8,10,{'rOFA','rFFA','rSTSF'},'/Applications/MathWorks/MATLAB Add-Ons/Collections/The Multivariate Granger Causality (MVGC) Toolbox' 
```

## License

GNU GENERAL PUBLIC LICENSE
**Free Software**
-------

[https://proactionlab.fpce.uc.pt](https://proactionlab.fpce.uc.pt)