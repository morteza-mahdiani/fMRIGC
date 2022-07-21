# Connectivity
_Investigating the connectivity of human brain's recorded signals_

The first version of the code was implemented by [Morteza Mahdiani](https://morteza-mahdiani.github.io/) and [Jon Walbrin](https://orcid.org/0000-0001-9740-4471) at [PROACTION Laboratory](https://proactionlab.fpce.uc.pt/).

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

- **ntrials**
	- number of trials
- **nobs**
	- number of observations per trial
- **regmode**
	- VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
- **icregmode**
	-	information criteria regression mode ('OLS', 'LWR' or empty for default)
- **morder**
	- model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
- **momax**
	- maximum model order for model order estimation
- **tstat**
	- statistical test for MVGC:  'chi2' for Geweke's chi2 test (default) or'F' for Granger's F-test
- **alpha**
 	- significance level for significance test
- **mhtc**
	- multiple hypothesis test correction (see routine 'significance')
- **seed**
	- random seed (0 for unseeded)

For example we can set them like this:

```bash
gc_instance = GCParameters(1,1100,'OLS', 'LWR', 'AIC', 20, 'F', 0.05, 'FDR', 0)
```

Secondly we need to create an instance for Connectivity class by calling its constructor as below:

```bash
c_instance = Connectivity(pathOfData, pathOfMasks, outputPath, fID, lID, ROIs, path_to_MVGC)
```

The parameters are:

- **pathOfData**
	- path to the nifti files of fMRI data
- **pathOfMasks**
	- path to the masks that we want to apply on fMRI data
- **outputPath**
	- path to a directory for storing the outputs
- **fID** 
	- the index of the first subject in the nifti files
- **lID**
	- the index of the last subject in the nifti files
- **ROIs**
	- the resions that we want to use for masking input data
- **path_to_MVGC**
	- path to the MVGC toolbox root in our local device

For example we can set them like this:

```bash
c_instance = Connectivity('/Documents/ResidualTimeCourse_THBFP_FIR/','/Documents/SubjReg_SearchSpaces_GM_ASMasked/','/Documents/out/', 8,10,{'rOFA','rFFA','rSTSF'},'/Applications/MathWorks/MATLAB Add-Ons/Collections/The Multivariate Granger Causality (MVGC) Toolbox' 
```

Finally you should use preprocess and GCM methods for the initialized instance of the Connectivity class! For preprocess function you should pass the number of trials as follow:

```bash
c_instance.preprocess(number_of_trials)
```

Use a number like:

```bash
c_instance.preprocess(1100)
```

and for GCM you should pass the instance of GCParameters class, actual order of the model and the path to the preprocessed data obtained from preprocess function.

```bash
c_instance.GCM(GCParameters_initialized_instance, actual_model_order, path_of_the_preprocessed_data)
```

 For example, you can use it as follow:

```bash
c_instance.GCM(gc_instance, 20, '/Documents/out/')
```

The 3D output matrix of Granger Causality analysis(region by region by number of subjects) after Z conversion will be stored in the 'GCMOutput' folder located in the output directory. Also, you will find the preprocessed data there and you should use them as when calling GCM function. You can visualize it with the following function:

```bash
c_instance.visualize(obj, path_to_data, mode, subjectID)
```

The parameters are:

- **path_to_data**
	- path to the 3D output matrix of connectivity function
- **mode**
	- the mode can be 'm' or 'sn' that represents whether you need to have the mean of all subjects results or just an specific subject result. If you set mode as 'sn', you should then provide the subject number. 
- **subjectID**
	- the ID of the subject. Set it when the mode is 'sn'.

For example, set parameters as below:

```bash
c_instance.visualize('/Users/Documents/out/GCMOutput/GC3DMat.mat','m') 
```

Or if you want to visualize the results of a specific subject, provide its ID like:

```bash
c_instance.visualize('/Users/Documents/out/GCMOutput/GC3DMat.mat','sn', 2) 
```

## To do

- I'm not sure it's very feasible to always write out .mat time course files during 'preprocessing'; for example, in future analyses, we'll want to run voxel-wise comparisons for 200+ parcels (i.e. nchoosek(200,2) = ~20000 mat files per subject! So I think it is important that we have a separate function that outputs the timecourse data ready for another function that runs MVGC wihtout needing to store the data

- After generating .mat time course files, re-running the GC script will default to taking all available .mat files, rather than those specified when re-running (e.g. if I first ran it with 20 subjects, then later specify only for 3 subjects, it runs for all 20 rather than the 3 I specify)

- My feeling is that perhaps this would be easiest to implement as a series of stand alone functions rather than as a single script (e.g. like the cosmo MVPA toolbox, and most other matlab toolboxes where each function is a separate .m file); I think the logic would certainly be a little more intuitive to MATLAB users and it's a little easier to inspect the outputs of each function rather than having everything wrapped up in this way

- GC3DMat - last dimension of 3D mat is greater than n subjects (perhaps because i ran it with 2 subjects missing (i.e. subj 4 and 7); is it better to pre-allocate the size of this, or include some other output (e.g. subject names) that allows for easy understanding of which matrix corresponds to a particular subject

- Most of the the z values are negative, try norminv(1-pval) to make the Z scores positive

- We're getting quite a few infinity values (due to pvalues of 0) which are inconveneiently outputted from the mvgc toolbox; can you see if there is a way of returning non-zero values for these (when using mvgc toolbox), so that we can calculate Z-scores from them? 

- can you remind me which axis of the output matrix is which (e.g. are vertical labels describing the causal influence of regions in the lower triangle of the matrix (and horizontal labels are causal for the upper triangle)?

- Perhaps a more intuitive naming of each function (e.g. connectivity is quite vague) 

## License

GNU GENERAL PUBLIC LICENSE
**Free Software**
-------

[https://proactionlab.fpce.uc.pt](https://proactionlab.fpce.uc.pt)
