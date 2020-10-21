# How to run the code

> version: CMSSW_10_2_18  
> description: To process PYTHIA level particles and do ZZ selection  

## Configure enviriment

On lxplus, do the following: 

```bash
cmsrel CMSSW_10_2_18
cd CMSSW_10_2_18/src
cmsenv
wget -O ${TMPDIR}/checkout.csh https://raw.githubusercontent.com/meng-xiao/PythiaAnalysis/master/checkout.csh
# wget -O ${TMPDIR}/checkout.csh https://raw.githubusercontent.com/qwezarty/PythiaAnalysis/master/checkout.csh
cd $CMSSW_BASE/src
chmod u+x ${TMPDIR}/checkout.csh
${TMPDIR}/checkout.csh
```

## Compile packages

```bash
cd $CMSSW_BASE/src
scram b -j8
```

## Run at clusters (Grid)

You need a permission / cert for reaching grid clusters. See [Getting a Personal Grid Certificate from CERN](https://uscms.org/uscms_at_work/physics/computing/getstarted/get_grid_cert.shtml) for more details.

```bash
# setup grid credentials
cd MiniAnalyzer/test
source grid.sh
# create jobs and submit
../scripts/batch_Condor.py samples_2017_MC.csv -i analyzer.py -o minloFilter
cd minloFilter && resubmit_Condor.csh
# to check how many jobs are still running or pending on Condor
condor_q
# once the jobs are done, from the same folder
checkProd.csh
```

See more details at [Submitting Jobs](https://github.com/CJLST/ZZAnalysis/wiki/SubmittingJobs).

## Dependencies
- [ZZAnalysis](https://github.com/CJLST/ZZAnalysis)
- [HiggsAnalysis-ZZMatrixElement](https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement)
- [MelaAnalytics](https://github.com/usarica/MelaAnalytics)
