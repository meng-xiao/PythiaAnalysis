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

## Run in clusters

```bash
cd MiniAnalyzer/test
source grid.sh
batch_Condor.py samples_2017_MC.csv -i analyzer_gen.py -o minloFilter
```
