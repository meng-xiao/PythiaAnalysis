#same usage as ZZAnalysis

#To process PYTHIA level particles and do ZZ selection


#Run in CMSSW_9_4_4

#On lxplus, do the following

cmsrel CMSSW_9_4_4
cd CMSSW_9_4_4/src
cmsenv
wget -O ${TMPDIR}/checkout.csh https://raw.githubusercontent.com/meng-xiao/PythiaAnalysis/master/checkout.csh

cd $CMSSW_BASE/src

chmod u+x ${TMPDIR}/checkout.csh

${TMPDIR}/checkout.csh

#Once the above step is done, compile the package

cd $CMSSW_BASE/src
scram b -j8

cd MiniAnalyzer/test
source grid.sh
batch_Condor.py samples_2017_MC.csv -i analyzer_gen.py -o minloFilter

