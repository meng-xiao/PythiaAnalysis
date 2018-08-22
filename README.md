same usage as ZZAnalysis
To process PYTHIA level particles and do ZZ selection


Run in CMSSW_9_4_4

wget -O ${TMPDIR}/checkout.csh https://raw.githubusercontent.com/meng-xiao/Test/master/checkout.csh

cd $CMSSW_BASE/src

cmsenv

chmod u+x ${TMPDIR}/checkout.csh

${TMPDIR}/checkout.csh
