#!/bin/tcsh -fe

# MELA Analytics
git clone https://github.com/usarica/MelaAnalytics.git
(cd MelaAnalytics; git checkout -b from-v19 v1.9)
# MELA
git clone https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement
(cd ZZMatrixElement; git checkout -b from-v222 v2.2.2)
# JHUGenMELA
git clone https://github.com/JHUGen/JHUGenMELA.git
(cd JHUGenMELA; git checkout -b from-v231 v2.3.1)
# Test
git clone https://github.com/meng-xiao/PythiaAnalysis.git Test
# git clone https://github.com/qwezarty/PythiaAnalysis.git Test

# compile JHUGenMELA
(					\
  cd ${CMSSW_BASE}/src/JHUGenMELA	;\
  ./setup.sh -j 8			;\
)


# replace ZZMatrixElement/MELA/setup.sh -j 8
(                                                                 \
  cd ${CMSSW_BASE}/src/ZZMatrixElement/MELA/COLLIER/             ;\
  set pkgname="collier-1.2.0"                                      ;\
  set pkgdir="COLLIER-1.2"                                       ;\
  set tarname=$pkgname".tar.gz"                                  ;\
  set tarweb="https://www.hepforge.org/archive/collier/"$tarname ;\
  set libname="libcollier.so"                                    ;\
  set tmpdir="colliertmp"                                        ;\
  wget $tarweb                                                   ;\
  mkdir $tmpdir                                                  ;\
  tar -xvzf $tarname -C $tmpdir                                  ;\
  rm $tarname                                                    ;\
  mv $tmpdir"/"$pkgdir"/src/"* ./                                ;\
  rm -rf $tmpdir                                                 ;\
  make                                                           ;\
  mv $libname "../data/"$SCRAM_ARCH"/"$libname                   ;\
)
(                                                                 \
  cd ${CMSSW_BASE}/src/ZZMatrixElement/MELA/fortran/             ;\
  make all                                                       ;\
  mv libjhugenmela.so ../data/${SCRAM_ARCH}/                     ;\
)

