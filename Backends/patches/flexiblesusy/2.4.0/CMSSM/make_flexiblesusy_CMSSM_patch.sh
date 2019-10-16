#!/usr/bin/env bash
cd ../../../../downloaded/
tar -xvf FlexibleSUSY-2.4.0.tar.gz
mkdir FlexibleSUSY-2.4.0/models
mkdir FlexibleSUSY-2.4.0/models/CMSSM
cp -r ../../Models/data/FlexibleSUSY/CMSSM/* ./FlexibleSUSY-2.4.0/models/CMSSM/
echo "diff -rupN FlexibleSUSY-2.4.0/models/CMSSM/CMSSM_slha_io.hpp ../installed/flexiblesusy/2.4.0/CMSSM/models/CMSSM/CMSSM_slha_io.hpp" >> patch_flexiblesusy_2.4.0_CMSSM.dif
diff -rupN FlexibleSUSY-2.4.0/models/CMSSM/CMSSM_slha_io.hpp ../installed/flexiblesusy/2.4.0/CMSSM/models/CMSSM/CMSSM_slha_io.hpp >> patch_flexiblesusy_2.4.0_CMSSM.dif
mv patch_flexiblesusy_2.4.0_CMSSM.dif ../patches/flexiblesusy/2.4.0/CMSSM
rm -r FlexibleSUSY-2.4.0
cd ../patches/flexiblesusy/2.4.0/CMSSM
