echo "Processing step 1 (GEN SIM DIGI RAW) MC samples to L1Ntuples with the timing bit, done with TDC threshold = 3*18.7 for investigations"
echo " "
cd ../../
scram b -j 8
cd HcalTrigger/Validation/
echo " "
echo "QCD"
cmsRun ntuple_maker_def_QCD_3xthreshold.py
mv L1Ntuple.root L1Ntuple_QCD_3xthreshold.root
echo " "
echo "Neutrino Gun"
cmsRun ntuple_maker_def_nugun.py
mv L1Ntuple.root L1Ntuple_NeutrinoGun.root
echo " "
echo "LLP mh=125 GeV, ctau=10 m"
cmsRun ntuple_maker_def_mh125_pl10000_3xthreshold.py
mv L1Ntuple.root L1Ntuple_mh125_mx50_pl10000_3xthreshold.root
echo " "
echo "LLP mh=125 GeV, ctau=1 m"
cmsRun ntuple_maker_def_mh125_pl1000_3xthreshold.py
mv L1Ntuple.root L1Ntuple_mh125_mx50_pl1000_3xthreshold.root
echo " "
echo "LLP mh=125 GeV, ctau=0.5 m"
cmsRun ntuple_maker_def_mh125_pl500_3xthreshold.py
mv L1Ntuple.root L1Ntuple_mh125_mx50_pl500_3xthreshold.root
echo " "
echo "LLP mh=1000 GeV, ctau=10 m"
cmsRun ntuple_maker_def_mh1000_pl10000_3xthreshold.py
mv L1Ntuple.root L1Ntuple_mh1000_mx450_pl10000_3xthreshold.root
echo " "
echo "LLP mh=1000 GeV, ctau=1 m"
cmsRun ntuple_maker_def_mh1000_pl1000_3xthreshold.py
mv L1Ntuple.root L1Ntuple_mh1000_mx450_pl1000_3xthreshold.root
echo " "
echo "LLP mh=1000 GeV, ctau=0.5 m"
cmsRun ntuple_maker_def_mh1000_pl500_3xthreshold.py
mv L1Ntuple.root L1Ntuple_mh1000_mx450_pl500_3xthreshold.root
echo " "
echo "LLP mh=350 GeV, ctau=10 m"
cmsRun ntuple_maker_def_mh350_pl10000_3xthreshold.py
mv L1Ntuple.root L1Ntuple_mh350_mx160_pl10000_3xthreshold.root
echo " "
echo "LLP mh=350 GeV, ctau=1 m"
cmsRun ntuple_maker_def_mh350_pl1000_3xthreshold.py
mv L1Ntuple.root L1Ntuple_mh350_mx160_pl1000_3xthreshold.root
echo " "
echo "LLP mh=350 GeV, ctau=0.5 m"
cmsRun ntuple_maker_def_mh350_pl500_3xthreshold.py
mv L1Ntuple.root L1Ntuple_mh350_mx160_pl500_3xthreshold.root
echo " "
echo "LLP mh=250 GeV, ctau=10 m"
cmsRun ntuple_maker_def_mh250_pl10000_3xthreshold.py
mv L1Ntuple.root L1Ntuple_mh250_mx120_pl10000_3xthreshold.root
echo " "
echo "LLP mh=250 GeV, ctau=1 m"
cmsRun ntuple_maker_def_mh250_pl1000_3xthreshold.py
mv L1Ntuple.root L1Ntuple_mh250_mx120_pl1000_3xthreshold.root
echo " "
echo "LLP mh=250 GeV, ctau=0.5 m"
cmsRun ntuple_maker_def_mh250_pl500_3xthreshold.py
mv L1Ntuple.root L1Ntuple_mh250_mx120_pl500_3xthreshold.root
