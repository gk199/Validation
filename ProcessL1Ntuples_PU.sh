echo "Processing step 1 (GEN SIM DIGI RAW) MC samples to L1Ntuples with the timing bit"
echo " "
cd ../../
scram b -j 8
cd HcalTrigger/Validation/
echo " "
echo "QCD"
cmsRun ntuple_maker_def_QCD.py
mv L1Ntuple.root L1Ntuple_QCD_PU.root
echo " "
echo "Neutrino Gun"
cmsRun ntuple_maker_def_NeutrinoGun.py
mv L1Ntuple.root L1Ntuple_NeutrinoGun_PU.root
echo " "
echo "LLP mh=125 GeV, ctau=30 m"
cmsRun ntuple_maker_def_mh125_pl30000.py
mv L1Ntuple.root L1Ntuple_mh125_mx50_pl30000_PU.root
echo " "
echo "LLP mh=125 GeV, ctau=3 m"
cmsRun ntuple_maker_def_mh125_pl3000.py
mv L1Ntuple.root L1Ntuple_mh125_mx50_pl3000_PU.root
echo " "
echo "LLP mh=250 GeV, ctau=10 m"
cmsRun ntuple_maker_def_mh250_pl10000.py
mv L1Ntuple.root L1Ntuple_mh250_mx120_pl10000_PU.root
echo " "
echo "LLP mh=250 GeV, ctau=1 m"
cmsRun ntuple_maker_def_mh250_pl1000.py
mv L1Ntuple.root L1Ntuple_mh250_mx120_pl1000_PU.root
echo " "
echo "LLP mh=250 GeV, ctau=0.5 m"
cmsRun ntuple_maker_def_mh250_pl500.py
mv L1Ntuple.root L1Ntuple_mh250_mx120_pl500_PU.root
echo " "
echo "LLP mh=350 GeV, ctau=10 m"
cmsRun ntuple_maker_def_mh350_pl10000.py
mv L1Ntuple.root L1Ntuple_mh350_mx160_pl10000_PU.root
echo " "
echo "LLP mh=350 GeV, ctau=1 m"
cmsRun ntuple_maker_def_mh350_pl1000.py
mv L1Ntuple.root L1Ntuple_mh350_mx160_pl1000_PU.root
echo " "
echo "LLP mh=350 GeV, ctau=0.5 m"
cmsRun ntuple_maker_def_mh350_pl500.py
mv L1Ntuple.root L1Ntuple_mh350_mx160_pl500_PU.root
echo " "
echo "LLP mh=1000 GeV, ctau=100 m"
cmsRun ntuple_maker_def_mh1000_pl100000.py
mv L1Ntuple.root L1Ntuple_mh1000_mx450_pl100000_PU.root
echo " "
echo "LLP mh=1000 GeV, ctau=10 m"
cmsRun ntuple_maker_def_mh1000_pl10000.py
mv L1Ntuple.root L1Ntuple_mh1000_mx450_pl10000_PU.root
