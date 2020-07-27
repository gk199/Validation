echo "Compiling rates.cxx, draw_timingbit.cxx files"
echo " "
scram b -j 8
echo " "
echo " "
echo "QCD"
rates.exe new /eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/TimingTrigger/QCD_Pt-15to3000_TuneCP5_Flat_14TeV_pythia8_HCAL/CRAB3_QCD_TDC_MC/200722_192212/0000/QCD
mv rates_new_cond.root rates_new_cond_QCD.root
echo " "
echo "Neutrino gun new and old conditions"
rates.exe new /eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/TimingTrigger/Neutrino_Pt-2to20_gun_HCAL/CRAB3_NuGun_TDC_MC/200722_142743/0000/
rates_original.exe def /eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/TimingTrigger/Neutrino_Pt-2to20_gun_HCAL/CRAB3_NuGun_TDC_MC/200722_142743/0000/
echo " "
echo "Rates plots from draw_rates.exe"
draw_rates.exe
mv plots/*Rates_emu.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/
echo " "
echo "LLP mh=125 GeV, ctau=3 m"
rates.exe new /eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/TimingTrigger/HTo2LongLivedTo4b_MH-125_MFF-50_CTau-3000mm_TuneCP5_14TeV_pythia8_HCAL/CRAB3_mh125_ctau3000_TDC_MC/200722_185526/0000
mv rates_new_cond.root rates_new_cond_LLP_pl3000.root
echo " "
echo "Making overlay plots from draw_timingbit.exe"
draw_timingbit_mh125.exe
mv plots/ADC50_3ns_4JetMultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/ADC50_3ns_no1HE4HB_4JetMultHBOverlay_mh125.pdf
mv plots/ADC50_3ns_4JetMultHEOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/ADC50_3ns_no1HE4HB_4JetMultHEOverlay_mh125.pdf
mv plots/ADC50_3ns_4JetMultHBHEOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/ADC50_3ns_no1HE4HB_4JetMultHBHEOverlay_mh125.pdf
echo " "
echo "LLP mh=250 GeV, ctau=0.5 m"
rates.exe new /eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/TimingTrigger/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-500mm_TuneCP5_14TeV_pythia8_HCAL/CRAB3_mh250_ctau500_TDC_MC/200722_193005/0000
mv rates_new_cond.root rates_new_cond_LLP_pl500.root
echo " "
echo "LLP mh=250 GeV, ctau=1 m"
rates.exe new /eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/TimingTrigger/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/CRAB3_mh250_ctau1000_TDC_MC/200722_185009/0000
mv rates_new_cond.root rates_new_cond_LLP_pl1000.root
echo " "
echo "LLP mh=250 GeV, ctau=10 m"
rates.exe new /eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/TimingTrigger/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-10000mm_TuneCP5_14TeV_pythia8_HCAL/CRAB3_mh250_ctau10000_TDC_MC/200722_190316/0000
mv rates_new_cond.root rates_new_cond_LLP_pl10000.root
echo " "
echo "Making overlay plots from draw_timingbit.exe"
draw_timingbit_mh250.exe
mv plots/ADC50_3ns_4JetMultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/ADC50_3ns_no1HE4HB_4JetMultHBOverlay_mh250.pdf
mv plots/ADC50_3ns_4JetMultHEOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/ADC50_3ns_no1HE4HB_4JetMultHEOverlay_mh250.pdf
mv plots/ADC50_3ns_4JetMultHBHEOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/ADC50_3ns_no1HE4HB_4JetMultHBHEOverlay_mh250.pdf
echo " "
echo " "
echo "LLP mh=350 GeV, ctau=1 m"
rates.exe new /eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/TimingTrigger/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/CRAB3_mh350_ctau1000_TDC_MC/200722_182856/0000
mv rates_new_cond.root rates_new_cond_LLP_pl1000.root
echo " "
echo "LLP mh=350 GeV, ctau=10 m"
rates.exe new /eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/TimingTrigger/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-10000mm_TuneCP5_14TeV_pythia8_HCAL/CRAB3_mh350_ctau10000_TDC_MC/200722_190619/0000
mv rates_new_cond.root rates_new_cond_LLP_pl10000.root
echo " "
echo "Making overlay plots from draw_timingbit.exe"
draw_timingbit_mh350.exe
mv plots/ADC50_3ns_4JetMultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/ADC50_3ns_no1HE4HB_4JetMultHBOverlay_mh350.pdf
mv plots/ADC50_3ns_4JetMultHEOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/ADC50_3ns_no1HE4HB_4JetMultHEOverlay_mh350.pdf
mv plots/ADC50_3ns_4JetMultHBHEOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/ADC50_3ns_no1HE4HB_4JetMultHBHEOverlay_mh350.pdf
echo " "
echo "LLP mh=1000 GeV, ctau=10 m"
rates.exe new /eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/TimingTrigger/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-10000mm_TuneCP5_14TeV_pythia8_HCAL/CRAB3_mh1000_ctau10000_TDC_MC/200722_191012/0000
mv rates_new_cond.root rates_new_cond_LLP_pl10000.root
echo " "
echo "LLP mh=1000 GeV, ctau=100 m"
rates.exe new /eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/TimingTrigger/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-100000mm_TuneCP5_14TeV_pythia8_HCAL/CRAB3_mh1000_ctau100000_TDC_MC/200722_191315/0000
mv rates_new_cond.root rates_new_cond_LLP_pl100000.root
echo " "
echo "Making overlay plots from draw_timingbit.exe"
draw_timingbit_mh1000.exe
mv plots/ADC50_3ns_4JetMultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/ADC50_3ns_no1HE4HB_4JetMultHBOverlay_mh1000.pdf
mv plots/ADC50_3ns_4JetMultHEOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/ADC50_3ns_no1HE4HB_4JetMultHEOverlay_mh1000.pdf
mv plots/ADC50_3ns_4JetMultHBHEOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/ADC50_3ns_no1HE4HB_4JetMultHBHEOverlay_mh1000.pdf
echo " "
echo " "
echo "Efficiency vs. rate plots"
RateEfficiencyPlots_LLP_ctau.exe
mv /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/EffRate* /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/
