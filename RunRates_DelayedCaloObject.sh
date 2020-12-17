echo "Compiling rates.cxx, draw_timingbit.cxx files"
echo " "
scram b -j 8
echo " "
echo " "
echo "LLP mh=125 GeV, ctau=3 m"
rates_delayed_cluster.exe new /eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/TimingTrigger/HTo2LongLivedTo4b_MH-125_MFF-50_CTau-3000mm_TuneCP5_14TeV_pythia8_HCAL/CRAB3_mh125_ctau3000_TDC_MC/200722_185526/0000 $1 $2 $3 $4 $5 $6
mv rates_new_cond.root rates_new_cond_LLP_pl3000.root
echo " "
echo "LLP mh=125 GeV, ctau=30 m"
rates_delayed_cluster.exe new /eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/TimingTrigger/HTo2LongLivedTo4b_MH-125_MFF-50_CTau-30000mm_TuneCP5_14TeV_pythia8_HCAL/CRAB3_mh125_ctau30000_TDC_MC/200728_205042/0000 $1 $2 $3 $4 $5 $6
mv rates_new_cond.root rates_new_cond_LLP_pl30000.root
echo " "
echo "QCD"
rates_delayed_cluster.exe new /eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/TimingTrigger/QCD_Pt-15to3000_TuneCP5_Flat_14TeV_pythia8_HCAL/CRAB3_QCD_TDC_MC/200722_192212/0000/QCDsmall $1 $2 $3 $4 $5 $6
mv rates_new_cond.root rates_new_cond_QCD.root
echo " "
echo "Neutrino gun new and old conditions"
rates_delayed_cluster.exe new /eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/TimingTrigger/Neutrino_Pt-2to20_gun_HCAL/CRAB3_NuGun_TDC_MC/200722_142743/0000/nugun $1 $2 $3 $4 $5 $6
rates_original.exe def /eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/TimingTrigger/Neutrino_Pt-2to20_gun_HCAL/CRAB3_NuGun_TDC_MC/200722_142743/0000/nugun
echo " "
echo " "
echo "Rates plots from draw_rates.exe"
draw_rates.exe
mv plots/*Rates_emu.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/DelayedCaloObject/
mv rates_new_cond.root rates_new_cond_110X_nugun.root
echo " "
echo "Making overlay plots from draw_timingbit.exe"
draw_timingbit_delayed_cluster_mh125.exe
mv rates_new_cond_LLP_pl30000.root rates_new_cond_LLP_mh125_pl30000.root
mv rates_new_cond_LLP_pl3000.root rates_new_cond_LLP_mh125_pl3000.root
mv plots/Delayed_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/DelayedCaloObject/Delayed_2x2_MultHBOverlay_mh125.pdf
mv plots/Prompt_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/DelayedCaloObject/Prompt_2x2_MultHBOverlay_mh125.pdf
mv plots/Prompt_Energy_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/DelayedCaloObject/Prompt_Energy_2x2_MultHBOverlay_mh125.pdf
mv plots/Mult_delayed_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/DelayedCaloObject/Mult_delayed_hitOverlay_mh125.pdf
mv plots/Mult_prompt_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/DelayedCaloObject/Mult_prompt_hitOverlay_mh125.pdf
mv plots/HTdistributionOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/DelayedCaloObject/HTdistributionOverlay_mh125.pdf
mv plots/HTdistribution_trigOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/DelayedCaloObject/HTdistribution_trig_Overlay_mh125.pdf
echo " "
echo "LLP mh=1000 GeV, ctau=10 m"
rates_delayed_cluster.exe new /eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/TimingTrigger/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-10000mm_TuneCP5_14TeV_pythia8_HCAL/CRAB3_mh1000_ctau10000_TDC_MC/200722_191012/0000 $1 $2 $3 $4 $5 $6
mv rates_new_cond.root rates_new_cond_LLP_pl10000.root
echo " "
echo "LLP mh=1000 GeV, ctau=100 m"
rates_delayed_cluster.exe new /eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/TimingTrigger/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-100000mm_TuneCP5_14TeV_pythia8_HCAL/CRAB3_mh1000_ctau100000_TDC_MC/200722_191315/0000 $1 $2 $3 $4 $5 $6
mv rates_new_cond.root rates_new_cond_LLP_pl100000.root
echo " "
echo "Making overlay plots from draw_timingbit.exe"
draw_timingbit_delayed_cluster_mh1000.exe
mv rates_new_cond_LLP_pl1000.root rates_new_cond_LLP_mh1000_pl1000.root
mv rates_new_cond_LLP_pl10000.root rates_new_cond_LLP_mh1000_pl10000.root

mv plots/Delayed_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/DelayedCaloObject/Delayed_2x2_MultHBOverlay_mh1000.pdf
mv plots/Prompt_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/DelayedCaloObject/Prompt_2x2_MultHBOverlay_mh1000.pdf
mv plots/Prompt_Energy_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/DelayedCaloObject/Prompt_Energy_2x2_MultHBOverlay_mh1000.pdf
mv plots/Mult_delayed_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/DelayedCaloObject/Mult_delayed_hitOverlay_mh1000.pdf
mv plots/Mult_prompt_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/DelayedCaloObject/Mult_prompt_hitOverlay_mh1000.pdf
mv plots/HTdistributionOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/DelayedCaloObject/HTdistributionOverlay_mh1000.pdf
mv plots/HTdistribution_trigOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/DelayedCaloObject/HTdistribution_trig_Overlay_mh1000.pdf
echo " "
echo "LLP mh=350 GeV, ctau=10 m"
rates_delayed_cluster.exe new /eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/TimingTrigger/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-10000mm_TuneCP5_14TeV_pythia8_HCAL/CRAB3_mh350_ctau10000_TDC_MC/200722_190619/0000 $1 $2 $3 $4 $5 $6
mv rates_new_cond.root rates_new_cond_LLP_pl10000.root
echo " "
echo "LLP mh=350 GeV, ctau=1 m"
rates_delayed_cluster.exe new /eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/TimingTrigger/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/CRAB3_mh350_ctau1000_TDC_MC/200722_182856/0000 $1 $2 $3 $4 $5 $6
mv rates_new_cond.root rates_new_cond_LLP_pl1000.root
echo " "
echo "LLP mh=350 GeV, ctau=0.5 m"
rates_delayed_cluster.exe new /eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/TimingTrigger/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-500mm_TuneCP5_14TeV_pythia8_HCAL/CRAB3_mh350_ctau500_TDC_MC/200728_204829/0000 $1 $2 $3 $4 $5 $6
mv rates_new_cond.root rates_new_cond_LLP_pl500.root
echo " "
echo "Making overlay plots from draw_timingbit.exe"
draw_timingbit_delayed_cluster_mh350.exe
mv rates_new_cond_LLP_pl500.root rates_new_cond_LLP_mh350_pl500.root
mv rates_new_cond_LLP_pl1000.root rates_new_cond_LLP_mh350_pl1000.root
mv rates_new_cond_LLP_pl10000.root rates_new_cond_LLP_mh350_pl10000.root
mv plots/Delayed_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/DelayedCaloObject/Delayed_2x2_MultHBOverlay_mh350.pdf
mv plots/Prompt_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/DelayedCaloObject/Prompt_2x2_MultHBOverlay_mh350.pdf
mv plots/Prompt_Energy_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/DelayedCaloObject/Prompt_Energy_2x2_MultHBOverlay_mh350.pdf
mv plots/Mult_delayed_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/DelayedCaloObject/Mult_delayed_hitOverlay_mh350.pdf
mv plots/Mult_prompt_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/DelayedCaloObject/Mult_prompt_hitOverlay_mh350.pdf
mv plots/HTdistributionOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/DelayedCaloObject/HTdistributionOverlay_mh350.pdf
mv plots/HTdistribution_trigOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/DelayedCaloObject/HTdistribution_trig_Overlay_mh350.pdf
echo " "
echo "LLP mh=250 GeV, ctau=0.5 m"
rates_delayed_cluster.exe new /eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/TimingTrigger/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-500mm_TuneCP5_14TeV_pythia8_HCAL/CRAB3_mh250_ctau500_TDC_MC/200722_193005/0000 $1 $2 $3 $4 $5 $6
mv rates_new_cond.root rates_new_cond_LLP_pl500.root
echo " "
echo "LLP mh=250 GeV, ctau=1 m"
rates_delayed_cluster.exe new /eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/TimingTrigger/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/CRAB3_mh250_ctau1000_TDC_MC/200722_185009/0000 $1 $2 $3 $4 $5 $6
mv rates_new_cond.root rates_new_cond_LLP_pl1000.root
echo " "
echo "LLP mh=250 GeV, ctau=10 m"
rates_delayed_cluster.exe new /eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/TimingTrigger/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-10000mm_TuneCP5_14TeV_pythia8_HCAL/CRAB3_mh250_ctau10000_TDC_MC/200722_190316/0000 $1 $2 $3 $4 $5 $6
mv rates_new_cond.root rates_new_cond_LLP_pl10000.root
echo " "
echo "Making overlay plots from draw_timingbit.exe"
draw_timingbit_delayed_cluster_mh250.exe
mv rates_new_cond_LLP_pl500.root rates_new_cond_LLP_mh250_pl500.root
mv rates_new_cond_LLP_pl1000.root rates_new_cond_LLP_mh250_pl1000.root
mv rates_new_cond_LLP_pl10000.root rates_new_cond_LLP_mh250_pl10000.root
mv plots/Delayed_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/DelayedCaloObject/Delayed_2x2_MultHBOverlay_mh250.pdf
mv plots/Prompt_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/DelayedCaloObject/Prompt_2x2_MultHBOverlay_mh250.pdf
mv plots/Prompt_Energy_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/DelayedCaloObject/Prompt_Energy_2x2_MultHBOverlay_mh250.pdf
mv plots/Mult_delayed_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/DelayedCaloObject/Mult_delayed_hitOverlay_mh250.pdf
mv plots/Mult_prompt_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/DelayedCaloObject/Mult_prompt_hitOverlay_mh250.pdf
mv plots/HTdistributionOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/DelayedCaloObject/HTdistributionOverlay_mh250.pdf
mv plots/HTdistribution_trigOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/DelayedCaloObject/HTdistribution_trig_Overlay_mh250.pdf
echo " "
echo " "
echo "Efficiency vs. rate plots"
RateEfficiencyPlots_LLP_ctau.exe
RateEfficiencyPlots_LLP_ctau_AddedEff.exe
Rate_HTbin.exe
Eff_ctau.exe
mv /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/EffRate* /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/DelayedCaloObject/
mv /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/Eff_ctau.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/DelayedCaloObject/
mv /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/EffHT_LLP* /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/DelayedCaloObject/
cp MultiplicityHits50ADC3ns_ht120_*.txt /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/DelayedCaloObject/
cp NuGunRates.txt /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/ADC50_3ns_no1HE4HB/DelayedCaloObject/
