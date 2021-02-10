echo "Compiling rates.cxx, draw_timingbit.cxx files"
echo " "
scram b -j 8
echo " "
echo " "
echo "LLP mh=125 GeV, ctau=3 m"
rates_delayed_cluster.exe new mh125__pl3000__ $1 $2 $3 $4 $5 $6
mv rates_new_cond.root rates_new_cond_LLP_pl3000.root
echo " "
echo "LLP mh=125 GeV, ctau=30 m"
rates_delayed_cluster.exe new mh125__pl30000_ $1 $2 $3 $4 $5 $6
mv rates_new_cond.root rates_new_cond_LLP_pl30000.root
echo " "
echo "QCD"
rates_delayed_cluster.exe new QCD $1 $2 $3 $4 $5 $6
mv rates_new_cond.root rates_new_cond_QCD.root
echo " "
echo "Neutrino gun new and old conditions"
rates_delayed_cluster.exe new RelValNuGun $1 $2 $3 $4 $5 $6
rates_original.exe def RelValNuGun
echo " "
echo " "
echo "Rates plots from draw_rates.exe"
draw_rates.exe
mv plots/*Rates_emu.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/
mv rates_new_cond.root rates_new_cond_110X_nugun.root
echo " "
echo "Making overlay plots from draw_timingbit.exe"
draw_timingbit_delayed_cluster_mh125.exe
mv rates_new_cond_LLP_pl30000.root rates_new_cond_LLP_mh125_pl30000.root
mv rates_new_cond_LLP_pl3000.root rates_new_cond_LLP_mh125_pl3000.root
mv plots/Delayed_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Delayed_2x2_MultHBOverlay_mh125.pdf
mv plots/Prompt_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Prompt_2x2_MultHBOverlay_mh125.pdf
mv plots/Prompt_Energy_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Prompt_Energy_2x2_MultHBOverlay_mh125.pdf
mv plots/Mult_delayed_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hitOverlay_mh125.pdf
mv plots/Mult_prompt_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_prompt_hitOverlay_mh125.pdf
mv plots/HTdistributionOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/HTdistributionOverlay_mh125.pdf
mv plots/HTdistribution_trigOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/HTdistribution_trig_Overlay_mh125.pdf
echo " "
echo " "
echo "LLP mh=250 GeV, ctau=10 m"
rates_delayed_cluster.exe new mh250__pl10000_ $1 $2 $3 $4 $5 $6
mv rates_new_cond.root rates_new_cond_LLP_pl10000.root
echo " "
echo "LLP mh=250 GeV, ctau=1 m"
rates_delayed_cluster.exe new mh250__pl1000__ $1 $2 $3 $4 $5 $6
mv rates_new_cond.root rates_new_cond_LLP_pl1000.root
echo " "
echo "LLP mh=250 GeV, ctau=0.5 m"
rates_delayed_cluster.exe new mh250__pl500___ $1 $2 $3 $4 $5 $6
mv rates_new_cond.root rates_new_cond_LLP_pl500.root
echo " "
echo "Making overlay plots from draw_timingbit.exe"
draw_timingbit_delayed_cluster_mh250.exe
mv rates_new_cond_LLP_pl500.root rates_new_cond_LLP_mh250_pl500.root
mv rates_new_cond_LLP_pl1000.root rates_new_cond_LLP_mh250_pl1000.root
mv rates_new_cond_LLP_pl10000.root rates_new_cond_LLP_mh250_pl10000.root
mv plots/Delayed_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Delayed_2x2_MultHBOverlay_mh250.pdf
mv plots/Prompt_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Prompt_2x2_MultHBOverlay_mh250.pdf
mv plots/Prompt_Energy_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Prompt_Energy_2x2_MultHBOverlay_mh250.pdf
mv plots/Mult_delayed_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hitOverlay_mh250.pdf
mv plots/Mult_prompt_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_prompt_hitOverlay_mh250.pdf
mv plots/HTdistributionOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/HTdistributionOverlay_mh250.pdf
mv plots/HTdistribution_trigOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/HTdistribution_trig_Overlay_mh250.pdf
echo " "
echo " "
echo "LLP mh=350 GeV, ctau=10 m"
rates_delayed_cluster.exe new mh350__pl10000_ $1 $2 $3 $4 $5 $6
mv rates_new_cond.root rates_new_cond_LLP_pl10000.root
echo " "
echo "LLP mh=350 GeV, ctau=1 m"
rates_delayed_cluster.exe new mh350__pl1000__ $1 $2 $3 $4 $5 $6
mv rates_new_cond.root rates_new_cond_LLP_pl1000.root
echo " "
echo "LLP mh=350 GeV, ctau=0.5 m"
rates_delayed_cluster.exe new mh350__pl500___ $1 $2 $3 $4 $5 $6
mv rates_new_cond.root rates_new_cond_LLP_pl500.root
echo " "
echo "Making overlay plots from draw_timingbit.exe"
draw_timingbit_delayed_cluster_mh350.exe
mv rates_new_cond_LLP_pl500.root rates_new_cond_LLP_mh350_pl500.root
mv rates_new_cond_LLP_pl1000.root rates_new_cond_LLP_mh350_pl1000.root
mv rates_new_cond_LLP_pl10000.root rates_new_cond_LLP_mh350_pl10000.root
mv plots/Delayed_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Delayed_2x2_MultHBOverlay_mh350.pdf
mv plots/Prompt_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Prompt_2x2_MultHBOverlay_mh350.pdf
mv plots/Prompt_Energy_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Prompt_Energy_2x2_MultHBOverlay_mh350.pdf
mv plots/Mult_delayed_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hitOverlay_mh350.pdf
mv plots/Mult_prompt_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_prompt_hitOverlay_mh350.pdf
mv plots/HTdistributionOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/HTdistributionOverlay_mh350.pdf
mv plots/HTdistribution_trigOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/HTdistribution_trig_Overlay_mh350.pdf
echo " "
echo "LLP mh=1000 GeV, ctau=10 m"
rates_delayed_cluster.exe new mh1000_pl10000_ $1 $2 $3 $4 $5 $6
mv rates_new_cond.root rates_new_cond_LLP_pl10000.root
echo " "
echo "LLP mh=1000 GeV, ctau=100 m"
rates_delayed_cluster.exe new mh1000_pl100000 $1 $2 $3 $4 $5 $6
mv rates_new_cond.root rates_new_cond_LLP_pl100000.root
echo " "
echo "Making overlay plots from draw_timingbit.exe"
draw_timingbit_delayed_cluster_mh1000.exe
mv rates_new_cond_LLP_pl10000.root rates_new_cond_LLP_mh1000_pl1000.root
mv rates_new_cond_LLP_pl100000.root rates_new_cond_LLP_mh1000_pl10000.root
mv plots/Delayed_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Delayed_2x2_MultHBOverlay_mh1000.pdf
mv plots/Prompt_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Prompt_2x2_MultHBOverlay_mh1000.pdf
mv plots/Prompt_Energy_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Prompt_Energy_2x2_MultHBOverlay_mh1000.pdf
mv plots/Mult_delayed_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hitOverlay_mh1000.pdf
mv plots/Mult_prompt_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_prompt_hitOverlay_mh1000.pdf
mv plots/HTdistributionOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/HTdistributionOverlay_mh1000.pdf
mv plots/HTdistribution_trigOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/HTdistribution_trig_Overlay_mh1000.pdf
echo " "
echo " "
echo "Efficiency vs. rate plots"
RateEfficiencyPlots_LLP_ctau.exe
RateEfficiencyPlots_LLP_ctau_AddedEff.exe
Rate_HTbin.exe
Eff_ctau.exe
mv /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/EffRate* /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/
mv /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/Eff_ctau.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/
mv /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/EffHT_LLP* /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/
cp MultiplicityHits50ADC3ns_ht120_*.txt /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/
cp NuGunRates.txt /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/
