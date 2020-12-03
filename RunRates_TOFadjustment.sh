echo " " 
echo "LLP mh=125 GeV, ctau=10 m"
rates_tof_adjustment.exe new GeV3_ADC50_ADC105_no1HE4HB/mh125__mx50__pl10000 $1 $2 $3 $4 $5 $6
mv rates_new_cond.root rates_new_cond_LLP_pl10000.root
echo " "
echo "LLP mh=125 GeV, ctau=1 m"
rates_tof_adjustment.exe new GeV3_ADC50_ADC105_no1HE4HB/mh125__mx50__pl1000_ $1 $2 $3 $4 $5 $6
mv rates_new_cond.root rates_new_cond_LLP_pl1000.root
echo " "
echo "LLP mh=125 GeV, ctau=0.5 m"
rates_tof_adjustment.exe new GeV3_ADC50_ADC105_no1HE4HB/mh125__mx50__pl500__ $1 $2 $3 $4 $5 $6
mv rates_new_cond.root rates_new_cond_LLP_pl500.root
echo " "
echo "QCD"
rates_tof_adjustment.exe new GeV3_ADC50_ADC105_no1HE4HB/QCD $1 $2 $3 $4 $5 $6
mv rates_new_cond.root rates_new_cond_QCD.root
echo " "
echo "Neutrino gun new and old conditions"
rates_tof_adjustment.exe new GeV3_ADC50_ADC105_no1HE4HB/NeutrinoGun $1 $2 $3 $4 $5 $6
rates_original.exe def GeV3_ADC50_ADC105_no1HE4HB/NeutrinoGun
echo " "
echo " "
echo "Rates plots from draw_rates.exe"
draw_rates.exe
mv plots/*Rates_emu.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/
mv rates_new_cond.root rates_new_cond_106X_nugun.root
echo " "
echo "Making overlay plots from draw_timingbit.exe"
draw_timingbit_delayed_cluster.exe
mv rates_new_cond_LLP_pl500.root rates_new_cond_LLP_mh125_pl500.root
mv rates_new_cond_LLP_pl1000.root rates_new_cond_LLP_mh125_pl1000.root
mv rates_new_cond_LLP_pl10000.root rates_new_cond_LLP_mh125_pl10000.root
mv plots/Delayed_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Delayed_2x2_MultHBOverlay_mh125.pdf
mv plots/Prompt_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Prompt_2x2_MultHBOverlay_mh125.pdf
mv plots/Prompt_Energy_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Prompt_Energy_2x2_MultHBOverlay_mh125.pdf
mv plots/mhit1Overlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/mhit1_mh125.pdf
mv plots/mhit2Overlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/mhit2_mh125.pdf
mv plots/mhit3Overlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/mhit3_mh125.pdf
# mv plots/DeltaR_L1_delayed_seedOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/DeltaR_L1_delayed_seedOverlay_mh125.pdf
# mv plots/DeltaR_L1_prompt_seedOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/DeltaR_L1_prompt_seedOverlay_mh125.pdf
# mv plots/DeltaR_L1_delayed_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/DeltaR_L1_delayed_hitOverlay_mh125.pdf
# mv plots/DeltaR_L1_prompt_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/DeltaR_L1_prompt_hitOverlay_mh125.pdf
mv plots/Mult_delayed_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Mult_delayed_hitOverlay_mh125.pdf
mv plots/Mult_prompt_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Mult_prompt_hitOverlay_mh125.pdf
mv plots/HTdistributionOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/HTdistributionOverlay_mh125.pdf
mv plots/HTdistribution_trigOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/HTdistribution_trig_Overlay_mh125.pdf
echo " "
echo "LLP mh=1000 GeV, ctau=10 m"
rates_tof_adjustment.exe new GeV3_ADC50_ADC105_no1HE4HB/mh1000_mx450_pl10000 $1 $2 $3 $4 $5 $6
mv rates_new_cond.root rates_new_cond_LLP_pl10000.root
echo " "
echo "LLP mh=1000 GeV, ctau=1 m"
rates_tof_adjustment.exe new GeV3_ADC50_ADC105_no1HE4HB/mh1000_mx450_pl1000_ $1 $2 $3 $4 $5 $6
mv rates_new_cond.root rates_new_cond_LLP_pl1000.root
echo " "
echo "LLP mh=1000 GeV, ctau=0.5 m"
rates_tof_adjustment.exe new GeV3_ADC50_ADC105_no1HE4HB/mh1000_mx450_pl500__ $1 $2 $3 $4 $5 $6
mv rates_new_cond.root rates_new_cond_LLP_pl500.root
echo " "
echo "Making overlay plots from draw_timingbit.exe"
draw_timingbit_delayed_cluster.exe
mv rates_new_cond_LLP_pl500.root rates_new_cond_LLP_mh1000_pl500.root
mv rates_new_cond_LLP_pl1000.root rates_new_cond_LLP_mh1000_pl1000.root
mv rates_new_cond_LLP_pl10000.root rates_new_cond_LLP_mh1000_pl10000.root
mv plots/Delayed_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Delayed_2x2_MultHBOverlay_mh1000.pdf
mv plots/Prompt_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Prompt_2x2_MultHBOverlay_mh1000.pdf
mv plots/Prompt_Energy_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Prompt_Energy_2x2_MultHBOverlay_mh1000.pdf
mv plots/mhit1Overlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/mhit1_mh1000.pdf
mv plots/mhit2Overlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/mhit2_mh1000.pdf
mv plots/mhit3Overlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/mhit3_mh1000.pdf
# mv plots/DeltaR_L1_delayed_seedOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/DeltaR_L1_delayed_seedOverlay_mh1000.pdf
# mv plots/DeltaR_L1_prompt_seedOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/DeltaR_L1_prompt_seedOverlay_mh1000.pdf
# mv plots/DeltaR_L1_delayed_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/DeltaR_L1_delayed_hitOverlay_mh1000.pdf
# mv plots/DeltaR_L1_prompt_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/DeltaR_L1_prompt_hitOverlay_mh1000.pdf
mv plots/Mult_delayed_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Mult_delayed_hitOverlay_mh1000.pdf
mv plots/Mult_prompt_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Mult_prompt_hitOverlay_mh1000.pdf
mv plots/HTdistributionOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/HTdistributionOverlay_mh1000.pdf
mv plots/HTdistribution_trigOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/HTdistribution_trig_Overlay_mh1000.pdf
echo " "
echo "LLP mh=350 GeV, ctau=10 m"
rates_tof_adjustment.exe new GeV3_ADC50_ADC105_no1HE4HB/mh350__mx160_pl10000 $1 $2 $3 $4 $5 $6
mv rates_new_cond.root rates_new_cond_LLP_pl10000.root
echo " "
echo "LLP mh=350 GeV, ctau=1 m"
rates_tof_adjustment.exe new GeV3_ADC50_ADC105_no1HE4HB/mh350__mx160_pl1000_ $1 $2 $3 $4 $5 $6
mv rates_new_cond.root rates_new_cond_LLP_pl1000.root
echo " "
echo "LLP mh=350 GeV, ctau=0.5 m"
rates_tof_adjustment.exe new GeV3_ADC50_ADC105_no1HE4HB/mh350__mx160_pl500__ $1 $2 $3 $4 $5 $6
mv rates_new_cond.root rates_new_cond_LLP_pl500.root
echo " "
echo "Making overlay plots from draw_timingbit.exe"
draw_timingbit_delayed_cluster.exe
mv rates_new_cond_LLP_pl500.root rates_new_cond_LLP_mh350_pl500.root
mv rates_new_cond_LLP_pl1000.root rates_new_cond_LLP_mh350_pl1000.root
mv rates_new_cond_LLP_pl10000.root rates_new_cond_LLP_mh350_pl10000.root
mv plots/Delayed_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Delayed_2x2_MultHBOverlay_mh350.pdf
mv plots/Prompt_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Prompt_2x2_MultHBOverlay_mh350.pdf
mv plots/Prompt_Energy_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Prompt_Energy_2x2_MultHBOverlay_mh350.pdf
mv plots/mhit1Overlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/mhit1_mh350.pdf
mv plots/mhit2Overlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/mhit2_mh350.pdf
mv plots/mhit3Overlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/mhit3_mh350.pdf
# mv plots/DeltaR_L1_delayed_seedOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/DeltaR_L1_delayed_seedOverlay_mh350.pdf
# mv plots/DeltaR_L1_prompt_seedOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/DeltaR_L1_prompt_seedOverlay_mh350.pdf
# mv plots/DeltaR_L1_delayed_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/DeltaR_L1_delayed_hitOverlay_mh350.pdf
# mv plots/DeltaR_L1_prompt_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/DeltaR_L1_prompt_hitOverlay_mh350.pdf
mv plots/Mult_delayed_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Mult_delayed_hitOverlay_mh350.pdf
mv plots/Mult_prompt_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Mult_prompt_hitOverlay_mh350.pdf
mv plots/HTdistributionOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/HTdistributionOverlay_mh350.pdf
mv plots/HTdistribution_trigOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/HTdistribution_trig_Overlay_mh350.pdf
echo " "
echo "LLP mh=250 GeV, ctau=1 m"
rates_tof_adjustment.exe new GeV3_ADC50_ADC105_no1HE4HB/mh250__mx120_pl1000_ $1 $2 $3 $4 $5 $6
echo " "
echo "Efficiency vs. rate plots"
RateEfficiencyPlots_LLP_ctau.exe
RateEfficiencyPlots_LLP_ctau_AddedEff.exe
Rate_HTbin.exe
Eff_ctau.exe
mv /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/EffRate* /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/
mv /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/Eff_ctau.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/
mv /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/EffHT_LLP* /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/
cp MultiplicityHits50ADC3ns_ht120_*.txt /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/
cp NuGunRates.txt /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/
