echo "Compiling rates.cxx, draw_timingbit.cxx files"
echo " "
scram b -j 8
echo " "
echo " "
echo "LLP mh=125 GeV, ctau=10 m"
rates_delayed_cluster.exe new GeV3_ADC50_ADC105_no1HE4HB/mh125__mx50__pl10000
mv rates_new_cond.root rates_new_cond_LLP_pl10000.root
echo " "
echo "LLP mh=125 GeV, ctau=1 m"
rates_delayed_cluster.exe new GeV3_ADC50_ADC105_no1HE4HB/mh125__mx50__pl1000_
mv rates_new_cond.root rates_new_cond_LLP_pl1000.root
echo " "
echo "LLP mh=125 GeV, ctau=0.5 m"
rates_delayed_cluster.exe new GeV3_ADC50_ADC105_no1HE4HB/mh125__mx50__pl500__
mv rates_new_cond.root rates_new_cond_LLP_pl500.root
echo " "
echo "QCD"
rates_delayed_cluster.exe new GeV3_ADC50_ADC105_no1HE4HB/QCD
mv rates_new_cond.root rates_new_cond_QCD.root
echo " "
echo "Neutrino gun new and old conditions"
rates_delayed_cluster.exe new GeV3_ADC50_ADC105_no1HE4HB/NeutrinoGun
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
mv rates_new_cond_LLP_pl1000.root rates_new_cond_LLP_mh125_pl1000.root
mv plots/Delayed_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Delayed_2x2_MultHBOverlay_mh125.pdf
mv plots/Delayed_6x6_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Delayed_6x6_MultHBOverlay_mh125.pdf
mv plots/Delayed_full_6x6_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Delayed_full_6x6_MultHBOverlay_mh125.pdf
mv plots/Delayed_2x2_0GeV_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Delayed_2x2_0GeV_MultHBOverlay_mh125.pdf
mv plots/Delayed_2x2_1GeV_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Delayed_2x2_1GeV_MultHBOverlay_mh125.pdf
mv plots/Delayed_2x2_2GeV_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Delayed_2x2_2GeV_MultHBOverlay_mh125.pdf
mv plots/Delayed_2x2_3GeV_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Delayed_2x2_3GeV_MultHBOverlay_mh125.pdf
mv plots/HTdistributionOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/HTdistributionOverlay_mh125.pdf
mv plots/HTdistribution_trigOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/HTdistribution_trig_Overlay_mh125.pdf
echo " "
echo "LLP mh=1000 GeV, ctau=10 m"
rates_delayed_cluster.exe new GeV3_ADC50_ADC105_no1HE4HB/mh1000_mx450_pl10000
mv rates_new_cond.root rates_new_cond_LLP_pl10000.root
echo " "
echo "LLP mh=1000 GeV, ctau=1 m"
rates_delayed_cluster.exe new GeV3_ADC50_ADC105_no1HE4HB/mh1000_mx450_pl1000_
mv rates_new_cond.root rates_new_cond_LLP_pl1000.root
echo " "
echo "LLP mh=1000 GeV, ctau=0.5 m"
rates_delayed_cluster.exe new GeV3_ADC50_ADC105_no1HE4HB/mh1000_mx450_pl500__
mv rates_new_cond.root rates_new_cond_LLP_pl500.root
echo " "
echo "Making overlay plots from draw_timingbit.exe"
draw_timingbit_delayed_cluster.exe
mv rates_new_cond_LLP_pl1000.root rates_new_cond_LLP_mh1000_pl1000.root
mv rates_new_cond_LLP_pl10000.root rates_new_cond_LLP_mh1000_pl10000.root
mv plots/Delayed_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Delayed_2x2_MultHBOverlay_mh1000.pdf
mv plots/Delayed_6x6_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Delayed_6x6_MultHBOverlay_mh1000.pdf
mv plots/Delayed_full_6x6_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Delayed_full_6x6_MultHBOverlay_mh1000.pdf
mv plots/Delayed_2x2_0GeV_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Delayed_2x2_0GeV_MultHBOverlay_mh1000.pdf
mv plots/Delayed_2x2_1GeV_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Delayed_2x2_1GeV_MultHBOverlay_mh1000.pdf
mv plots/Delayed_2x2_2GeV_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Delayed_2x2_2GeV_MultHBOverlay_mh1000.pdf
mv plots/Delayed_2x2_3GeV_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Delayed_2x2_3GeV_MultHBOverlay_mh1000.pdf
mv plots/HTdistributionOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/HTdistributionOverlay_mh1000.pdf
mv plots/HTdistribution_trigOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/HTdistribution_trig_Overlay_mh1000.pdf
echo " "
echo "LLP mh=350 GeV, ctau=10 m"
rates_delayed_cluster.exe new GeV3_ADC50_ADC105_no1HE4HB/mh350__mx160_pl10000
mv rates_new_cond.root rates_new_cond_LLP_pl10000.root
echo " "
echo "LLP mh=350 GeV, ctau=1 m"
rates_delayed_cluster.exe new GeV3_ADC50_ADC105_no1HE4HB/mh350__mx160_pl1000_
mv rates_new_cond.root rates_new_cond_LLP_pl1000.root
echo " "
echo "LLP mh=350 GeV, ctau=0.5 m"
rates_delayed_cluster.exe new GeV3_ADC50_ADC105_no1HE4HB/mh350__mx160_pl500__
mv rates_new_cond.root rates_new_cond_LLP_pl500.root
echo " "
echo "Making overlay plots from draw_timingbit.exe"
draw_timingbit_delayed_cluster.exe
mv rates_new_cond_LLP_pl1000.root rates_new_cond_LLP_mh350_pl1000.root
mv plots/Delayed_2x2_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Delayed_2x2_MultHBOverlay_mh350.pdf
mv plots/Delayed_6x6_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Delayed_6x6_MultHBOverlay_mh350.pdf
mv plots/Delayed_full_6x6_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Delayed_full_6x6_MultHBOverlay_mh350.pdf
mv plots/Delayed_2x2_0GeV_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Delayed_2x2_0GeV_MultHBOverlay_mh350.pdf
mv plots/Delayed_2x2_1GeV_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Delayed_2x2_1GeV_MultHBOverlay_mh350.pdf
mv plots/Delayed_2x2_2GeV_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Delayed_2x2_2GeV_MultHBOverlay_mh350.pdf
mv plots/Delayed_2x2_3GeV_MultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/Delayed_2x2_3GeV_MultHBOverlay_mh350.pdf
mv plots/HTdistributionOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/HTdistributionOverlay_mh350.pdf
mv plots/HTdistribution_trigOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/HTdistribution_trig_Overlay_mh350.pdf
echo " "
echo "LLP mh=250 GeV, ctau=1 m"
rates_delayed_cluster.exe new GeV3_ADC50_ADC105_no1HE4HB/mh250__mx120_pl1000_
echo " "
echo "Efficiency vs. rate plots"
RateEfficiencyPlots_LLP_ctau.exe
Rate_HTbin.exe
Eff_ctau.exe
mv /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/EffRate* /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/
mv /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/Eff_ctau.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/
mv /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/EffHT_LLP* /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/
cp MultiplicityHits50ADC3ns_ht120_*.txt /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/
cp NuGunRates.txt /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/DelayedCaloObject/
