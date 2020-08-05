echo "Compiling rates.cxx, draw_timingbit.cxx files"
echo " "
scram b -j 8
echo " "
echo " "
echo "LLP mh=125 GeV, ctau=10 m"
rates.exe new GeV3_ADC50_ADC105_no1HE4HB/mh125__mx50__pl10000
mv rates_new_cond.root rates_new_cond_LLP_pl10000.root
echo " "
echo "LLP mh=125 GeV, ctau=1 m"
rates.exe new GeV3_ADC50_ADC105_no1HE4HB/mh125__mx50__pl1000_
mv rates_new_cond.root rates_new_cond_LLP_pl1000.root
echo " "
echo "LLP mh=125 GeV, ctau=0.5 m"
rates.exe new GeV3_ADC50_ADC105_no1HE4HB/mh125__mx50__pl500__
mv rates_new_cond.root rates_new_cond_LLP_pl500.root
echo " "
echo "QCD"
rates.exe new GeV3_ADC50_ADC105_no1HE4HB/QCD
mv rates_new_cond.root rates_new_cond_QCD.root
echo " "
echo "Neutrino gun new and old conditions"
rates.exe new GeV3_ADC50_ADC105_no1HE4HB/NeutrinoGun
rates_original.exe def GeV3_ADC50_ADC105_no1HE4HB/NeutrinoGun
echo " "
echo " "
echo "Rates plots from draw_rates.exe"
draw_rates.exe
mv plots/*Rates_emu.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/
mv rates_new_cond.root rates_new_cond_106X_nugun.root
echo " "
echo "Making overlay plots from draw_timingbit.exe"
draw_timingbit.exe
mv rates_new_cond_LLP_pl1000.root rates_new_cond_LLP_mh125_pl1000.root
mv plots/ADC50_3ns_4JetMultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/ADC50_3ns_no1HE4HB_4JetMultHBOverlay_mh125.pdf
mv plots/ADC50_3ns_4JetMultHEOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/ADC50_3ns_no1HE4HB_4JetMultHEOverlay_mh125.pdf
mv plots/ADC50_3ns_4JetMultHBHEOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/ADC50_3ns_no1HE4HB_4JetMultHBHEOverlay_mh125.pdf
echo " "
echo "LLP mh=1000 GeV, ctau=10 m"
rates.exe new GeV3_ADC50_ADC105_no1HE4HB/mh1000_mx450_pl10000
mv rates_new_cond.root rates_new_cond_LLP_pl10000.root
echo " "
echo "LLP mh=1000 GeV, ctau=1 m"
rates.exe new GeV3_ADC50_ADC105_no1HE4HB/mh1000_mx450_pl1000_
mv rates_new_cond.root rates_new_cond_LLP_pl1000.root
echo " "
echo "LLP mh=1000 GeV, ctau=0.5 m"
rates.exe new GeV3_ADC50_ADC105_no1HE4HB/mh1000_mx450_pl500__
mv rates_new_cond.root rates_new_cond_LLP_pl500.root
echo " "
echo "Making overlay plots from draw_timingbit.exe"
draw_timingbit.exe
mv rates_new_cond_LLP_pl1000.root rates_new_cond_LLP_mh1000_pl1000.root
mv rates_new_cond_LLP_pl10000.root rates_new_cond_LLP_mh1000_pl10000.root
mv plots/ADC50_3ns_4JetMultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/ADC50_3ns_no1HE4HB_4JetMultHBOverlay_mh1000.pdf
mv plots/ADC50_3ns_4JetMultHEOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/ADC50_3ns_no1HE4HB_4JetMultHEOverlay_mh1000.pdf
mv plots/ADC50_3ns_4JetMultHBHEOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/ADC50_3ns_no1HE4HB_4JetMultHBHEOverlay_mh1000.pdf
echo " "
echo "LLP mh=350 GeV, ctau=10 m"
rates.exe new GeV3_ADC50_ADC105_no1HE4HB/mh350__mx160_pl10000
mv rates_new_cond.root rates_new_cond_LLP_pl10000.root
echo " "
echo "LLP mh=350 GeV, ctau=1 m"
rates.exe new GeV3_ADC50_ADC105_no1HE4HB/mh350__mx160_pl1000_
mv rates_new_cond.root rates_new_cond_LLP_pl1000.root
echo " "
echo "LLP mh=350 GeV, ctau=0.5 m"
rates.exe new GeV3_ADC50_ADC105_no1HE4HB/mh350__mx160_pl500__
mv rates_new_cond.root rates_new_cond_LLP_pl500.root
echo " "
echo "Making overlay plots from draw_timingbit.exe"
draw_timingbit.exe
mv rates_new_cond_LLP_pl1000.root rates_new_cond_LLP_mh350_pl1000.root
mv plots/ADC50_3ns_4JetMultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/ADC50_3ns_no1HE4HB_4JetMultHBOverlay_mh350.pdf
mv plots/ADC50_3ns_4JetMultHEOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/ADC50_3ns_no1HE4HB_4JetMultHEOverlay_mh350.pdf
mv plots/ADC50_3ns_4JetMultHBHEOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/ADC50_3ns_no1HE4HB_4JetMultHBHEOverlay_mh350.pdf
echo " "
echo "LLP mh=250 GeV, ctau=1 m"
rates.exe new GeV3_ADC50_ADC105_no1HE4HB/mh250__mx120_pl1000_
mv rates_new_cond.root rates_new_cond_LLP_mh250_pl1000.root
echo " "
echo "Making overlay plots from draw_1mtimingbit.exe"
draw_1mtimingbit.exe
mv plots/ADC50_3ns_4JetMultHBOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/ADC50_3ns_no1HE4HB_4JetMultHBOverlay_1m.pdf
mv plots/ADC50_3ns_4JetMultHEOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/ADC50_3ns_no1HE4HB_4JetMultHEOverlay_1m.pdf
mv plots/ADC50_3ns_4JetMultHBHEOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/ADC50_3ns_no1HE4HB_4JetMultHBHEOverlay_1m.pdf
echo " "
echo "Efficiency vs. rate plots"
RateEfficiencyPlots_LLP_ctau.exe
Rate_HTbin.exe
Eff_ctau.exe
mv /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/EffRate* /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/
mv /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/Eff_ctau.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/
mv /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/EffHT_LLP* /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/
cp MultiplicityHits50ADC3ns_ht120_*.txt /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/
cp NuGunRates.txt /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/ADC50_3ns_no1HE4HB/
