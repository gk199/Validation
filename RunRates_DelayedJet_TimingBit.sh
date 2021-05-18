echo "Compiling rates.cxx, draw_timingbit.cxx files"
echo " "
scram b -j 8
echo " "
echo " "
echo "LLP mh=125 GeV, ctau=3 m"
rates_LLPflag.exe new TimingBit/mh125__pl3000__PU $1 $2 $3
mv rates_new_cond.root rates_new_cond_LLP_pl3000.root
echo " "
echo "LLP mh=125 GeV, ctau=30 m"
rates_LLPflag.exe new TimingBit/mh125__pl30000_PU $1 $2 $3
mv rates_new_cond.root rates_new_cond_LLP_pl30000.root
echo " "
echo "QCD"
rates_LLPflag.exe new TimingBit/QCD_PU $1 $2 $3
mv rates_new_cond.root rates_new_cond_QCD.root
echo " "
echo "Neutrino gun new and old conditions"
rates_LLPflag.exe new TimingBit/RelValNuGun $1 $2 $3
rates_original.exe def TimingBit/RelValNuGun
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
mv plots/Mult_delayed_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hitOverlay_mh125.pdf
mv plots/Mult_delayed_hit_jetETOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hit_jetETOverlay_mh125.pdf
mv plots/Mult_delayed_hit_promptVOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hit_promptVOverlay_mh125.pdf
mv plots/Mult_delayed_hit_jetETpromptVOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hit_jetETpromptVOverlay_mh125.pdf
mv plots/Mult_delayed_hitOverlay_Log.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hitOverlay_Log_mh125.pdf
mv plots/Mult_delayed_hit_jetETOverlay_Log.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hit_jetETOverlay_Log_mh125.pdf
mv plots/Mult_delayed_hit_promptVOverlay_Log.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hit_promptVOverlay_Log_mh125.pdf
mv plots/Mult_delayed_hit_jetETpromptVOverlay_Log.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hit_jetETpromptVOverlay_Log_mh125.pdf
mv plots/Mult_prompt_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_prompt_hitOverlay_mh125.pdf
mv plots/HTdistributionOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/HTdistributionOverlay_mh125.pdf
mv plots/HTdistribution_trigOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/HTdistribution_trig_Overlay_mh125.pdf
echo " "
echo " "
echo "LLP mh=250 GeV, ctau=0.5 m"
rates_LLPflag.exe new TimingBit/mh250__pl500___PU $1 $2 $3
mv rates_new_cond.root rates_new_cond_LLP_pl500.root
echo " "
echo "LLP mh=250 GeV, ctau=1 m"
rates_LLPflag.exe new TimingBit/mh250__pl1000__PU $1 $2 $3
mv rates_new_cond.root rates_new_cond_LLP_pl1000.root
echo " "
echo "LLP mh=250 GeV, ctau=10 m"
rates_LLPflag.exe new TimingBit/mh250__pl10000_PU $1 $2 $3
mv rates_new_cond.root rates_new_cond_LLP_pl10000.root
echo " "
echo " "
echo "Making overlay plots from draw_timingbit.exe"
draw_timingbit_delayed_cluster_mh250.exe
mv rates_new_cond_LLP_pl500.root rates_new_cond_LLP_mh250_pl500.root
mv rates_new_cond_LLP_pl1000.root rates_new_cond_LLP_mh250_pl1000.root
mv rates_new_cond_LLP_pl10000.root rates_new_cond_LLP_mh250_pl10000.root
mv plots/Mult_delayed_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hitOverlay_mh250.pdf
mv plots/Mult_delayed_hit_jetETOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hit_jetETOverlay_mh250.pdf
mv plots/Mult_delayed_hit_promptVOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hit_promptVOverlay_mh250.pdf
mv plots/Mult_delayed_hit_jetETpromptVOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hit_jetETpromptVOverlay_mh250.pdf
mv plots/Mult_delayed_hitOverlay_Log.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hitOverlay_Log_mh250.pdf
mv plots/Mult_delayed_hit_jetETOverlay_Log.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hit_jetETOverlay_Log_mh250.pdf
mv plots/Mult_delayed_hit_promptVOverlay_Log.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hit_promptVOverlay_Log_mh250.pdf
mv plots/Mult_delayed_hit_jetETpromptVOverlay_Log.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hit_jetETpromptVOverlay_Log_mh250.pdf
mv plots/Mult_prompt_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_prompt_hitOverlay_mh250.pdf
mv plots/HTdistributionOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/HTdistributionOverlay_mh250.pdf
mv plots/HTdistribution_trigOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/HTdistribution_trig_Overlay_mh250.pdf
echo " "
echo " "
echo "LLP mh=350 GeV, ctau=0.5 m"
rates_LLPflag.exe new TimingBit/mh350__pl500___PU $1 $2 $3
mv rates_new_cond.root rates_new_cond_LLP_pl500.root
echo " "
echo "LLP mh=350 GeV, ctau=1 m"
rates_LLPflag.exe new TimingBit/mh350__pl1000__PU $1 $2 $3
mv rates_new_cond.root rates_new_cond_LLP_pl1000.root
echo " "
echo "LLP mh=350 GeV, ctau=10 m"
rates_LLPflag.exe new TimingBit/mh350__pl10000_PU $1 $2 $3
mv rates_new_cond.root rates_new_cond_LLP_pl10000.root
echo " "
echo " "
echo "Making overlay plots from draw_timingbit.exe"
draw_timingbit_delayed_cluster_mh350.exe
mv rates_new_cond_LLP_pl500.root rates_new_cond_LLP_mh350_pl500.root
mv rates_new_cond_LLP_pl1000.root rates_new_cond_LLP_mh350_pl1000.root
mv rates_new_cond_LLP_pl10000.root rates_new_cond_LLP_mh350_pl10000.root
mv plots/Mult_delayed_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hitOverlay_mh350.pdf
mv plots/Mult_delayed_hit_jetETOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hit_jetETOverlay_mh350.pdf
mv plots/Mult_delayed_hit_promptVOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hit_promptVOverlay_mh350.pdf
mv plots/Mult_delayed_hit_jetETpromptVOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hit_jetETpromptVOverlay_mh350.pdf
mv plots/Mult_delayed_hitOverlay_Log.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hitOverlay_Log_mh350.pdf
mv plots/Mult_delayed_hit_jetETOverlay_Log.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hit_jetETOverlay_Log_mh350.pdf
mv plots/Mult_delayed_hit_promptVOverlay_Log.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hit_promptVOverlay_Log_mh350.pdf
mv plots/Mult_delayed_hit_jetETpromptVOverlay_Log.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hit_jetETpromptVOverlay_Log_mh350.pdf
mv plots/Mult_prompt_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_prompt_hitOverlay_mh350.pdf
mv plots/HTdistributionOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/HTdistributionOverlay_mh350.pdf
mv plots/HTdistribution_trigOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/HTdistribution_trig_Overlay_mh350.pdf
echo " "
echo "LLP mh=1000 GeV, ctau=10 m"
rates_LLPflag.exe new TimingBit/mh1000_pl10000_PU $1 $2 $3
mv rates_new_cond.root rates_new_cond_LLP_pl10000.root
echo " "
echo "LLP mh=1000 GeV, ctau=100 m"
rates_LLPflag.exe new TimingBit/mh1000_pl100000PU $1 $2 $3
mv rates_new_cond.root rates_new_cond_LLP_pl100000.root
echo " "
echo "Making overlay plots from draw_timingbit.exe"
draw_timingbit_delayed_cluster_mh1000.exe
mv rates_new_cond_LLP_pl10000.root rates_new_cond_LLP_mh1000_pl10000.root
mv rates_new_cond_LLP_pl100000.root rates_new_cond_LLP_mh1000_pl100000.root
mv plots/Mult_delayed_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hitOverlay_mh1000.pdf
mv plots/Mult_delayed_hit_jetETOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hit_jetETOverlay_mh1000.pdf
mv plots/Mult_delayed_hit_promptVOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hit_promptVOverlay_mh1000.pdf
mv plots/Mult_delayed_hit_jetETpromptVOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hit_jetETpromptVOverlay_mh1000.pdf
mv plots/Mult_delayed_hitOverlay_Log.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hitOverlay_Log_mh1000.pdf
mv plots/Mult_delayed_hit_jetETOverlay_Log.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hit_jetETOverlay_Log_mh1000.pdf
mv plots/Mult_delayed_hit_promptVOverlay_Log.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hit_promptVOverlay_Log_mh1000.pdf
mv plots/Mult_delayed_hit_jetETpromptVOverlay_Log.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_delayed_hit_jetETpromptVOverlay_Log_mh1000.pdf
mv plots/Mult_prompt_hitOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/Mult_prompt_hitOverlay_mh1000.pdf
mv plots/HTdistributionOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/HTdistributionOverlay_mh1000.pdf
mv plots/HTdistribution_trigOverlay.pdf /eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/HTdistribution_trig_Overlay_mh1000.pdf
echo " "
