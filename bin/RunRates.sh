echo "Compiling rates.cxx, rates_TPin4groups.cxx, and draw_rates.cxx files"
echo " "
scram b -j 8
echo " "
echo " "
echo "Starting with running rates.cxx where all HCAL TPs in a DR cone of L1 jet are considered OR HCAL TP double counting is prevented -- change in code!"
echo "LLP mh=1000 GeV, ctau=1 m"
rates.exe new ../mh1000_pl1000/
mv rates_new_cond.root rates_new_cond_pl1000_4L1Jets.root
echo " "
echo "LLP mh=1000 GeV, ctau=0.5 m"
rates.exe new ../mh1000_pl500/
mv rates_new_cond.root rates_new_cond_pl500_4L1Jets.root
echo " "
echo "LLP mh=1000 GeV, ctau=0.5 m, no PU"
rates.exe new ../mh1000_pl500_noPU/
mv rates_new_cond.root rates_new_cond_pl500_noPU_4L1Jets.root
echo " "
echo "QCD"
rates.exe new ../QCD/
mv rates_new_cond.root rates_new_cond_QCD_4L1Jets.root
echo " "
echo "Neutrino gun new conditions"
rates.exe new ../NeutrinoGun/
mv rates_new_cond.root rates_new_cond_nugunTDC_4L1Jets.root
echo " "
echo "Making overlay plots from draw_rates.exe"
draw_rates.exe
echo "Moving output plots to a new directory"
mv plots/* plots_DRaroundTPandL1_mh1000/
echo " "
echo " "
echo "LLP mh=350 GeV, ctau=1 m"
rates.exe new ../mh350_pl1000_mx160/
mv rates_new_cond.root rates_new_cond_pl1000_4L1Jets.root
echo " "
echo "LLP mh=350 GeV, ctau=0.5 m"
rates.exe new ../mh350_pl500_mx160/
mv rates_new_cond.root rates_new_cond_pl500_4L1Jets.root
echo " "
echo "Making overlay plots from draw_rates.exe"
draw_rates.exe
echo "Moving output plots to a new directory"
mv plots/* plots_DRaroundTPandL1_mh350/
echo "Making overlay ROC curves"
SignalBkgEff_DelayedHitFrac_Plot.exe
echo " "
