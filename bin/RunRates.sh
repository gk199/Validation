echo "Compiling rates.cxx and draw_rates.cxx files"
echo " "
scram b -j 8
echo " "
echo "LLP mh=1000 GeV, ctau=10 m"
rates.exe new ../mh1000_pl10000/
mv rates_new_cond.root rates_new_cond_pl10000_4L1Jets.root
echo " "
echo "LLP mh=1000 GeV, ctau=1 m"
rates.exe new ../mh1000_pl1000/
mv rates_new_cond.root rates_new_cond_pl1000_4L1Jets.root
echo " "
echo "LLP mh=1000 GeV, ctau=0.5 m"
rates.exe new ../mh1000_pl500/
mv rates_new_cond.root rates_new_cond_pl500_4L1Jets.root
echo " "
echo "QCD"
rates.exe new ../QCD/
mv rates_new_cond.root rates_new_cond_QCD_4L1Jets.root
echo " "
echo "Neutrino gun"
rates.exe new ../NeutrinoGun/
mv rates_new_cond.root rates_new_cond_nugunTDC_4L1Jets.root
echo " "
echo "Making overlay plots from draw_rates.exe"
draw_rates.exe
