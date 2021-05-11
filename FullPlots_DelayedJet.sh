echo "Delayed jet trigger with no jet ET requirement"
echo " " 
./RunRates_DelayedJet_PU.sh 4 4 0
hadd -f rates_new_cond_4combinedPU_jet0.root rates_new_cond_LLP_mh125_pl3000.root rates_new_cond_LLP_mh250_pl10000.root rates_new_cond_LLP_mh350_pl10000.root rates_new_cond_LLP_mh1000_pl10000.root

echo "Jet ET efficiency plots"
echo " "
cd L1plots_eff_rates
./JetEfficiency_Plots.sh
root <<EOF
.L efficiency_jetPt_combined.C++
.x efficiency_jetPt_combined.C
.q
EOF

echo "Delayed jet trigger with jet ET requirement"
echo " " 
cd ..
./RunRates_DelayedJet_PU.sh 4 4 40
hadd -f rates_new_cond_4combinedPU_jet40.root rates_new_cond_LLP_mh125_pl3000.root rates_new_cond_LLP_mh250_pl10000.root rates_new_cond_LLP_mh350_pl10000.root rates_new_cond_LLP_mh1000_pl10000.root

echo "Efficiency plots vs event HT, LLP ctau, and TOF delay"
echo " " 
cd L1plots_eff_rates
./Efficiency_Plots.sh
./LLPEfficiency_Plots.sh
./TOFEfficiency_Plots.sh

root <<EOF
.L efficiency_ctau_combined.C++
.x efficiency_ctau_combined.C
.q
EOF

root <<EOF
.L efficiency_ctau_zoom_combined.C++
.x efficiency_ctau_zoom_combined.C
.q
EOF

root <<EOF
.L efficiency_TOF_combined.C++
.x efficiency_TOF_combined.C
.q
EOF
