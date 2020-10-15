echo "Compiling all files with scram b"
scram b -j 8
for TDC_HB in 3 4
do
for TDC_HE in 3 4
do
for GeV_HB in 2 3
do
for GeV_HE in 1 2
do
for prompt_TP in 3 5
do
for prompt_2x2 in 4 6
do
echo " "
echo "Running rates delayed calo object with parameters of: ROC_delayed-${TDC_HB}ns${GeV_HB}GeVHB-${TDC_HE}ns${GeV_HE}GeVHE_promptTP-${prompt_TP}GeV_prompt2x2-${prompt_2x2}GeV"
./RunRates_DelayedCaloObject.sh $TDC_HB $TDC_HE $GeV_HB $GeV_HE ${prompt_TP} ${prompt_2x2}
mkdir ROC_delayed-${TDC_HB}ns${GeV_HB}GeVHB-${TDC_HE}ns${GeV_HE}GeVHE_promptTP-${prompt_TP}GeV_prompt2x2-${prompt_2x2}GeV
ROC_signal_background.exe
mv plots/ROC_120timing* ROC_delayed-${TDC_HB}ns${GeV_HB}GeVHB-${TDC_HE}ns${GeV_HE}GeVHE_promptTP-${prompt_TP}GeV_prompt2x2-${prompt_2x2}GeV/
done
done
done
done
done
done
