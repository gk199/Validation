# HCAL Trigger Studies for LLPs

## Setup
These scripts use the L1Ntuple framework, which should be set up as described here:

[https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TStage2Instructions#Environment_Setup_with_Integrati](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TStage2Instructions#Environment_Setup_with_Integrati)

For depth (timing, energy) information, some branches have not yet been rebased in 112X. To get the information, in a 110X environment that has been set up with either of:
```
git cms-merge-topic --unsafe lwang046:upgradeHcalTPs-timeingbit-110X
git cms-merge-topic --unsafe georgia14:upgradeHcalTPs-l1t-110X
```
the [linked files](https://github.com/lwang046/cmssw/commit/a06ab4c691a9a3f640bdeefe99d096d2c6e9a9ab) can be copied over to the 112X environment, and all compiled.

After setting up the L1Ntuple environment, issue the following:
```
git clone git@github.com:cms-hcal-trigger/Validation.git HcalTrigger/Validation
```

To add gen particle information (parton position, momentum, etc) in `L1Trigger/L1TNtuples/plugins/L1GenTreeProducer.cc` add:
```
l1GenData_->partVx.push_back(p.vertex().X());	
l1GenData_->partVy.push_back(p.vertex().Y());
l1GenData_->partVz.push_back(p.vertex().Z());
l1GenData_->partPx.push_back(p.px());
l1GenData_->partPy.push_back(p.py());
l1GenData_->partPz.push_back(p.pz());
l1GenData_->partHardProcess.push_back(p.isHardProcess());
```
and in `L1Trigger/L1TNtuples/interface/L1AnalysisGeneratorDataFormat.h` add:
```
partVx.resize(0);
partVy.resize(0);
partVz.resize(0);
partPx.resize(0);
partPy.resize(0);
partPz.resize(0);
partHardProcess.resize(0);
std::vector<float> partVx;
std::vector<float> partVy;
std::vector<float> partVz;
std::vector<float> partPx;
std::vector<float> partPy;
std::vector<float> partPz;
std::vector<int> partHardProcess;
```
Make sure this is done before processing L1Ntuples, otherwise the needed info is not saved in them.

## L1Ntuples
L1Ntuples are made from GEN-SIM-DIGI-RAW MC samples, as made [here](https://github.com/gk199/MonteCarlo_PrivateProduction/tree/master/LLP_TDC). To produce L1Ntuples, edit the input file name in `ntuple_maker_def_*.py`, and simply do `cmsRun ntuple_maker_def_*.py`. To make all the needed L1Ntuples at once, the bash script `ProcessL1Ntuples.sh` is used with `./ProcessL1Ntuples.sh`. 

## Trigger Emulation -- Rates and Efficiencies
HCAL timing and energy values are avaliable at each cell, allowing for the formation of a delayed jet trigger based on clusters of delayed cells. Delayed jets can be rejected based on multiple prompt towers nearby. `rates_delayed_cluster.cxx` simulates the trigger algorithm (only using events where the LLP is expected to intersect the calorimeter, our triggerability restriction). Then plots are made and output to the [EOS webpage](https://gkopp.web.cern.ch/gkopp/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/). This is all run with 
```
./RunRates_DelayedCaloObject.sh 4 4 3 3 2.5 4
./RunRates_DelayedCaloObject_PU.sh 4 4 3 3 2.5 4
```
where the arguments are HB TDC, HE TDC, HB energy (GeV), HE energy (GeV), prompt TP energy veto, and prompt 2x2 energy veto. Currently 2 cells in a 4x4 region make it "delayed" and a jet is rejected if there are any prompt high energy regions in it. These last two parameters are set in `rates_delayed_cluster.cxx`.

A few things to check in the code:
* Number of cells to form a delayed seed (set by `delayed_4x4_variable`)
* Delayed seeds required (set by `delayed_calo_objects`)
* Prompt veto imposed (set by `prompt_energy`)
* Triggerability restrictions (set by `triggerableJets`)
* Jet pT cut (set by `l1emu_->jetEt[jetIt] < `)
* 1 bit per tower, or 4 bits per 2x2 (set by `TPdelayed`)

In the current format, the first two parameters (TDC values) are set from scanning over the background QCD sample, and determining what TDC value 90% (for example) of the background is below. This is done in:
```
rates_bkgTDCdist.exe new QCD 2 2
```
where the last two parameters are the HB and HE energy values. This then makes multiple txt files listing the ieta, and for each depth, what TDC value most of the background is below.

An updated trigger algorithm using a per tower prompt veto is in `rates_delayed_jet.cxx`. This requires 2 LLP flagged TT per jet, where a LLP flagged TT has a delayed cell, and is vetoed by a prompt cell in the tower. The only arguments needed are the HB and HE GeV (per cell), and the triggered jet ET requirements:
```
./RunRates_DelayedJet_PU.sh 4 4 40
./RunRates_DelayedJet_TimingBit.sh 40
```
The timing bit file bases the trigger on the set timing bit (6 bits = depth, 10, 10, 01, 01, veto) in the L1Ntuple. Bit masks are used in `bin/rates_LLPflag.cxx` to access the depth and timing bits. This is set in:
```
../../SimCalorimetry/HcalTrigPrimAlgos/src/HcalTimingBit.cc
../../SimCalorimetry/HcalTrigPrimAlgos/interface/HcalTimingBit.h
../../DataFormats/HcalDigi/src/HcalUpgradeTriggerPrimitiveSample.cc
../../DataFormats/HcalDigi/interface/HcalUpgradeTriggerPrimitiveSample.h 
```

## Efficiency Plots
Efficiency plots are in the directory `L1plots_eff_rates`. For running over individual samples, the bash scripts can be used
```
./Efficiency_Plots.sh
./LLPEfficiency_Plots.sh
./TOFEfficiency_Plots.sh
./ADCEfficiency_Plots.sh
# note -- for jet eff plots, need to rerun DelayedJet script with jet energy requirement of 0, details below
./JetEfficiency_Plots.sh
```
and for combined samples (`hadd rates_new_cond_*combined.root` for the mh125 3m, mh250 10m, mh350 10m, and mh1000 10m samples), the plotting scripts can be run in ROOT:
```
root
.L efficiency_ctau_combined.C++
.x efficiency_ctau_combined.C
.q

root
.L efficiency_ctau_zoom_combined.C++
.x efficiency_ctau_zoom_combined.C
.q

root
.L efficiency_TOF_combined.C++
.x efficiency_TOF_combined.C
.q

root
.L efficiency_ADC_combined.C++
.x efficiency_ADC_combined.C
.q

# note -- for jet eff plots, need to rerun DelayedJet script with jet energy requirement of 0, details below
root
.L efficiency_jetPt_combined.C++
.x efficiency_jetPt_combined.C
.q
```
Each of these outputs plots to the EOS directory as well.

Full order of how to run:
```
./RunRates_DelayedJet_PU.sh 4 4 0
hadd -f rates_new_cond_4combinedPU_jet0.root rates_new_cond_LLP_mh125_pl3000.root rates_new_cond_LLP_mh250_pl10000.root rates_new_cond_LLP_mh350_pl10000.root rates_new_cond_LLP_mh1000_pl10000.root

cd L1plots_eff_rates
./JetEfficiency_Plots.sh
root
.L efficiency_jetPt_combined.C++
.x efficiency_jetPt_combined.C
.q

cd ..

./RunRates_DelayedJet_PU.sh 4 4 40
hadd -f rates_new_cond_4combinedPU_jet40.root rates_new_cond_LLP_mh125_pl3000.root rates_new_cond_LLP_mh250_pl10000.root rates_new_cond_LLP_mh350_pl10000.root rates_new_cond_LLP_mh1000_pl10000.root

cd L1plots_eff_rates
./Efficiency_Plots.sh
./LLPEfficiency_Plots.sh
./TOFEfficiency_Plots.sh
./ADCEfficiency_Plots.sh

root
.L efficiency_ctau_combined.C++
.x efficiency_ctau_combined.C
.q

root
.L efficiency_ctau_zoom_combined.C++
.x efficiency_ctau_zoom_combined.C
.q

root
.L efficiency_TOF_combined.C++
.x efficiency_TOF_combined.C
.q

root
.L efficiency_ADC_combined.C++
.x efficiency_ADC_combined.C
.q
```
All of this is done in the script `FullPlots_DelayedJet.sh`, and this takes about 1 hour to run.

# HCAL Trigger studies for LLPs -- Depth Trigger
Towards a new L1 seed to trigger on LLP signatures with HCAL using H/E + depth 

Start with L1Ntuple framework, and setup the L1T environment as here:

https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TStage2Instructions#Environment_Setup_with_Integrati

In 11_0_X:
```
#git cms-merge-topic --unsafe georgia14:upgradeHcalTPs-l1t-110X
#git cms-merge-topic --unsafe lwang046:upgradeHcalTPs-timeingbit-106X 
git cms-merge-topic --unsafe lwang046:upgradeHcalTPs-timeingbit-110X
git clone git@github.com:cms-hcal-trigger/validation.git HcalTrigger/Validation 
```

Example of energy depth profiles of HCAL TPs is shown here: 

https://github.com/gk199/Validation/blob/TimingAndDepth/bin/rates.cxx#L1168-L1189 

H/E on L1Jets can be applied as here:

https://github.com/gk199/Validation/blob/HoE_RatesWork/bin/rates.cxx#L373-L434
