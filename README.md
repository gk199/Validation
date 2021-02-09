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
./RunRates_DelayedCaloObject.sh 4 4 2 2 3 4
```
where the arguments are HB TDC, HE TDC, HB energy (GeV), HE energy (GeV), prompt TP energy veto, and prompt 2x2 energy veto. Currently 2 cells in a 4x4 region make it "delayed" and a jet is rejected if there are greater than or equal to 2 prompt high energy towers in it. These last two parameters are set in `rates_delayed_cluster.cxx`.

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
