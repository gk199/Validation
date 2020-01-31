# Setup
These scripts use the L1Ntuple framework, which should be set up as described here: [L1Ntuple environment set up instructions](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TStage2Instructions#Environment_Setup_with_Integrati).

To include the timing and depth branch, after compiling the L1Ntuple environment, do
```
git cms-merge-topic --unsafe georgia14:upgradeHcalTPs-l1t-106X
```
After compiling with the L1Ntuple environment and the upgrade timing and depth branch, issue the following:
```
git clone git@github.com:cms-hcal-trigger/Validation.git HcalTrigger/Validation
```
To have timing and energy information in the HCAL depth layers, note that L1Ntuples will need to be processed in this directory with the inclusion of Georgia's upgrade HCAL TPs branch. L1Ntuples are made with `ntuple_maker_def.py`.

# HCAL Rates and Multiplicity Studies
The script that submit CRAB jobs is called `submit_jobs.py`. Its required arguments are a good run lumimask, a dataset name, the new HcalL1TriggerObjects tag, and the storage site for the output. For example:
```
./submit_jobs.py -l lumimask_302472.json -d /ZeroBias/Run2017D-v1/RAW -t HcalL1TriggerObjects_2017Plan1_v13.0 -o T2_CH_CERN
```
To run on a LLP MC sample, list the sample in `bin/submit_jobs.py`, and then run with (no dataset, lumimask required):
```
./scripts/submit_jobs.py -o T2_CH_CERN -n NO_EXEC
```
This makes CMSSW config files (only need default ones for the current studies, new is used for checking HCAL conditions update). To make ntuples:
```
cmsRun ntuple_maker_def.py
mv L1Ntuple.root def_dir/L1Ntuple_def.root
```
To run on another MC sample, can simply edit its name in the `ntuple_maker_def.py' file listed under fileNames. Then compile and run the rates and plotting macros on this, referencing the correct directory where the L1Ntuple was moved to:
```
USER_CXXFLAGS="-Wno-error=unused-but-set-variable" scram b
rates.exe def def_dir/
draw_rates.exe
```
Follow this twiki on how to use the Rates scripts for more details:
[HCAL Days L1 Rates workshop](https://twiki.cern.ch/twiki/bin/view/Sandbox/L1TriggerAtHCALdays2019#HCAL_conditions_impact_at_L1_rat)

## Multiplicity Studies with Timing and Energy Information
With the inclusion of Georgia's upgrade HCAL TPs branch, HCAL timing (TDC) and energy (ACD) information are avaliable at every HCAL depth layer. Multiplicity calculations are done in `rates.cxx` for global counts and regional counts. The regional multiplicity counts are done in the region of a L1 Jet, by matching HCAL TPs within deltaR < 0.5 of the leading L1 Jet. `draw_rates.cxx` makes plots of multiplicity counts (global and regional) with a timescan, overlaying for QCD and multiple LLPs, and plots average energy and timing information per depth level, along with the ratio of energy deposited in 2 HCAL layers to total.