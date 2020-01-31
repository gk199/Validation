# Setup
These scripts use the L1Ntuple framework, which should be set up as described here: [L1Ntuple environment set up instructions](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TStage2Instructions#Environment_Setup_with_Integrati).

After setting up the L1Ntuple environment and compiling, issue the following:
```
git clone git@github.com:cms-hcal-trigger/Validation.git HcalTrigger/Validation
```

# HCAL conditions validation 
Scripts for HCAL radiation damage correction validation.

The script that submit CRAB jobs is called `submit_jobs.py`. Its required arguments are a dataset name (-d) and the storage site for the output (-o). For example:
```
./scripts/submit_jobs.py -d /RelValNuGun/CMSSW_10_6_1_patch1-PU_106X_mcRun3_2021_realistic_v3_rsb-v2/GEN-SIM-DIGI-RAW -o T2_CH_CERN -n NO_EXEC
```

# HCAL Rates and H/E Studies
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

## H/E Studies
`bin/rates.cxx` has the code to study the effects of adding a H/E requirement on the jet (single, double, triple, quad) rates, and can be compared against the default rates from `rates_original.cxx`. H/E code is from Matthew Citron, and can be calculated from 1x1, 3x3, or 9x9 regions of the HCAL and ECAL. Plots of 1x1 and 3x3 H/E are made, and overlayed in `draw_rates.cxx'.

## Rate Studies
In `rates.cxx' jet, egamma, tau, energy, and MET rates are calculated. New conditions rates can be calculated by adding restrictions (as is done with jet rates and H/E) and then plotting the current (default) and new rates in `draw_rates.cxx'. In addition, the ratio between new and current (default) rates will be plotted.

