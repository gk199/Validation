# Setup
These scripts use the L1Ntuple framework, which should be set up as described here: [L1Ntuple environment set up instructions](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TStage2Instructions#Environment_Setup_with_Integrati).

After setting up the L1Ntuple environment, issue the following:
```
git clone git@github.com:cms-hcal-trigger/Validation.git HcalTrigger/Validation
```

# HCAL conditions validation 
Scripts for HCAL radiation damage correction validation.

The script that submit CRAB jobs is called `submit_jobs.py`. Its required arguments are a dataset name (-d) and the storage site for the output (-o). For example:
```
./scripts/submit_jobs.py -d /RelValNuGun/CMSSW_10_6_1_patch1-PU_106X_mcRun3_2021_realistic_v3_rsb-v2/GEN-SIM-DIGI-RAW -o T2_CH_CERN -n NO_EXEC
```

Follow this twiki on how to use the Rates scripts:
[HCAL Days L1 Rates workshop](https://twiki.cern.ch/twiki/bin/view/Sandbox/L1TriggerAtHCALdays2019#HCAL_conditions_impact_at_L1_rat)

`bin/rates_withHoE.cxx` has the code to study the effects of adding a H/E requirement on the jet (single, double, triple, quad) rates, and can be compared against the default rates from `rates_original.cxx`.

