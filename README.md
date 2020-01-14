# Setup
These scripts use the L1Ntuple framework, which should be set up as described here:

[https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TStage2Instructions#Environment_Setup_with_Integrati](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TStage2Instructions#Environment_Setup_with_Integrati)

To include the timing and depth branch, after compiling the L1Ntuple environment, do
```
git cms-merge-topic --unsafe georgia14:upgradeHcalTPs-l1t-106X
```
After setting up the L1Ntuple environment, issue the following:
```
git clone git@github.com:cms-hcal-trigger/Validation.git HcalTrigger/Validation
```

# HCAL conditions validation 
Scripts for HCAL radiation damage correction validation.

The script that submit CRAB jobs is called `submit_jobs.py`. Its required arguments are a good run lumimask, a dataset name, the new HcalL1TriggerObjects tag, and the storage site for the output. For example:
```
./submit_jobs.py -l lumimask_302472.json -d /ZeroBias/Run2017D-v1/RAW -t HcalL1TriggerObjects_2017Plan1_v13.0 -o T2_CH_CERN
```
Follow this twiki on how to use the Rates scripts:
https://twiki.cern.ch/twiki/bin/view/Sandbox/L1TriggerAtHCALdays2019#HCAL_conditions_impact_at_L1_rat

To run on a LLP MC sample, list the sample in submit_jobs.py, and then run with (no dataset, lumimask required):
```
./scripts/submit_jobs.py -o T2_CH_CERN -n NO_EXEC
```
Instructions from the HCAL Days L1 Rates workshop are here: [https://twiki.cern.ch/twiki/bin/view/Sandbox/L1TriggerAtHCALdays2019]. 