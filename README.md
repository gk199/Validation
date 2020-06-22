# Setup
These scripts use the L1Ntuple framework, which should be set up as described here: [L1Ntuple environment set up instructions](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TStage2Instructions#Environment_Setup_with_Integrati).

To include the timing and depth branch with the timing bit set (calculates the multiplicity per TP, based on cells above a time and energy value), after compiling the L1Ntuple environment, do (first option in 10_6_X)
```
git cms-merge-topic --unsafe lwang046:upgradeHcalTPs-timeingbit-106X 
```
After compiling with the L1Ntuple environment and the upgrade timing and depth branch, issue the following:
```
git clone git@github.com:cms-hcal-trigger/Validation.git HcalTrigger/Validation
```
To have timing and energy information in the HCAL depth layers, note that L1Ntuples will need to be processed in this directory with the inclusion of Georgia's / Long's upgrade HCAL TPs branch. L1Ntuples are made with `ntuple_maker_def.py`.

After the NTuple root files are made in this branch, there will be an extra timingbit branch in the root file that gives you the multiplicity of hits. This is calculated based on number of cells in a TP above time and energy thresholds, which are set in:
```
https://github.com/lwang046/cmssw/blob/32a8f7cd2f9dc50fef0ef1e48fc82883ff93bace/SimCalorimetry/HcalTrigPrimAlgos/src/HcalTimingBit.cc#L19-L24
cmssw/SimCalorimetry/HcalTrigPrimAlgos/src/HcalTimingBit.cc
```
The loop and cuts can be modified to check different multiplicity definitions.

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
With the inclusion of Georgia's upgrade HCAL TPs branch, HCAL timing (TDC) and energy (ACD) information are avaliable at every HCAL depth layer. Multiplicity calculations are done in `rates.cxx` for global counts and regional counts. The regional multiplicity counts are done either: in the region of a L1 Jet, by matching HCAL TPs within deltaR < 0.5 of the leading L1 Jet; or by associating each HCAL TP to the closest L1 jet, and then optionally adding a DR restriction as well. `draw_rates.cxx` makes plots of multiplicity counts (global and regional) with a timescan, overlaying for QCD and multiple LLPs, and plots average energy and timing information per depth level, along with the ratio of energy deposited in 2 HCAL layers to total.

If `rates.cxx` is edited (changing multiplicity counting, DR values, energy thresholds, multiplicity thresholds), the new rate NTuples (neutrino gun, QCD, LLP ctau = 10000, 1000, 500 mm) and resulting overlay plots can all be made by running
```
./RunRates.sh
```

## TMVA
TMVA is the [Toolkit for Multivariate Data Analysis with ROOT](https://root.cern.ch/root/html/guides/tmva/TMVAUsersGuide.pdf) and uses ML techniques to optimize signal and background separation. This is run with `mvaAnalysisTemplate_multiplicity.py' with
```
python mvaAnalysisTemplate_multiplicity.py
root
TMVA::TMVAGui("output.root")
root -l -b -q TMVAClassificationApplication.C
```
Which opens the GUI to see the classifier cut efficiencies. Data from ROOT files is read in to the analysis template with 'background_filename' and 'signal_filename'. TMVAClassificationApplication.C outputs a single root file, TMVApp.root with a discriminator distribution for the input file, using the trained classifier and weights files in dataset/weights.