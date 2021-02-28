// 2020: edited by Gillian Kopp for HCAL L1 LLP trigger studies, based on using a timing bit set by TP cell multiplicity
// Script for calculating rate histograms
// Originally from Aaron Bundock
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TChain.h"
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include "L1Trigger/L1TNtuples/interface/L1AnalysisEventDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoVertexDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisCaloTPDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisGeneratorDataFormat.h"
//#include "FastSimDataFormats/NuclearInteractions/interface/FSimDisplacedVertex.h"

/* TODO: put errors in rates...
creates the the rates and distributions for l1 trigger objects
How to use:
1. input the number of bunches in the run (~line 35)
2. change the variables "newConditionsNtuples" and "oldConditionsNtuples" to ntuple paths
3. If good run JSON is not applied during ntuple production, modify isGoodLumiSection()

Optionally, if you want to rescale to a given instantaneous luminosity:
1. input the instantaneous luminosity of the run (~line 32) [only if we scale to 2016 nominal]
2. select whether you rescale to L=1.5e34 (~line606??...) generally have it setup to rescale
nb: for 2&3 I have provided the info in runInfoForRates.txt
*/

// configurable parameters
double numBunch = 2556; //1537; //the number of bunches colliding for the run of interest, from Run2 settings https://home.cern/news/news/accelerators/lhc-report-full-house-lhc
double runLum = 0.02; // 0.44: 275783  0.58:  276363 //luminosity of the run of interest (*10^34)
double expectedLum = 1.15; //expected luminosity of 2016 runs (*10^34)

void rates_delayed_cluster(bool newConditions, const std::string& inputFileDirectory, double GeV_HB_variable, double GeV_HE_variable);

int main(int argc, char *argv[])
{
  bool newConditions = true;
  std::string ntuplePath("");
  double GeV_HB;
  double GeV_HE;

  if (argc != 5) {
    std::cout << "Usage: rates.exe [new/def] [path to ntuples]\n"
	      << "[new/def] indicates new or default (existing) conditions"
	      << "then values for HB and HE cell GeV requirements" << std::endl;
    exit(1);
  }
  else {
    std::string par1(argv[1]);
    std::transform(par1.begin(), par1.end(), par1.begin(), ::tolower);
    if(par1.compare("new") == 0) newConditions = true;
    else if(par1.compare("def") == 0) newConditions = false;
    else {
      std::cout << "First parameter must be \"new\" or \"def\"" << std::endl;
      exit(1);
    }
    ntuplePath = argv[2];
    GeV_HB = atof(argv[3]);
    GeV_HE = atof(argv[4]);
  }

  rates_delayed_cluster(newConditions, ntuplePath, GeV_HB, GeV_HE);

  return 0;
}

// only need to edit this section if good run JSON
// is not used during ntuple production
bool isGoodLumiSection(int lumiBlock)
{
  if (lumiBlock >= 1
      || lumiBlock <= 10000) {
    return true;
  }

  return false;
}

// functions to calculate eta from ieta, phi from iphi, delta eta, delta phi, and deltaR. Code from https://github.com/gk199/cms-hcal-debug/blob/PulseShape/plugins/HcalCompareUpgradeChains.cc#L894-L954
double etaVal(int ieta) { // calculate eta given ieta
  double etavl;
  if (ieta <= -24){
    etavl = .1695*ieta + 1.9931;
  }
  else if (ieta <= -1){
    etavl = .0875*ieta + .0489;
  }
  else if (ieta < 24){
    etavl = .0875*ieta - .0489;
  }
  else {
    etavl = .1695*ieta - 1.9931;
  }
  return etavl;
}

double phiVal(int iphi) { // calculate phi given iphi
  double phiBins=72.;
  double phivl;
  phivl=double(iphi)*(2.*TMath::Pi()/phiBins);
  if (iphi > 36) phivl -= 2.*TMath::Pi();
  return phivl;
}

void rates_delayed_cluster(bool newConditions, const std::string& inputFileDirectory, double GeV_HB_variable, double GeV_HE_variable){
  
  bool hwOn = true;   //are we using data from hardware? (upgrade trigger had to be running!!!)
  bool emuOn = true;  //are we using data from emulator?

  if (hwOn==false && emuOn==false){
    std::cout << "exiting as neither hardware or emulator selected" << std::endl;
    return;
  }

  std::string inputFile(inputFileDirectory);
  inputFile += "/L1Ntuple_*.root";
  std::string outputDirectory = "emu";  //***runNumber, triggerType, version, hw/emu/both***MAKE SURE IT EXISTS
  std::string outputFilename = "rates_def.root";
  if(newConditions) outputFilename = "rates_new_cond.root";
  TFile* kk = TFile::Open( outputFilename.c_str() , "recreate");

  // make trees
  std::cout << "Loading up the TChain..." << std::endl;
  TChain * treeL1emu = new TChain("l1UpgradeEmuTree/L1UpgradeTree");
  if (emuOn){
    treeL1emu->Add(inputFile.c_str());
  }
  TChain * treeL1hw = new TChain("l1UpgradeTree/L1UpgradeTree");
  if (hwOn){
    treeL1hw->Add(inputFile.c_str());
  }
  TChain * eventTree = new TChain("l1EventTree/L1EventTree");
  eventTree->Add(inputFile.c_str());
  TChain * treeL1TPemu = new TChain("l1CaloTowerEmuTree/L1CaloTowerTree");
  if (emuOn){
    treeL1TPemu->Add(inputFile.c_str());
  }
  TChain * treeL1TPhw = new TChain("l1CaloTowerTree/L1CaloTowerTree");
  if (hwOn){
    treeL1TPhw->Add(inputFile.c_str());
  }
  TChain * treeL1CaloTPemu = new TChain("l1CaloTowerEmuTree/L1CaloTowerTree");
  if (emuOn){
    treeL1CaloTPemu->Add(inputFile.c_str());
  }
  TChain * genTree = new TChain("l1GeneratorTree/L1GenTree");
  genTree->Add(inputFile.c_str());

  // In case you want to include PU info
  // TChain * vtxTree = new TChain("l1RecoTree/RecoTree");
  // if(binByPileUp){
  //   vtxTree->Add(inputFile.c_str());
  // }

  // These are the branches on the above trees, so can use multiple branches on same tree
  L1Analysis::L1AnalysisL1UpgradeDataFormat    *l1emu_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  treeL1emu->SetBranchAddress("L1Upgrade", &l1emu_);
  L1Analysis::L1AnalysisL1UpgradeDataFormat    *l1hw_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  treeL1hw->SetBranchAddress("L1Upgrade", &l1hw_);
  L1Analysis::L1AnalysisEventDataFormat    *event_ = new L1Analysis::L1AnalysisEventDataFormat();
  eventTree->SetBranchAddress("Event", &event_);
  L1Analysis::L1AnalysisCaloTPDataFormat    *l1TPemu_ = new L1Analysis::L1AnalysisCaloTPDataFormat();
  treeL1TPemu->SetBranchAddress("CaloTP", &l1TPemu_);
  L1Analysis::L1AnalysisCaloTPDataFormat    *l1TPhw_ = new L1Analysis::L1AnalysisCaloTPDataFormat();
  treeL1TPhw->SetBranchAddress("CaloTP", &l1TPhw_);
  L1Analysis::L1AnalysisCaloTPDataFormat *l1CaloTPemu_ = new L1Analysis::L1AnalysisCaloTPDataFormat();
  treeL1CaloTPemu->SetBranchAddress("CaloTP", &l1CaloTPemu_);
  L1Analysis::L1AnalysisGeneratorDataFormat    *generator_ = new L1Analysis::L1AnalysisGeneratorDataFormat();
  genTree->SetBranchAddress("Generator", &generator_);
  // L1Analysis::L1AnalysisRecoVertexDataFormat    *vtx_ = new L1Analysis::L1AnalysisRecoVertexDataFormat();
  // vtxTree->SetBranchAddress("Vertex", &vtx_); 

  // get number of entries
  Long64_t nentries;
  if (emuOn) nentries = treeL1emu->GetEntries();
  else nentries = treeL1hw->GetEntries();
  int goodLumiEventCount = 0;

  std::string outputTxtFilename = "output_rates/" + outputDirectory + "/extraInfo.txt";
  std::ofstream myfile; // save info about the run, including rates for a given lumi section, and number of events we used.
  myfile.open(outputTxtFilename.c_str());
  eventTree->GetEntry(0);
  myfile << "run number = " << event_->run << std::endl;

  // set parameters for histograms
  // jet bins
  int nJetBins = 400;
  float jetLo = 0.;
  float jetHi = 400.;
  float jetBinWidth = (jetHi-jetLo)/nJetBins;

  // EG bins
  int nEgBins = 300;
  float egLo = 0.;
  float egHi = 300.;
  float egBinWidth = (egHi-egLo)/nEgBins;

  // tau bins
  int nTauBins = 300;
  float tauLo = 0.;
  float tauHi = 300.;
  float tauBinWidth = (tauHi-tauLo)/nTauBins;

  // htSum bins
  int nHtSumBins = 1200;
  float htSumLo = 0.;
  float htSumHi = 1200.;
  float htSumBinWidth = (htSumHi-htSumLo)/nHtSumBins;

  // mhtSum bins
  int nMhtSumBins = 300;
  float mhtSumLo = 0.;
  float mhtSumHi = 300.;
  float mhtSumBinWidth = (mhtSumHi-mhtSumLo)/nMhtSumBins;

  // etSum bins
  int nEtSumBins = 600;
  float etSumLo = 0.;
  float etSumHi = 600.;
  float etSumBinWidth = (etSumHi-etSumLo)/nEtSumBins;

  // metSum bins
  int nMetSumBins = 300;
  float metSumLo = 0.;
  float metSumHi = 300.;
  float metSumBinWidth = (metSumHi-metSumLo)/nMetSumBins;

  // metHFSum bins
  int nMetHFSumBins = 300;
  float metHFSumLo = 0.;
  float metHFSumHi = 300.;
  float metHFSumBinWidth = (metHFSumHi-metHFSumLo)/nMetHFSumBins;

  // tp bins
  int nTpBins = 100;
  float tpLo = 0.;
  float tpHi = 100.;

  std::string axR = ";Threshold E_{T} (GeV);rate (Hz)";
  std::string axD = ";E_{T} (GeV);events/bin";

  //make histos
  TH1F* singleJetRates_emu = new TH1F("singleJetRates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* doubleJetRates_emu = new TH1F("doubleJetRates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* tripleJetRates_emu = new TH1F("tripleJetRates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* quadJetRates_emu = new TH1F("quadJetRates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* singleEgRates_emu = new TH1F("singleEgRates_emu", axR.c_str(), nEgBins, egLo, egHi);
  TH1F* doubleEgRates_emu = new TH1F("doubleEgRates_emu", axR.c_str(), nEgBins, egLo, egHi);
  TH1F* singleTauRates_emu = new TH1F("singleTauRates_emu", axR.c_str(), nTauBins, tauLo, tauHi);
  TH1F* doubleTauRates_emu = new TH1F("doubleTauRates_emu", axR.c_str(), nTauBins, tauLo, tauHi);
  TH1F* singleISOEgRates_emu = new TH1F("singleISOEgRates_emu", axR.c_str(), nEgBins, egLo, egHi);
  TH1F* doubleISOEgRates_emu = new TH1F("doubleISOEgRates_emu", axR.c_str(), nEgBins, egLo, egHi);
  TH1F* singleISOTauRates_emu = new TH1F("singleISOTauRates_emu", axR.c_str(), nTauBins, tauLo, tauHi);
  TH1F* doubleISOTauRates_emu = new TH1F("doubleISOTauRates_emu", axR.c_str(), nTauBins, tauLo, tauHi);
  TH1F* htSumRates_emu = new TH1F("htSumRates_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
  TH1F* mhtSumRates_emu = new TH1F("mhtSumRates_emu",axR.c_str(), nMhtSumBins, mhtSumLo, mhtSumHi);
  TH1F* etSumRates_emu = new TH1F("etSumRates_emu",axR.c_str(), nEtSumBins, etSumLo, etSumHi);
  TH1F* metSumRates_emu = new TH1F("metSumRates_emu",axR.c_str(), nMetSumBins, metSumLo, metSumHi); 
  TH1F* metHFSumRates_emu = new TH1F("metHFSumRates_emu",axR.c_str(), nMetHFSumBins, metHFSumLo, metHFSumHi); 
  
  TH1F* singleJetRates_hw = new TH1F("singleJetRates_hw", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* doubleJetRates_hw = new TH1F("doubleJetRates_hw", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* tripleJetRates_hw = new TH1F("tripleJetRates_hw", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* quadJetRates_hw = new TH1F("quadJetRates_hw", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* singleEgRates_hw = new TH1F("singleEgRates_hw", axR.c_str(), nEgBins, egLo, egHi);
  TH1F* doubleEgRates_hw = new TH1F("doubleEgRates_hw", axR.c_str(), nEgBins, egLo, egHi);
  TH1F* singleTauRates_hw = new TH1F("singleTauRates_hw", axR.c_str(), nTauBins, tauLo, tauHi);
  TH1F* doubleTauRates_hw = new TH1F("doubleTauRates_hw", axR.c_str(), nTauBins, tauLo, tauHi);
  TH1F* singleISOEgRates_hw = new TH1F("singleISOEgRates_hw", axR.c_str(), nEgBins, egLo, egHi);
  TH1F* doubleISOEgRates_hw = new TH1F("doubleISOEgRates_hw", axR.c_str(), nEgBins, egLo, egHi);
  TH1F* singleISOTauRates_hw = new TH1F("singleISOTauRates_hw", axR.c_str(), nTauBins, tauLo, tauHi);
  TH1F* doubleISOTauRates_hw = new TH1F("doubleISOTauRates_hw", axR.c_str(), nTauBins, tauLo, tauHi);
  TH1F* htSumRates_hw = new TH1F("htSumRates_hw",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
  TH1F* mhtSumRates_hw = new TH1F("mhtSumRates_hw",axR.c_str(), nMhtSumBins, mhtSumLo, mhtSumHi);
  TH1F* etSumRates_hw = new TH1F("etSumRates_hw",axR.c_str(), nEtSumBins, etSumLo, etSumHi);
  TH1F* metSumRates_hw = new TH1F("metSumRates_hw",axR.c_str(), nMetHFSumBins, metHFSumLo, metHFSumHi); 
  TH1F* metHFSumRates_hw = new TH1F("metHFSumRates_hw",axR.c_str(), nMetHFSumBins, metHFSumLo, metHFSumHi); 

  TH1F* hcalTP_emu = new TH1F("hcalTP_emu", ";TP E_{T}; # Entries", nTpBins, tpLo, tpHi);
  TH1F* ecalTP_emu = new TH1F("ecalTP_emu", ";TP E_{T}; # Entries", nTpBins, tpLo, tpHi);

  TH1F* hcalTP_hw = new TH1F("hcalTP_hw", ";TP E_{T}; # Entries", nTpBins, tpLo, tpHi);
  TH1F* ecalTP_hw = new TH1F("ecalTP_hw", ";TP E_{T}; # Entries", nTpBins, tpLo, tpHi);

  TH1F * htSumDistribution = new TH1F("htSumDistribution","htSum Distribution;L1_HT (GeV);Number of Events",200,0,2000); // x bins, x low, x up

  int eta_depth_tdc[30][8][51] = {{{0}}}; // increment this array element each time a cell above energy thresholds has this TDC value (=2x simulated TDC in ns to get QIE11 value)                          

  /////////////////////////////////
  // loop through all the entries//
  /////////////////////////////////
  for (Long64_t jentry=0; jentry<nentries; jentry++){
    if((jentry%10000)==0) std::cout << "Done " << jentry  << " events of " << nentries << std::endl;
    eventTree->GetEntry(jentry);
    //skip the corresponding event
    if (!isGoodLumiSection(event_->lumi)) continue;
    goodLumiEventCount++;

    //do routine for L1 emulator quantites
    if (emuOn){

      // get entries from the trees
      treeL1CaloTPemu->GetEntry(jentry);
      genTree->GetEntry(jentry); // used for gen particle matching later 
      treeL1TPemu->GetEntry(jentry);
      double tpEt(0.);
      
      for(int i=0; i < l1TPemu_->nHCALTP; i++){
	tpEt = l1TPemu_->hcalTPet[i];
	hcalTP_emu->Fill(tpEt);
      }
      for(int i=0; i < l1TPemu_->nECALTP; i++){
	tpEt = l1TPemu_->ecalTPet[i];
	ecalTP_emu->Fill(tpEt);
      }

      treeL1emu->GetEntry(jentry);
      // get jetEt*, egEt*, tauEt, htSum, mhtSum, etSum, metSum
      // ALL EMU OBJECTS HAVE BX=0...
      double jetEt_1 = 0;
      double jetEt_2 = 0;
      double jetEt_3 = 0;
      double jetEt_4 = 0;
      if (l1emu_->nJets>0) jetEt_1 = l1emu_->jetEt[0];
      if (l1emu_->nJets>1) jetEt_2 = l1emu_->jetEt[1];
      if (l1emu_->nJets>2) jetEt_3 = l1emu_->jetEt[2];
      if (l1emu_->nJets>3) jetEt_4 = l1emu_->jetEt[3];

      double egEt_1 = 0;
      double egEt_2 = 0;
      //EG pt's are not given in descending order...bx?
      for (UInt_t c=0; c<l1emu_->nEGs; c++){
        if (l1emu_->egEt[c] > egEt_1){
          egEt_2 = egEt_1;
          egEt_1 = l1emu_->egEt[c];
        }
        else if (l1emu_->egEt[c] <= egEt_1 && l1emu_->egEt[c] > egEt_2){
          egEt_2 = l1emu_->egEt[c];
        }
      }

      double tauEt_1 = 0;
      double tauEt_2 = 0;
      //tau pt's are not given in descending order
      for (UInt_t c=0; c<l1emu_->nTaus; c++){
        if (l1emu_->tauEt[c] > tauEt_1){
          tauEt_2 = tauEt_1;
          tauEt_1 = l1emu_->tauEt[c];
        }
        else if (l1emu_->tauEt[c] <= tauEt_1 && l1emu_->tauEt[c] > tauEt_2){
          tauEt_2 = l1emu_->tauEt[c];
        }
      }

      double egISOEt_1 = 0;
      double egISOEt_2 = 0;
      //EG pt's are not given in descending order...bx?
      for (UInt_t c=0; c<l1emu_->nEGs; c++){
        if (l1emu_->egEt[c] > egISOEt_1 && l1emu_->egIso[c]==1){
          egISOEt_2 = egISOEt_1;
          egISOEt_1 = l1emu_->egEt[c];
        }
        else if (l1emu_->egEt[c] <= egISOEt_1 && l1emu_->egEt[c] > egISOEt_2 && l1emu_->egIso[c]==1){
          egISOEt_2 = l1emu_->egEt[c];
        }
      }

      double tauISOEt_1 = 0;
      double tauISOEt_2 = 0;
      //tau pt's are not given in descending order
      for (UInt_t c=0; c<l1emu_->nTaus; c++){
        if (l1emu_->tauEt[c] > tauISOEt_1 && l1emu_->tauIso[c]>0){
          tauISOEt_2 = tauISOEt_1;
          tauISOEt_1 = l1emu_->tauEt[c];
        }
        else if (l1emu_->tauEt[c] <= tauISOEt_1 && l1emu_->tauEt[c] > tauISOEt_2 && l1emu_->tauIso[c]>0){
          tauISOEt_2 = l1emu_->tauEt[c];
        }
      }

      double htSum(0.0);
      double mhtSum(0.0);
      double etSum(0.0);
      double metSum(0.0);
      double metHFSum(0.0);
      for (unsigned int c=0; c<l1emu_->nSums; c++){
          if( l1emu_->sumBx[c] != 0 ) continue;
          if( l1emu_->sumType[c] == L1Analysis::kTotalEt ) etSum = l1emu_->sumEt[c];
          if( l1emu_->sumType[c] == L1Analysis::kTotalHt ) htSum = l1emu_->sumEt[c];
          if( l1emu_->sumType[c] == L1Analysis::kMissingEt ) metSum = l1emu_->sumEt[c];
	  if( l1emu_->sumType[c] == L1Analysis::kMissingEtHF ) metHFSum = l1emu_->sumEt[c];
          if( l1emu_->sumType[c] == L1Analysis::kMissingHt ) mhtSum = l1emu_->sumEt[c];
      }

      // for each bin fill according to whether our object has a larger corresponding energy
      for(int bin=0; bin<nJetBins; bin++){
        if( (jetEt_1) >= jetLo + (bin*jetBinWidth) ) singleJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
      } 

      for(int bin=0; bin<nJetBins; bin++){
        if( (jetEt_2) >= jetLo + (bin*jetBinWidth) ) doubleJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
      }  

      for(int bin=0; bin<nJetBins; bin++){
        if( (jetEt_3) >= jetLo + (bin*jetBinWidth) ) tripleJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
      }  

      for(int bin=0; bin<nJetBins; bin++){
        if( (jetEt_4) >= jetLo + (bin*jetBinWidth) ) quadJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
      }  
             
      for(int bin=0; bin<nEgBins; bin++){
        if( (egEt_1) >= egLo + (bin*egBinWidth) ) singleEgRates_emu->Fill(egLo+(bin*egBinWidth));  //GeV
      } 

      for(int bin=0; bin<nEgBins; bin++){
        if( (egEt_2) >= egLo + (bin*egBinWidth) ) doubleEgRates_emu->Fill(egLo+(bin*egBinWidth));  //GeV
      }  

      for(int bin=0; bin<nTauBins; bin++){
        if( (tauEt_1) >= tauLo + (bin*tauBinWidth) ) singleTauRates_emu->Fill(tauLo+(bin*tauBinWidth));  //GeV
      }

      for(int bin=0; bin<nTauBins; bin++){
        if( (tauEt_2) >= tauLo + (bin*tauBinWidth) ) doubleTauRates_emu->Fill(tauLo+(bin*tauBinWidth));  //GeV
      } 

      for(int bin=0; bin<nEgBins; bin++){
        if( (egISOEt_1) >= egLo + (bin*egBinWidth) ) singleISOEgRates_emu->Fill(egLo+(bin*egBinWidth));  //GeV
      } 

      for(int bin=0; bin<nEgBins; bin++){
        if( (egISOEt_2) >= egLo + (bin*egBinWidth) ) doubleISOEgRates_emu->Fill(egLo+(bin*egBinWidth));  //GeV
      }  

      for(int bin=0; bin<nTauBins; bin++){
        if( (tauISOEt_1) >= tauLo + (bin*tauBinWidth) ) singleISOTauRates_emu->Fill(tauLo+(bin*tauBinWidth));  //GeV
      }

      for(int bin=0; bin<nTauBins; bin++){
        if( (tauISOEt_2) >= tauLo + (bin*tauBinWidth) ) doubleISOTauRates_emu->Fill(tauLo+(bin*tauBinWidth));  //GeV
      } 

      for(int bin=0; bin<nHtSumBins; bin++){
        if( (htSum) >= htSumLo+(bin*htSumBinWidth)) htSumRates_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV
      }

      for(int bin=0; bin<nMhtSumBins; bin++){
        if( (mhtSum) >= mhtSumLo+(bin*mhtSumBinWidth) ) mhtSumRates_emu->Fill(mhtSumLo+(bin*mhtSumBinWidth)); //GeV           
      }

      for(int bin=0; bin<nEtSumBins; bin++){
        if( (etSum) >= etSumLo+(bin*etSumBinWidth) ) etSumRates_emu->Fill(etSumLo+(bin*etSumBinWidth)); //GeV           
      }

      for(int bin=0; bin<nMetSumBins; bin++){
        if( (metSum) >= metSumLo+(bin*metSumBinWidth) ) metSumRates_emu->Fill(metSumLo+(bin*metSumBinWidth)); //GeV           
      }
      for(int bin=0; bin<nMetHFSumBins; bin++){
        if( (metHFSum) >= metHFSumLo+(bin*metHFSumBinWidth) ) metHFSumRates_emu->Fill(metHFSumLo+(bin*metHFSumBinWidth)); //GeV           
      }


      //////////////////////////////////////
      ////////// HCAL TP Loop //////////
      //////////////////////////////////////
      double nCaloTPemu = l1CaloTPemu_->nHCALTP; // number of TPs varies from 400-1400 per event, approximately Gaussian
      for (int HcalTPIt = 0; HcalTPIt < nCaloTPemu; HcalTPIt++){ // loop over HCAL TPs
	int tpEtaemu = l1CaloTPemu_->hcalTPieta[HcalTPIt]; // ieta
	//	double tpPhiemu = l1CaloTPemu_->hcalTPCaliphi[HcalTPIt]; // iphi

	if (abs(tpEtaemu) > 29 ) continue; // 16 just HB. don't consider HCAL TPs outside of HB or HE -- note, 29 should be summed with 28 still presumably 
	// End of HB is eta=1.479m, end of HE is eta=3 from here http://www.hephy.at/user/friedl/diss/html/node8.html.

	if (abs(tpEtaemu) > 16) {
	  if ( l1CaloTPemu_->hcalTPDepth1[HcalTPIt] >= GeV_HE_variable && l1CaloTPemu_->hcalTPtiming1[HcalTPIt] >= 0) eta_depth_tdc[abs(tpEtaemu)][1][static_cast<int>(l1CaloTPemu_->hcalTPtiming1[HcalTPIt] * 2 + 0.5)] += 1;
	  if ( l1CaloTPemu_->hcalTPDepth2[HcalTPIt] >= GeV_HE_variable && l1CaloTPemu_->hcalTPtiming2[HcalTPIt] >= 0) eta_depth_tdc[abs(tpEtaemu)][2][static_cast<int>(l1CaloTPemu_->hcalTPtiming2[HcalTPIt] * 2 + 0.5)] += 1;
	  if ( l1CaloTPemu_->hcalTPDepth3[HcalTPIt] >= GeV_HE_variable && l1CaloTPemu_->hcalTPtiming3[HcalTPIt] >= 0) eta_depth_tdc[abs(tpEtaemu)][3][static_cast<int>(l1CaloTPemu_->hcalTPtiming3[HcalTPIt] * 2 + 0.5)] += 1;
	  if ( l1CaloTPemu_->hcalTPDepth4[HcalTPIt] >= GeV_HE_variable && l1CaloTPemu_->hcalTPtiming4[HcalTPIt] >= 0) eta_depth_tdc[abs(tpEtaemu)][4][static_cast<int>(l1CaloTPemu_->hcalTPtiming4[HcalTPIt] * 2 + 0.5)] += 1;
	  if ( l1CaloTPemu_->hcalTPDepth5[HcalTPIt] >= GeV_HE_variable && l1CaloTPemu_->hcalTPtiming5[HcalTPIt] >= 0) eta_depth_tdc[abs(tpEtaemu)][5][static_cast<int>(l1CaloTPemu_->hcalTPtiming5[HcalTPIt] * 2 + 0.5)] += 1;
	  if ( l1CaloTPemu_->hcalTPDepth6[HcalTPIt] >= GeV_HE_variable && l1CaloTPemu_->hcalTPtiming6[HcalTPIt] >= 0) eta_depth_tdc[abs(tpEtaemu)][6][static_cast<int>(l1CaloTPemu_->hcalTPtiming6[HcalTPIt] * 2 + 0.5)] += 1;
	  if ( l1CaloTPemu_->hcalTPDepth7[HcalTPIt] >= GeV_HE_variable && l1CaloTPemu_->hcalTPtiming7[HcalTPIt] >= 0) eta_depth_tdc[abs(tpEtaemu)][7][static_cast<int>(l1CaloTPemu_->hcalTPtiming7[HcalTPIt] * 2 + 0.5)] += 1;
	}

	if (abs(tpEtaemu) <= 16) { 
	  if ( l1CaloTPemu_->hcalTPDepth1[HcalTPIt] >= GeV_HB_variable && l1CaloTPemu_->hcalTPtiming1[HcalTPIt] >= 0) eta_depth_tdc[abs(tpEtaemu)][1][static_cast<int>(l1CaloTPemu_->hcalTPtiming1[HcalTPIt] * 2 + 0.5)] += 1;
	  if ( l1CaloTPemu_->hcalTPDepth2[HcalTPIt] >= GeV_HB_variable && l1CaloTPemu_->hcalTPtiming2[HcalTPIt] >= 0) eta_depth_tdc[abs(tpEtaemu)][2][static_cast<int>(l1CaloTPemu_->hcalTPtiming2[HcalTPIt] * 2 + 0.5)] += 1;
	  if ( l1CaloTPemu_->hcalTPDepth3[HcalTPIt] >= GeV_HB_variable && l1CaloTPemu_->hcalTPtiming3[HcalTPIt] >= 0) eta_depth_tdc[abs(tpEtaemu)][3][static_cast<int>(l1CaloTPemu_->hcalTPtiming3[HcalTPIt] * 2 + 0.5)] += 1;
	  if ( l1CaloTPemu_->hcalTPDepth4[HcalTPIt] >= GeV_HB_variable && l1CaloTPemu_->hcalTPtiming4[HcalTPIt] >= 0) eta_depth_tdc[abs(tpEtaemu)][4][static_cast<int>(l1CaloTPemu_->hcalTPtiming4[HcalTPIt] * 2 + 0.5)] += 1;
	}
	
      }

      htSumDistribution->Fill(htSum);

    }// closes if 'emuOn' is true


    //do routine for L1 hardware quantities
    if (hwOn){

    }// closes if 'hwOn' is true

  }// closes loop through events


  //  TFile g( outputFilename.c_str() , "new");
  kk->cd();
  // normalisation factor for rate histograms (11kHz is the orbit frequency)
  double norm = 11246*(numBunch/goodLumiEventCount); // no lumi rescale
  //  double norm = 11246*(numBunch/goodLumiEventCount)*(expectedLum/runLum); //scale to nominal lumi

  if (emuOn){
    singleJetRates_emu->Scale(norm);
    doubleJetRates_emu->Scale(norm);
    tripleJetRates_emu->Scale(norm);
    quadJetRates_emu->Scale(norm);
    singleEgRates_emu->Scale(norm);
    doubleEgRates_emu->Scale(norm);
    singleTauRates_emu->Scale(norm);
    doubleTauRates_emu->Scale(norm);
    singleISOEgRates_emu->Scale(norm);
    doubleISOEgRates_emu->Scale(norm);
    singleISOTauRates_emu->Scale(norm);
    doubleISOTauRates_emu->Scale(norm);
    htSumRates_emu->Scale(norm);
    mhtSumRates_emu->Scale(norm);
    etSumRates_emu->Scale(norm);
    metSumRates_emu->Scale(norm);
    metHFSumRates_emu->Scale(norm);

    //set the errors for the rates
    //want error -> error * sqrt(norm) ?

    hcalTP_emu->Write();
    ecalTP_emu->Write();
    singleJetRates_emu->Write();
    doubleJetRates_emu->Write();
    tripleJetRates_emu->Write();
    quadJetRates_emu->Write();
    singleEgRates_emu->Write();
    doubleEgRates_emu->Write();
    singleTauRates_emu->Write();
    doubleTauRates_emu->Write();
    singleISOEgRates_emu->Write();
    doubleISOEgRates_emu->Write();
    singleISOTauRates_emu->Write();
    doubleISOTauRates_emu->Write();
    htSumRates_emu->Write();
    mhtSumRates_emu->Write();
    etSumRates_emu->Write();
    metSumRates_emu->Write();
    metHFSumRates_emu->Write();

    htSumDistribution->Write();
  }

  if (hwOn){

    singleJetRates_hw->Scale(norm);
    doubleJetRates_hw->Scale(norm);
    tripleJetRates_hw->Scale(norm);
    quadJetRates_hw->Scale(norm);
    singleEgRates_hw->Scale(norm);
    doubleEgRates_hw->Scale(norm);
    singleTauRates_hw->Scale(norm);
    doubleTauRates_hw->Scale(norm);
    singleISOEgRates_hw->Scale(norm);
    doubleISOEgRates_hw->Scale(norm);
    singleISOTauRates_hw->Scale(norm);
    doubleISOTauRates_hw->Scale(norm);
    htSumRates_hw->Scale(norm);
    mhtSumRates_hw->Scale(norm);
    etSumRates_hw->Scale(norm);
    metSumRates_hw->Scale(norm);
    metHFSumRates_hw->Scale(norm);

    hcalTP_hw->Write();
    ecalTP_hw->Write();
    singleJetRates_hw->Write();
    doubleJetRates_hw->Write();
    tripleJetRates_hw->Write();
    quadJetRates_hw->Write();
    singleEgRates_hw->Write();
    doubleEgRates_hw->Write();
    singleTauRates_hw->Write();
    doubleTauRates_hw->Write();
    singleISOEgRates_hw->Write();
    doubleISOEgRates_hw->Write();
    singleISOTauRates_hw->Write();
    doubleISOTauRates_hw->Write();
    htSumRates_hw->Write();
    mhtSumRates_hw->Write();
    etSumRates_hw->Write();
    metSumRates_hw->Write();
    metHFSumRates_hw->Write();
  }

  // calculate partial sum / cumulative sum from the TDC distributions
  double sum[30][8] = {{0}}; // indexed by eta, depth
  double partial_sum[30][8][51] = {{{0}}}; // counts partial sum up to each tdc value
  double cumulative_frac[30][8][51] = {{{0}}};

  for (int eta = 0; eta < 30; eta++) {
    for (int depth = 0; depth < 8; depth++) {
      for (int tdc = 0; tdc<=50; tdc++) {
	sum[eta][depth] += eta_depth_tdc[eta][depth][tdc];
	partial_sum[eta][depth][tdc] = sum[eta][depth];
      }
    }
  }
  std::cout << partial_sum[0][1][50] << std::endl;

  int bkg50[30][8] = {{0}};
  int bkg60[30][8] = {{0}};
  int bkg70[30][8] = {{0}};
  int bkg80[30][8] = {{0}};
  int bkg80_delayed[30][8] = {{0}};
  int bkg90[30][8] = {{0}};
  int bkg90_delayed[30][8] = {{0}};
  int bkg95[30][8] = {{0}};
  for (int eta = 0; eta < 30; eta++) {
    for (int depth = 0; depth < 8; depth++) {
      for (int tdc = 0; tdc<=50; tdc++) {
	cumulative_frac[eta][depth][tdc] = (double) ( partial_sum[eta][depth][tdc] / partial_sum[eta][depth][50] );
	if (cumulative_frac[eta][depth][tdc] > 0.5 && bkg50[eta][depth] == 0 ) bkg50[eta][depth] = tdc;
	if (cumulative_frac[eta][depth][tdc] > 0.6 && bkg60[eta][depth] == 0 ) bkg60[eta][depth] = tdc;
	if (cumulative_frac[eta][depth][tdc] > 0.7 && bkg70[eta][depth] == 0 ) bkg70[eta][depth] = tdc;
	if (cumulative_frac[eta][depth][tdc] > 0.8 && bkg80[eta][depth] == 0 ) bkg80[eta][depth] = tdc;
        if (cumulative_frac[eta][depth][tdc] > 0.8 && bkg80_delayed[eta][depth] == 0 ) bkg80_delayed[eta][depth] = tdc+2;
	if (cumulative_frac[eta][depth][tdc] > 0.9 && bkg90[eta][depth] == 0 ) bkg90[eta][depth] = tdc;
        if (cumulative_frac[eta][depth][tdc] > 0.9 && bkg90_delayed[eta][depth] == 0 ) bkg90_delayed[eta][depth] = tdc+2;
	//        if (eta == 15 && depth == 1) std::cout << cumulative_frac[eta][depth][tdc] << " depth 1 at TDC = " << tdc << " with ps = " << partial_sum[eta][depth][tdc] << std::endl;
	//        if (eta == 15 && depth == 2) std::cout << cumulative_frac[eta][depth][tdc] << " depth 2 at TDC = " << tdc << " with ps = " << partial_sum[eta][depth][tdc] << std::endl;
	//        if (eta == 15 && depth == 3) std::cout << cumulative_frac[eta][depth][tdc] << " depth 3 at TDC = " << tdc << " with ps = " << partial_sum[eta][depth][tdc] << std::endl;
	//	if (eta == 15 && depth == 4) std::cout << cumulative_frac[eta][depth][tdc] << " depth 4 at TDC = " << tdc << " with ps = " << partial_sum[eta][depth][tdc] << std::endl;
	if (cumulative_frac[eta][depth][tdc] > 0.95 && bkg95[eta][depth] == 0 ) bkg95[eta][depth] = tdc; // tdc value where 95% bkg below
      }
    }
  }

  // for determining avg over HB, HE depths
  double partialHB_sum[8][51] = {{0}};
  double partialHB_sum_smoothed[17][8][51] = {{{0}}};
  double partialHE_sum[8][51] = {{0}};
  double cumulativeHB_frac[8][51] = {{0}};
  double cumulativeHB_smoothed_frac[17][8][51] = {{{0}}};
  double cumulativeHE_frac[8][51] = {{0}};
  for (int eta = 0; eta < 30; eta++) {
    for (int depth = 0; depth < 8; depth++) {
      for (int tdc = 0; tdc<=50; tdc++) {
        if ( (eta < 16) || (eta == 16 && depth < 4) ) {
	  if (eta%2 == 1 && (eta < 13 || depth < 4) ) {
	    partialHB_sum_smoothed[eta][depth][tdc] += partial_sum[eta][depth][tdc] + partial_sum[eta+1][depth][tdc]; // combine depth 4 in 2x2 groups, 13+14+15
	    partialHB_sum_smoothed[eta+1][depth][tdc] += partial_sum[eta][depth][tdc] + partial_sum[eta+1][depth][tdc];
	  }
	  if (eta == 13 && depth == 4) {
	    partialHB_sum_smoothed[eta][depth][tdc] += partial_sum[eta][depth][tdc] + partial_sum[eta+1][depth][tdc] + partial_sum[eta+2][depth][tdc];
            partialHB_sum_smoothed[eta+1][depth][tdc] += partial_sum[eta][depth][tdc] + partial_sum[eta+1][depth][tdc] + partial_sum[eta+2][depth][tdc];
            partialHB_sum_smoothed[eta+2][depth][tdc] += partial_sum[eta][depth][tdc] + partial_sum[eta+1][depth][tdc] + partial_sum[eta+2][depth][tdc];
	  }
	  partialHB_sum[depth][tdc] += partial_sum[eta][depth][tdc];
	}
        if ( (eta > 16) || (eta == 16 && depth == 4) ) partialHE_sum[depth][tdc] += partial_sum[eta][depth][tdc];
      }
    }
  }
  int bkg90_depthHE[8] = {{0}};
  int bkg90_depthHB[8] = {{0}};
  int bkg90_depthHB_smoothed[17][8] = {{0}};
  for (int depth = 0; depth < 8; depth++) {
    for (int tdc = 0; tdc<=50; tdc++) {
      cumulativeHB_frac[depth][tdc] = (double) ( partialHB_sum[depth][tdc] / partialHB_sum[depth][50] );
      for (int eta = 1; eta < 17; eta++) cumulativeHB_smoothed_frac[eta][depth][tdc] = (double) ( partialHB_sum_smoothed[eta][depth][tdc] / partialHB_sum_smoothed[eta][depth][50] );
      cumulativeHE_frac[depth][tdc] = (double) ( partialHE_sum[depth][tdc] / partialHE_sum[depth][50] );
      if (cumulativeHB_frac[depth][tdc] > 0.9 && bkg90_depthHB[depth] == 0 ) bkg90_depthHB[depth] = tdc;
      if (cumulativeHE_frac[depth][tdc] > 0.9 && bkg90_depthHE[depth] == 0 ) bkg90_depthHE[depth] = tdc;
      for (int eta = 1; eta < 17; eta++) if (cumulativeHB_smoothed_frac[eta][depth][tdc] > 0.9 && bkg90_depthHB_smoothed[eta][depth] == 0) bkg90_depthHB_smoothed[eta][depth] = tdc;
    }
    std::cout << "Average values per depth, HB and HE: " << std::endl;
    std::cout << "HE Depth = " << depth << " and TDC value = " << bkg90_depthHE[depth] << std::endl;
    std::cout << "HB Depth = " << depth << " and TDC value = " << bkg90_depthHB[depth] <<std::endl;
  }

  //  for (int tdc = 0; tdc <= 50; tdc ++) std::cout << partial_sum[1][1][tdc] / partial_sum[1][1][50] << " partial sum tdc / partial sum 50 , and then the bkg95 = " << bkg95[1][1] << std::endl;

  // background efficiencies 
  if (inputFile.substr(0,3) == "QCD" ) {
    std::ofstream TDCdistribution_Background;
    std::ofstream TDCdistribution_Background95;
    std::ofstream TDCdistribution_Background90;
    std::ofstream TDCdistribution_Background90_delayed;
    std::ofstream TDCdistribution_Background90smooth;
    std::ofstream TDCdistribution_Background80;
    std::ofstream TDCdistribution_Background80_delayed;
    std::ofstream TDCdistribution_Background70;
    std::ofstream TDCdistribution_Background60;
    std::ofstream TDCdistribution_Background50;

    TDCdistribution_Background.open(Form("TDCdistribution_Background_%s.txt",inputFile.substr(0,3).c_str()));
    TDCdistribution_Background95.open(Form("TDCdistribution_Background95_%s.txt",inputFile.substr(0,3).c_str()));
    TDCdistribution_Background90.open(Form("TDCdistribution_Background90_%s.txt",inputFile.substr(0,3).c_str()));
    TDCdistribution_Background90_delayed.open(Form("TDCdistribution_Background90_delayed_%s.txt",inputFile.substr(0,3).c_str()));
    TDCdistribution_Background80.open(Form("TDCdistribution_Background80_%s.txt",inputFile.substr(0,3).c_str()));
    TDCdistribution_Background80_delayed.open(Form("TDCdistribution_Background80_delayed_%s.txt",inputFile.substr(0,3).c_str()));
    TDCdistribution_Background70.open(Form("TDCdistribution_Background70_%s.txt",inputFile.substr(0,3).c_str()));
    TDCdistribution_Background60.open(Form("TDCdistribution_Background60_%s.txt",inputFile.substr(0,3).c_str()));
    TDCdistribution_Background50.open(Form("TDCdistribution_Background50_%s.txt",inputFile.substr(0,3).c_str()));
    TDCdistribution_Background90smooth.open(Form("TDCdistribution_Background90_%s_smoothed.txt",inputFile.substr(0,3).c_str()));

    for (int eta = 1; eta < 30; eta++) {
      for (int depth = 1; depth < 8; depth++) {
	TDCdistribution_Background << "ieta = " << eta << ", depth = " << depth << ",      TDC50 = " << bkg50[eta][depth] << ", TDC60 = " << bkg60[eta][depth] << ", TDC70 = " << bkg70[eta][depth] << ", TDC80 = " << bkg80[eta][depth] << ", TDC90 = " << bkg90[eta][depth] << ", TDC95 = " << bkg95[eta][depth] << std::endl; 
	//        TDCdistribution_Background95 << "ieta = " << eta << ", depth = " << depth << ", TDC95 = " << bkg95[eta][depth] << std::endl;
      }

      TDCdistribution_Background95 << eta << ", " << bkg95[eta][1] << ", " << bkg95[eta][2] << ", " << bkg95[eta][3] << ", " << bkg95[eta][4] << ", " << bkg95[eta][5] << ", " << bkg95[eta][6] << ", " << bkg95[eta][7] << std::endl;
      TDCdistribution_Background90 << eta << ", " << bkg90[eta][1] << ", " << bkg90[eta][2] << ", " << bkg90[eta][3] << ", " << bkg90[eta][4] << ", " << bkg90[eta][5] << ", " << bkg90[eta][6] << ", " << bkg90[eta][7] << std::endl;
      TDCdistribution_Background90_delayed << eta << ", " << bkg90_delayed[eta][1] << ", " << bkg90_delayed[eta][2] << ", " << bkg90_delayed[eta][3] << ", " << bkg90_delayed[eta][4] << ", " << bkg90_delayed[eta][5] << ", " << bkg90_delayed[eta][6] << ", " << bkg90_delayed[eta][7] << std::endl;
      if (eta < 16) TDCdistribution_Background90smooth << eta << ", " << bkg90_depthHB_smoothed[eta][1] << ", " << bkg90_depthHB_smoothed[eta][2] << ", " << bkg90_depthHB_smoothed[eta][3] << ", " << bkg90_depthHB_smoothed[eta][4] << ", " << bkg90_depthHB_smoothed[eta][5] << ", " << bkg90_depthHB_smoothed[eta][6] << ", " << bkg90_depthHB_smoothed[eta][7] << std::endl;
      if (eta == 16) TDCdistribution_Background90smooth << eta << ", " << bkg90_depthHB_smoothed[eta][1] << ", " << bkg90_depthHB_smoothed[eta][2] << ", " << bkg90_depthHB_smoothed[eta][3] << ", " << bkg90[eta][4] << ", " << bkg90[eta][5] << ", " << bkg90[eta][6] << ", " << bkg90[eta][7] << std::endl;
      if (eta > 16) TDCdistribution_Background90smooth << eta << ", " << bkg90[eta][1] << ", " << bkg90[eta][2] << ", " << bkg90[eta][3] << ", " << bkg90[eta][4] << ", " << bkg90[eta][5] << ", " << bkg90[eta][6] << ", " << bkg90[eta][7] << std::endl;
      TDCdistribution_Background80 << eta << ", " << bkg80[eta][1] << ", " << bkg80[eta][2] << ", " << bkg80[eta][3] << ", " << bkg80[eta][4] << ", " << bkg80[eta][5] << ", " << bkg80[eta][6] << ", " << bkg80[eta][7] << std::endl;
      TDCdistribution_Background80_delayed << eta << ", " << bkg80_delayed[eta][1] << ", " << bkg80_delayed[eta][2] << ", " << bkg80_delayed[eta][3] << ", " << bkg80_delayed[eta][4] << ", " << bkg80_delayed[eta][5] << ", " << bkg80_delayed[eta][6] << ", " << bkg80_delayed[eta][7] << std::endl;
      TDCdistribution_Background70 << eta << ", " << bkg70[eta][1] << ", " << bkg70[eta][2] << ", " << bkg70[eta][3] << ", " << bkg70[eta][4] << ", " << bkg70[eta][5] << ", " << bkg70[eta][6] << ", " << bkg70[eta][7] << std::endl;
      TDCdistribution_Background60 << eta << ", " << bkg60[eta][1] << ", " << bkg60[eta][2] << ", " << bkg60[eta][3] << ", " << bkg60[eta][4] << ", " << bkg60[eta][5] << ", " << bkg60[eta][6] << ", " << bkg60[eta][7] << std::endl;
      TDCdistribution_Background50 << eta << ", " << bkg50[eta][1] << ", " << bkg50[eta][2] << ", " << bkg50[eta][3] << ", " << bkg50[eta][4] << ", " << bkg50[eta][5] << ", " << bkg50[eta][6] << ", " << bkg50[eta][7] << std::endl;
    }
    TDCdistribution_Background.close();
    TDCdistribution_Background95.close();
    TDCdistribution_Background90.close();
    TDCdistribution_Background90_delayed.close();
    TDCdistribution_Background90smooth.close();
    TDCdistribution_Background80.close();
    TDCdistribution_Background80_delayed.close();
    TDCdistribution_Background70.close();
    TDCdistribution_Background60.close();
    TDCdistribution_Background50.close();
  }

  else {
    std::ofstream TDCdistribution_Background;
    std::ofstream TDCdistribution_Background90;
    TDCdistribution_Background.open(Form("TDCdistribution_Background_%s.txt",inputFile.substr(0,15).c_str()));
    TDCdistribution_Background90.open(Form("TDCdistribution_Background90_%s.txt",inputFile.substr(0,15).c_str()));

    for (int eta = 1; eta < 30; eta++) {
      for (int depth = 1; depth < 8; depth++) {
        TDCdistribution_Background << "ieta = " << eta << ", depth = " << depth << ",      TDC50 = " << bkg50[eta][depth] << ", TDC60 = " << bkg60[eta][depth] << ", TDC70 = " << bkg70[eta][depth] << ", TDC80 = " << bkg80[eta][depth] << ", TDC90 = " << bkg90[eta][depth] << ", TDC95 = " << bkg95[eta][depth] << std::endl;
        TDCdistribution_Background90 << "ieta = " << eta << ", depth = " << depth << ", TDC90 = " << bkg90[eta][depth] << std::endl;
      }
    }
    TDCdistribution_Background.close();
    TDCdistribution_Background90.close();
  }

  myfile << "using the following ntuple: " << inputFile << std::endl;
  myfile << "number of colliding bunches = " << numBunch << std::endl;
  myfile << "run luminosity = " << runLum << std::endl;
  myfile << "expected luminosity = " << expectedLum << std::endl;
  myfile << "norm factor used = " << norm << std::endl;
  myfile << "number of good events = " << goodLumiEventCount << std::endl;
  myfile.close(); 
}//closes the function 'rates'
