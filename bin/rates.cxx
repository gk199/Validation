// Script for calculating rate histograms
// Originally from Aaron Bundock
// Edited by Gillian Kopp for multiplicity studies for LLP L1 trigger using HCAL depth and timing (2020)
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include <iostream>
#include <fstream>
#include <string>
#include "TLorentzVector.h"
#include "Math/LorentzVector.h" 
#include "L1Trigger/L1TNtuples/interface/L1AnalysisEventDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoVertexDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisCaloTPDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1CaloTowerDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisGeneratorDataFormat.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

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
double numBunch = 2556; //1537; //the number of bunches colliding for the run of interest -- comparing rates to Run II conditions, so use the number of bunches at end of Run II
double runLum = 0.02; // 0.44: 275783  0.58:  276363 //luminosity of the run of interest (*10^34)
double expectedLum = 1.15; //expected luminosity of 2016 runs (*10^34)

void rates(bool newConditions, const std::string& inputFileDirectory);

int main(int argc, char *argv[])
{
  bool newConditions = true;
  std::string ntuplePath("");

  if (argc != 3) {
    std::cout << "Usage: rates.exe [new/def] [path to ntuples]\n"
	      << "[new/def] indicates new or default (existing) conditions" << std::endl;
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
  }

  rates(newConditions, ntuplePath);

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
double deltaPhi(double phi1, double phi2) { 
  double result = phi1 - phi2;
  if(fabs(result) > 9999) return result;
  while (result > TMath::Pi()) result -= 2*TMath::Pi();
  while (result <= -TMath::Pi()) result += 2*TMath::Pi();
  return result;
}

double deltaR(double eta1, double phi1, double eta2, double phi2) {
  double deta = eta1 - eta2;
  double dphi = deltaPhi(phi1, phi2);
  return sqrt(deta*deta + dphi*dphi);
}

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

// sorting Lorentz vector
struct ptsort: public std::binary_function<LorentzVector, LorentzVector, bool> 
{
  bool operator () (const LorentzVector & x, const LorentzVector & y) 
  { return  ( x.pt() > y.pt() ) ; }
};

// this is a correction function for the eta phi value of the intersection of a particle not resulting from the IP       
// used since the b quark is matched to the TP based on eta phi, but the b quark results from a LLP decay so has a different vertex 
// from HcalCompareUpgradeChains intersect function https://github.com/gk199/cms-hcal-debug/blob/PulseShape/plugins/HcalCompareUpgradeChains.cc#L961
std::vector<double> intersect(double vx, double vy,double vz, double px, double py, double pz) {
  double lightSpeed = 29979245800;
  double radius = 179; // 130 for calorimeters (ECAL + HCAL)
  double length = 388; // 300 for calorimeters (ECAL + HCAL)
  double energy = sqrt(px*px + py*py + pz*pz);
  // First work out intersection with cylinder (barrel)        
  double a = (px*px + py*py)*lightSpeed*lightSpeed/(energy*energy);
  double b = 2*(vx*px + vy*py)*lightSpeed/energy;
  double c = (vx*vx + vy*vy) - radius*radius;
  double sqrt_disc = sqrt(b*b - 4*a*c);
  double tCircle1 = (-b + sqrt_disc)/(2*a);
  double tCircle2 = (-b - sqrt_disc)/(2*a);
  // If intersection in the past it doesn't count         
  if (tCircle1 < 0) tCircle1 = 1E9;
  if (tCircle2 < 0) tCircle2 = 1E9;
  // If the intsersection occurs outside the barrel length it doesn't count                       
  double zPosCircle1 = tCircle1*(pz/energy)*lightSpeed + vz;
  double zPosCircle2 = tCircle2*(pz/energy)*lightSpeed + vz;
  if (zPosCircle1 > length) tCircle1 = 1E9;
  if (zPosCircle2 > length) tCircle2 = 1E9;
  // Now work out if it intersects the endcap                      
  double tPlane1 = (length-vz)*energy/(pz*lightSpeed);
  double tPlane2 = (-length-vz)*energy/(pz*lightSpeed);
  // If intersection in the past it doesn't count                     
  if (tPlane1 < 0) tPlane1 = 1E9;
  if (tPlane2 < 0) tPlane2 = 1E9;
  double xPosPlane1 = tPlane1*(px/energy)*lightSpeed + vx;
  double yPosPlane1 = tPlane1*(py/energy)*lightSpeed + vy;
  double xPosPlane2 = tPlane2*(px/energy)*lightSpeed + vx;
  double yPosPlane2 = tPlane2*(py/energy)*lightSpeed + vy;
  // If the intsersection occurs outside the endcap radius it doesn't count     
  if (sqrt(xPosPlane1*xPosPlane1 + yPosPlane1*yPosPlane1) > radius) tPlane1 = 1E9;
  if (sqrt(xPosPlane2*xPosPlane2+yPosPlane2*yPosPlane2) > radius) tPlane2 = 1E9;
  // Find the first intersection                          
  double tInter = std::min({tCircle1,tCircle2,tPlane1,tPlane2});
  // Return 1000,1000 if not intersection with barrel or endcap             
  std::vector<double> etaphi;
  if (tInter > 1E6)
    {
      etaphi.push_back(1000);
      etaphi.push_back(1000);
      return etaphi;
    }
  // Find position of intersection                          
  double xPos = tInter*(px/energy)*lightSpeed + vx;
  double yPos = tInter*(py/energy)*lightSpeed + vy;
  double zPos = tInter*(pz/energy)*lightSpeed + vz;
  // Find eta/phi of intersection                          
  double phi = atan2(yPos,xPos); // return the arc tan in radians                                                                                                                               
  double theta = acos(zPos/sqrt(xPos*xPos + yPos*yPos + zPos*zPos));
  double eta = -log(tan(theta/2.));
  etaphi.push_back(eta);
  etaphi.push_back(phi);
  return etaphi;
}
  

/*
// matching L1 jets with b quark decays from LLP
bool isMatched(L1JetObject &myL1Jet) {
}

// seeing if LLP decays within HCAL volume
bool isWithinHCALVolume(L1Jet &myJet) {
}
*/

void rates(bool newConditions, const std::string& inputFileDirectory){
  
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
  // if (kk!=0){
  //   cout << "TERMINATE: not going to overwrite file " << outputFilename << endl;
  //   return;
  // }

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
  TChain * treeL1Towemu = new TChain("l1CaloTowerEmuTree/L1CaloTowerTree");
  if (emuOn){
    treeL1Towemu->Add(inputFile.c_str());
  }
  TChain * treeL1CaloTPemu = new TChain("l1CaloTowerEmuTree/L1CaloTowerTree");
  if (emuOn){
    treeL1CaloTPemu->Add(inputFile.c_str());
  }
  TChain * treeL1CaloTPhw = new TChain("l1CaloTowerTree/L1CaloTowerTree");
  if (hwOn){
    treeL1CaloTPhw->Add(inputFile.c_str());
  }
  TChain * eventTree = new TChain("l1EventTree/L1EventTree");
  eventTree->Add(inputFile.c_str());
  TChain * genTree = new TChain("l1GeneratorTree/L1GenTree");
  genTree->Add(inputFile.c_str());

  // In case you want to include PU info
  // TChain * vtxTree = new TChain("l1RecoTree/RecoTree");
  // if(binByPileUp){
  //   vtxTree->Add(inputFile.c_str());
  // }


  L1Analysis::L1AnalysisL1UpgradeDataFormat    *l1emu_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  treeL1emu->SetBranchAddress("L1Upgrade", &l1emu_);
  L1Analysis::L1AnalysisL1UpgradeDataFormat    *l1hw_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  treeL1hw->SetBranchAddress("L1Upgrade", &l1hw_);
  L1Analysis::L1AnalysisL1CaloTowerDataFormat     *l1Towemu_ = new L1Analysis::L1AnalysisL1CaloTowerDataFormat();
  treeL1Towemu->SetBranchAddress("L1CaloTower", &l1Towemu_);
  L1Analysis::L1AnalysisCaloTPDataFormat     *l1CaloTPemu_ = new L1Analysis::L1AnalysisCaloTPDataFormat();
  treeL1CaloTPemu->SetBranchAddress("CaloTP", &l1CaloTPemu_);
  L1Analysis::L1AnalysisEventDataFormat    *event_ = new L1Analysis::L1AnalysisEventDataFormat();
  eventTree->SetBranchAddress("Event", &event_);
  L1Analysis::L1AnalysisCaloTPDataFormat    *l1CaloTPhw_ = new L1Analysis::L1AnalysisCaloTPDataFormat();
  treeL1CaloTPhw->SetBranchAddress("CaloTP", &l1CaloTPhw_);
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
  int nHtSumBins = 1200; // 600
  float htSumLo = 0.;
  float htSumHi = 1200.;
  float htSumBinWidth = (htSumHi-htSumLo)/nHtSumBins;

  // mhtSum bins
  int nMhtSumBins = 300;
  float mhtSumLo = 0.;
  float mhtSumHi = 300.;
  float mhtSumBinWidth = (mhtSumHi-mhtSumLo)/nMhtSumBins;

  // etSum bins
  int nEtSumBins = 1200; // 600
  float etSumLo = 0.;
  float etSumHi = 1200.;
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
  int nTpBins = 130;
  float tpLo = 0.;
  float tpHi = 130.;

  std::string axR = ";Threshold E_{T} (GeV);rate (Hz)";
  std::string axD = ";E_{T} (GeV);events/bin";
  std::string mult = ";Hit Multiplicity;Number of Entries";

  //make histos
  TH1F* singleJetRates_emu = new TH1F("singleJetRates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* doubleJetRates_emu = new TH1F("doubleJetRates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* tripleJetRates_emu = new TH1F("tripleJetRates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* quadJetRates_emu = new TH1F("quadJetRates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* singleJetGlobalRates_emu = new TH1F("singleJetGlobalRates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* doubleJetGlobalRates_emu = new TH1F("doubleJetGlobalRates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* tripleJetGlobalRates_emu = new TH1F("tripleJetGlobalRates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* quadJetGlobalRates_emu = new TH1F("quadJetGlobalRates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* singleEgRates_emu = new TH1F("singleEgRates_emu", axR.c_str(), nEgBins, egLo, egHi);
  TH1F* doubleEgRates_emu = new TH1F("doubleEgRates_emu", axR.c_str(), nEgBins, egLo, egHi);
  TH1F* singleTauRates_emu = new TH1F("singleTauRates_emu", axR.c_str(), nTauBins, tauLo, tauHi);
  TH1F* doubleTauRates_emu = new TH1F("doubleTauRates_emu", axR.c_str(), nTauBins, tauLo, tauHi);
  TH1F* singleISOEgRates_emu = new TH1F("singleISOEgRates_emu", axR.c_str(), nEgBins, egLo, egHi);
  TH1F* doubleISOEgRates_emu = new TH1F("doubleISOEgRates_emu", axR.c_str(), nEgBins, egLo, egHi);
  TH1F* singleISOTauRates_emu = new TH1F("singleISOTauRates_emu", axR.c_str(), nTauBins, tauLo, tauHi);
  TH1F* doubleISOTauRates_emu = new TH1F("doubleISOTauRates_emu", axR.c_str(), nTauBins, tauLo, tauHi);
  TH1F* htSumRates_emu = new TH1F("htSumRates_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
  TH1F* htSumGlobalRates_emu = new TH1F("htSumGlobalRates_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
  TH1F* mhtSumRates_emu = new TH1F("mhtSumRates_emu",axR.c_str(), nMhtSumBins, mhtSumLo, mhtSumHi);
  TH1F* etSumRates_emu = new TH1F("etSumRates_emu",axR.c_str(), nEtSumBins, etSumLo, etSumHi);
  TH1F* etSumGlobalRates_emu = new TH1F("etSumGlobalRates_emu",axR.c_str(), nEtSumBins, etSumLo, etSumHi);
  TH1F* metSumRates_emu = new TH1F("metSumRates_emu",axR.c_str(), nMetSumBins, metSumLo, metSumHi); 
  TH1F* metHFSumRates_emu = new TH1F("metHFSumRates_emu",axR.c_str(), nMetHFSumBins, metHFSumLo, metHFSumHi); 
  
  TH1F* singleJetRates_hw = new TH1F("singleJetRates_hw", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* doubleJetRates_hw = new TH1F("doubleJetRates_hw", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* tripleJetRates_hw = new TH1F("tripleJetRates_hw", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* quadJetRates_hw = new TH1F("quadJetRates_hw", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* singleJetGlobalRates_hw = new TH1F("singleJetGlobalRates_hw", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* doubleJetGlobalRates_hw = new TH1F("doubleJetGlobalRates_hw", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* tripleJetGlobalRates_hw = new TH1F("tripleJetGlobalRates_hw", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* quadJetGlobalRates_hw = new TH1F("quadJetGlobalRates_hw", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* singleEgRates_hw = new TH1F("singleEgRates_hw", axR.c_str(), nEgBins, egLo, egHi);
  TH1F* doubleEgRates_hw = new TH1F("doubleEgRates_hw", axR.c_str(), nEgBins, egLo, egHi);
  TH1F* singleTauRates_hw = new TH1F("singleTauRates_hw", axR.c_str(), nTauBins, tauLo, tauHi);
  TH1F* doubleTauRates_hw = new TH1F("doubleTauRates_hw", axR.c_str(), nTauBins, tauLo, tauHi);
  TH1F* singleISOEgRates_hw = new TH1F("singleISOEgRates_hw", axR.c_str(), nEgBins, egLo, egHi);
  TH1F* doubleISOEgRates_hw = new TH1F("doubleISOEgRates_hw", axR.c_str(), nEgBins, egLo, egHi);
  TH1F* singleISOTauRates_hw = new TH1F("singleISOTauRates_hw", axR.c_str(), nTauBins, tauLo, tauHi);
  TH1F* doubleISOTauRates_hw = new TH1F("doubleISOTauRates_hw", axR.c_str(), nTauBins, tauLo, tauHi);
  TH1F* htSumRates_hw = new TH1F("htSumRates_hw",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
  TH1F* htSumGlobalRates_hw = new TH1F("htSumGlobalRates_hw",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
  TH1F* mhtSumRates_hw = new TH1F("mhtSumRates_hw",axR.c_str(), nMhtSumBins, mhtSumLo, mhtSumHi);
  TH1F* etSumRates_hw = new TH1F("etSumRates_hw",axR.c_str(), nEtSumBins, etSumLo, etSumHi);
  TH1F* etSumGlobalRates_hw = new TH1F("etSumGlobalRates_hw",axR.c_str(), nEtSumBins, etSumLo, etSumHi);
  TH1F* metSumRates_hw = new TH1F("metSumRates_hw",axR.c_str(), nMetHFSumBins, metHFSumLo, metHFSumHi); 
  TH1F* metHFSumRates_hw = new TH1F("metHFSumRates_hw",axR.c_str(), nMetHFSumBins, metHFSumLo, metHFSumHi); 

  TH1F* hcalTP_emu = new TH1F("hcalTP_emu", "HCAL TP ET;TP E_{T}; # Entries", nTpBins, tpLo, tpHi);
  TH1F* ecalTP_emu = new TH1F("ecalTP_emu", "ECAL TP ET;TP E_{T}; # Entries", nTpBins, tpLo, tpHi);

  TH1F* hcalTP_hw = new TH1F("hcalTP_hw", ";TP E_{T}; # Entries", nTpBins, tpLo, tpHi);
  TH1F* ecalTP_hw = new TH1F("ecalTP_hw", ";TP E_{T}; # Entries", nTpBins, tpLo, tpHi);

  TH1D * hJetEt = new TH1D("jetET","L1 Jet ET;E_{T};# Entries",100,0,1000);
  TH1D * hJetEt_1 = new TH1D("jetET_1","1st L1 Jet ET;E_{T};# Entries",100,0,1000);
  TH1D * hJetEt_2 = new TH1D("jetET_2","2nd L1 Jet ET;E_{T};# Entries",100,0,1000);
  TH1D * hJetEt_3 = new TH1D("jetET_3","3rd L1 Jet ET;E_{T};# Entries",100,0,1000);
  TH1D * hJetEt_4 = new TH1D("jetET_4","4th L1 Jet ET;E_{T};# Entries",100,0,1000);

  // 3 GeV energy cuts, scanning time cuts
  // inclusive
  TH1F * dt3GeV1nsMult_emu = new TH1F("dt3GeV1nsMult_emu","Multiplicity of 1ns delayed cells above 3 GeV (inclusive);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV2nsMult_emu = new TH1F("dt3GeV2nsMult_emu","Multiplicity of 2ns delayed cells above 3 GeV (inclusive);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV3nsMult_emu = new TH1F("dt3GeV3nsMult_emu","Multiplicity of 3ns delayed cells above 3 GeV (inclusive);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV4nsMult_emu = new TH1F("dt3GeV4nsMult_emu","Multiplicity of 4ns delayed cells above 3 GeV (inclusive);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV5nsMult_emu = new TH1F("dt3GeV5nsMult_emu","Multiplicity of 5ns delayed cells above 3 GeV (inclusive);Hit Multiplicity;Number of Entries",120,0,120);
  // HE
  TH1F * dt3GeV1nsHEMult_emu = new TH1F("dt3GeV1nsHEMult_emu","Multiplicity of 1ns delayed cells above 3 GeV (HE);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV2nsHEMult_emu = new TH1F("dt3GeV2nsHEMult_emu","Multiplicity of 2ns delayed cells above 3 GeV (HE);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV3nsHEMult_emu = new TH1F("dt3GeV3nsHEMult_emu","Multiplicity of 3ns delayed cells above 3 GeV (HE);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV4nsHEMult_emu = new TH1F("dt3GeV4nsHEMult_emu","Multiplicity of 4ns delayed cells above 3 GeV (HE);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV5nsHEMult_emu = new TH1F("dt3GeV5nsHEMult_emu","Multiplicity of 5ns delayed cells above 3 GeV (HE);Hit Multiplicity;Number of Entries",120,0,120);
  // HB
  TH1F * dt3GeV1nsHBMult_emu = new TH1F("dt3GeV1nsHBMult_emu","Multiplicity of 1ns delayed cells above 3 GeV (HB);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV2nsHBMult_emu = new TH1F("dt3GeV2nsHBMult_emu","Multiplicity of 2ns delayed cells above 3 GeV (HB);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV3nsHBMult_emu = new TH1F("dt3GeV3nsHBMult_emu","Multiplicity of 3ns delayed cells above 3 GeV (HB);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV4nsHBMult_emu = new TH1F("dt3GeV4nsHBMult_emu","Multiplicity of 4ns delayed cells above 3 GeV (HB);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV5nsHBMult_emu = new TH1F("dt3GeV5nsHBMult_emu","Multiplicity of 5ns delayed cells above 3 GeV (HB);Hit Multiplicity;Number of Entries",120,0,120);
  // for HCAL TP matched with L1 Jet
  // inclusive
  TH1F * dt3GeV1nsJetMult_emu = new TH1F("dt3GeV1nsJetMult_emu","Multiplicity of 1ns delayed cells above 3 GeV (inclusive, match TP with L1Jet);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV2nsJetMult_emu = new TH1F("dt3GeV2nsJetMult_emu","Multiplicity of 2ns delayed cells above 3 GeV (inclusive, match TP with L1Jet);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV3nsJetMult_emu = new TH1F("dt3GeV3nsJetMult_emu","Multiplicity of 3ns delayed cells above 3 GeV (inclusive, match TP with L1Jet);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV4nsJetMult_emu = new TH1F("dt3GeV4nsJetMult_emu","Multiplicity of 4ns delayed cells above 3 GeV (inclusive, match TP with L1Jet);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV5nsJetMult_emu = new TH1F("dt3GeV5nsJetMult_emu","Multiplicity of 5ns delayed cells above 3 GeV (inclusive, match TP with L1Jet);Hit Multiplicity;Number of Entries",120,0,120);
  // HE
  TH1F * dt3GeV1nsHEJetMult_emu = new TH1F("dt3GeV1nsHEJetMult_emu","Multiplicity of 1ns delayed cells above 3 GeV (HE, match TP with L1Jet);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV2nsHEJetMult_emu = new TH1F("dt3GeV2nsHEJetMult_emu","Multiplicity of 2ns delayed cells above 3 GeV (HE, match TP with L1Jet);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV3nsHEJetMult_emu = new TH1F("dt3GeV3nsHEJetMult_emu","Multiplicity of 3ns delayed cells above 3 GeV (HE, match TP with L1Jet);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV4nsHEJetMult_emu = new TH1F("dt3GeV4nsHEJetMult_emu","Multiplicity of 4ns delayed cells above 3 GeV (HE, match TP with L1Jet);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV5nsHEJetMult_emu = new TH1F("dt3GeV5nsHEJetMult_emu","Multiplicity of 5ns delayed cells above 3 GeV (HE, match TP with L1Jet);Hit Multiplicity;Number of Entries",120,0,120);
  // HB
  TH1F * dt3GeV1nsHBJetMult_emu = new TH1F("dt3GeV1nsHBJetMult_emu","Multiplicity of 1ns delayed cells above 3 GeV (HB, match TP with L1Jet);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV2nsHBJetMult_emu = new TH1F("dt3GeV2nsHBJetMult_emu","Multiplicity of 2ns delayed cells above 3 GeV (HB, match TP with L1Jet);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV3nsHBJetMult_emu = new TH1F("dt3GeV3nsHBJetMult_emu","Multiplicity of 3ns delayed cells above 3 GeV (HB, match TP with L1Jet);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV4nsHBJetMult_emu = new TH1F("dt3GeV4nsHBJetMult_emu","Multiplicity of 4ns delayed cells above 3 GeV (HB, match TP with L1Jet);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV5nsHBJetMult_emu = new TH1F("dt3GeV5nsHBJetMult_emu","Multiplicity of 5ns delayed cells above 3 GeV (HB, match TP with L1Jet);Hit Multiplicity;Number of Entries",120,0,120);
  // include all HCAL TPs included in DR cone of a L1 jet
  TH1F * dt3GeV3nsHBJet1Mult_emu = new TH1F("dt3GeV3nsHBJet1Mult_emu","Multiplicity of 3ns delayed cells above 3 GeV (HB, all TPs within DR<2 1st L1Jet);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV3nsHBJet2Mult_emu = new TH1F("dt3GeV3nsHBJet2Mult_emu","Multiplicity of 3ns delayed cells above 3 GeV (HB, all TPs within DR<2 2nd L1Jet);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV3nsHBJet3Mult_emu = new TH1F("dt3GeV3nsHBJet3Mult_emu","Multiplicity of 3ns delayed cells above 3 GeV (HB, all TPs within DR<2 3rd L1Jet);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeV3nsHBJet4Mult_emu = new TH1F("dt3GeV3nsHBJet4Mult_emu","Multiplicity of 3ns delayed cells above 3 GeV (HB, all TPs within DR<2 4th L1Jet);Hit Multiplicity;Number of Entries",120,0,120);

  // 2 GeV energy cuts, scanning time cut
  // inclusive
  TH1F * dt2GeV1nsMult_emu = new TH1F("dt2GeV1nsMult_emu","Multiplicity of 1ns delayed cells above 2 GeV (inclusive);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV2nsMult_emu = new TH1F("dt2GeV2nsMult_emu","Multiplicity of 2ns delayed cells above 2 GeV (inclusive);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV3nsMult_emu = new TH1F("dt2GeV3nsMult_emu","Multiplicity of 3ns delayed cells above 2 GeV (inclusive);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV4nsMult_emu = new TH1F("dt2GeV4nsMult_emu","Multiplicity of 4ns delayed cells above 2 GeV (inclusive);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV5nsMult_emu = new TH1F("dt2GeV5nsMult_emu","Multiplicity of 5ns delayed cells above 2 GeV (inclusive);Hit Multiplicity;Number of Entries",200,0,200);
  // HE
  TH1F * dt2GeV1nsHEMult_emu = new TH1F("dt2GeV1nsHEMult_emu","Multiplicity of 1ns delayed cells above 2 GeV (HE);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV2nsHEMult_emu = new TH1F("dt2GeV2nsHEMult_emu","Multiplicity of 2ns delayed cells above 2 GeV (HE);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV3nsHEMult_emu = new TH1F("dt2GeV3nsHEMult_emu","Multiplicity of 3ns delayed cells above 2 GeV (HE);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV4nsHEMult_emu = new TH1F("dt2GeV4nsHEMult_emu","Multiplicity of 4ns delayed cells above 2 GeV (HE);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV5nsHEMult_emu = new TH1F("dt2GeV5nsHEMult_emu","Multiplicity of 5ns delayed cells above 2 GeV (HE);Hit Multiplicity;Number of Entries",200,0,200);
  // HB
  TH1F * dt2GeV1nsHBMult_emu = new TH1F("dt2GeV1nsHBMult_emu","Multiplicity of 1ns delayed cells above 2 GeV (HB);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV2nsHBMult_emu = new TH1F("dt2GeV2nsHBMult_emu","Multiplicity of 2ns delayed cells above 2 GeV (HB);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV3nsHBMult_emu = new TH1F("dt2GeV3nsHBMult_emu","Multiplicity of 3ns delayed cells above 2 GeV (HB);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV4nsHBMult_emu = new TH1F("dt2GeV4nsHBMult_emu","Multiplicity of 4ns delayed cells above 2 GeV (HB);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV5nsHBMult_emu = new TH1F("dt2GeV5nsHBMult_emu","Multiplicity of 5ns delayed cells above 2 GeV (HB);Hit Multiplicity;Number of Entries",200,0,200);
  // for HCAL TP matched with L1 Jet
  // inclusive
  TH1F * dt2GeV1nsJetMult_emu = new TH1F("dt2GeV1nsJetMult_emu","Multiplicity of 1ns delayed cells above 2 GeV (inclusive, match TP with L1Jet);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV2nsJetMult_emu = new TH1F("dt2GeV2nsJetMult_emu","Multiplicity of 2ns delayed cells above 2 GeV (inclusive, match TP with L1Jet);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV3nsJetMult_emu = new TH1F("dt2GeV3nsJetMult_emu","Multiplicity of 3ns delayed cells above 2 GeV (inclusive, match TP with L1Jet);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV4nsJetMult_emu = new TH1F("dt2GeV4nsJetMult_emu","Multiplicity of 4ns delayed cells above 2 GeV (inclusive, match TP with L1Jet);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV5nsJetMult_emu = new TH1F("dt2GeV5nsJetMult_emu","Multiplicity of 5ns delayed cells above 2 GeV (inclusive, match TP with L1Jet);Hit Multiplicity;Number of Entries",200,0,200);
  // HE
  TH1F * dt2GeV1nsHEJetMult_emu = new TH1F("dt2GeV1nsHEJetMult_emu","Multiplicity of 1ns delayed cells above 2 GeV (HE, match TP with L1Jet);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV2nsHEJetMult_emu = new TH1F("dt2GeV2nsHEJetMult_emu","Multiplicity of 2ns delayed cells above 2 GeV (HE, match TP with L1Jet);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV3nsHEJetMult_emu = new TH1F("dt2GeV3nsHEJetMult_emu","Multiplicity of 3ns delayed cells above 2 GeV (HE, match TP with L1Jet);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV4nsHEJetMult_emu = new TH1F("dt2GeV4nsHEJetMult_emu","Multiplicity of 4ns delayed cells above 2 GeV (HE, match TP with L1Jet);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV5nsHEJetMult_emu = new TH1F("dt2GeV5nsHEJetMult_emu","Multiplicity of 5ns delayed cells above 2 GeV (HE, match TP with L1Jet);Hit Multiplicity;Number of Entries",200,0,200);
  // HB
  TH1F * dt2GeV1nsHBJetMult_emu = new TH1F("dt2GeV1nsHBJetMult_emu","Multiplicity of 1ns delayed cells above 2 GeV (HB, match TP with L1Jet);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV2nsHBJetMult_emu = new TH1F("dt2GeV2nsHBJetMult_emu","Multiplicity of 2ns delayed cells above 2 GeV (HB, match TP with L1Jet);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV3nsHBJetMult_emu = new TH1F("dt2GeV3nsHBJetMult_emu","Multiplicity of 3ns delayed cells above 2 GeV (HB, match TP with L1Jet);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV4nsHBJetMult_emu = new TH1F("dt2GeV4nsHBJetMult_emu","Multiplicity of 4ns delayed cells above 2 GeV (HB, match TP with L1Jet);Hit Multiplicity;Number of Entries",200,0,200);
  TH1F * dt2GeV5nsHBJetMult_emu = new TH1F("dt2GeV5nsHBJetMult_emu","Multiplicity of 5ns delayed cells above 2 GeV (HB, match TP with L1Jet);Hit Multiplicity;Number of Entries",200,0,200);
  // 1 GeV energy cuts, scanning time cuts
  // inclusive
  TH1F * dt1GeV1nsMult_emu = new TH1F("dt1GeV1nsMult_emu","Multiplicity of 1ns delayed cells above 1 GeV (inclusive);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV2nsMult_emu = new TH1F("dt1GeV2nsMult_emu","Multiplicity of 2ns delayed cells above 1 GeV (inclusive);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV3nsMult_emu = new TH1F("dt1GeV3nsMult_emu","Multiplicity of 3ns delayed cells above 1 GeV (inclusive);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV4nsMult_emu = new TH1F("dt1GeV4nsMult_emu","Multiplicity of 4ns delayed cells above 1 GeV (inclusive);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV5nsMult_emu = new TH1F("dt1GeV5nsMult_emu","Multiplicity of 5ns delayed cells above 1 GeV (inclusive);Hit Multiplicity;Number of Entries",400,0,400);
  // HE
  TH1F * dt1GeV1nsHEMult_emu = new TH1F("dt1GeV1nsHEMult_emu","Multiplicity of 1ns delayed cells above 1 GeV (HE);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV2nsHEMult_emu = new TH1F("dt1GeV2nsHEMult_emu","Multiplicity of 2ns delayed cells above 1 GeV (HE);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV3nsHEMult_emu = new TH1F("dt1GeV3nsHEMult_emu","Multiplicity of 3ns delayed cells above 1 GeV (HE);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV4nsHEMult_emu = new TH1F("dt1GeV4nsHEMult_emu","Multiplicity of 4ns delayed cells above 1 GeV (HE);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV5nsHEMult_emu = new TH1F("dt1GeV5nsHEMult_emu","Multiplicity of 5ns delayed cells above 1 GeV (HE);Hit Multiplicity;Number of Entries",400,0,400);
  // HB
  TH1F * dt1GeV1nsHBMult_emu = new TH1F("dt1GeV1nsHBMult_emu","Multiplicity of 1ns delayed cells above 1 GeV (HB);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV2nsHBMult_emu = new TH1F("dt1GeV2nsHBMult_emu","Multiplicity of 2ns delayed cells above 1 GeV (HB);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV3nsHBMult_emu = new TH1F("dt1GeV3nsHBMult_emu","Multiplicity of 3ns delayed cells above 1 GeV (HB);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV4nsHBMult_emu = new TH1F("dt1GeV4nsHBMult_emu","Multiplicity of 4ns delayed cells above 1 GeV (HB);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV5nsHBMult_emu = new TH1F("dt1GeV5nsHBMult_emu","Multiplicity of 5ns delayed cells above 1 GeV (HB);Hit Multiplicity;Number of Entries",400,0,400);
  // for HCAL TP matched with L1 Jet
  // inclusive
  TH1F * dt1GeV1nsJetMult_emu = new TH1F("dt1GeV1nsJetMult_emu","Multiplicity of 1ns delayed cells above 1 GeV (inclusive, match TP with L1Jet);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV2nsJetMult_emu = new TH1F("dt1GeV2nsJetMult_emu","Multiplicity of 2ns delayed cells above 1 GeV (inclusive, match TP with L1Jet);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV3nsJetMult_emu = new TH1F("dt1GeV3nsJetMult_emu","Multiplicity of 3ns delayed cells above 1 GeV (inclusive, match TP with L1Jet);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV4nsJetMult_emu = new TH1F("dt1GeV4nsJetMult_emu","Multiplicity of 4ns delayed cells above 1 GeV (inclusive, match TP with L1Jet);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV5nsJetMult_emu = new TH1F("dt1GeV5nsJetMult_emu","Multiplicity of 5ns delayed cells above 1 GeV (inclusive, match TP with L1Jet);Hit Multiplicity;Number of Entries",400,0,400);
  // HE
  TH1F * dt1GeV1nsHEJetMult_emu = new TH1F("dt1GeV1nsHEJetMult_emu","Multiplicity of 1ns delayed cells above 1 GeV (HE, match TP with L1Jet);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV2nsHEJetMult_emu = new TH1F("dt1GeV2nsHEJetMult_emu","Multiplicity of 2ns delayed cells above 1 GeV (HE, match TP with L1Jet);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV3nsHEJetMult_emu = new TH1F("dt1GeV3nsHEJetMult_emu","Multiplicity of 3ns delayed cells above 1 GeV (HE, match TP with L1Jet);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV4nsHEJetMult_emu = new TH1F("dt1GeV4nsHEJetMult_emu","Multiplicity of 4ns delayed cells above 1 GeV (HE, match TP with L1Jet);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV5nsHEJetMult_emu = new TH1F("dt1GeV5nsHEJetMult_emu","Multiplicity of 5ns delayed cells above 1 GeV (HE, match TP with L1Jet);Hit Multiplicity;Number of Entries",400,0,400);
  // HB
  TH1F * dt1GeV1nsHBJetMult_emu = new TH1F("dt1GeV1nsHBJetMult_emu","Multiplicity of 1ns delayed cells above 1 GeV (HB, match TP with L1Jet);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV2nsHBJetMult_emu = new TH1F("dt1GeV2nsHBJetMult_emu","Multiplicity of 2ns delayed cells above 1 GeV (HB, match TP with L1Jet);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV3nsHBJetMult_emu = new TH1F("dt1GeV3nsHBJetMult_emu","Multiplicity of 3ns delayed cells above 1 GeV (HB, match TP with L1Jet);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV4nsHBJetMult_emu = new TH1F("dt1GeV4nsHBJetMult_emu","Multiplicity of 4ns delayed cells above 1 GeV (HB, match TP with L1Jet);Hit Multiplicity;Number of Entries",400,0,400);
  TH1F * dt1GeV5nsHBJetMult_emu = new TH1F("dt1GeV5nsHBJetMult_emu","Multiplicity of 5ns delayed cells above 1 GeV (HB, match TP with L1Jet);Hit Multiplicity;Number of Entries",400,0,400);
  // HB by ieta regions to scan energy cuts. These are in barrel region, considering caloTowers of 4x4 in ieta, iphi, and using all same ieta values, with 2pi of iphi (4 total regions)
  TH1F * dt1GeVcaloT1Mult_emu = new TH1F("dt1GeVcaloT1Mult_emu","Multiplicity of cells above 1 GeV, CaloTower 1 (iEta 1-4);Hit Multiplicity;Number of Entries",120,0,120); 
  TH1F * dt1GeVcaloT2Mult_emu = new TH1F("dt1GeVcaloT2Mult_emu","Multiplicity of cells above 1 GeV, CaloTower 2 (iEta 5-8);Hit Multiplicity;Number of Entries",120,0,120); 
  TH1F * dt1GeVcaloT3Mult_emu = new TH1F("dt1GeVcaloT3Mult_emu","Multiplicity of cells above 1 GeV, CaloTower 3 (iEta 9-12);Hit Multiplicity;Number of Entries",120,0,120); 
  TH1F * dt1GeVcaloT4Mult_emu = new TH1F("dt1GeVcaloT4Mult_emu","Multiplicity of cells above 1 GeV, CaloTower 4 (iEta 13-16);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt2GeVcaloT1Mult_emu = new TH1F("dt2GeVcaloT1Mult_emu","Multiplicity of cells above 2 GeV, CaloTower 1 (iEta 1-4);Hit Multiplicity;Number of Entries",120,0,120);  
  TH1F * dt2GeVcaloT2Mult_emu = new TH1F("dt2GeVcaloT2Mult_emu","Multiplicity of cells above 2 GeV, CaloTower 2 (iEta 5-8);Hit Multiplicity;Number of Entries",120,0,120);  
  TH1F * dt2GeVcaloT3Mult_emu = new TH1F("dt2GeVcaloT3Mult_emu","Multiplicity of cells above 2 GeV, CaloTower 3 (iEta 9-12);Hit Multiplicity;Number of Entries",120,0,120); 
  TH1F * dt2GeVcaloT4Mult_emu = new TH1F("dt2GeVcaloT4Mult_emu","Multiplicity of cells above 2 GeV, CaloTower 4 (iEta 13-16);Hit Multiplicity;Number of Entries",120,0,120);
  TH1F * dt3GeVcaloT1Mult_emu = new TH1F("dt3GeVcaloT1Mult_emu","Multiplicity of cells above 3 GeV, CaloTower 1 (iEta 1-4);Hit Multiplicity;Number of Entries",120,0,120); 
  TH1F * dt3GeVcaloT2Mult_emu = new TH1F("dt3GeVcaloT2Mult_emu","Multiplicity of cells above 3 GeV, CaloTower 2 (iEta 5-8);Hit Multiplicity;Number of Entries",120,0,120); 
  TH1F * dt3GeVcaloT3Mult_emu = new TH1F("dt3GeVcaloT3Mult_emu","Multiplicity of cells above 3 GeV, CaloTower 3 (iEta 9-12);Hit Multiplicity;Number of Entries",120,0,120);  
  TH1F * dt3GeVcaloT4Mult_emu = new TH1F("dt3GeVcaloT4Mult_emu","Multiplicity of cells above 3 GeV, CaloTower 4 (iEta 13-16);Hit Multiplicity;Number of Entries",120,0,120); 
  // making TH2F for the energy depth plots
  TH2F * Energy_Depth = new TH2F("Energy_Depth", "TP Energy Fraction vs. Depth;HCAL Depth;Energy Fraction", 8, -0.5, 7.5, 60, 0, 1.2);
  TH2F * Timing_Depth = new TH2F("Timing_Depth", "TP Timing Value vs. Depth;HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5, 60, 0, 30);
  TH2F * Energy_DepthHE = new TH2F("Energy_DepthHE", "TP Energy Fraction vs. Depth in HE;HCAL Depth;Energy Fraction", 8, -0.5, 7.5, 60, 0, 1.2);
  TH2F * Timing_DepthHE = new TH2F("Timing_DepthHE", "TP Timing Value vs. Depth in HE;HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5, 60, 0, 30);
  TH2F * Energy_DepthHB = new TH2F("Energy_DepthHB", "TP Energy Fraction vs. Depth in HB;HCAL Depth;Energy Fraction", 8, -0.5, 7.5, 60, 0, 1.2);
  TH2F * Timing_DepthHB = new TH2F("Timing_DepthHB", "TP Timing Value vs. Depth in HB;HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5, 60, 0, 30);
  // making TH2F for the energy depth plots for high energy TPs
  TH2F * Energy_Depth_HighE = new TH2F("Energy_Depth_HighE", "TP Energy Fraction vs. Depth for TP E_{T} > 10 GeV;HCAL Depth;Energy Fraction", 8, -0.5, 7.5, 60, 0, 1.2);
  TH2F * Energy_DepthHE_HighE = new TH2F("Energy_DepthHE_HighE", "TP Energy Fraction vs. Depth in HE for TP E_{T} > 10 GeV;HCAL Depth;Energy Fraction", 8, -0.5, 7.5, 60, 0, 1.2);
  TH2F * Energy_DepthHB_HighE = new TH2F("Energy_DepthHB_HighE", "TP Energy Fraction vs. Depth in HB for TP E_{T} > 10 GeV;HCAL Depth;Energy Fraction", 8, -0.5, 7.5, 60, 0, 1.2);
  // TH2F for energy depth, where matched to jets
  TH2F * Energy_Depth_Jets = new TH2F("Energy_Depth_Jets", "TP Energy Fraction vs. Depth (TP matched w/ L1 Jets);HCAL Depth;Energy Fraction", 8, -0.5, 7.5, 60, 0, 1.2);
  TH2F * Timing_Depth_Jets = new TH2F("Timing_Depth_Jets", "TP Timing Value vs. Depth (TP matched w/ L1 Jets);HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5, 60, 0, 30);
  TH2F * Energy_DepthHE_Jets = new TH2F("Energy_DepthHE_Jets", "TP Energy Fraction vs. Depth in HE (TP matched w/ L1 Jets);HCAL Depth;Energy Fraction", 8, -0.5, 7.5, 60, 0, 1.2);
  TH2F * Timing_DepthHE_Jets = new TH2F("Timing_DepthHE_Jets", "TP Timing Value vs. Depth in HE (TP matched w/ L1 Jets);HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5, 60, 0, 30);
  TH2F * Energy_DepthHB_Jets = new TH2F("Energy_DepthHB_Jets", "TP Energy Fraction vs. Depth in HB (TP matched w/ L1 Jets);HCAL Depth;Energy Fraction", 8, -0.5, 7.5, 60, 0, 1.2);
  TH2F * Timing_DepthHB_Jets = new TH2F("Timing_DepthHB_Jets", "TP Timing Value vs. Depth in HB (TP matched w/ L1 Jets);HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5, 60, 0, 30);
  // TH2F for energy depth, where matched to jets for high energy TPs
  TH2F * Energy_Depth_Jets_HighE = new TH2F("Energy_Depth_Jets_HighE", "TP Energy Fraction vs. Depth for TP E_{T} > 10 GeV (TP matched w/ L1 Jets);HCAL Depth;Energy Fraction", 8, -0.5, 7.5, 60, 0, 1.2);
  TH2F * Energy_DepthHE_Jets_HighE = new TH2F("Energy_DepthHE_Jets_HighE", "TP Energy Fraction vs. Depth in HE for TP E_{T} > 10 GeV (TP matched w/ L1 Jets);HCAL Depth;Energy Fraction", 8, -0.5, 7.5, 60, 0, 1.2);
  TH2F * Energy_DepthHB_Jets_HighE = new TH2F("Energy_DepthHB_Jets_HighE", "TP Energy Fraction vs. Depth in HB for TP E_{T} > 10 GeV (TP matched w/ L1 Jets);HCAL Depth;Energy Fraction", 8, -0.5, 7.5, 60, 0, 1.2);
  // making TH1D for the ProfileX() from the TH2F of energy_depth or timing_depth plots
  TH1D * Energy_Depth_avg = new TH1D("Energy_Depth_avg", "TP Avg Energy Fraction vs. Depth;HCAL Depth;Energy Fraction", 8, -0.5, 7.5);
  TH1D * Timing_Depth_avg = new TH1D("Timing_Depth_avg", "TP Avg Timing Value vs. Depth;HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5);
  TH1D * Energy_DepthHE_avg = new TH1D("Energy_DepthHE_avg", "TP Avg Energy Fraction vs. Depth in HE;HCAL Depth;Energy Fraction", 8, -0.5, 7.5);
  TH1D * Timing_DepthHE_avg = new TH1D("Timing_DepthHE_avg", "TP Avg Timing Value vs. Depth in HE;HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5);
  TH1D * Energy_DepthHB_avg = new TH1D("Energy_DepthHB_avg", "TP Avg Energy Fraction vs. Depth in HB;HCAL Depth;Energy Fraction", 8, -0.5, 7.5);
  TH1D * Timing_DepthHB_avg = new TH1D("Timing_DepthHB_avg", "TP Avg Timing Value vs. Depth in HB;HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5);
  // TH1D for the ProfileX(), where matched to jets
  TH1D * Energy_Depth_avg_Jets = new TH1D("Energy_Depth_avg_Jets", "TP Avg Energy Fraction vs. Depth (TP matched w/ L1 Jets);HCAL Depth;Energy Fraction", 8, -0.5, 7.5);
  TH1D * Timing_Depth_avg_Jets = new TH1D("Timing_Depth_avg_Jets", "TP Avg Timing Value vs. Depth (TP matched w/ L1 Jets);HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5);
  TH1D * Energy_DepthHE_avg_Jets = new TH1D("Energy_DepthHE_avg_Jets", "TP Avg Energy Fraction vs. Depth in HE (TP matched w/ L1 Jets);HCAL Depth;Energy Fraction", 8, -0.5, 7.5);
  TH1D * Timing_DepthHE_avg_Jets = new TH1D("Timing_DepthHE_avg_Jets", "TP Avg Timing Value vs. Depth in HE (TP matched w/ L1 Jets);HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5);
  TH1D * Energy_DepthHB_avg_Jets = new TH1D("Energy_DepthHB_avg_Jets", "TP Avg Energy Fraction vs. Depth in HB (TP matched w/ L1 Jets);HCAL Depth;Energy Fraction", 8, -0.5, 7.5);
  TH1D * Timing_DepthHB_avg_Jets = new TH1D("Timing_DepthHB_avg_Jets", "TP Avg Timing Value vs. Depth in HB (TP matched w/ L1 Jets);HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5);
  // Ratio of energy in HCAL depth layers
  TH1F * Ratio_Depth = new TH1F("Ratio_Depth", "Ratio of First 2 HCAL Layers to E_{T};Ratio;Number of Events", 50,0,1);
  TH1F * Ratio_DepthHE = new TH1F("Ratio_DepthHE", "Ratio of First 2 HCAL Layers to E_{T} in HE;Ratio;Number of Events", 50,0,1);
  TH1F * Ratio_DepthHB = new TH1F("Ratio_DepthHB", "Ratio of First 2 HCAL Layers to E_{T} in HB;Ratio;Number of Events", 50,0,1);
  TH1F * Ratio_Depth_Jets = new TH1F("Ratio_Depth_Jets", "Ratio of First 2 HCAL Layers to E_{T}, matched w/Jets;Ratio;Number of Events", 50,0,1);
  TH1F * Ratio_DepthHE_Jets = new TH1F("Ratio_DepthHE_Jets", "Ratio of First 2 HCAL Layers to E_{T} in HE, matched w/Jets;Ratio;Number of Events", 50,0,1);
  TH1F * Ratio_DepthHB_Jets = new TH1F("Ratio_DepthHB_Jets", "Ratio of First 2 HCAL Layers to E_{T} in HB, matched w/Jets;Ratio;Number of Events", 50,0,1);
  // delta R plot for HCAL TP max energy near L1 jet
  TH1F * DeltaR_TP_L1Jet_1 = new TH1F("DeltaR_TP_L1Jet_1", "DeltaR Between First L1 Jet and Closest HCAL TP; DeltaR;Number of Events",50,0,0.5);
  TH1F * DeltaR_TP_L1Jet_2 = new TH1F("DeltaR_TP_L1Jet_2", "DeltaR Between Second L1 Jet and Closest HCAL TP; DeltaR;Number of Events",50,0,0.5);
  TH1F * DeltaR_TP_L1Jet_3 = new TH1F("DeltaR_TP_L1Jet_3", "DeltaR Between Third L1 Jet and Closest HCAL TP; DeltaR;Number of Events",50,0,0.5);
  TH1F * DeltaR_TP_L1Jet_4 = new TH1F("DeltaR_TP_L1Jet_4", "DeltaR Between Fourth L1 Jet and Closest HCAL TP; DeltaR;Number of Events",50,0,0.5);
  TH1F * DeltaR_partons_L1Jet_1 = new TH1F("DeltaR_partons_L1Jet_1", "DeltaR Between First L1 Jet and All Partons; DeltaR;Number of Events",50,0,5);
  TH1F * DeltaR_partons_L1Jet_2 = new TH1F("DeltaR_partons_L1Jet_2", "DeltaR Between Second L1 Jet and All Partons; DeltaR;Number of Events",50,0,5);
  TH1F * DeltaR_partons_L1Jet_3 = new TH1F("DeltaR_partons_L1Jet_3", "DeltaR Between Third L1 Jet and All Partons; DeltaR;Number of Events",50,0,5);
  TH1F * DeltaR_partons_L1Jet_4 = new TH1F("DeltaR_partons_L1Jet_4", "DeltaR Between Fourth L1 Jet and All Partons; DeltaR;Number of Events",50,0,5);
  TH1F * DeltaR_partonsHCAL_L1Jet_1 = new TH1F("DeltaR_partonsHCAL_L1Jet_1", "DeltaR Between First L1 Jet and All Partons with HCAL volume vertex; DeltaR;Number of Events",50,0,5);
  TH1F * DeltaR_partonsHCAL_L1Jet_2 = new TH1F("DeltaR_partonsHCAL_L1Jet_2", "DeltaR Between Second L1 Jet and All Partons with HCAL volume vertex; DeltaR;Number of Events",50,0,5);
  TH1F * DeltaR_partonsHCAL_L1Jet_3 = new TH1F("DeltaR_partonsHCAL_L1Jet_3", "DeltaR Between Third L1 Jet and All Partons with HCAL volume vertex; DeltaR;Number of Events",50,0,5);
  TH1F * DeltaR_partonsHCAL_L1Jet_4 = new TH1F("DeltaR_partonsHCAL_L1Jet_4", "DeltaR Between Fourth L1 Jet and All Partons with HCAL volume vertex; DeltaR;Number of Events",50,0,5);
  TH1F * DeltaR_TPall_L1Jet_1 = new TH1F("DeltaR_TPall_L1Jet_1", "DeltaR Between First L1 Jet and All HCAL TPs; DeltaR;Number of Events",50,0,5);
  TH1F * DeltaR_TPall_L1Jet_2 = new TH1F("DeltaR_TPall_L1Jet_2", "DeltaR Between Second L1 Jet and All HCAL TPs; DeltaR;Number of Events",50,0,5);
  TH1F * DeltaR_TPall_L1Jet_3 = new TH1F("DeltaR_TPall_L1Jet_3", "DeltaR Between Third L1 Jet and All HCAL TPs; DeltaR;Number of Events",50,0,5);
  TH1F * DeltaR_TPall_L1Jet_4 = new TH1F("DeltaR_TPall_L1Jet_4", "DeltaR Between Fourth L1 Jet and All HCAL TPs; DeltaR;Number of Events",50,0,5);
  TH1F * DeltaR_TPet_L1Jet_1 = new TH1F("DeltaR_TPet_L1Jet_1", "DeltaR Between First L1 Jet and HCAL TPs above energy, time cuts; DeltaR;Number of Events",50,0,5);
  TH1F * DeltaR_TPet_L1Jet_2 = new TH1F("DeltaR_TPet_L1Jet_2", "DeltaR Between Second L1 Jet and HCAL TPs above energy, time cuts; DeltaR;Number of Events",50,0,5);
  TH1F * DeltaR_TPet_L1Jet_3 = new TH1F("DeltaR_TPet_L1Jet_3", "DeltaR Between Third L1 Jet and HCAL TPs above energy, time cuts; DeltaR;Number of Events",50,0,5);
  TH1F * DeltaR_TPet_L1Jet_4 = new TH1F("DeltaR_TPet_L1Jet_4", "DeltaR Between Fourth L1 Jet and HCAL TPs above energy, time cuts; DeltaR;Number of Events",50,0,5);
  // delta R plots between each L1 jet
  TH1F * DeltaR_L1Jets_1_2 = new TH1F("DeltaR_L1Jets_1_2", "DeltaR Between 1st and 2nd L1 Jets;DeltaR;Number of Events",90,0,6);
  TH1F * DeltaR_L1Jets_1_3 = new TH1F("DeltaR_L1Jets_1_3", "DeltaR Between 1st and 3rd L1 Jets;DeltaR;Number of Events",90,0,6);
  TH1F * DeltaR_L1Jets_1_4 = new TH1F("DeltaR_L1Jets_1_4", "DeltaR Between 1st and 4th L1 Jets;DeltaR;Number of Events",90,0,6);
  TH1F * DeltaR_L1Jets_2_3 = new TH1F("DeltaR_L1Jets_2_3", "DeltaR Between 2nd and 3rd L1 Jets;DeltaR;Number of Events",90,0,6);
  TH1F * DeltaR_L1Jets_2_4 = new TH1F("DeltaR_L1Jets_2_4", "DeltaR Between 2nd and 4th L1 Jets;DeltaR;Number of Events",90,0,6);
  TH1F * DeltaR_L1Jets_3_4 = new TH1F("DeltaR_L1Jets_3_4", "DeltaR Between 3rd and 4th L1 Jets;DeltaR;Number of Events",90,0,6);
  // timing values for center of barrel
  TH1F * centralTiming = new TH1F("centralTiming", "Time of arrival - TOF (central barrel iEta);Time (ns);Number of Events",50,-10,40);
  // eta, ieta, phi, iphi plots to check 
  TH1F * JetiEta_1 = new TH1F("JetiEta_1", "iEta position of First L1Jet;iEta;Number of Events",80,-40,40);
  TH1F * JetiPhi_1 = new TH1F("JetiPhi_1", "iPhi position of First L1Jet;iPhi;Number of Events",74,0,74);
  TH1F * JetEta_1 = new TH1F("JetEta_1", "Eta position of First L1Jet;Eta;Number of Events",100,-4,4);
  TH1F * JetEta_2 = new TH1F("JetEta_2", "Eta position of Second L1Jet;Eta;Number of Events",100,-4,4);
  TH1F * JetEta_3 = new TH1F("JetEta_3", "Eta position of Third L1Jet;Eta;Number of Events",100,-4,4);
  TH1F * JetEta_4 = new TH1F("JetEta_4", "Eta position of Fourth L1Jet;Eta;Number of Events",100,-4,4);
  TH1F * JetPhi_1 = new TH1F("JetPhi_1", "Phi position of First L1Jet;Phi;Number of Events",77,-3.32,3.32);
  TH1F * JetPhi_2 = new TH1F("JetPhi_2", "Phi position of Second L1Jet;Phi;Number of Events",77,-3.32,3.32);
  TH1F * JetPhi_3 = new TH1F("JetPhi_3", "Phi position of Third L1Jet;Phi;Number of Events",77,-3.32,3.32);
  TH1F * JetPhi_4 = new TH1F("JetPhi_4", "Phi position of Fourth L1Jet;Phi;Number of Events",77,-3.32,3.32);

  TH1F * HCALTPEta = new TH1F("HCALTPEta", "Eta position of HCALTP;Eta;Number of Events",100,-4,4);
  TH1F * HCALTPiEta = new TH1F("HCALTPiEta", "iEta position of HCALTP;iEta;Number of Events",80,-40,40);
  TH1F * HCALTPPhi = new TH1F("HCALTPPhi", "Phi position of HCALTP;Phi;Number of Events",77,-3.32,3.32);
  TH1F * HCALTPiPhi = new TH1F("HCALTPiPhi", "iPhi position of HCALTP;iPhi;Number of Events",74,-1,73);

  // TGraph for eta phi of HCAL TPs and L1 jets
  TGraph * etaphiTP = new TGraph();
  TGraph * etaphiJet = new TGraph();
  
  // Create a TTree Object so multiplicities can be sent to TMVA analyzer
  TTree *tree = new TTree("MultForTMVA","MultForTMVA");
  //  Double_t mult3GeV3nsHB_Jets(0);
  //  Int_t event;
  Float_t mult_Jet1, mult_Global, ET_Jet1(0);
  Int_t event;
  tree->Branch("mult_Jet1",&mult_Jet1,"mult_Jet1/F");
  tree->Branch("mult_Global",&mult_Global,"mult_Global/F");
  tree->Branch("ET_Jet1",&ET_Jet1,"ET_Jet1/F");
  tree->Branch("event",&event,"event/I");
  //  tree->Branch("EventBranch", "Event", &event);
  //  tree->Branch("LocalMultBranch", "LocalMult", &mult3GeV3nsHB_Jets);

  // counting LLP efficiencies
  double totalJets(0), passedMultJets(0);
  double totalGlobal(0), passedMultGlobal(0);
  // change these variables to change the multiplicity thresholds that are set
  uint GeV3ns3Global_threshold = 3;
  uint GeV3ns3Jet_threshold = 3;
  double DR_threshold = 0.5;

  /////////////////////////////////
  // loop through all the entries//
  /////////////////////////////////
  for (Long64_t jentry=0; jentry<nentries; jentry++){
    if((jentry%10000)==0) std::cout << "Done " << jentry  << " events of " << nentries << std::endl;

    //lumi break clause
    eventTree->GetEntry(jentry);
    //skip the corresponding event
    if (!isGoodLumiSection(event_->lumi)) continue;
    goodLumiEventCount++;

    //do routine for L1 emulator quantites
    if (emuOn){

      treeL1Towemu->GetEntry(jentry);
      treeL1CaloTPemu->GetEntry(jentry);
      treeL1emu->GetEntry(jentry);
      genTree->GetEntry(jentry); // used for gen particle matching later 

      double tpEt(0.);      
      for(int i=0; i < l1CaloTPemu_->nHCALTP; i++){
	tpEt = l1CaloTPemu_->hcalTPet[i];
	hcalTP_emu->Fill(tpEt);
      }
      for(int i=0; i < l1CaloTPemu_->nECALTP; i++){
	tpEt = l1CaloTPemu_->ecalTPet[i];
	ecalTP_emu->Fill(tpEt);
      }

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

      int seedTowerIEta(-1);
      int seedTowerIPhi(-1);
      int nDepth(-1);
      double nCaloTPemu(0), tpEtaemu(0), tpPhiemu(0), tpEtemu(0);
      uint nJetemu(0);
      nCaloTPemu = l1CaloTPemu_->nHCALTP;
      nJetemu = l1emu_->nJets;

      // hcalTPdepth and hcalTPtiming will store the timing and depth variables from the 7 HCAL layers
      double hcalTPdepth[7] = {0};
      double hcalTPtiming[7] = {0};
      std::map<const TString, std::vector<double> > TimingVariablesAllJets;
      std::map<const TString, std::vector<double> > DepthVariablesAllJets;
      // multiplicity for all HCAL TPs in an entry
      double mult3GeV1ns(0), mult3GeV2ns(0), mult3GeV3ns(0), mult3GeV4ns(0), mult3GeV5ns(0);
      double mult3GeV1nsHE(0), mult3GeV2nsHE(0), mult3GeV3nsHE(0), mult3GeV4nsHE(0), mult3GeV5nsHE(0);
      double mult3GeV1nsHB(0), mult3GeV2nsHB(0), mult3GeV3nsHB(0), mult3GeV4nsHB(0), mult3GeV5nsHB(0);
      double mult2GeV1ns(0), mult2GeV2ns(0), mult2GeV3ns(0), mult2GeV4ns(0), mult2GeV5ns(0);
      double mult2GeV1nsHE(0), mult2GeV2nsHE(0), mult2GeV3nsHE(0), mult2GeV4nsHE(0), mult2GeV5nsHE(0);
      double mult2GeV1nsHB(0), mult2GeV2nsHB(0), mult2GeV3nsHB(0), mult2GeV4nsHB(0), mult2GeV5nsHB(0);
      double mult1GeV1ns(0), mult1GeV2ns(0), mult1GeV3ns(0), mult1GeV4ns(0), mult1GeV5ns(0);
      double mult1GeV1nsHE(0), mult1GeV2nsHE(0), mult1GeV3nsHE(0), mult1GeV4nsHE(0), mult1GeV5nsHE(0);
      double mult1GeV1nsHB(0), mult1GeV2nsHB(0), mult1GeV3nsHB(0), mult1GeV4nsHB(0), mult1GeV5nsHB(0);
      // multiplicity for ieta regions of caloTowers (4x4 ieta iphi) in the barrel regions                    
      double mult1GeVcaloT1(0), mult2GeVcaloT1(0), mult3GeVcaloT1(0); // abs(ieta) between 1-4
      double mult1GeVcaloT2(0), mult2GeVcaloT2(0), mult3GeVcaloT2(0); // abs(ieta) between 5-8
      double mult1GeVcaloT3(0), mult2GeVcaloT3(0), mult3GeVcaloT3(0); // abs(ieta) between 9-12
      double mult1GeVcaloT4(0), mult2GeVcaloT4(0), mult3GeVcaloT4(0); // abs(ieta) between 13-16
      // multiplicity for when HCAL TP is matched with Jets
      double mult3GeV1ns_Jets(0), mult3GeV2ns_Jets(0), mult3GeV3ns_Jets(0), mult3GeV4ns_Jets(0), mult3GeV5ns_Jets(0);
      double mult3GeV1nsHE_Jets(0), mult3GeV2nsHE_Jets(0), mult3GeV3nsHE_Jets(0), mult3GeV4nsHE_Jets(0), mult3GeV5nsHE_Jets(0);
      double mult3GeV1nsHB_Jets(0), mult3GeV2nsHB_Jets(0), mult3GeV3nsHB_Jets(0), mult3GeV4nsHB_Jets(0), mult3GeV5nsHB_Jets(0);
      double mult2GeV1ns_Jets(0), mult2GeV2ns_Jets(0), mult2GeV3ns_Jets(0), mult2GeV4ns_Jets(0), mult2GeV5ns_Jets(0);
      double mult2GeV1nsHE_Jets(0), mult2GeV2nsHE_Jets(0), mult2GeV3nsHE_Jets(0), mult2GeV4nsHE_Jets(0), mult2GeV5nsHE_Jets(0);
      double mult2GeV1nsHB_Jets(0), mult2GeV2nsHB_Jets(0), mult2GeV3nsHB_Jets(0), mult2GeV4nsHB_Jets(0), mult2GeV5nsHB_Jets(0);
      double mult1GeV1ns_Jets(0), mult1GeV2ns_Jets(0), mult1GeV3ns_Jets(0), mult1GeV4ns_Jets(0), mult1GeV5ns_Jets(0);
      double mult1GeV1nsHE_Jets(0), mult1GeV2nsHE_Jets(0), mult1GeV3nsHE_Jets(0), mult1GeV4nsHE_Jets(0), mult1GeV5nsHE_Jets(0);
      double mult1GeV1nsHB_Jets(0), mult1GeV2nsHB_Jets(0), mult1GeV3nsHB_Jets(0), mult1GeV4nsHB_Jets(0), mult1GeV5nsHB_Jets(0);
      double mult3GeV3nsHB_Jet0(0), mult3GeV3nsHB_Jet1(0), mult3GeV3nsHB_Jet2(0), mult3GeV3nsHB_Jet3(0);
      double JetEta1(0), JetEta2(0), JetEta3(0), JetEta4(0);
      double JetPhi1(0), JetPhi2(0), JetPhi3(0), JetPhi4(0);


      // loop over L1 jets, and only do first four (4 highest energy L1 jets from 4 leptons)
      for(uint jetIt=0; jetIt < nJetemu && jetIt < 4; jetIt++){
	hJetEt->Fill(l1emu_->jetEt[jetIt]); // these are already in order of highest E_T
	if ((jetIt == 0) && (l1emu_->jetEt[jetIt] < 1000) ) ET_Jet1 = l1emu_->jetEt[jetIt];
        if (jetIt == 0 ) hJetEt_1->Fill(l1emu_->jetEt[jetIt]);
        if (jetIt == 1 ) hJetEt_2->Fill(l1emu_->jetEt[jetIt]);
        if (jetIt == 2 ) hJetEt_3->Fill(l1emu_->jetEt[jetIt]);
        if (jetIt == 3 ) hJetEt_4->Fill(l1emu_->jetEt[jetIt]);

	seedTowerIPhi = l1emu_->jetTowerIPhi[jetIt];
	seedTowerIEta = l1emu_->jetTowerIEta[jetIt];

	double Jet_eta;
	double Jet_phi;
	Jet_eta = l1emu_->jetEta[jetIt]; // etaVal(seedTowerIEta);
	Jet_phi = l1emu_->jetPhi[jetIt];  //phiVal(seedTowerIPhi);

	// positions and deltaR between each jet to understand spatial information
	if (jetIt == 0 ) {
	  JetiEta_1->Fill(seedTowerIEta);
	  JetEta_1->Fill(Jet_eta);
	  JetiPhi_1->Fill(seedTowerIPhi);
	  JetPhi_1->Fill(Jet_phi);
	  JetEta1 = Jet_eta;
	  JetPhi1 = Jet_phi;
	}
	if (jetIt == 1) {
	  JetEta_2->Fill(Jet_eta);
	  JetPhi_2->Fill(Jet_phi);
	  JetEta2= Jet_eta;
          JetPhi2= Jet_phi;
	}
	if (jetIt == 2){
          JetEta_3->Fill(Jet_eta);
          JetPhi_3->Fill(Jet_phi);
	  JetEta3= Jet_eta;
          JetPhi3= Jet_phi;
	}
	if (jetIt == 3){
          JetEta_4->Fill(Jet_eta);
          JetPhi_4->Fill(Jet_phi);
	  JetEta4= Jet_eta;
          JetPhi4= Jet_phi;
	}

	if(jentry == 1) { // eta phi position of 4 L1 jets
	  if (l1emu_->jetEt[jetIt] < 20) continue; // only consider jets that are greater than 20 GeV
	  etaphiJet->SetMarkerColor(2); // red marker color for jets
	  etaphiJet->SetPoint(jetIt,Jet_eta,Jet_phi);
	}
	double min_DeltaR = 100; // used for storing min deltaR value between a L1 jet and a HCAL TP
	if (l1emu_->jetEt[jetIt] < 20 ) continue; // require jet is greater than 20 GeV to attempt matching to HCAL TP

	// ************* GEN PARTICLE MATCHING CODE ******************                                                         
	// already got the entries from genTree at start of event loop                  
	double min_deltaR_partonGen_L1jet = 1000;
	double min_deltaR_partonGenHCAL_L1jet = 1000; // min dR defined per each L1 jet
	for (int partonN = 0; partonN < generator_->nPart; partonN ++) {
	  // reset values for eta, phi of partons (w/ and w/o requirement of vertex within HCAL volume), and for dR between parton and L1 jet on a per-parton basis
	  double partonEta = 1000;
	  double partonPhi = 1000;
	  double partonEtaHCAL = 1000;
	  double partonPhiHCAL = 1000;
	  double deltaR_partonGen_L1jet = 1000;
	  double deltaR_partonGenHCAL_L1jet = 1000;
	  if (generator_->partHardProcess[partonN] == 0 ) continue; // require hard process - this is what LLP results from 
	  if ( (abs(generator_->partId[partonN]) >= 1 && abs(generator_->partId[partonN]) <=5 ) || (abs(generator_->partId[partonN]) == 21) ) { // only consider generator particles that are partons (quarks, gluons) from the LLP decay. Top quark is unstable 
	    // detector cuts for the HB and HE regions
	    if ((generator_->partPt[partonN] > 20) && (abs(intersect(generator_->partVx[partonN],generator_->partVy[partonN],generator_->partVz[partonN], generator_->partPx[partonN],generator_->partPy[partonN],generator_->partPz[partonN])[0]) < 3) ) { // require parton pT > 20 GeV to be considered for gen matching and in the HE HB region
	      //	      std::cout << generator_->partId[partonN] << std::endl;
	      partonEta = intersect(generator_->partVx[partonN],generator_->partVy[partonN],generator_->partVz[partonN], generator_->partPx[partonN],generator_->partPy[partonN],generator_->partPz[partonN])[0];
	      partonPhi = intersect(generator_->partVx[partonN],generator_->partVy[partonN],generator_->partVz[partonN], generator_->partPx[partonN],generator_->partPy[partonN],generator_->partPz[partonN])[1];
	      // see if LLP decayed in HCAL volume
	      // requiring quark gen vertex to be in the HCAL volume
	      // HE region 3.88 - 5.68 m and out to 2.95 m radius
	      // HB region z below 3.88m, radius 1.79 - 2.95m
	      double radius = sqrt(generator_->partVx[partonN]*generator_->partVx[partonN] + generator_->partVy[partonN]*generator_->partVy[partonN]);
	      if ( abs(generator_->partVz[partonN]) < 568 && (radius < 295) && ( (abs(generator_->partVz[partonN]) > 388) || radius > 179 ) ) {
		partonEtaHCAL = intersect(generator_->partVx[partonN],generator_->partVy[partonN],generator_->partVz[partonN], generator_->partPx[partonN],generator_->partPy[partonN],generator_->partPz[partonN])[0];
		partonPhiHCAL = intersect(generator_->partVx[partonN],generator_->partVy[partonN],generator_->partVz[partonN], generator_->partPx[partonN],generator_->partPy[partonN],generator_->partPz[partonN])[1];
	      }
	    } // end parton Pt cut loop 
	  } // end loop restricting to quarks or gluons 
	  // calculate the delta R between the L1 jet and the generator parton that passed all cuts
	  deltaR_partonGen_L1jet = deltaR(Jet_eta,Jet_phi,partonEta,partonPhi);
	  if (deltaR_partonGen_L1jet < min_deltaR_partonGen_L1jet ) min_deltaR_partonGen_L1jet = deltaR_partonGen_L1jet; // find the min dR between each L1 jet and the gen partons. Use this for gen matching criterion
	  deltaR_partonGenHCAL_L1jet = deltaR(Jet_eta,Jet_phi,partonEtaHCAL,partonPhiHCAL);
          if (deltaR_partonGenHCAL_L1jet < min_deltaR_partonGenHCAL_L1jet ) min_deltaR_partonGenHCAL_L1jet = deltaR_partonGenHCAL_L1jet; // find the min dR between each L1 jet and the gen partons. Use this for gen matching criterion, when LLP decays in HCAL volume

	  if (jetIt == 0) DeltaR_partons_L1Jet_1->Fill(deltaR_partonGen_L1jet); // deltaR plot of distance between L1 jet 1 and all partons passing cuts
	  if (jetIt == 1) DeltaR_partons_L1Jet_2->Fill(deltaR_partonGen_L1jet);
	  if (jetIt == 2) DeltaR_partons_L1Jet_3->Fill(deltaR_partonGen_L1jet);
	  if (jetIt == 3) DeltaR_partons_L1Jet_4->Fill(deltaR_partonGen_L1jet);
	  if (jetIt == 0) DeltaR_partonsHCAL_L1Jet_1->Fill(deltaR_partonGenHCAL_L1jet); // deltaR plot of distance between L1 jet 1 and all partons passing cuts that have a vertex in the HCAL volume
          if (jetIt == 1) DeltaR_partonsHCAL_L1Jet_2->Fill(deltaR_partonGenHCAL_L1jet);
          if (jetIt == 2) DeltaR_partonsHCAL_L1Jet_3->Fill(deltaR_partonGenHCAL_L1jet);
          if (jetIt == 3) DeltaR_partonsHCAL_L1Jet_4->Fill(deltaR_partonGenHCAL_L1jet);
	} // end parton loop

	// DeltaR cut between L1 jet and generator partons -- this is where GEN MATCHING is enforced
	if (min_deltaR_partonGen_L1jet > 0.3) continue; 

	// multiplicity plots including HCAL TPs in a DeltaR region of the L1 jets
	// loop over HCAL TPs
	for (int HcalTPIt = 0; HcalTPIt < nCaloTPemu; HcalTPIt++){
	  tpEtaemu = l1CaloTPemu_->hcalTPieta[HcalTPIt]; // use for HB HE restrictions                                 
	  tpPhiemu = l1CaloTPemu_->hcalTPCaliphi[HcalTPIt];  // hcalTPiphi[HcalTPIt]; // use for deltaR calculation
	  tpEtemu = l1CaloTPemu_->hcalTPet[HcalTPIt]; // used for energy normalization in the energy depth plots    
	  nDepth = l1CaloTPemu_->hcalTPnDepths[HcalTPIt]; // how many HCAL depth layers
	  if (nDepth == 0) continue; // skip events
	  // Energy deposited in each depth layer for every HCAL TP (4 in HB, 7 in HE)  
	  hcalTPdepth[0] = l1CaloTPemu_->hcalTPDepth1[HcalTPIt];
	  hcalTPdepth[1] = l1CaloTPemu_->hcalTPDepth2[HcalTPIt];
	  hcalTPdepth[2] = l1CaloTPemu_->hcalTPDepth3[HcalTPIt];
	  hcalTPdepth[3] = l1CaloTPemu_->hcalTPDepth4[HcalTPIt];
	  hcalTPdepth[4] = l1CaloTPemu_->hcalTPDepth5[HcalTPIt];
	  hcalTPdepth[5] = l1CaloTPemu_->hcalTPDepth6[HcalTPIt];
	  hcalTPdepth[6] = l1CaloTPemu_->hcalTPDepth7[HcalTPIt];
	  // timing info for each layer, in 25 ns with resolution 0.5 ns 
	  hcalTPtiming[0] = l1CaloTPemu_->hcalTPtiming1[HcalTPIt];
	  hcalTPtiming[1] = l1CaloTPemu_->hcalTPtiming2[HcalTPIt];
	  hcalTPtiming[2] = l1CaloTPemu_->hcalTPtiming3[HcalTPIt];
	  hcalTPtiming[3] = l1CaloTPemu_->hcalTPtiming4[HcalTPIt];
	  hcalTPtiming[4] = l1CaloTPemu_->hcalTPtiming5[HcalTPIt];
	  hcalTPtiming[5] = l1CaloTPemu_->hcalTPtiming6[HcalTPIt];
	  hcalTPtiming[6] = l1CaloTPemu_->hcalTPtiming7[HcalTPIt];

	  double TP_eta;
	  double TP_phi;
	  TP_eta = etaVal(tpEtaemu);
	  TP_phi = phiVal(tpPhiemu);
	  // minimum deltaR between HCAL TP and each L1 jet
	  if (deltaR(Jet_eta,Jet_phi,TP_eta,TP_phi) < min_DeltaR ) min_DeltaR = deltaR(Jet_eta,Jet_phi,TP_eta,TP_phi);
	  // fill DR plots between each L1 jet and every HCAL TP
	  if (jetIt == 0) DeltaR_TPall_L1Jet_1->Fill(deltaR(Jet_eta,Jet_phi,TP_eta,TP_phi));
	  if (jetIt == 1) DeltaR_TPall_L1Jet_2->Fill(deltaR(Jet_eta,Jet_phi,TP_eta,TP_phi));
	  if (jetIt == 2) DeltaR_TPall_L1Jet_3->Fill(deltaR(Jet_eta,Jet_phi,TP_eta,TP_phi));
	  if (jetIt == 3) DeltaR_TPall_L1Jet_4->Fill(deltaR(Jet_eta,Jet_phi,TP_eta,TP_phi));
	  // fill DR plots between each L1 jet and HCAL TPs with at least one layer above energy and timing cuts
	  if ( (hcalTPtiming[0]>3 && hcalTPdepth[0]>3) || (hcalTPtiming[1]>3 && hcalTPdepth[1]>3) || (hcalTPtiming[2]>3 && hcalTPdepth[2]>3) || (hcalTPtiming[3]>3 && hcalTPdepth[3]>3) || (hcalTPtiming[4]>3 && hcalTPdepth[4]>3) || (hcalTPtiming[5]>3 && hcalTPdepth[5]>3) || (hcalTPtiming[6]>3 && hcalTPdepth[6]>3) ) {
	    if (jetIt == 0) DeltaR_TPet_L1Jet_1->Fill(deltaR(Jet_eta,Jet_phi,TP_eta,TP_phi));
	    if (jetIt == 1) DeltaR_TPet_L1Jet_2->Fill(deltaR(Jet_eta,Jet_phi,TP_eta,TP_phi));
	    if (jetIt == 2) DeltaR_TPet_L1Jet_3->Fill(deltaR(Jet_eta,Jet_phi,TP_eta,TP_phi));
	    if (jetIt == 3) DeltaR_TPet_L1Jet_4->Fill(deltaR(Jet_eta,Jet_phi,TP_eta,TP_phi));
	  }

	  if (deltaR(Jet_eta,Jet_phi,TP_eta,TP_phi) > DR_threshold ) continue; // DR between L1 jet and HCAL TP for the matching in a DR cone

	  // Ratio of energy in first two HCAL layers to all HCAL layers. Only consider for high energy TPs > 10 GeV
	  if ( tpEtemu > 10 ) {
	    Ratio_Depth_Jets->Fill( (hcalTPdepth[0]+hcalTPdepth[1]) / tpEtemu);
	    if (abs(tpEtaemu) < 16) {
	      Ratio_DepthHB_Jets->Fill( (hcalTPdepth[0]+hcalTPdepth[1]) / tpEtemu);
	    }
	    if ((abs(tpEtaemu) > 16) && (abs(tpEtaemu) < 29) ) {
	      Ratio_DepthHE_Jets->Fill( (hcalTPdepth[0]+hcalTPdepth[1]) / tpEtemu);
	    }
	  }

	  // loop over HCAL depths for the HCAL TP
	  // filling fractional energy deposit and avg time plots for each of 7 HCAL depths. Done for all HCAL TPs within DR 0.5 of the L1 Jet
	  for (int depthIt = 0; depthIt < nDepth-1; depthIt++){
	    Energy_Depth_Jets->Fill(depthIt+1,hcalTPdepth[depthIt]/tpEtemu); // normalized by total energy in event so is fractional energy in each layer
	    Timing_Depth_Jets->Fill(depthIt+1,hcalTPtiming[depthIt]); // raw timing value in each layer
	    if (tpEtemu > 10 ) {
	      Energy_Depth_Jets_HighE->Fill(depthIt+1,hcalTPdepth[depthIt]/tpEtemu); // restricting to high energy HCAL TPs
	    }
	    if (abs(tpEtaemu) < 16) {
	      Energy_DepthHB_Jets->Fill(depthIt+1,hcalTPdepth[depthIt]/tpEtemu);
	      Timing_DepthHB_Jets->Fill(depthIt+1,hcalTPtiming[depthIt]);
	      if (tpEtemu > 10 ) {
		Energy_DepthHB_Jets_HighE->Fill(depthIt+1,hcalTPdepth[depthIt]/tpEtemu);
	      }
	    }
	    if (abs(tpEtaemu) > 16 && abs(tpEtaemu) < 29 ) {
	      Energy_DepthHE_Jets->Fill(depthIt+1,hcalTPdepth[depthIt]/tpEtemu);
	      Timing_DepthHE_Jets->Fill(depthIt+1,hcalTPtiming[depthIt]);
	      if (tpEtemu > 10 ) {
		Energy_DepthHE_Jets_HighE->Fill(depthIt+1,hcalTPdepth[depthIt]/tpEtemu);
	      }
	    }

	    // multiplicity counter now
	    if ( (hcalTPdepth[depthIt] > 3) && (hcalTPtiming[depthIt] > 3) && (abs(tpEtaemu) < 16) ) { // 3 GeV 3ns in HB region                                             
	      if (jetIt == 0) mult3GeV3nsHB_Jet0 += 1;
	      if (jetIt == 1) mult3GeV3nsHB_Jet1 += 1;
	      if (jetIt == 2) mult3GeV3nsHB_Jet2 += 1;
	      if (jetIt == 3) mult3GeV3nsHB_Jet3 += 1;
	    } // close HB region loop   
	    // count multiplicity of layers given a timing and energy threshold   
	    // 3 GeV energy cut
	    if (hcalTPdepth[depthIt] > 3){
	      if (hcalTPtiming[depthIt] > 1) mult3GeV1ns_Jets += 1;
	      if (hcalTPtiming[depthIt] > 2) mult3GeV2ns_Jets += 1;
	      if (hcalTPtiming[depthIt] > 3) mult3GeV3ns_Jets += 1;
	      if (hcalTPtiming[depthIt] > 4) mult3GeV4ns_Jets += 1;
	      if (hcalTPtiming[depthIt] > 5) mult3GeV5ns_Jets += 1;
	      // 3 GeV HB regions
	      if ( (hcalTPtiming[depthIt] > 1) && (abs(tpEtaemu) < 16) ) {
		mult3GeV1nsHB_Jets += 1;
		if (hcalTPtiming[depthIt] > 2) mult3GeV2nsHB_Jets += 1;
		if (hcalTPtiming[depthIt] > 3) mult3GeV3nsHB_Jets += 1;
		if (hcalTPtiming[depthIt] > 4) mult3GeV4nsHB_Jets += 1;
		if (hcalTPtiming[depthIt] > 5) mult3GeV5nsHB_Jets += 1;
	      }
	      // 3 GeV HE regions
	      if ( (hcalTPtiming[depthIt] > 1) && (abs(tpEtaemu) > 16) && (abs(tpEtaemu) < 29) ){
		mult3GeV1nsHE_Jets += 1;
		if (hcalTPtiming[depthIt] > 2) mult3GeV2nsHE_Jets += 1;
		if (hcalTPtiming[depthIt] > 3) mult3GeV3nsHE_Jets += 1;
		if (hcalTPtiming[depthIt] > 4) mult3GeV4nsHE_Jets += 1;
		if (hcalTPtiming[depthIt] > 5) mult3GeV5nsHE_Jets += 1;
	      }
	    } // closing 3 GeV energy cut loop
	    // 2 GeV energy cut
	    if (hcalTPdepth[depthIt] > 2){
	      if (hcalTPtiming[depthIt] > 1) mult2GeV1ns_Jets += 1;
	      if (hcalTPtiming[depthIt] > 2) mult2GeV2ns_Jets += 1;
	      if (hcalTPtiming[depthIt] > 3) mult2GeV3ns_Jets += 1;
	      if (hcalTPtiming[depthIt] > 4) mult2GeV4ns_Jets += 1;
	      if (hcalTPtiming[depthIt] > 5) mult2GeV5ns_Jets += 1;
	      // 2 GeV HB regions                                
	      if ( (hcalTPtiming[depthIt] > 1) && (abs(tpEtaemu) < 16) ){
		mult2GeV1nsHB_Jets += 1;
		if (hcalTPtiming[depthIt] > 2) mult2GeV2nsHB_Jets += 1;
		if (hcalTPtiming[depthIt] > 3) mult2GeV3nsHB_Jets += 1;
		if (hcalTPtiming[depthIt] > 4) mult2GeV4nsHB_Jets += 1;
		if (hcalTPtiming[depthIt] > 5) mult2GeV5nsHB_Jets += 1;
	      }
	      // 2 GeV HE regions
	      if ( (hcalTPtiming[depthIt] > 1) && (abs(tpEtaemu) > 16) && (abs(tpEtaemu) < 29) ){
		mult2GeV1nsHE_Jets += 1;
		if (hcalTPtiming[depthIt] > 2) mult2GeV2nsHE_Jets += 1;
		if (hcalTPtiming[depthIt] > 3) mult2GeV3nsHE_Jets += 1;
		if (hcalTPtiming[depthIt] > 4) mult2GeV4nsHE_Jets += 1;
		if (hcalTPtiming[depthIt] > 5) mult2GeV5nsHE_Jets += 1;
	      }
	    } // closing 2 GeV energy cut loop
	    // 1 GeV energy cut
	    if (hcalTPdepth[depthIt] > 1){
	      if (hcalTPtiming[depthIt] > 1) mult1GeV1ns_Jets += 1;
	      if (hcalTPtiming[depthIt] > 2) mult1GeV2ns_Jets += 1;
	      if (hcalTPtiming[depthIt] > 3) mult1GeV3ns_Jets += 1;
	      if (hcalTPtiming[depthIt] > 4) mult1GeV4ns_Jets += 1;
	      if (hcalTPtiming[depthIt] > 5) mult1GeV5ns_Jets += 1;
	      // 1 GeV HB regions                                                      
	      if ( (hcalTPtiming[depthIt] > 1) && (abs(tpEtaemu) < 16) ){
		mult1GeV1nsHB_Jets += 1;
		if (hcalTPtiming[depthIt] > 2) mult1GeV2nsHB_Jets += 1;
		if (hcalTPtiming[depthIt] > 3) mult1GeV3nsHB_Jets += 1;
		if (hcalTPtiming[depthIt] > 4) mult1GeV4nsHB_Jets += 1;
		if (hcalTPtiming[depthIt] > 5) mult1GeV5nsHB_Jets += 1;
	      }
	      // 1 GeV HE regions
	      if ( (hcalTPtiming[depthIt] > 1) && (abs(tpEtaemu) > 16) && (abs(tpEtaemu) < 29) ){
		mult1GeV1nsHE_Jets += 1;
		if (hcalTPtiming[depthIt] > 2) mult1GeV2nsHE_Jets += 1;
		if (hcalTPtiming[depthIt] > 3) mult1GeV3nsHE_Jets += 1;
		if (hcalTPtiming[depthIt] > 4) mult1GeV4nsHE_Jets += 1;
		if (hcalTPtiming[depthIt] > 5) mult1GeV5nsHE_Jets += 1;
	      }
	      if (abs(tpEtaemu) < 5 ) mult1GeVcaloT1 += 1;
	      if ((abs(tpEtaemu) > 4) && (abs(tpEtaemu) < 9) ) mult1GeVcaloT2 += 1;
	      if ((abs(tpEtaemu) > 8) && (abs(tpEtaemu) < 13) ) mult1GeVcaloT3 += 1;
	      if ((abs(tpEtaemu) > 12) && (abs(tpEtaemu) < 17) ) mult1GeVcaloT4 += 1;
	    } // closing 1 GeV energy cut loop
	  } // closing HCAL depth loop
	} // closing HCAL TP loop
	// minimum DeltaR between L1 jet and HCAL TP on a per event basis
	if (jetIt == 0) DeltaR_TP_L1Jet_1->Fill(min_DeltaR); // (deltaR(Jet_eta,Jet_phi,TP_eta,TP_phi)); // (min_DeltaR);
	if (jetIt == 1) DeltaR_TP_L1Jet_2->Fill(min_DeltaR); // (deltaR(Jet_eta,Jet_phi,TP_eta,TP_phi)); // (min_DeltaR);
	if (jetIt == 2) DeltaR_TP_L1Jet_3->Fill(min_DeltaR); // (deltaR(Jet_eta,Jet_phi,TP_eta,TP_phi)); // (min_DeltaR);                
	if (jetIt == 3) DeltaR_TP_L1Jet_4->Fill(min_DeltaR); // (deltaR(Jet_eta,Jet_phi,TP_eta,TP_phi)); // (min_DeltaR);

	if (jetIt == 0) dt3GeV3nsHBJet1Mult_emu->Fill(mult3GeV3nsHB_Jet0);
	if (jetIt == 1)	dt3GeV3nsHBJet2Mult_emu->Fill(mult3GeV3nsHB_Jet1);
	if (jetIt == 2) dt3GeV3nsHBJet3Mult_emu->Fill(mult3GeV3nsHB_Jet2);
	if (jetIt == 3) dt3GeV3nsHBJet4Mult_emu->Fill(mult3GeV3nsHB_Jet3);
      } // closing L1 Jets loop
      // fill DeltaR histograms for DR between each L1 Jet to understand spatial distribution
      DeltaR_L1Jets_1_2->Fill(deltaR(JetEta1,JetPhi1,JetEta2,JetPhi2));
      DeltaR_L1Jets_1_3->Fill(deltaR(JetEta1,JetPhi1,JetEta3,JetPhi3));
      DeltaR_L1Jets_1_4->Fill(deltaR(JetEta1,JetPhi1,JetEta4,JetPhi4));
      DeltaR_L1Jets_2_3->Fill(deltaR(JetEta2,JetPhi2,JetEta3,JetPhi3));
      DeltaR_L1Jets_2_4->Fill(deltaR(JetEta2,JetPhi2,JetEta4,JetPhi4));
      DeltaR_L1Jets_3_4->Fill(deltaR(JetEta3,JetPhi3,JetEta4,JetPhi4));

      // HCAL TP loop
      // Used for plots with all HCAL TPs (not matched to L1 jets), and when each HCAL TP is associated to one L1 Jet + apply DR restrictions (in L1 jet loop inside of HCAL TP loop)
      for (int HcalTPIt = 0; HcalTPIt < nCaloTPemu; HcalTPIt++){
	tpEtaemu = l1CaloTPemu_->hcalTPieta[HcalTPIt]; // use for HB HE restrictions                                 
	tpPhiemu = l1CaloTPemu_->hcalTPCaliphi[HcalTPIt];  // hcalTPiphi[HcalTPIt]; // use for deltaR calculation
	tpEtemu = l1CaloTPemu_->hcalTPet[HcalTPIt]; // used for energy normalization in the energy depth plots    
	nDepth = l1CaloTPemu_->hcalTPnDepths[HcalTPIt]; // how many HCAL depth layers to go over
	
       	if (nDepth == 0) continue; // skipping events where depth = 0, since here timing = -1 and energy = 0 (invalid event)

	// Energy deposited in each depth layer for every HCAL TP (4 in HB, 7 in HE)  
	hcalTPdepth[0] = l1CaloTPemu_->hcalTPDepth1[HcalTPIt];
	hcalTPdepth[1] = l1CaloTPemu_->hcalTPDepth2[HcalTPIt];
	hcalTPdepth[2] = l1CaloTPemu_->hcalTPDepth3[HcalTPIt];
	hcalTPdepth[3] = l1CaloTPemu_->hcalTPDepth4[HcalTPIt];
	hcalTPdepth[4] = l1CaloTPemu_->hcalTPDepth5[HcalTPIt];
	hcalTPdepth[5] = l1CaloTPemu_->hcalTPDepth6[HcalTPIt];
	hcalTPdepth[6] = l1CaloTPemu_->hcalTPDepth7[HcalTPIt];
	// timing info for each layer, in 25 ns with resolution 0.5 ns 
	hcalTPtiming[0] = l1CaloTPemu_->hcalTPtiming1[HcalTPIt];
	hcalTPtiming[1] = l1CaloTPemu_->hcalTPtiming2[HcalTPIt];
	hcalTPtiming[2] = l1CaloTPemu_->hcalTPtiming3[HcalTPIt];
	hcalTPtiming[3] = l1CaloTPemu_->hcalTPtiming4[HcalTPIt];
	hcalTPtiming[4] = l1CaloTPemu_->hcalTPtiming5[HcalTPIt];
	hcalTPtiming[5] = l1CaloTPemu_->hcalTPtiming6[HcalTPIt];
	hcalTPtiming[6] = l1CaloTPemu_->hcalTPtiming7[HcalTPIt];

	for (int i = 0; i < 4; i++ ) {
	  if ( (abs(tpEtaemu) == 1) && (hcalTPtiming[i] > -0.5) ) {
	    centralTiming->Fill(hcalTPtiming[i]-5.5); // filling distribution of time of hits (corrected by expected TOF) in the center of the barrel (ieta = 1)
	  }
	}

        // filling energy and time plots for each of 7 HCAL depths  
	for (int i = 0; i < 7; i++){
	  Energy_Depth->Fill(i+1,hcalTPdepth[i]/tpEtemu); // normalized by total energy in event so is fractional energy in each layer      
	  Timing_Depth->Fill(i+1,hcalTPtiming[i]); // raw timing value in each layer  
          if (tpEtemu > 10 ) {
            Energy_Depth_HighE->Fill(i+1,hcalTPdepth[i]/tpEtemu); // energy depth for high energy HCAL TPs
          }
	  if (abs(tpEtaemu) < 16) {
	    Energy_DepthHB->Fill(i+1,hcalTPdepth[i]/tpEtemu);
	    Timing_DepthHB->Fill(i+1,hcalTPtiming[i]);
	    if (tpEtemu > 10 ) {
	      Energy_DepthHB_HighE->Fill(i+1,hcalTPdepth[i]/tpEtemu);
	    }
	  }
	  if (abs(tpEtaemu) > 16 && abs(tpEtaemu) < 29 ) {
	    Energy_DepthHE->Fill(i+1,hcalTPdepth[i]/tpEtemu);
	    Timing_DepthHE->Fill(i+1,hcalTPtiming[i]);
	    if (tpEtemu > 10 ) {
	      Energy_DepthHE_HighE->Fill(i+1,hcalTPdepth[i]/tpEtemu);
	    }
	  }
	}

	// ratio of energy in first two HCAL layers to total energy in HCAL, only for high energy TPs
	if ( tpEtemu > 10 ) {
	  Ratio_Depth->Fill( (hcalTPdepth[0]+hcalTPdepth[1]) / tpEtemu);
	  if (abs(tpEtaemu) < 16) {
	    Ratio_DepthHB->Fill( (hcalTPdepth[0]+hcalTPdepth[1]) / tpEtemu);
	  }
	  if (abs(tpEtaemu) > 16 && abs(tpEtaemu) < 29 ) {
	    Ratio_DepthHE->Fill( (hcalTPdepth[0]+hcalTPdepth[1]) / tpEtemu);
	  }
	}

	// loop over HCAL depths for every HCAL TP
        for (int depthIt = 0; depthIt < nDepth-1; depthIt++){
          // count multiplicity of layers given a timing and energy threshold   
          // 3 GeV energy cut
          if (hcalTPdepth[depthIt] > 3){
	    if (hcalTPtiming[depthIt] > 1) mult3GeV1ns += 1;
	    if (hcalTPtiming[depthIt] > 2) mult3GeV2ns += 1;
	    if (hcalTPtiming[depthIt] > 3) mult3GeV3ns += 1;
	    if (hcalTPtiming[depthIt] > 4) mult3GeV4ns += 1;
	    if (hcalTPtiming[depthIt] > 5) mult3GeV5ns += 1;
	    // 3 GeV HB regions
	    if ( (hcalTPtiming[depthIt] > 1) && (abs(tpEtaemu) < 16) ) {
	      mult3GeV1nsHB += 1;
	      if (hcalTPtiming[depthIt] > 2) mult3GeV2nsHB += 1;
	      if (hcalTPtiming[depthIt] > 3) mult3GeV3nsHB += 1;
	      if (hcalTPtiming[depthIt] > 4) mult3GeV4nsHB += 1;
	      if (hcalTPtiming[depthIt] > 5) mult3GeV5nsHB += 1;
	    }
	    // 3 GeV HE regions
	    if ( (hcalTPtiming[depthIt] > 1) && (abs(tpEtaemu) > 16) && (abs(tpEtaemu) < 29) ){
	      mult3GeV1nsHE += 1;
	      if (hcalTPtiming[depthIt] > 2) mult3GeV2nsHE += 1;
	      if (hcalTPtiming[depthIt] > 3) mult3GeV3nsHE += 1;
	      if (hcalTPtiming[depthIt] > 4) mult3GeV4nsHE += 1;
	      if (hcalTPtiming[depthIt] > 5) mult3GeV5nsHE += 1;
	    }
	    if (abs(tpEtaemu) < 5 ) mult3GeVcaloT1 += 1;
	    if ((abs(tpEtaemu) > 4) && (abs(tpEtaemu) < 9) ) mult3GeVcaloT2 += 1;
	    if ((abs(tpEtaemu) > 8) && (abs(tpEtaemu) < 13) ) mult3GeVcaloT3 += 1;
	    if ((abs(tpEtaemu) > 12) && (abs(tpEtaemu) < 17) ) mult3GeVcaloT4 += 1;
	  } // closing 3 GeV energy cut loop
	  // 2 GeV energy cut
	  if (hcalTPdepth[depthIt] > 2){
	    if (hcalTPtiming[depthIt] > 1) mult2GeV1ns += 1;
	    if (hcalTPtiming[depthIt] > 2) mult2GeV2ns += 1;
	    if (hcalTPtiming[depthIt] > 3) mult2GeV3ns += 1;
	    if (hcalTPtiming[depthIt] > 4) mult2GeV4ns += 1;
	    if (hcalTPtiming[depthIt] > 5) mult2GeV5ns += 1;
	    // 2 GeV HB regions                                
	    if ( (hcalTPtiming[depthIt] > 1) && (abs(tpEtaemu) < 16) ){
	      mult2GeV1nsHB += 1;
	      if (hcalTPtiming[depthIt] > 2) mult2GeV2nsHB += 1;
	      if (hcalTPtiming[depthIt] > 3) mult2GeV3nsHB += 1;
	      if (hcalTPtiming[depthIt] > 4) mult2GeV4nsHB += 1;
	      if (hcalTPtiming[depthIt] > 5) mult2GeV5nsHB += 1;
	    }
	    // 2 GeV HE regions
	    if ( (hcalTPtiming[depthIt] > 1) && (abs(tpEtaemu) > 16) && (abs(tpEtaemu) < 29) ){
	      mult2GeV1nsHE += 1;
	      if (hcalTPtiming[depthIt] > 2) mult2GeV2nsHE += 1;
	      if (hcalTPtiming[depthIt] > 3) mult2GeV3nsHE += 1;
	      if (hcalTPtiming[depthIt] > 4) mult2GeV4nsHE += 1;
	      if (hcalTPtiming[depthIt] > 5) mult2GeV5nsHE += 1;
	    }
	    if (abs(tpEtaemu) < 5 ) mult2GeVcaloT1 += 1;
	    if ((abs(tpEtaemu) > 4) && (abs(tpEtaemu) < 9) ) mult2GeVcaloT2 += 1;
	    if ((abs(tpEtaemu) > 8) && (abs(tpEtaemu) < 13) ) mult2GeVcaloT3 += 1;
	    if ((abs(tpEtaemu) > 12) && (abs(tpEtaemu) < 17) ) mult2GeVcaloT4 += 1;
	  } // closing 2 GeV energy cut loop
	  // 1 GeV energy cut
	  if (hcalTPdepth[depthIt] > 1){
	    if (hcalTPtiming[depthIt] > 1) mult1GeV1ns += 1;
	    if (hcalTPtiming[depthIt] > 2) mult1GeV2ns += 1;
	    if (hcalTPtiming[depthIt] > 3) mult1GeV3ns += 1;
	    if (hcalTPtiming[depthIt] > 4) mult1GeV4ns += 1;
	    if (hcalTPtiming[depthIt] > 5) mult1GeV5ns += 1;
	    // 1 GeV HB regions                                                      
	    if ( (hcalTPtiming[depthIt] > 1) && (abs(tpEtaemu) < 16) ){
	      mult1GeV1nsHB += 1;
	      if (hcalTPtiming[depthIt] > 2) mult1GeV2nsHB += 1;
	      if (hcalTPtiming[depthIt] > 3) mult1GeV3nsHB += 1;
	      if (hcalTPtiming[depthIt] > 4) mult1GeV4nsHB += 1;
	      if (hcalTPtiming[depthIt] > 5) mult1GeV5nsHB += 1;
	    }
	    // 1 GeV HE regions
	    if ( (hcalTPtiming[depthIt] > 1) && (abs(tpEtaemu) > 16) && (abs(tpEtaemu) < 29) ){
	      mult1GeV1nsHE += 1;
	      if (hcalTPtiming[depthIt] > 2) mult1GeV2nsHE += 1;
	      if (hcalTPtiming[depthIt] > 3) mult1GeV3nsHE += 1;
	      if (hcalTPtiming[depthIt] > 4) mult1GeV4nsHE += 1;
	      if (hcalTPtiming[depthIt] > 5) mult1GeV5nsHE += 1;
	    }
	    if (abs(tpEtaemu) < 5 ) mult1GeVcaloT1 += 1;
            if ((abs(tpEtaemu) > 4) && (abs(tpEtaemu) < 9) ) mult1GeVcaloT2 += 1;
            if ((abs(tpEtaemu) > 8) && (abs(tpEtaemu) < 13) ) mult1GeVcaloT3 += 1;
	    if ((abs(tpEtaemu) > 12) && (abs(tpEtaemu) < 17) ) mult1GeVcaloT4 += 1;
	  } // closing 1 GeV energy cut loop
	} // closing HCAL depths loop
	// here have calculated multiplicity for a single HCAL TP in an event

	// convert HCAL TP ieta, iphi to eta phi so that deltaR can be used for L1 jet matching. From https://github.com/gk199/cms-hcal-debug/blob/PulseShape/plugins/HcalCompareUpgradeChains.cc#L915-L956
	double TP_eta;
	double TP_phi;
	TP_eta = etaVal(tpEtaemu);
	TP_phi = phiVal(tpPhiemu);
	// eta, ieta, phi, iphi for HCAL TPs (not yet matched to L1 Jets)
	HCALTPEta->Fill(TP_eta);
	HCALTPiEta->Fill(tpEtaemu);
	HCALTPPhi->Fill(TP_phi);
	HCALTPiPhi->Fill(tpPhiemu);
	if(jentry == 1) { // eta phi positions for HCAL TPs
	  etaphiTP->SetMarkerColor(4); // blue marker color for HCAL TPs
	  etaphiTP->SetPoint(HcalTPIt+4,TP_eta,TP_phi);
       	}
      } // closing HCAL TP loop
    

      totalGlobal += 1;
      totalJets += 1;
      if (mult3GeV3nsHB > GeV3ns3Global_threshold ) passedMultGlobal += 1;
      if (mult3GeV3nsHB_Jets > GeV3ns3Jet_threshold ) passedMultJets += 1;

      // after HCAL depth and HCAL TP loops fill the histograms with multiplicity variables. The multiplicity counter is reset on each loop iteration 
      // 3 GeV histograms
      dt3GeV1nsMult_emu->Fill(mult3GeV1ns);
      dt3GeV2nsMult_emu->Fill(mult3GeV2ns);
      dt3GeV3nsMult_emu->Fill(mult3GeV3ns);
      dt3GeV4nsMult_emu->Fill(mult3GeV4ns);
      dt3GeV5nsMult_emu->Fill(mult3GeV5ns);
      dt3GeV1nsHEMult_emu->Fill(mult3GeV1nsHE);
      dt3GeV2nsHEMult_emu->Fill(mult3GeV2nsHE);
      dt3GeV3nsHEMult_emu->Fill(mult3GeV3nsHE);
      dt3GeV4nsHEMult_emu->Fill(mult3GeV4nsHE);
      dt3GeV5nsHEMult_emu->Fill(mult3GeV5nsHE);
      dt3GeV1nsHBMult_emu->Fill(mult3GeV1nsHB);
      dt3GeV2nsHBMult_emu->Fill(mult3GeV2nsHB);
      dt3GeV3nsHBMult_emu->Fill(mult3GeV3nsHB);
      dt3GeV4nsHBMult_emu->Fill(mult3GeV4nsHB);
      dt3GeV5nsHBMult_emu->Fill(mult3GeV5nsHB);
      // 2 GeV histograms
      dt2GeV1nsMult_emu->Fill(mult2GeV1ns);
      dt2GeV2nsMult_emu->Fill(mult2GeV2ns);
      dt2GeV3nsMult_emu->Fill(mult2GeV3ns);
      dt2GeV4nsMult_emu->Fill(mult2GeV4ns);
      dt2GeV5nsMult_emu->Fill(mult2GeV5ns);
      dt2GeV1nsHEMult_emu->Fill(mult2GeV1nsHE);
      dt2GeV2nsHEMult_emu->Fill(mult2GeV2nsHE);
      dt2GeV3nsHEMult_emu->Fill(mult2GeV3nsHE);
      dt2GeV4nsHEMult_emu->Fill(mult2GeV4nsHE);
      dt2GeV5nsHEMult_emu->Fill(mult2GeV5nsHE);
      dt2GeV1nsHBMult_emu->Fill(mult2GeV1nsHB);
      dt2GeV2nsHBMult_emu->Fill(mult2GeV2nsHB);
      dt2GeV3nsHBMult_emu->Fill(mult2GeV3nsHB);
      dt2GeV4nsHBMult_emu->Fill(mult2GeV4nsHB);
      dt2GeV5nsHBMult_emu->Fill(mult2GeV5nsHB);
      // 1 GeV histograms
      dt1GeV1nsMult_emu->Fill(mult1GeV1ns);
      dt1GeV2nsMult_emu->Fill(mult1GeV2ns);
      dt1GeV3nsMult_emu->Fill(mult1GeV3ns);
      dt1GeV4nsMult_emu->Fill(mult1GeV4ns);
      dt1GeV5nsMult_emu->Fill(mult1GeV5ns);
      dt1GeV1nsHEMult_emu->Fill(mult1GeV1nsHE);
      dt1GeV2nsHEMult_emu->Fill(mult1GeV2nsHE);
      dt1GeV3nsHEMult_emu->Fill(mult1GeV3nsHE);
      dt1GeV4nsHEMult_emu->Fill(mult1GeV4nsHE);
      dt1GeV5nsHEMult_emu->Fill(mult1GeV5nsHE);
      dt1GeV1nsHBMult_emu->Fill(mult1GeV1nsHB);
      dt1GeV2nsHBMult_emu->Fill(mult1GeV2nsHB);
      dt1GeV3nsHBMult_emu->Fill(mult1GeV3nsHB);
      dt1GeV4nsHBMult_emu->Fill(mult1GeV4nsHB);
      dt1GeV5nsHBMult_emu->Fill(mult1GeV5nsHB);
      // filling the ieta energy region scan histograms
      dt1GeVcaloT1Mult_emu->Fill(mult1GeVcaloT1);
      dt1GeVcaloT2Mult_emu->Fill(mult1GeVcaloT2);
      dt1GeVcaloT3Mult_emu->Fill(mult1GeVcaloT3);
      dt1GeVcaloT4Mult_emu->Fill(mult1GeVcaloT4);
      dt2GeVcaloT1Mult_emu->Fill(mult2GeVcaloT1);
      dt2GeVcaloT2Mult_emu->Fill(mult2GeVcaloT2);
      dt2GeVcaloT3Mult_emu->Fill(mult2GeVcaloT3);
      dt2GeVcaloT4Mult_emu->Fill(mult2GeVcaloT4);
      dt3GeVcaloT1Mult_emu->Fill(mult3GeVcaloT1);
      dt3GeVcaloT2Mult_emu->Fill(mult3GeVcaloT2);
      dt3GeVcaloT3Mult_emu->Fill(mult3GeVcaloT3);
      dt3GeVcaloT4Mult_emu->Fill(mult3GeVcaloT4);

      // after HCAL depth loop and L1 Jet loop fill histograms with multiplicity variables. Multiplicity counter reset on each loop iteration. These are for where HCAL TP is matched to the L1 Jet
      // 3 GeV histograms
      dt3GeV1nsJetMult_emu->Fill(mult3GeV1ns_Jets);
      dt3GeV2nsJetMult_emu->Fill(mult3GeV2ns_Jets);
      dt3GeV3nsJetMult_emu->Fill(mult3GeV3ns_Jets);
      dt3GeV4nsJetMult_emu->Fill(mult3GeV4ns_Jets);
      dt3GeV5nsJetMult_emu->Fill(mult3GeV5ns_Jets);
      dt3GeV1nsHEJetMult_emu->Fill(mult3GeV1nsHE_Jets);
      dt3GeV2nsHEJetMult_emu->Fill(mult3GeV2nsHE_Jets);
      dt3GeV3nsHEJetMult_emu->Fill(mult3GeV3nsHE_Jets);
      dt3GeV4nsHEJetMult_emu->Fill(mult3GeV4nsHE_Jets);
      dt3GeV5nsHEJetMult_emu->Fill(mult3GeV5nsHE_Jets);
      dt3GeV1nsHBJetMult_emu->Fill(mult3GeV1nsHB_Jets);
      dt3GeV2nsHBJetMult_emu->Fill(mult3GeV2nsHB_Jets);
      dt3GeV3nsHBJetMult_emu->Fill(mult3GeV3nsHB_Jets);
      dt3GeV4nsHBJetMult_emu->Fill(mult3GeV4nsHB_Jets);
      dt3GeV5nsHBJetMult_emu->Fill(mult3GeV5nsHB_Jets);
      // 2 GeV histograms
      dt2GeV1nsJetMult_emu->Fill(mult2GeV1ns_Jets);
      dt2GeV2nsJetMult_emu->Fill(mult2GeV2ns_Jets);
      dt2GeV3nsJetMult_emu->Fill(mult2GeV3ns_Jets);
      dt2GeV4nsJetMult_emu->Fill(mult2GeV4ns_Jets);
      dt2GeV5nsJetMult_emu->Fill(mult2GeV5ns_Jets);
      dt2GeV1nsHEJetMult_emu->Fill(mult2GeV1nsHE_Jets);
      dt2GeV2nsHEJetMult_emu->Fill(mult2GeV2nsHE_Jets);
      dt2GeV3nsHEJetMult_emu->Fill(mult2GeV3nsHE_Jets);
      dt2GeV4nsHEJetMult_emu->Fill(mult2GeV4nsHE_Jets);
      dt2GeV5nsHEJetMult_emu->Fill(mult2GeV5nsHE_Jets);
      dt2GeV1nsHBJetMult_emu->Fill(mult2GeV1nsHB_Jets);
      dt2GeV2nsHBJetMult_emu->Fill(mult2GeV2nsHB_Jets);
      dt2GeV3nsHBJetMult_emu->Fill(mult2GeV3nsHB_Jets);
      dt2GeV4nsHBJetMult_emu->Fill(mult2GeV4nsHB_Jets);
      dt2GeV5nsHBJetMult_emu->Fill(mult2GeV5nsHB_Jets);
      // 1 GeV histograms
      dt1GeV1nsJetMult_emu->Fill(mult1GeV1ns_Jets);
      dt1GeV2nsJetMult_emu->Fill(mult1GeV2ns_Jets);
      dt1GeV3nsJetMult_emu->Fill(mult1GeV3ns_Jets);
      dt1GeV4nsJetMult_emu->Fill(mult1GeV4ns_Jets);
      dt1GeV5nsJetMult_emu->Fill(mult1GeV5ns_Jets);
      dt1GeV1nsHEJetMult_emu->Fill(mult1GeV1nsHE_Jets);
      dt1GeV2nsHEJetMult_emu->Fill(mult1GeV2nsHE_Jets);
      dt1GeV3nsHEJetMult_emu->Fill(mult1GeV3nsHE_Jets);
      dt1GeV4nsHEJetMult_emu->Fill(mult1GeV4nsHE_Jets);
      dt1GeV5nsHEJetMult_emu->Fill(mult1GeV5nsHE_Jets);
      dt1GeV1nsHBJetMult_emu->Fill(mult1GeV1nsHB_Jets);
      dt1GeV2nsHBJetMult_emu->Fill(mult1GeV2nsHB_Jets);
      dt1GeV3nsHBJetMult_emu->Fill(mult1GeV3nsHB_Jets);
      dt1GeV4nsHBJetMult_emu->Fill(mult1GeV4nsHB_Jets);
      dt1GeV5nsHBJetMult_emu->Fill(mult1GeV5nsHB_Jets);

      // setting the things for the tree used with TMVA, ET_Jet1 has already been set at the top of the event and jet loop
      event = jentry;
      mult_Jet1= mult3GeV3nsHB_Jet0;
      mult_Global = mult3GeV3nsHB;
      tree->Fill();


      // for each bin fill according to whether our object has a larger corresponding energy
      // Global. nJetBins = 400, jetBinWidht = 1 so this is jetLo + 1, jetLo + 2...etc
      for(int bin=0; bin<nJetBins; bin++){
        if( ((jetEt_1) >= (jetLo + bin*jetBinWidth)) && (mult3GeV3nsHB > GeV3ns3Global_threshold) ) singleJetGlobalRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
        if( ((jetEt_2) >= (jetLo + bin*jetBinWidth)) && (mult3GeV3nsHB > GeV3ns3Global_threshold) ) doubleJetGlobalRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
	if( ((jetEt_3) >= (jetLo + bin*jetBinWidth)) && (mult3GeV3nsHB > GeV3ns3Global_threshold) ) tripleJetGlobalRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
	if( ((jetEt_4) >= (jetLo + bin*jetBinWidth)) && (mult3GeV3nsHB > GeV3ns3Global_threshold) ) quadJetGlobalRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
	// matched to L1 Jets, cut mult3GeV3nsHB_Jets>1 for 1 L1jet matched, mult3GeV3nsHB_Jets>2 for 4 L1 jets matched
        if( ((jetEt_1) >= (jetLo + bin*jetBinWidth)) && (mult3GeV3nsHB_Jets > GeV3ns3Jet_threshold) ) singleJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV     
	if( ((jetEt_2) >= (jetLo + bin*jetBinWidth)) && (mult3GeV3nsHB_Jets > GeV3ns3Jet_threshold) ) doubleJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV  
	if( ((jetEt_3) >= (jetLo + bin*jetBinWidth)) && (mult3GeV3nsHB_Jets > GeV3ns3Jet_threshold) ) tripleJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV        
	if( ((jetEt_4) >= (jetLo + bin*jetBinWidth)) && (mult3GeV3nsHB_Jets > GeV3ns3Jet_threshold) ) quadJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
      }

      // eGamma rates             
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
        if( ((htSum) >= htSumLo+(bin*htSumBinWidth)) && (mult3GeV3nsHB_Jets > GeV3ns3Jet_threshold) ) htSumRates_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV
        if( ((htSum) >= htSumLo+(bin*htSumBinWidth)) && (mult3GeV3nsHB > GeV3ns3Global_threshold) ) htSumGlobalRates_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV  
      }

      for(int bin=0; bin<nMhtSumBins; bin++){
        if( (mhtSum) >= mhtSumLo+(bin*mhtSumBinWidth) ) mhtSumRates_emu->Fill(mhtSumLo+(bin*mhtSumBinWidth)); //GeV           
      }

      for(int bin=0; bin<nEtSumBins; bin++){
        if( ((etSum) >= etSumLo+(bin*etSumBinWidth)) && (mult3GeV3nsHB_Jets > GeV3ns3Jet_threshold) ) etSumRates_emu->Fill(etSumLo+(bin*etSumBinWidth)); //GeV           
        if( ((etSum) >= etSumLo+(bin*etSumBinWidth)) && (mult3GeV3nsHB > GeV3ns3Global_threshold) ) etSumGlobalRates_emu->Fill(etSumLo+(bin*etSumBinWidth)); //GeV  
      }

      for(int bin=0; bin<nMetSumBins; bin++){
        if( (metSum) >= metSumLo+(bin*metSumBinWidth) ) metSumRates_emu->Fill(metSumLo+(bin*metSumBinWidth)); //GeV           
      }
      for(int bin=0; bin<nMetHFSumBins; bin++){
        if( (metHFSum) >= metHFSumLo+(bin*metHFSumBinWidth) ) metHFSumRates_emu->Fill(metHFSumLo+(bin*metHFSumBinWidth)); //GeV           
      }
    }// closes if 'emuOn' is true


    //do routine for L1 hardware quantities
    if (hwOn){

      treeL1CaloTPhw->GetEntry(jentry);
      double tpEt(0.);
      
      for(int i=0; i < l1CaloTPhw_->nHCALTP; i++){
	tpEt = l1CaloTPhw_->hcalTPet[i];
	hcalTP_hw->Fill(tpEt);
      }
      for(int i=0; i < l1CaloTPhw_->nECALTP; i++){
	tpEt = l1CaloTPhw_->ecalTPet[i];
	ecalTP_hw->Fill(tpEt);
      }

      treeL1hw->GetEntry(jentry);
      // get jetEt*, egEt*, tauEt, htSum, mhtSum, etSum, metSum
      // ***INCLUDES NON_ZERO bx*** can't just read values off
      double jetEt_1 = 0;
      double jetEt_2 = 0;
      double jetEt_3 = 0;
      double jetEt_4 = 0;
      for (UInt_t c=0; c<l1hw_->nJets; c++){
        if (l1hw_->jetBx[c]==0 && l1hw_->jetEt[c] > jetEt_1){
          jetEt_4 = jetEt_3;
          jetEt_3 = jetEt_2;
          jetEt_2 = jetEt_1;
          jetEt_1 = l1hw_->jetEt[c];
        }
        else if (l1hw_->jetBx[c]==0 && l1hw_->jetEt[c] <= jetEt_1 && l1hw_->jetEt[c] > jetEt_2){
          jetEt_4 = jetEt_3;
          jetEt_3 = jetEt_2;      
          jetEt_2 = l1hw_->jetEt[c];
        }
        else if (l1hw_->jetBx[c]==0 && l1hw_->jetEt[c] <= jetEt_2 && l1hw_->jetEt[c] > jetEt_3){
          jetEt_4 = jetEt_3;     
          jetEt_3 = l1hw_->jetEt[c];
        }
        else if (l1hw_->jetBx[c]==0 && l1hw_->jetEt[c] <= jetEt_3 && l1hw_->jetEt[c] > jetEt_4){   
          jetEt_4 = l1hw_->jetEt[c];
        }
      }

      double egEt_1 = 0;
      double egEt_2 = 0;
      for (UInt_t c=0; c<l1hw_->nEGs; c++){
        if (l1hw_->egBx[c]==0 && l1hw_->egEt[c] > egEt_1){
          egEt_2 = egEt_1;
          egEt_1 = l1hw_->egEt[c];
        }
        else if (l1hw_->egBx[c]==0 && l1hw_->egEt[c] <= egEt_1 && l1hw_->egEt[c] > egEt_2){
          egEt_2 = l1hw_->egEt[c];
        }
      }

      double tauEt_1 = 0;
      double tauEt_2 = 0;
      //tau pt's are not given in descending order
      for (UInt_t c=0; c<l1hw_->nTaus; c++){
        if (l1hw_->tauBx[c]==0 && l1hw_->tauEt[c] > tauEt_1){
          tauEt_1 = l1hw_->tauEt[c];
        }
        else if (l1hw_->tauBx[c]==0 && l1hw_->tauEt[c] <= tauEt_1 && l1hw_->tauEt[c] > tauEt_2){
          tauEt_2 = l1hw_->tauEt[c];
        }
      }

      double egISOEt_1 = 0;
      double egISOEt_2 = 0;
      //EG pt's are not given in descending order...bx?
      for (UInt_t c=0; c<l1hw_->nEGs; c++){
        if (l1hw_->egBx[c]==0 && l1hw_->egEt[c] > egISOEt_1 && l1hw_->egIso[c]==1){
          egISOEt_2 = egISOEt_1;
          egISOEt_1 = l1hw_->egEt[c];
        }
        else if (l1hw_->egBx[c]==0 && l1hw_->egEt[c] <= egISOEt_1 && l1hw_->egEt[c] > egISOEt_2 && l1hw_->egIso[c]==1){
          egISOEt_2 = l1hw_->egEt[c];
        }
      }

      double tauISOEt_1 = 0;
      double tauISOEt_2 = 0;
      //tau pt's are not given in descending order
      for (UInt_t c=0; c<l1hw_->nTaus; c++){
        if (l1hw_->tauBx[c]==0 && l1hw_->tauEt[c] > tauISOEt_1 && l1hw_->tauIso[c]>0){
          tauISOEt_2 = tauISOEt_1;
          tauISOEt_1 = l1hw_->tauEt[c];
        }
        else if (l1hw_->tauBx[c]==0 && l1hw_->tauEt[c] <= tauISOEt_1 && l1hw_->tauEt[c] > tauISOEt_2 && l1hw_->tauIso[c]>0){
          tauISOEt_2 = l1hw_->tauEt[c];
        }
      }

      double htSum = 0;
      double mhtSum = 0;
      double etSum = 0;
      double metSum = 0;
      double metHFSum = 0;
      // HW includes -2,-1,0,1,2 bx info (hence the different numbers, could cause a seg fault if this changes)
      for (unsigned int c=0; c<l1hw_->nSums; c++){
          if( l1hw_->sumBx[c] != 0 ) continue;
          if( l1hw_->sumType[c] == L1Analysis::kTotalEt ) etSum = l1hw_->sumEt[c];
          if( l1hw_->sumType[c] == L1Analysis::kTotalHt ) htSum = l1hw_->sumEt[c];
          if( l1hw_->sumType[c] == L1Analysis::kMissingEt ) metSum = l1hw_->sumEt[c];
	  if( l1hw_->sumType[c] == L1Analysis::kMissingEtHF ) metHFSum = l1hw_->sumEt[c];
          if( l1hw_->sumType[c] == L1Analysis::kMissingHt ) mhtSum = l1hw_->sumEt[c];
      }

      // for each bin fill according to whether our object has a larger corresponding energy
      for(int bin=0; bin<nJetBins; bin++){
        if( (jetEt_1) >= jetLo + (bin*jetBinWidth) ) singleJetRates_hw->Fill(jetLo+(bin*jetBinWidth));  //GeV
      } 

      for(int bin=0; bin<nJetBins; bin++){
        if( (jetEt_2) >= jetLo + (bin*jetBinWidth) ) doubleJetRates_hw->Fill(jetLo+(bin*jetBinWidth));  //GeV
      }  

      for(int bin=0; bin<nJetBins; bin++){
        if( (jetEt_3) >= jetLo + (bin*jetBinWidth) ) tripleJetRates_hw->Fill(jetLo+(bin*jetBinWidth));  //GeV
      }  

      for(int bin=0; bin<nJetBins; bin++){
        if( (jetEt_4) >= jetLo + (bin*jetBinWidth) ) quadJetRates_hw->Fill(jetLo+(bin*jetBinWidth));  //GeV
      }  

      for(int bin=0; bin<nJetBins; bin++){
        if( (jetEt_1) >= jetLo + (bin*jetBinWidth) ) singleJetGlobalRates_hw->Fill(jetLo+(bin*jetBinWidth));  //GeV
      }
      for(int bin=0; bin<nJetBins; bin++){
        if( (jetEt_2) >= jetLo + (bin*jetBinWidth) ) doubleJetGlobalRates_hw->Fill(jetLo+(bin*jetBinWidth));  //GeV
      }
      for(int bin=0; bin<nJetBins; bin++){
        if( (jetEt_3) >= jetLo + (bin*jetBinWidth) ) tripleJetGlobalRates_hw->Fill(jetLo+(bin*jetBinWidth));  //GeV
      }
      for(int bin=0; bin<nJetBins; bin++){
        if( (jetEt_4) >= jetLo + (bin*jetBinWidth) ) quadJetGlobalRates_hw->Fill(jetLo+(bin*jetBinWidth));  //GeV
      }
             
      for(int bin=0; bin<nEgBins; bin++){
        if( (egEt_1) >= egLo + (bin*egBinWidth) ) singleEgRates_hw->Fill(egLo+(bin*egBinWidth));  //GeV
      } 

      for(int bin=0; bin<nEgBins; bin++){
        if( (egEt_2) >= egLo + (bin*egBinWidth) ) doubleEgRates_hw->Fill(egLo+(bin*egBinWidth));  //GeV
      }  

      for(int bin=0; bin<nTauBins; bin++){
        if( (tauEt_1) >= tauLo + (bin*tauBinWidth) ) singleTauRates_hw->Fill(tauLo+(bin*tauBinWidth));  //GeV
      } 

      for(int bin=0; bin<nTauBins; bin++){
        if( (tauEt_2) >= tauLo + (bin*tauBinWidth) ) doubleTauRates_hw->Fill(tauLo+(bin*tauBinWidth));  //GeV
      } 

      for(int bin=0; bin<nEgBins; bin++){
        if( (egISOEt_1) >= egLo + (bin*egBinWidth) ) singleISOEgRates_hw->Fill(egLo+(bin*egBinWidth));  //GeV
      } 

      for(int bin=0; bin<nEgBins; bin++){
        if( (egISOEt_2) >= egLo + (bin*egBinWidth) ) doubleISOEgRates_hw->Fill(egLo+(bin*egBinWidth));  //GeV
      }  

      for(int bin=0; bin<nTauBins; bin++){
        if( (tauISOEt_1) >= tauLo + (bin*tauBinWidth) ) singleISOTauRates_hw->Fill(tauLo+(bin*tauBinWidth));  //GeV
      }

      for(int bin=0; bin<nTauBins; bin++){
        if( (tauISOEt_2) >= tauLo + (bin*tauBinWidth) ) doubleISOTauRates_hw->Fill(tauLo+(bin*tauBinWidth));  //GeV
      } 

      for(int bin=0; bin<nHtSumBins; bin++){
        if( (htSum) >= htSumLo+(bin*htSumBinWidth) ) htSumRates_hw->Fill(htSumLo+(bin*htSumBinWidth)); //GeV
        if( (htSum) >= htSumLo+(bin*htSumBinWidth) ) htSumGlobalRates_hw->Fill(htSumLo+(bin*htSumBinWidth)); //GeV  
      }

      for(int bin=0; bin<nMhtSumBins; bin++){
        if( (mhtSum) >= mhtSumLo+(bin*mhtSumBinWidth) ) mhtSumRates_hw->Fill(mhtSumLo+(bin*mhtSumBinWidth)); //GeV           
      }

      for(int bin=0; bin<nEtSumBins; bin++){
        if( (etSum) >= etSumLo+(bin*etSumBinWidth) ) etSumRates_hw->Fill(etSumLo+(bin*etSumBinWidth)); //GeV         
	if( (etSum) >= etSumLo+(bin*etSumBinWidth) ) etSumGlobalRates_hw->Fill(etSumLo+(bin*etSumBinWidth)); //GeV  
      }

      for(int bin=0; bin<nMetSumBins; bin++){
        if( (metSum) >= metSumLo+(bin*metSumBinWidth) ) metSumRates_hw->Fill(metSumLo+(bin*metSumBinWidth)); //GeV           
      } 
      for(int bin=0; bin<nMetHFSumBins; bin++){
        if( (metHFSum) >= metHFSumLo+(bin*metHFSumBinWidth) ) metHFSumRates_hw->Fill(metHFSumLo+(bin*metHFSumBinWidth)); //GeV           
      } 
    }// closes if 'hwOn' is true
  }// closes loop through events

  std::cout << "LLP passed jet matched multiplicity cut of mult3GeV3nsHB_Jets>" << GeV3ns3Jet_threshold << " / total LLPs = " << passedMultJets/totalJets << std::endl;
  std::cout << "LLP passed global multiplicity cut of mult3GeV3nsHB>" << GeV3ns3Global_threshold << " / total LLPs = " << passedMultGlobal/totalGlobal << std::endl;

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
    singleJetGlobalRates_emu->Scale(norm);
    doubleJetGlobalRates_emu->Scale(norm);
    tripleJetGlobalRates_emu->Scale(norm);
    quadJetGlobalRates_emu->Scale(norm);
    singleEgRates_emu->Scale(norm);
    doubleEgRates_emu->Scale(norm);
    singleTauRates_emu->Scale(norm);
    doubleTauRates_emu->Scale(norm);
    singleISOEgRates_emu->Scale(norm);
    doubleISOEgRates_emu->Scale(norm);
    singleISOTauRates_emu->Scale(norm);
    doubleISOTauRates_emu->Scale(norm);
    htSumRates_emu->Scale(norm);
    htSumGlobalRates_emu->Scale(norm);
    mhtSumRates_emu->Scale(norm);
    etSumRates_emu->Scale(norm);
    etSumGlobalRates_emu->Scale(norm);
    metSumRates_emu->Scale(norm);
    metHFSumRates_emu->Scale(norm);

    // TH1D as a ProfileX of the depth TH2F for timing and energy
    Energy_Depth_avg = Energy_Depth->ProfileX();
    Energy_DepthHE_avg = Energy_DepthHE->ProfileX();
    Energy_DepthHB_avg = Energy_DepthHB->ProfileX();
    Timing_Depth_avg = Timing_Depth->ProfileX();
    Timing_DepthHE_avg = Timing_DepthHE->ProfileX();
    Timing_DepthHB_avg = Timing_DepthHB->ProfileX();
    // for the ones matched with L1 jets
    Energy_Depth_avg_Jets = Energy_Depth_Jets->ProfileX();
    Energy_DepthHE_avg_Jets = Energy_DepthHE_Jets->ProfileX();
    Energy_DepthHB_avg_Jets = Energy_DepthHB_Jets->ProfileX();
    Timing_Depth_avg_Jets = Timing_Depth_Jets->ProfileX();
    Timing_DepthHE_avg_Jets = Timing_DepthHE_Jets->ProfileX();
    Timing_DepthHB_avg_Jets = Timing_DepthHB_Jets->ProfileX();

    //set the errors for the rates
    //want error -> error * sqrt(norm) ?

    if (inputFile.substr(0,11) == "../Neutrino" ) {
      int SJ_60GeV_bin = singleJetGlobalRates_emu->GetXaxis()->FindBin(60); // get x value of bin of interest                      
      int SJet60GeV = singleJetGlobalRates_emu->GetBinContent(SJ_60GeV_bin); // get rate value for single jet at 60 GeV ET threshold
      int HT_120GeV_bin = htSumGlobalRates_emu->GetXaxis()->FindBin(120); // get x value of bin of interest
      int htSum120GeV = htSumGlobalRates_emu->GetBinContent(HT_120GeV_bin); // get rate value for ht sum rate at 120 GeV ET threshold
      int HT_350GeV_bin = htSumGlobalRates_emu->GetXaxis()->FindBin(350); // get x value of bin of interest         
      int htSum350GeV = htSumGlobalRates_emu->GetBinContent(HT_350GeV_bin); // get rate value for ht sum rate at 350 GeV ET threshold
      std::cout << "For global multiplicity threshold of >" << GeV3ns3Global_threshold << " the single jet rate = " << SJet60GeV << " and htSum=120 rate = " << htSum120GeV <<  " and htSum=350 rate = " << htSum350GeV << std::endl;

      int SJ_60GeV_bin_l = singleJetRates_emu->GetXaxis()->FindBin(60); // get x value of bin of interest                          
      int SJet60GeV_l = singleJetRates_emu->GetBinContent(SJ_60GeV_bin_l); // get rate value for single jet at 60 GeV ET threshold   
      int HT_120GeV_bin_l = htSumRates_emu->GetXaxis()->FindBin(120); // get x value of bin of interest            
      int htSum120GeV_l = htSumRates_emu->GetBinContent(HT_120GeV_bin_l); // get rate value for ht sum rate at 120 GeV ET threshold
      int HT_350GeV_bin_l = htSumRates_emu->GetXaxis()->FindBin(350); // get x value of bin of interest                            
      int htSum350GeV_l = htSumRates_emu->GetBinContent(HT_350GeV_bin_l); // get rate value for ht sum rate at 350 GeV ET threshold  
      std::cout <<  "For L1 jet matched multiplicity threshold of > " << GeV3ns3Jet_threshold << " the single jet rate = " << SJet60GeV_l << " and htSum=120 rate = " << htSum120GeV_l << " and htSum=350 rate = " << htSum350GeV_l << std::endl;
      std::cout << inputFile << std::endl;
    }

    hJetEt->Write();
    hJetEt_1->Write();
    hJetEt_2->Write();
    hJetEt_3->Write();
    hJetEt_4->Write();
    hcalTP_emu->Write();
    ecalTP_emu->Write();
    singleJetRates_emu->Write();
    doubleJetRates_emu->Write();
    tripleJetRates_emu->Write();
    quadJetRates_emu->Write();
    singleJetGlobalRates_emu->Write();
    doubleJetGlobalRates_emu->Write();
    tripleJetGlobalRates_emu->Write();
    quadJetGlobalRates_emu->Write();
    singleEgRates_emu->Write();
    doubleEgRates_emu->Write();
    singleTauRates_emu->Write();
    doubleTauRates_emu->Write();
    singleISOEgRates_emu->Write();
    doubleISOEgRates_emu->Write();
    singleISOTauRates_emu->Write();
    doubleISOTauRates_emu->Write();
    htSumRates_emu->Write();
    htSumGlobalRates_emu->Write();
    mhtSumRates_emu->Write();
    etSumRates_emu->Write();
    etSumGlobalRates_emu->Write();
    metSumRates_emu->Write();
    metHFSumRates_emu->Write();
    // 3 GeV
    dt3GeV1nsMult_emu->Write();
    dt3GeV2nsMult_emu->Write();
    dt3GeV3nsMult_emu->Write();
    dt3GeV4nsMult_emu->Write();
    dt3GeV5nsMult_emu->Write();
    dt3GeV1nsHEMult_emu->Write();
    dt3GeV2nsHEMult_emu->Write();
    dt3GeV3nsHEMult_emu->Write();
    dt3GeV4nsHEMult_emu->Write();
    dt3GeV5nsHEMult_emu->Write();
    dt3GeV1nsHBMult_emu->Write();
    dt3GeV2nsHBMult_emu->Write();
    dt3GeV3nsHBMult_emu->Write();
    dt3GeV4nsHBMult_emu->Write();
    dt3GeV5nsHBMult_emu->Write();
    // 2 GeV
    dt2GeV1nsMult_emu->Write();
    dt2GeV2nsMult_emu->Write();
    dt2GeV3nsMult_emu->Write();
    dt2GeV4nsMult_emu->Write();
    dt2GeV5nsMult_emu->Write();
    dt2GeV1nsHEMult_emu->Write();
    dt2GeV2nsHEMult_emu->Write();
    dt2GeV3nsHEMult_emu->Write();
    dt2GeV4nsHEMult_emu->Write();
    dt2GeV5nsHEMult_emu->Write();
    dt2GeV1nsHBMult_emu->Write();
    dt2GeV2nsHBMult_emu->Write();
    dt2GeV3nsHBMult_emu->Write();
    dt2GeV4nsHBMult_emu->Write();
    dt2GeV5nsHBMult_emu->Write();
    // 1 GeV
    dt1GeV1nsMult_emu->Write();
    dt1GeV2nsMult_emu->Write();
    dt1GeV3nsMult_emu->Write();
    dt1GeV4nsMult_emu->Write();
    dt1GeV5nsMult_emu->Write();
    dt1GeV1nsHEMult_emu->Write();
    dt1GeV2nsHEMult_emu->Write();
    dt1GeV3nsHEMult_emu->Write();
    dt1GeV4nsHEMult_emu->Write();
    dt1GeV5nsHEMult_emu->Write();
    dt1GeV1nsHBMult_emu->Write();
    dt1GeV2nsHBMult_emu->Write();
    dt1GeV3nsHBMult_emu->Write();
    dt1GeV4nsHBMult_emu->Write();
    dt1GeV5nsHBMult_emu->Write();
    // ieta energy scan
    dt1GeVcaloT1Mult_emu->Write();
    dt1GeVcaloT2Mult_emu->Write();
    dt1GeVcaloT3Mult_emu->Write();
    dt1GeVcaloT4Mult_emu->Write();
    dt2GeVcaloT1Mult_emu->Write();
    dt2GeVcaloT2Mult_emu->Write();
    dt2GeVcaloT3Mult_emu->Write();
    dt2GeVcaloT4Mult_emu->Write();
    dt3GeVcaloT1Mult_emu->Write();
    dt3GeVcaloT2Mult_emu->Write();
    dt3GeVcaloT3Mult_emu->Write();
    dt3GeVcaloT4Mult_emu->Write();

    // and for the multiplicity counter for HCAL TPs matched with L1 Jets
    // 3 GeV
    dt3GeV1nsJetMult_emu->Write();
    dt3GeV2nsJetMult_emu->Write();
    dt3GeV3nsJetMult_emu->Write();
    dt3GeV4nsJetMult_emu->Write();
    dt3GeV5nsJetMult_emu->Write();
    dt3GeV1nsHEJetMult_emu->Write();
    dt3GeV2nsHEJetMult_emu->Write();
    dt3GeV3nsHEJetMult_emu->Write();
    dt3GeV4nsHEJetMult_emu->Write();
    dt3GeV5nsHEJetMult_emu->Write();
    dt3GeV1nsHBJetMult_emu->Write();
    dt3GeV2nsHBJetMult_emu->Write();
    dt3GeV3nsHBJetMult_emu->Write();
    dt3GeV4nsHBJetMult_emu->Write();
    dt3GeV5nsHBJetMult_emu->Write();

    dt3GeV3nsHBJet1Mult_emu->Write();
    dt3GeV3nsHBJet2Mult_emu->Write();
    dt3GeV3nsHBJet3Mult_emu->Write();
    dt3GeV3nsHBJet4Mult_emu->Write();

    // 2 GeV
    dt2GeV1nsJetMult_emu->Write();
    dt2GeV2nsJetMult_emu->Write();
    dt2GeV3nsJetMult_emu->Write();
    dt2GeV4nsJetMult_emu->Write();
    dt2GeV5nsJetMult_emu->Write();
    dt2GeV1nsHEJetMult_emu->Write();
    dt2GeV2nsHEJetMult_emu->Write();
    dt2GeV3nsHEJetMult_emu->Write();
    dt2GeV4nsHEJetMult_emu->Write();
    dt2GeV5nsHEJetMult_emu->Write();
    dt2GeV1nsHBJetMult_emu->Write();
    dt2GeV2nsHBJetMult_emu->Write();
    dt2GeV3nsHBJetMult_emu->Write();
    dt2GeV4nsHBJetMult_emu->Write();
    dt2GeV5nsHBJetMult_emu->Write();
    // 1 GeV      
    dt1GeV1nsJetMult_emu->Write();
    dt1GeV2nsJetMult_emu->Write();
    dt1GeV3nsJetMult_emu->Write();
    dt1GeV4nsJetMult_emu->Write();
    dt1GeV5nsJetMult_emu->Write();
    dt1GeV1nsHEJetMult_emu->Write();
    dt1GeV2nsHEJetMult_emu->Write();
    dt1GeV3nsHEJetMult_emu->Write();
    dt1GeV4nsHEJetMult_emu->Write();
    dt1GeV5nsHEJetMult_emu->Write();
    dt1GeV1nsHBJetMult_emu->Write();
    dt1GeV2nsHBJetMult_emu->Write();
    dt1GeV3nsHBJetMult_emu->Write();
    dt1GeV4nsHBJetMult_emu->Write();
    dt1GeV5nsHBJetMult_emu->Write();

    Energy_Depth->Write();
    Timing_Depth->Write();
    Energy_DepthHB->Write();
    Timing_DepthHB->Write();
    Energy_DepthHE->Write();
    Timing_DepthHE->Write();

    Energy_Depth_HighE->Write();
    Energy_DepthHB_HighE->Write();
    Energy_DepthHE_HighE->Write();

    Energy_Depth_Jets->Write();
    Timing_Depth_Jets->Write();
    Energy_DepthHB_Jets->Write();
    Timing_DepthHB_Jets->Write();
    Energy_DepthHE_Jets->Write();
    Timing_DepthHE_Jets->Write();

    Energy_Depth_Jets_HighE->Write();
    Energy_DepthHB_Jets_HighE->Write();
    Energy_DepthHE_Jets_HighE->Write();

    Energy_Depth_avg->Write();
    Energy_DepthHE_avg->Write();
    Energy_DepthHB_avg->Write();
    Timing_Depth_avg->Write();
    Timing_DepthHE_avg->Write();
    Timing_DepthHB_avg->Write();

    Energy_Depth_avg_Jets->Write();
    Energy_DepthHE_avg_Jets->Write();
    Energy_DepthHB_avg_Jets->Write();
    Timing_Depth_avg_Jets->Write();
    Timing_DepthHE_avg_Jets->Write();
    Timing_DepthHB_avg_Jets->Write();

    Ratio_Depth->Write();
    Ratio_DepthHE->Write();
    Ratio_DepthHB->Write();
    Ratio_Depth_Jets->Write();
    Ratio_DepthHE_Jets->Write();
    Ratio_DepthHB_Jets->Write();

    DeltaR_TP_L1Jet_1->Write();
    DeltaR_TP_L1Jet_2->Write();
    DeltaR_TP_L1Jet_3->Write();
    DeltaR_TP_L1Jet_4->Write();
    DeltaR_partons_L1Jet_1->Write();
    DeltaR_partons_L1Jet_2->Write();
    DeltaR_partons_L1Jet_3->Write();
    DeltaR_partons_L1Jet_4->Write();
    DeltaR_partonsHCAL_L1Jet_1->Write();
    DeltaR_partonsHCAL_L1Jet_2->Write();
    DeltaR_partonsHCAL_L1Jet_3->Write();
    DeltaR_partonsHCAL_L1Jet_4->Write();
    DeltaR_TPall_L1Jet_1->Write();
    DeltaR_TPall_L1Jet_2->Write();
    DeltaR_TPall_L1Jet_3->Write();
    DeltaR_TPall_L1Jet_4->Write();
    DeltaR_TPet_L1Jet_1->Write();
    DeltaR_TPet_L1Jet_2->Write();
    DeltaR_TPet_L1Jet_3->Write();
    DeltaR_TPet_L1Jet_4->Write();
    DeltaR_L1Jets_1_2->Write();
    DeltaR_L1Jets_1_3->Write();
    DeltaR_L1Jets_1_4->Write();
    DeltaR_L1Jets_2_3->Write();
    DeltaR_L1Jets_2_4->Write();
    DeltaR_L1Jets_3_4->Write();

    centralTiming->Write();

    JetiEta_1->Write();
    JetiPhi_1->Write();

    JetEta_1->Write();
    JetEta_2->Write();
    JetEta_3->Write();
    JetEta_4->Write();
    JetPhi_1->Write();
    JetPhi_2->Write();
    JetPhi_3->Write();
    JetPhi_4->Write();

    HCALTPEta->Write();
    HCALTPiEta->Write();
    HCALTPPhi->Write();
    HCALTPiPhi->Write();

    etaphiJet->GetYaxis()->SetRangeUser(-3.5,3.5);
    etaphiJet->GetXaxis()->SetLimits(-3.5,3.5);
    etaphiJet->Draw("AP*");
    etaphiJet->Write();
    etaphiTP->GetYaxis()->SetRangeUser(-3.5,3.5);
    etaphiTP->GetXaxis()->SetLimits(-3.5,3.5);
    etaphiTP->Draw("AP*");
    etaphiTP->Write();

    tree->Write();
    //    treeL1CaloTPemu->Write();
    //    kk->Write();
  }

  if (hwOn){

    singleJetRates_hw->Scale(norm);
    doubleJetRates_hw->Scale(norm);
    tripleJetRates_hw->Scale(norm);
    quadJetRates_hw->Scale(norm);
    singleJetGlobalRates_hw->Scale(norm);
    doubleJetGlobalRates_hw->Scale(norm);
    tripleJetGlobalRates_hw->Scale(norm);
    quadJetGlobalRates_hw->Scale(norm);
    singleEgRates_hw->Scale(norm);
    doubleEgRates_hw->Scale(norm);
    singleTauRates_hw->Scale(norm);
    doubleTauRates_hw->Scale(norm);
    singleISOEgRates_hw->Scale(norm);
    doubleISOEgRates_hw->Scale(norm);
    singleISOTauRates_hw->Scale(norm);
    doubleISOTauRates_hw->Scale(norm);
    htSumRates_hw->Scale(norm);
    htSumGlobalRates_hw->Scale(norm);
    mhtSumRates_hw->Scale(norm);
    etSumRates_hw->Scale(norm);
    etSumGlobalRates_hw->Scale(norm);
    metSumRates_hw->Scale(norm);
    metHFSumRates_hw->Scale(norm);

    hcalTP_hw->Write();
    ecalTP_hw->Write();
    singleJetRates_hw->Write();
    doubleJetRates_hw->Write();
    tripleJetRates_hw->Write();
    quadJetRates_hw->Write();
    singleJetGlobalRates_hw->Write();
    doubleJetGlobalRates_hw->Write();
    tripleJetGlobalRates_hw->Write();
    quadJetGlobalRates_hw->Write();
    singleEgRates_hw->Write();
    doubleEgRates_hw->Write();
    singleTauRates_hw->Write();
    doubleTauRates_hw->Write();
    singleISOEgRates_hw->Write();
    doubleISOEgRates_hw->Write();
    singleISOTauRates_hw->Write();
    doubleISOTauRates_hw->Write();
    htSumRates_hw->Write();
    htSumGlobalRates_hw->Write();
    mhtSumRates_hw->Write();
    etSumRates_hw->Write();
    etSumGlobalRates_hw->Write();
    metSumRates_hw->Write();
    metHFSumRates_hw->Write();
  }
  myfile << "using the following ntuple: " << inputFile << std::endl;
  myfile << "number of colliding bunches = " << numBunch << std::endl;
  myfile << "run luminosity = " << runLum << std::endl;
  myfile << "expected luminosity = " << expectedLum << std::endl;
  myfile << "norm factor used = " << norm << std::endl;
  myfile << "number of good events = " << goodLumiEventCount << std::endl;
  myfile.close(); 
}//closes the function 'rates'
