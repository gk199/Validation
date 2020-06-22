// Script for calculating rate histograms
// Originally from Aaron Bundock
// Edited by Gillian Kopp for multiplicity studies for LLP L1 trigger using HCAL depth and timing (2020)
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
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
  

void rates(bool newConditions, const std::string& inputFileDirectory){
  
  bool hwOn = true;   //are we using data from hardware? (upgrade trigger had to be running!!!)
  bool emuOn = true;  //are we using data from emulator?

  if (hwOn==false && emuOn==false)
    {
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

  std::ofstream EffRate_file;
  EffRate_file.open("Eff_Rates.txt", std::ios_base::app); // append instead of overwrite

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

  std::string axR = ";Threshold E_{T} (GeV);rate (Hz, unnormalized)";
  std::string axD = ";E_{T} (GeV);events/bin";
  std::string mult = ";Hit Multiplicity;Fraction of Entries (normalized)";

  //make histos
  TH1F* singleJetRates_emu = new TH1F("singleJetRates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* doubleJetRates_emu = new TH1F("doubleJetRates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* tripleJetRates_emu = new TH1F("tripleJetRates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* quadJetRates_emu = new TH1F("quadJetRates_emu", axR.c_str(), nJetBins, jetLo, jetHi);

  TH1F* singleJet_th1_Rates_emu = new TH1F("singleJet_th1_Rates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* singleJet_th2_Rates_emu = new TH1F("singleJet_th2_Rates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* singleJet_th3_Rates_emu = new TH1F("singleJet_th3_Rates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* singleJet_th4_Rates_emu = new TH1F("singleJet_th4_Rates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* singleJet_th5_Rates_emu = new TH1F("singleJet_th5_Rates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* singleJet_th6_Rates_emu = new TH1F("singleJet_th6_Rates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* singleJet_th7_Rates_emu = new TH1F("singleJet_th7_Rates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* quadJet_th1_Rates_emu = new TH1F("quadJet_th1_Rates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* quadJet_th2_Rates_emu = new TH1F("quadJet_th2_Rates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* quadJet_th3_Rates_emu = new TH1F("quadJet_th3_Rates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* quadJet_th4_Rates_emu = new TH1F("quadJet_th4_Rates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* quadJet_th5_Rates_emu = new TH1F("quadJet_th5_Rates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* quadJet_th6_Rates_emu = new TH1F("quadJet_th6_Rates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* quadJet_th7_Rates_emu = new TH1F("quadJet_th7_Rates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* htSum_th1_Rates_emu = new TH1F("htSum_th1_Rates_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
  TH1F* htSum_th2_Rates_emu = new TH1F("htSum_th2_Rates_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
  TH1F* htSum_th3_Rates_emu = new TH1F("htSum_th3_Rates_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
  TH1F* htSum_th4_Rates_emu = new TH1F("htSum_th4_Rates_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
  TH1F* htSum_th5_Rates_emu = new TH1F("htSum_th5_Rates_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
  TH1F* htSum_th6_Rates_emu = new TH1F("htSum_th6_Rates_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
  TH1F* htSum_th7_Rates_emu = new TH1F("htSum_th7_Rates_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
  // for mult and frac around one jet
  std::map<int, TH1F*> htSum_timing_1GeV_Rates_emu, htSum_timing_2GeV_Rates_emu, htSum_timing_3GeV_Rates_emu, htSum_timing_4GeV_Rates_emu;
  for (int TDCns=0; TDCns < 8; TDCns++) {
    htSum_timing_1GeV_Rates_emu[TDCns] = new TH1F(Form("htSum_timing_1GeV%dns_Rates_emu",TDCns),axR.c_str(), nHtSumBins, htSumLo, htSumHi);
    htSum_timing_2GeV_Rates_emu[TDCns] = new TH1F(Form("htSum_timing_2GeV%dns_Rates_emu",TDCns),axR.c_str(), nHtSumBins, htSumLo, htSumHi);
    htSum_timing_3GeV_Rates_emu[TDCns] = new TH1F(Form("htSum_timing_3GeV%dns_Rates_emu",TDCns),axR.c_str(), nHtSumBins, htSumLo, htSumHi);
    htSum_timing_4GeV_Rates_emu[TDCns] = new TH1F(Form("htSum_timing_4GeV%dns_Rates_emu",TDCns),axR.c_str(), nHtSumBins, htSumLo, htSumHi);
  }

  TH1F* singleJetGlobalRates_emu = new TH1F("singleJetGlobalRates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* doubleJetGlobalRates_emu = new TH1F("doubleJetGlobalRates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* tripleJetGlobalRates_emu = new TH1F("tripleJetGlobalRates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* quadJetGlobalRates_emu = new TH1F("quadJetGlobalRates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* singleJetOriginalRates_emu = new TH1F("singleJetOriginalRates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* quadJetOriginalRates_emu = new TH1F("quadJetOriginalRates_emu", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* singleEgRates_emu = new TH1F("singleEgRates_emu", axR.c_str(), nEgBins, egLo, egHi);
  TH1F* doubleEgRates_emu = new TH1F("doubleEgRates_emu", axR.c_str(), nEgBins, egLo, egHi);
  TH1F* singleTauRates_emu = new TH1F("singleTauRates_emu", axR.c_str(), nTauBins, tauLo, tauHi);
  TH1F* doubleTauRates_emu = new TH1F("doubleTauRates_emu", axR.c_str(), nTauBins, tauLo, tauHi);
  TH1F* singleISOEgRates_emu = new TH1F("singleISOEgRates_emu", axR.c_str(), nEgBins, egLo, egHi);
  TH1F* doubleISOEgRates_emu = new TH1F("doubleISOEgRates_emu", axR.c_str(), nEgBins, egLo, egHi);
  TH1F* singleISOTauRates_emu = new TH1F("singleISOTauRates_emu", axR.c_str(), nTauBins, tauLo, tauHi);
  TH1F* doubleISOTauRates_emu = new TH1F("doubleISOTauRates_emu", axR.c_str(), nTauBins, tauLo, tauHi);
  TH1F* htSum1jetRates_emu = new TH1F("htSum1jetRates_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
  TH1F* htSum4jetRates_emu = new TH1F("htSum4jetRates_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
  TH1F* htSumGlobalRates_emu = new TH1F("htSumGlobalRates_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
  TH1F* htSumOriginalRates_emu = new TH1F("htSumOriginalRates_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
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
  TH1F* singleJetOriginalRates_hw = new TH1F("singleJetOriginalRates_hw", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* quadJetOriginalRates_hw = new TH1F("quadJetOriginalRates_hw", axR.c_str(), nJetBins, jetLo, jetHi);
  TH1F* singleEgRates_hw = new TH1F("singleEgRates_hw", axR.c_str(), nEgBins, egLo, egHi);
  TH1F* doubleEgRates_hw = new TH1F("doubleEgRates_hw", axR.c_str(), nEgBins, egLo, egHi);
  TH1F* singleTauRates_hw = new TH1F("singleTauRates_hw", axR.c_str(), nTauBins, tauLo, tauHi);
  TH1F* doubleTauRates_hw = new TH1F("doubleTauRates_hw", axR.c_str(), nTauBins, tauLo, tauHi);
  TH1F* singleISOEgRates_hw = new TH1F("singleISOEgRates_hw", axR.c_str(), nEgBins, egLo, egHi);
  TH1F* doubleISOEgRates_hw = new TH1F("doubleISOEgRates_hw", axR.c_str(), nEgBins, egLo, egHi);
  TH1F* singleISOTauRates_hw = new TH1F("singleISOTauRates_hw", axR.c_str(), nTauBins, tauLo, tauHi);
  TH1F* doubleISOTauRates_hw = new TH1F("doubleISOTauRates_hw", axR.c_str(), nTauBins, tauLo, tauHi);
  TH1F* htSum1jetRates_hw = new TH1F("htSum1jetRates_hw",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
  TH1F* htSum4jetRates_hw = new TH1F("htSum4jetRates_hw",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
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

  // Delayed hit fraction 
  std::map<int, TH1F*> DelayedHitFraction1GeV, DelayedHitFraction2GeV, DelayedHitFraction3GeV, DelayedHitFraction4GeV;
  for (int TDCns=0; TDCns < 21; TDCns++) {
    // 1D histograms plotting the fraction of hits over energy and scanned time threshold
    DelayedHitFraction1GeV[TDCns] = new TH1F(Form("DelayedHitFraction1GeV%dns",TDCns),Form("Fraction of Hits that are > 1GeV and also > %dns;Fraction of Hits;Number of Events",TDCns),22,0,1.1);
    DelayedHitFraction2GeV[TDCns] = new TH1F(Form("DelayedHitFraction2GeV%dns",TDCns),Form("Fraction of Hits that are > 2GeV and also > %dns;Fraction of Hits;Number of Events",TDCns),22,0,1.1);
    DelayedHitFraction3GeV[TDCns] = new TH1F(Form("DelayedHitFraction3GeV%dns",TDCns),Form("Fraction of Hits that are > 3GeV and also > %dns;Fraction of Hits;Number of Events",TDCns),22,0,1.1);
    DelayedHitFraction4GeV[TDCns] = new TH1F(Form("DelayedHitFraction4GeV%dns",TDCns),Form("Fraction of Hits that are > 4GeV and also > %dns;Fraction of Hits;Number of Events",TDCns),22,0,1.1);
  }
  // 2D histograms plotting the fraction of hits over energy and scanned time threshold
  TH2F * DelayedHit2D_Fraction1GeV = new TH2F("DelayedHit2D_Fraction1GeV","Fraction of Delayed Hits vs. Delayed Timing Threshold above 1GeV;Delayed Timing Threshold;Fraction of Delayed Hits",30,0,30,50,0,1);
  TH2F * DelayedHit2D_Fraction2GeV = new TH2F("DelayedHit2D_Fraction2GeV","Fraction of Delayed Hits vs. Delayed Timing Threshold above 2GeV;Delayed Timing Threshold;Fraction of Delayed Hits",30,0,30,50,0,1);
  TH2F * DelayedHit2D_Fraction3GeV = new TH2F("DelayedHit2D_Fraction3GeV","Fraction of Delayed Hits vs. Delayed Timing Threshold above 3GeV;Delayed Timing Threshold;Fraction of Delayed Hits",30,0,30,50,0,1);
  TH2F * DelayedHit2D_Fraction4GeV = new TH2F("DelayedHit2D_Fraction4GeV","Fraction of Delayed Hits vs. Delayed Timing Threshold above 4GeV;Delayed Timing Threshold;Fraction of Delayed Hits",30,0,30,50,0,1);
  // 2D histograms plotting the number of hits over energy and time threshold
  TH2F * DelayedHit2D_Number1GeV = new TH2F("DelayedHit2D_Number1GeV","Number of Delayed Hits vs. Delayed Timing Threshold above 1GeV;Delayed Timing Threshold;Number of Delayed Hits",30,0,30,50,0,50);
  TH2F * DelayedHit2D_Number2GeV = new TH2F("DelayedHit2D_Number2GeV","Number of Delayed Hits vs. Delayed Timing Threshold above 2GeV;Delayed Timing Threshold;Number of Delayed Hits",30,0,30,50,0,50);
  TH2F * DelayedHit2D_Number3GeV = new TH2F("DelayedHit2D_Number3GeV","Number of Delayed Hits vs. Delayed Timing Threshold above 3GeV;Delayed Timing Threshold;Number of Delayed Hits",30,0,30,50,0,50);
  TH2F * DelayedHit2D_Number4GeV = new TH2F("DelayedHit2D_Number4GeV","Number of Delayed Hits vs. Delayed Timing Threshold above 4GeV;Delayed Timing Threshold;Number of Delayed Hits",30,0,30,50,0,50);


  // Number of TPs above energy threshold as a function of timing values
  //  TH2F * NumberTPtiming_3GeV = new TH2F("NumberTPtiming_3GeV","Number of TPs Above 3GeV vs. TDC Values;TDC Value;Number of TPs",20,0,20,100,0,100);
  std::map<int, TH1F*> NumberTPtiming, NumberTPavgtiming_jet1, NumberTPavgtiming_jet2, NumberTPavgtiming_jet3, NumberTPavgtiming_jet4, NumberTPavgtiming_jetMax, NumberTPtotalhits_jet1, NumberTPtotalhits_jet2, NumberTPtotalhits_jet3, NumberTPtotalhits_jet4, NumberTPtiming_depth1, NumberTPtiming_depth2, NumberTPtiming_depth3, NumberTPtiming_depth4;
  for (int GeV=0; GeV<6; GeV++) {
    NumberTPtiming[GeV] = new TH1F(Form("NumberTPtiming_%dGeV",GeV),Form("Number of Cells Above %dGeV vs. Timing Value;Timing (ns);Number of Cells",GeV),25,0,25);
    NumberTPavgtiming_jet1[GeV] = new TH1F(Form("NumberTPavgtiming_jet1_%dGeV",GeV),Form("Number of Events Above %dGeV vs. Average Timing Value;Average Timing (ns);Number of Events",GeV),25,0,25);
    NumberTPavgtiming_jet2[GeV] = new TH1F(Form("NumberTPavgtiming_jet2_%dGeV",GeV),Form("Number of Events Above %dGeV vs. Average Timing Value;Average Timing (ns);Number of Events",GeV),25,0,25);
    NumberTPavgtiming_jet3[GeV] = new TH1F(Form("NumberTPavgtiming_jet3_%dGeV",GeV),Form("Number of Events Above %dGeV vs. Average Timing Value;Average Timing (ns);Number of Events",GeV),25,0,25);
    NumberTPavgtiming_jet4[GeV] = new TH1F(Form("NumberTPavgtiming_jet4_%dGeV",GeV),Form("Number of Events Above %dGeV vs. Average Timing Value;Average Timing (ns);Number of Events",GeV),25,0,25);
    NumberTPavgtiming_jetMax[GeV] = new TH1F(Form("NumberTPavgtiming_jetMax_%dGeV",GeV),Form("Number of Events Above %dGeV vs. Average Timing Value;Average Timing (ns);Number of Events",GeV),25,0,25);
    NumberTPtotalhits_jet1[GeV] = new TH1F(Form("NumberTPtotalhits_jet1_%dGeV",GeV),Form("Number of Events Above %dGeV vs. Average Timing Value;Number hits;Number of Events",GeV),25,0,25);
    NumberTPtotalhits_jet2[GeV] = new TH1F(Form("NumberTPtotalhits_jet2_%dGeV",GeV),Form("Number of Events Above %dGeV vs. Average Timing Value;Number hits;Number of Events",GeV),25,0,25);
    NumberTPtotalhits_jet3[GeV] = new TH1F(Form("NumberTPtotalhits_jet3_%dGeV",GeV),Form("Number of Events Above %dGeV vs. Average Timing Value;Number hits;Number of Events",GeV),25,0,25);
    NumberTPtotalhits_jet4[GeV] = new TH1F(Form("NumberTPtotalhits_jet4_%dGeV",GeV),Form("Number of Events Above %dGeV vs. Average Timing Value;Number hits;Number of Events",GeV),25,0,25);
    NumberTPtiming_depth1[GeV] = new TH1F(Form("NumberTPtiming_%dGeV_depth1",GeV),Form("Number of Cells Above %dGeV vs. Timing Value;Timing (ns);Number of Cells",GeV),25,0,25);
    NumberTPtiming_depth2[GeV] = new TH1F(Form("NumberTPtiming_%dGeV_depth2",GeV),Form("Number of Cells Above %dGeV vs. Timing Value;Timing (ns);Number of Cells",GeV),25,0,25);
    NumberTPtiming_depth3[GeV] = new TH1F(Form("NumberTPtiming_%dGeV_depth3",GeV),Form("Number of Cells Above %dGeV vs. Timing Value;Timing (ns);Number of Cells",GeV),25,0,25);
    NumberTPtiming_depth4[GeV] = new TH1F(Form("NumberTPtiming_%dGeV_depth4",GeV),Form("Number of Cells Above %dGeV vs. Timing Value;Timing (ns);Number of Cells",GeV),25,0,25);
  }

  TH2F * NumberTPtiming_energy = new TH2F("NumberTPtiming_energy","Number of Cells vs. Energy vs. Timing Value;Timing (ns);Cell Energy Value",30,0,30,45,0,45);
  TH2F * NumberEvents_Fraction_Mult = new TH2F("NumberEvents_Fraction_Mult","Number of Events vs. Delayed Fraction vs. Delayed Hit Count;Delayed Fraction;Number of Delayed Hits",50,0,1,40,0,40);

  // Fraction of TPs above energy threshold as a function of timing values
  //  TH2F * FractionTPtiming_3GeV = new TH2F("FractionTPtiming_3GeV","Fraction of TPs Above 3GeV vs. TDC Values;TDC Value;Fraction of TPs",20,0,20,20,0,20);

  // 3 GeV energy cuts, scanning time cuts
  // inclusive
  std::map<int, TH1F*> dt3GeVMult_emu, dt3GeVHEMult_emu, dt3GeVHBMult_emu, dt2GeVMult_emu, dt2GeVHEMult_emu, dt2GeVHBMult_emu, dt1GeVMult_emu, dt1GeVHEMult_emu, dt1GeVHBMult_emu;
  std::map<int, TH1F*> dt0GeVHBJetMult_emu, dt1GeVJetMult_emu, dt1GeVHEJetMult_emu, dt1GeVHBJetMult_emu, dt2GeVJetMult_emu, dt2GeVHEJetMult_emu, dt2GeVHBJetMult_emu, dt3GeVJetMult_emu, dt3GeVHEJetMult_emu, dt3GeVHBJetMult_emu;
  for (int TDCns=1; TDCns <= 5; TDCns++) {
    // inclusive, no L1 jet matching
    // 1GeV
    dt1GeVMult_emu[TDCns] = new TH1F(Form("dt1GeV%dnsMult_emu",TDCns),Form("Multiplicity of %dns delayed cells above 1 GeV (inclusive);Hit Multiplicity;Fraction of Entries (normalized)",TDCns),400,0,400);
    dt1GeVHEMult_emu[TDCns] = new TH1F(Form("dt1GeV%dnsHEMult_emu",TDCns),Form("Multiplicity of %dns delayed cells above 1 GeV (HE);Hit Multiplicity;Fraction of Entries (normalized)",TDCns),400,0,400);
    dt1GeVHBMult_emu[TDCns] = new TH1F(Form("dt1GeV%dnsHBMult_emu",TDCns),Form("Multiplicity of %dns delayed cells above 1 GeV (HB);Hit Multiplicity;Fraction of Entries (normalized)",TDCns),400,0,400);
    // 2GeV
    dt2GeVMult_emu[TDCns] = new TH1F(Form("dt2GeV%dnsMult_emu",TDCns),Form("Multiplicity of %dns delayed cells above 2 GeV (inclusive);Hit Multiplicity;Fraction of Entries (normalized)",TDCns),200,0,200);
    dt2GeVHEMult_emu[TDCns] = new TH1F(Form("dt2GeV%dnsHEMult_emu",TDCns),Form("Multiplicity of %dns delayed cells above 2 GeV (HE);Hit Multiplicity;Fraction of Entries (normalized)",TDCns),200,0,200);
    dt2GeVHBMult_emu[TDCns] = new TH1F(Form("dt2GeV%dnsHBMult_emu",TDCns),Form("Multiplicity of %dns delayed cells above 2 GeV (HB);Hit Multiplicity;Fraction of Entries (normalized)",TDCns),200,0,200);
    // 3GeV
    dt3GeVMult_emu[TDCns] = new TH1F(Form("dt3GeV%dnsMult_emu",TDCns),Form("Multiplicity of %dns delayed cells above 3 GeV (inclusive);Hit Multiplicity;Fraction of Entries (normalized)",TDCns),120,0,120);
    dt3GeVHEMult_emu[TDCns] = new TH1F(Form("dt3GeV%dnsHEMult_emu",TDCns),Form("Multiplicity of %dns delayed cells above 3 GeV (HE);Hit Multiplicity;Fraction of Entries (normalized)",TDCns),120,0,120);
    dt3GeVHBMult_emu[TDCns] = new TH1F(Form("dt3GeV%dnsHBMult_emu",TDCns),Form("Multiplicity of %dns delayed cells above 3 GeV (HB);Hit Multiplicity;Fraction of Entries (normalized)",TDCns),120,0,120);

    // for HCAL TP matched with L1 jet
    // 0GeV
    dt0GeVHBJetMult_emu[TDCns] = new TH1F(Form("dt0GeV%dnsHBJetMult_emu",TDCns),Form("Multiplicity of %dns delayed cells above 0 GeV (HB, match TP with L1Jet);Hit Multiplicity;Fraction of Entries (normalized)",TDCns),120,0,120); // HB
    // 1GeV
    dt1GeVJetMult_emu[TDCns] = new TH1F(Form("dt1GeV%dnsJetMult_emu",TDCns),Form("Multiplicity of %dns delayed cells above 1 GeV (inclusive, match TP with L1Jet);Hit Multiplicity;Fraction of Entries (normalized)",TDCns),400,0,400); // inclusive
    dt1GeVHEJetMult_emu[TDCns] = new TH1F(Form("dt1GeV%dnsHEJetMult_emu",TDCns),Form("Multiplicity of %dns delayed cells above 1 GeV (HE, match TP with L1Jet);Hit Multiplicity;Fraction of Entries (normalized)",TDCns),400,0,400); // HE
    dt1GeVHBJetMult_emu[TDCns] = new TH1F(Form("dt1GeV%dnsHBJetMult_emu",TDCns),Form("Multiplicity of %dns delayed cells above 1 GeV (HB, match TP with L1Jet);Hit Multiplicity;Fraction of Entries (normalized)",TDCns),400,0,400); // HB  
    // 2GeV
    dt2GeVJetMult_emu[TDCns] = new TH1F(Form("dt2GeV%dnsJetMult_emu",TDCns),Form("Multiplicity of %dns delayed cells above 2 GeV (inclusive, match TP with L1Jet);Hit Multiplicity;Fraction of Entries (normalized)",TDCns),200,0,200); // inclusive 
    dt2GeVHEJetMult_emu[TDCns] = new TH1F(Form("dt2GeV%dnsHEJetMult_emu",TDCns),Form("Multiplicity of %dns delayed cells above 2 GeV (HE, match TP with L1Jet);Hit Multiplicity;Fraction of Entries (normalized)",TDCns),200,0,120); // HE
    dt2GeVHBJetMult_emu[TDCns] = new TH1F(Form("dt2GeV%dnsHBJetMult_emu",TDCns),Form("Multiplicity of %dns delayed cells above 2 GeV (HB, match TP with L1Jet);Hit Multiplicity;Fraction of Entries (normalized)",TDCns),200,0,200); // HB  
    // 3GeV
    dt3GeVJetMult_emu[TDCns] = new TH1F(Form("dt3GeV%dnsJetMult_emu",TDCns),Form("Multiplicity of %dns delayed cells above 3 GeV (inclusive, match TP with L1Jet);Hit Multiplicity;Fraction of Entries (normalized)",TDCns),120,0,120); // inclusive
    dt3GeVHEJetMult_emu[TDCns] = new TH1F(Form("dt3GeV%dnsHEJetMult_emu",TDCns),Form("Multiplicity of %dns delayed cells above 3 GeV (HE, match TP with L1Jet);Hit Multiplicity;Fraction of Entries (normalized)",TDCns),120,0,120); // HE
    dt3GeVHBJetMult_emu[TDCns] = new TH1F(Form("dt3GeV%dnsHBJetMult_emu",TDCns),Form("Multiplicity of %dns delayed cells above 3 GeV (HB, match TP with L1Jet);Hit Multiplicity;Fraction of Entries (normalized)",TDCns),120,0,120); // HB
  }

  // include all HCAL TPs included in DR cone of a L1 jet
  TH1F * dt3GeV3nsHBQuadJetMult_emu = new TH1F("dt3GeV3nsHBQuadJetMult_emu","Quad Jet Sum Multiplicity of 3ns delayed cells above 3 GeV (HB, dR(TP,L1Jet)<0.5);Hit Multiplicity;Fraction of Entries (normalized)",120,0,120);
  TH1F * dt3GeV3nsHBTripleJetMult_emu = new TH1F("dt3GeV3nsHBTripleJetMult_emu","Triple Jet Sum Multiplicity of 3ns delayed cells above 3 GeV (HB, dR(TP,L1Jet)<0.5);Hit Multiplicity;Fraction of Entries (normalized)",120,0,120);
  TH1F * dt3GeV3nsHBDoubleJetMult_emu = new TH1F("dt3GeV3nsHBDoubleJetMult_emu","Double Jet Sum Multiplicity of 3ns delayed cells above 3 GeV (HB, dR(TP,L1Jet)<0.5);Hit Multiplicity;Fraction of Entries (normalized)",120,0,120);
  TH1F * dt3GeV3nsHBSingleJetMult_emu = new TH1F("dt3GeV3nsHBSingleJetMult_emu","Single Jet Sum Multiplicity of 3ns delayed cells above 3 GeV (HB, dR(TP,L1Jet)<0.5);Hit Multiplicity;Fraction of Entries (normalized)",120,0,120);
  TH1F * dt3GeV3nsHBJet1Mult_emu = new TH1F("dt3GeV3nsHBJet1Mult_emu","Multiplicity of 3ns delayed cells above 3 GeV (HB, all TPs within DR<0.5 1st L1Jet);Hit Multiplicity;Fraction of Entries (normalized)",120,0,120);
  TH1F * dt3GeV3nsHBJet2Mult_emu = new TH1F("dt3GeV3nsHBJet2Mult_emu","Multiplicity of 3ns delayed cells above 3 GeV (HB, all TPs within DR<0.5 2nd L1Jet);Hit Multiplicity;Fraction of Entries (normalized)",120,0,120);
  TH1F * dt3GeV3nsHBJet3Mult_emu = new TH1F("dt3GeV3nsHBJet3Mult_emu","Multiplicity of 3ns delayed cells above 3 GeV (HB, all TPs within DR<0.5 3rd L1Jet);Hit Multiplicity;Fraction of Entries (normalized)",120,0,120);
  TH1F * dt3GeV3nsHBJet4Mult_emu = new TH1F("dt3GeV3nsHBJet4Mult_emu","Multiplicity of 3ns delayed cells above 3 GeV (HB, all TPs within DR<0.5 4th L1Jet);Hit Multiplicity;Fraction of Entries (normalized)",120,0,120);
  TH1F * dt3GeV2nsHBQuadJet_depth1_Mult_emu = new TH1F("dt3GeV2nsHBQuadJet_depth1_Mult_emu","Quad Jet Multiplicity of 2ns delayed cells above 3 GeV (HB, dR(TP,L1Jet)<0.5, depth 1);Hit Multiplicity;Fraction of Entries (normalized)",120,0,120);
  TH1F * dt3GeV2nsHBQuadJet_depth2_Mult_emu = new TH1F("dt3GeV2nsHBQuadJet_depth2_Mult_emu","Quad Jet Multiplicity of 2ns delayed cells above 3 GeV (HB, dR(TP,L1Jet)<0.5, depth 2);Hit Multiplicity;Fraction of Entries (normalized)",120,0,120);
  TH1F * dt3GeV2nsHBQuadJet_depth3_Mult_emu = new TH1F("dt3GeV2nsHBQuadJet_depth3_Mult_emu","Quad Jet Multiplicity of 2ns delayed cells above 3 GeV (HB, dR(TP,L1Jet)<0.5, depth 3);Hit Multiplicity;Fraction of Entries (normalized)",120,0,120);
  TH1F * dt3GeV2nsHBQuadJet_depth4_Mult_emu = new TH1F("dt3GeV2nsHBQuadJet_depth4_Mult_emu","Quad Jet Multiplicity of 2ns delayed cells above 3 GeV (HB, dR(TP,L1Jet)<0.5, depth 4);Hit Multiplicity;Fraction of Entries (normalized)",120,0,120);
  TH1F * dt3GeV2nsHBHEQuadJet_depth1_Mult_emu = new TH1F("dt3GeV2nsHBHEQuadJet_depth1_Mult_emu","Quad Jet Multiplicity of 2ns delayed cells above 3 GeV (HBHE, dR(TP,L1Jet)<0.5, depth 1);Hit Multiplicity;Fraction of Entries (normalized)",120,0,120);
  TH1F * dt3GeV2nsHBHEQuadJet_depth2_Mult_emu = new TH1F("dt3GeV2nsHBHEQuadJet_depth2_Mult_emu","Quad Jet Multiplicity of 2ns delayed cells above 3 GeV (HBHE, dR(TP,L1Jet)<0.5, depth 2);Hit Multiplicity;Fraction of Entries (normalized)",120,0,120);
  TH1F * dt3GeV2nsHBHEQuadJet_depth3_Mult_emu = new TH1F("dt3GeV2nsHBHEQuadJet_depth3_Mult_emu","Quad Jet Multiplicity of 2ns delayed cells above 3 GeV (HBHE, dR(TP,L1Jet)<0.5, depth 3);Hit Multiplicity;Fraction of Entries (normalized)",120,0,120);
  TH1F * dt3GeV2nsHBHEQuadJet_depth4_Mult_emu = new TH1F("dt3GeV2nsHBHEQuadJet_depth4_Mult_emu","Quad Jet Multiplicity of 2ns delayed cells above 3 GeV (HBHE, dR(TP,L1Jet)<0.5, depth 4);Hit Multiplicity;Fraction of Entries (normalized)",120,0,120);

  // include all HCAL TPs included in DR cone of a L1 jet
  TH1F * dt0GeV5nsHBQuadJetMult_emu = new TH1F("dt0GeV5nsHBQuadJetMult_emu","Quad Jet Sum Multiplicity of 5ns delayed cells above 0 GeV (HB, dR(TP,L1Jet)<0.5);Hit Multiplicity;Fraction of Entries (normalized)",120,0,120);
  TH1F * dt0GeV5nsHBTripleJetMult_emu = new TH1F("dt0GeV5nsHBTripleJetMult_emu","Triple Jet Sum Multiplicity of 5ns delayed cells above 0 GeV (HB, dR(TP,L1Jet)<0.5);Hit Multiplicity;Fraction of Entries (normalized)",120,0,120);
  TH1F * dt0GeV5nsHBDoubleJetMult_emu = new TH1F("dt0GeV5nsHBDoubleJetMult_emu","Double Jet Sum Multiplicity of 5ns delayed cells above 0 GeV (HB, dR(TP,L1Jet)<0.5);Hit Multiplicity;Fraction of Entries (normalized)",120,0,120);
  TH1F * dt0GeV5nsHBSingleJetMult_emu = new TH1F("dt0GeV5nsHBSingleJetMult_emu","Single Jet Sum Multiplicity of 5ns delayed cells above 0 GeV (HB, dR(TP,L1Jet)<0.5);Hit Multiplicity;Fraction of Entries (normalized)",120,0,120);
  TH1F * dt0GeV5nsHBJet1Mult_emu = new TH1F("dt0GeV5nsHBJet1Mult_emu","Multiplicity of 5ns delayed cells above 0 GeV (HB, all TPs within DR<2 1st L1Jet);Hit Multiplicity;Fraction of Entries (normalized)",120,0,120);
  TH1F * dt0GeV5nsHBJet2Mult_emu = new TH1F("dt0GeV5nsHBJet2Mult_emu","Multiplicity of 5ns delayed cells above 0 GeV (HB, all TPs within DR<2 2nd L1Jet);Hit Multiplicity;Fraction of Entries (normalized)",120,0,120);
  TH1F * dt0GeV5nsHBJet3Mult_emu = new TH1F("dt0GeV5nsHBJet3Mult_emu","Multiplicity of 5ns delayed cells above 0 GeV (HB, all TPs within DR<2 3rd L1Jet);Hit Multiplicity;Fraction of Entries (normalized)",120,0,120);
  TH1F * dt0GeV5nsHBJet4Mult_emu = new TH1F("dt0GeV5nsHBJet4Mult_emu","Multiplicity of 5ns delayed cells above 0 GeV (HB, all TPs within DR<2 4th L1Jet);Hit Multiplicity;Fraction of Entries (normalized)",120,0,120);
  TH1F * DepthVariable = new TH1F("DepthVariable","Depth 4 / Depth 1;Hit Multiplicity; Fraction of Entries (normalized)",50,0,50);

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
  TH2F * Energy_DepthHE_Jets_1ns = new TH2F("Energy_DepthHE_Jets_1ns", "TP Energy Fraction vs. Depth in HE (TP matched w/ L1 Jets, hits over 1ns);HCAL Depth;Energy Fraction", 8, -0.5, 7.5, 60, 0, 20);
  TH2F * Energy_DepthHE_Jets_2ns = new TH2F("Energy_DepthHE_Jets_2ns", "TP Energy Fraction vs. Depth in HE (TP matched w/ L1 Jets, hits over 2ns);HCAL Depth;Energy Fraction", 8, -0.5, 7.5, 60, 0, 20);
  TH2F * Timing_DepthHE_Jets = new TH2F("Timing_DepthHE_Jets", "TP Timing Value vs. Depth in HE (TP matched w/ L1 Jets);HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5, 60, 0, 30);
  TH2F * Energy_DepthHB_Jets = new TH2F("Energy_DepthHB_Jets", "TP Energy Fraction vs. Depth in HB (TP matched w/ L1 Jets);HCAL Depth;Energy Fraction", 8, -0.5, 7.5, 60, 0, 1.2);
  TH2F * Energy_DepthHB_Jets_1ns = new TH2F("Energy_DepthHB_Jets_1ns", "TP Energy Fraction vs. Depth in HB (TP matched w/ L1 Jets, hits over 1ns);HCAL Depth;Energy Fraction", 8, -0.5, 7.5, 60, 0, 20);
  TH2F * Energy_DepthHB_Jets_2ns = new TH2F("Energy_DepthHB_Jets_2ns", "TP Energy Fraction vs. Depth in HB (TP matched w/ L1 Jets, hits over 2ns);HCAL Depth;Energy Fraction", 8, -0.5, 7.5, 60, 0, 20);
  TH2F * Timing_DepthHB_Jets = new TH2F("Timing_DepthHB_Jets", "TP Timing Value vs. Depth in HB (TP matched w/ L1 Jets);HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5, 60, 0, 30);
  // TH2F for energy depth, where matched to jets for high energy TPs
  TH2F * Energy_Depth_Jets_HighE = new TH2F("Energy_Depth_Jets_HighE", "TP Energy Fraction vs. Depth for TP E_{T} > 10 GeV (TP matched w/ L1 Jets);HCAL Depth;Energy Fraction", 8, -0.5, 7.5, 60, 0, 1.2);
  TH2F * Energy_DepthHE_Jets_HighE = new TH2F("Energy_DepthHE_Jets_HighE", "TP Energy Fraction vs. Depth in HE for TP E_{T} > 10 GeV (TP matched w/ L1 Jets);HCAL Depth;Energy Fraction", 8, -0.5, 7.5, 60, 0, 1.2);
  TH2F * Energy_DepthHB_Jets_HighE = new TH2F("Energy_DepthHB_Jets_HighE", "TP Energy Fraction vs. Depth in HB for TP E_{T} > 10 GeV (TP matched w/ L1 Jets);HCAL Depth;Energy Fraction", 8, -0.5, 7.5, 60, 0, 1.2); 
  TH2F * VolTiming_Depth_Jets = new TH2F("VolTiming_Depth_Jets", "TP Timing Value vs. Depth in HCAL Volume (TP matched w/ L1 Jets);HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5, 60, 0, 30);
  TH2F * VolTiming_DepthHB_Jets = new TH2F("VolTiming_DepthHB_Jets", "TP Timing Value vs. Depth in HCAL Volume, HB (TP matched w/ L1 Jets);HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5, 60, 0, 30);
  TH2F * VolTiming_DepthHE_Jets = new TH2F("VolTiming_DepthHE_Jets", "TP Timing Value vs. Depth in HCAL Volume, HE (TP matched w/ L1 Jets);HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5, 60, 0, 30);

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
  TH1D * Energy_DepthHE_avg_Jets_1ns = new TH1D("Energy_DepthHE_avg_Jets_1ns", "TP Avg Energy Fraction vs. Depth in HE (TP matched w/ L1 Jets, hits over 1ns);HCAL Depth;Energy Fraction", 8, -0.5, 7.5);
  TH1D * Energy_DepthHE_avg_Jets_2ns = new TH1D("Energy_DepthHE_avg_Jets_2ns", "TP Avg Energy Fraction vs. Depth in HE (TP matched w/ L1 Jets, hits over 2ns);HCAL Depth;Energy Fraction", 8, -0.5, 7.5);
  TH1D * Timing_DepthHE_avg_Jets = new TH1D("Timing_DepthHE_avg_Jets", "TP Avg Timing Value vs. Depth in HE (TP matched w/ L1 Jets);HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5);
  TH1D * Energy_DepthHB_avg_Jets = new TH1D("Energy_DepthHB_avg_Jets", "TP Avg Energy Fraction vs. Depth in HB (TP matched w/ L1 Jets);HCAL Depth;Energy Fraction", 8, -0.5, 7.5);
  TH1D * Energy_DepthHB_avg_Jets_1ns = new TH1D("Energy_DepthHB_avg_Jets_1ns", "TP Avg Energy Fraction vs. Depth in HB (TP matched w/ L1 Jets, hits over 1ns);HCAL Depth;Energy Fraction", 8, -0.5, 7.5);
  TH1D * Energy_DepthHB_avg_Jets_2ns = new TH1D("Energy_DepthHB_avg_Jets_2ns", "TP Avg Energy Fraction vs. Depth in HB (TP matched w/ L1 Jets, hits over 2ns);HCAL Depth;Energy Fraction", 8, -0.5, 7.5);
  TH1D * Timing_DepthHB_avg_Jets = new TH1D("Timing_DepthHB_avg_Jets", "TP Avg Timing Value vs. Depth in HB (TP matched w/ L1 Jets);HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5);
  TH1D * VolTiming_Depth_avg_Jets = new TH1D("VolTiming_Depth_avg_Jets", "TP Avg Timing Value vs. Depth  in HCAL Volume (TP matched w/ L1 Jets);HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5);
  TH1D * VolTiming_DepthHE_avg_Jets = new TH1D("VolTiming_DepthHE_avg_Jets", "TP Avg Timing Value vs. Depth in HCAL Volume, HE (TP matched w/ L1 Jets);HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5);
  TH1D * VolTiming_DepthHB_avg_Jets = new TH1D("VolTiming_DepthHB_avg_Jets", "TP Avg Timing Value vs. Depth in HCAL Volume, HB (TP matched w/ L1 Jets);HCAL Depth;Timing Value (ns)", 8, -0.5, 7.5);

  // Ratio of energy in HCAL depth layers
  TH1F * Ratio_Depth = new TH1F("Ratio_Depth", "Ratio of 3,4 HCAL Layers to E_{T};Ratio;Number of Events", 50,0,1);
  TH1F * Ratio_DepthHE = new TH1F("Ratio_DepthHE", "Ratio of 3,4 HCAL Layers to E_{T} in HE;Ratio;Number of Events", 50,0,1);
  TH1F * Ratio_DepthHB = new TH1F("Ratio_DepthHB", "Ratio of 3,4 HCAL Layers to E_{T} in HB;Ratio;Number of Events", 50,0,1);
  TH1F * Ratio_Depth_Jets = new TH1F("Ratio_Depth_Jets", "Ratio of 3,4 HCAL Layers to E_{T}, matched w/Jets;Ratio;Number of Events", 50,0,1);
  TH1F * Ratio_DepthHE_Jets = new TH1F("Ratio_DepthHE_Jets", "Ratio of 3,4 HCAL Layers to E_{T} in HE, matched w/Jets;Ratio;Number of Events", 50,0,1);
  TH1F * Ratio_DepthHB_Jets = new TH1F("Ratio_DepthHB_Jets", "Ratio of 3,4 HCAL Layers to E_{T} in HB, matched w/Jets;Ratio;Number of Events", 50,0,1);

  // delta R plot for HCAL TP max energy near L1 jet
  TH1F * DeltaR_TP_L1Jet_1 = new TH1F("DeltaR_TP_L1Jet_1", "DeltaR Between First L1 Jet and Closest HCAL TP; DeltaR;Number of Events",50,0,0.5);
  TH1F * DeltaR_TP_L1Jet_2 = new TH1F("DeltaR_TP_L1Jet_2", "DeltaR Between Second L1 Jet and Closest HCAL TP; DeltaR;Number of Events",50,0,0.5);
  TH1F * DeltaR_TP_L1Jet_3 = new TH1F("DeltaR_TP_L1Jet_3", "DeltaR Between Third L1 Jet and Closest HCAL TP; DeltaR;Number of Events",50,0,0.5);
  TH1F * DeltaR_TP_L1Jet_4 = new TH1F("DeltaR_TP_L1Jet_4", "DeltaR Between Fourth L1 Jet and Closest HCAL TP; DeltaR;Number of Events",50,0,0.5);
  TH1F * DeltaR_partons_L1Jet_1 = new TH1F("DeltaR_partons_L1Jet_1", "DeltaR Between First L1 Jet and Closest Parton; DeltaR;Number of Events",50,0,5);
  TH1F * DeltaR_partons_L1Jet_2 = new TH1F("DeltaR_partons_L1Jet_2", "DeltaR Between Second L1 Jet and Closest Parton; DeltaR;Number of Events",50,0,5);
  TH1F * DeltaR_partons_L1Jet_3 = new TH1F("DeltaR_partons_L1Jet_3", "DeltaR Between Third L1 Jet and Closest Parton; DeltaR;Number of Events",50,0,5);
  TH1F * DeltaR_partons_L1Jet_4 = new TH1F("DeltaR_partons_L1Jet_4", "DeltaR Between Fourth L1 Jet and Closest Parton; DeltaR;Number of Events",50,0,5);
  TH1F * DeltaR_partonsHCAL_L1Jet_1 = new TH1F("DeltaR_partonsHCAL_L1Jet_1", "DeltaR Between First L1 Jet and Closest Parton with HCAL volume vertex; DeltaR;Number of Events",50,0,5);
  TH1F * DeltaR_partonsHCAL_L1Jet_2 = new TH1F("DeltaR_partonsHCAL_L1Jet_2", "DeltaR Between Second L1 Jet and Closest Parton with HCAL volume vertex; DeltaR;Number of Events",50,0,5);
  TH1F * DeltaR_partonsHCAL_L1Jet_3 = new TH1F("DeltaR_partonsHCAL_L1Jet_3", "DeltaR Between Third L1 Jet and Closest Parton with HCAL volume vertex; DeltaR;Number of Events",50,0,5);
  TH1F * DeltaR_partonsHCAL_L1Jet_4 = new TH1F("DeltaR_partonsHCAL_L1Jet_4", "DeltaR Between Fourth L1 Jet and Closest Parton with HCAL volume vertex; DeltaR;Number of Events",50,0,5);
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
  // plot how many LLPs decay in HCAL
  TH1F * LLPdecayDetAcceptance = new TH1F("LLPdecayDetAcceptance", "LLPs decaying within detector acceptance; Number of LLPs;Number of Events",5,0,5);
  TH2F * LLPdecayRadiusDetAcceptance = new TH2F("LLPdecayRadiusDetAcceptance", "Decay Radius for LLPs within detector acceptance;Decay Position (z); Decay Radius (cm)",100,0,600,50,0,300);
  TH1F * LLPdecayXyzDetAcceptance = new TH1F("LLPdecayXyzDetAcceptance", "LLPs decaying within detector acceptance;Decay Radius (x,y,z);Number of Events",100,0,600);
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

  TH1F * htSumDistribution = new TH1F("htSumDistribution","htSum Distribution;HT Sum;Number of Events",100,0,1000);

  // TGraph for eta phi of HCAL TPs and L1 jets
  TGraph * etaphiTP = new TGraph();
  TGraph * etaphiJet = new TGraph();
  
  // Create a TTree Object so multiplicities can be sent to TMVA analyzer
  TTree *tree = new TTree("MultForTMVA","MultForTMVA");
  //  Int_t event;
  //  Float_t mult_Jet1, mult_Jet2, mult_Jet3, mult_Jet4;
  //  Float_t event_htSum, AllHits2GeV, AllHits2GeV_jet2, AllHits2GeV_jet3, AllHits2GeV_jet4;
  //  Float_t DelayedHits2GeV, DelayedHits2GeV_jet2, DelayedHits2GeV_jet3, DelayedHits2GeV_jet4;
  //  Float_t Jet1eta, Jet2eta, Jet3eta, Jet4eta;
  //  Float_t  ET_Jet1(0); 
  //  Int_t event;
  Float_t x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20;
  tree->Branch("x0",&x0,"x0/F"); // mult_Jet1
  tree->Branch("x1",&x1,"x1/F"); // mult_Jet2
  tree->Branch("x2",&x2,"x2/F"); // mult_Jet3
  tree->Branch("x3",&x3,"x3/F"); // mult_Jet4
  tree->Branch("x4",&x4,"x4/F"); // event_htSum
  tree->Branch("x5",&x5,"x5/F"); // AllHits2GeV
  tree->Branch("x6",&x6,"x6/F"); // AllHits2GeV_jet2
  tree->Branch("x7",&x7,"x7/F"); // AllHits2GeV_jet3
  tree->Branch("x8",&x8,"x8/F"); // AllHits2GeV_jet4
  tree->Branch("x9",&x9,"x9/F"); // DelayedHits2GeV
  tree->Branch("x10",&x10,"x10/F"); // DelayedHits2GeV_jet2
  tree->Branch("x11",&x11,"x11/F"); // DelayedHits2GeV_jet3
  tree->Branch("x12",&x12,"x12/F"); // DelayedHits2GeV_jet4
  tree->Branch("x13",&x13,"x13/F"); // Jet1eta
  tree->Branch("x14",&x14,"x14/F"); // Jet2eta
  tree->Branch("x15",&x15,"x15/F"); // Jet3eta
  tree->Branch("x16",&x16,"x16/F"); // Jet4eta
  tree->Branch("x17",&x17,"x17/F"); // jetEt_1
  tree->Branch("x18",&x18,"x18/F"); // event
  tree->Branch("x19",&x19,"x19/F"); // 3GeV 2ns HE
  tree->Branch("x20",&x20,"x20/F"); // 3GeV 3ns HE
  /*
  tree->Branch("mult_Jet1",&mult_Jet1,"mult_Jet1/F");
  tree->Branch("mult_Jet2",&mult_Jet2,"mult_Jet2/F");
  tree->Branch("mult_Jet3",&mult_Jet3,"mult_Jet3/F");
  tree->Branch("mult_Jet4",&mult_Jet4,"mult_Jet4/F");
  tree->Branch("event_htSum",&event_htSum,"event_htSum/F");
  tree->Branch("AllHits2GeV",&AllHits2GeV,"AllHits2GeV/F");
  tree->Branch("AllHits2GeV_jet2",&AllHits2GeV_jet2,"AllHits2GeV_jet2/F");
  tree->Branch("AllHits2GeV_jet3",&AllHits2GeV_jet3,"AllHits2GeV_jet3/F");
  tree->Branch("AllHits2GeV_jet4",&AllHits2GeV_jet4,"AllHits2GeV_jet4/F");
  tree->Branch("DelayedHits2GeV",&DelayedHits2GeV,"DelayedHits2GeV/F");
  tree->Branch("DelayedHits2GeV_jet2",&DelayedHits2GeV_jet2,"DelayedHits2GeV_jet2/F");
  tree->Branch("DelayedHits2GeV_jet3",&DelayedHits2GeV_jet3,"DelayedHits2GeV_jet3/F");
  tree->Branch("DelayedHits2GeV_jet4",&DelayedHits2GeV_jet4,"DelayedHits2GeV_jet4/F");
  tree->Branch("Jet1eta",&Jet1eta,"Jet1eta/F");
  tree->Branch("Jet2eta",&Jet2eta,"Jet2eta/F");
  tree->Branch("Jet3eta",&Jet3eta,"Jet3eta/F");
  tree->Branch("Jet4eta",&Jet4eta,"Jet4eta/F");
  tree->Branch("ET_Jet1",&ET_Jet1,"ET_Jet1/F");
  tree->Branch("event",&event,"event/I");
  */
  // counting LLP efficiencies -- need to reset on a per event basis
  double passedAvgTimeCut(0);
  double totalJets(0), passedMultJets(0), passedMultJets_120(0), passedMultJets_350(0);
  double passedMultJets3GeV3_1(0), passedMultJets3GeV3_2(0), passedMultJets3GeV3_3(0), passedMultJets3GeV3_4(0), passedMultJets3GeV3_5(0), passedMultJets3GeV3_6(0), passedMultJets3GeV3_7(0);
  double passedMultJets3GeV3_ht120_1(0), passedMultJets3GeV3_ht120_2(0), passedMultJets3GeV3_ht120_3(0), passedMultJets3GeV3_ht120_4(0), passedMultJets3GeV3_ht120_5(0), passedMultJets3GeV3_ht120_6(0), passedMultJets3GeV3_ht120_7(0);
  double totalGlobal(0), passedMultGlobal(0), passedMultGlobal_120(0), passedMultGlobal_350(0), passedHtSum120(0), passedHtSum350(0), passedHtSum360(0);
  double passedDelayedHitFraction1GeV_ht120[11][11][25] = {{{0}}}; // [frac delayed][number delayed][ns delayed]
  double passedDelayedHitFraction2GeV_ht120[11][11][25] = {{{0}}};
  double passedDelayedHitFraction3GeV_ht120[11][11][25] = {{{0}}};
  double passedDelayedHitFraction4GeV_ht120[11][11][25] = {{{0}}};
  // change these variables to change the multiplicity thresholds that are set
  uint GeV3ns3Global_threshold = 3;
  //  uint GeV3ns0Jet_threshold = 3;
  uint GeV3ns3Jet_threshold = 2;
  double DR_threshold = 0.5;
  double frac_delayed_scan[11] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
  double frac_delayed = 0.6;
  uint frac_delayed_x10 = 6;
  uint min_num_delayed_scan[11] = {0,1,2,3,4,5,6,7,8,9,10};
  uint min_num_delayed = 6;
  uint min_num_delayed_jetsum = 0;
  int min_num_jets_passed = 1;

  /////////////////////////////////
  // loop through all the entries//
  /////////////////////////////////
  for (Long64_t jentry=0; jentry<nentries; jentry++){
    //    std::cout << "new event" << std::endl;
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

      double mult3GeV[6] = {0};
      double mult3GeVHE[6] = {0};
      double mult3GeVHB[6] = {0};
      double mult2GeV[6] = {0};
      double mult2GeVHE[6] = {0};
      double mult2GeVHB[6] = {0};
      double mult1GeV[6] = {0};
      double mult1GeVHE[6] = {0};
      double mult1GeVHB[6] = {0};

      // multiplicity for ieta regions of caloTowers (4x4 ieta iphi) in the barrel regions                    
      double mult1GeVcaloT1(0), mult2GeVcaloT1(0), mult3GeVcaloT1(0); // abs(ieta) between 1-4
      double mult1GeVcaloT2(0), mult2GeVcaloT2(0), mult3GeVcaloT2(0); // abs(ieta) between 5-8
      double mult1GeVcaloT3(0), mult2GeVcaloT3(0), mult3GeVcaloT3(0); // abs(ieta) between 9-12
      double mult1GeVcaloT4(0), mult2GeVcaloT4(0), mult3GeVcaloT4(0); // abs(ieta) between 13-16

      // multiplicity for when HCAL TP is matched with Jets
      double mult3GeV2ns_Jets_depth1(0), mult3GeV2ns_Jets_depth2(0), mult3GeV2ns_Jets_depth3(0), mult3GeV2ns_Jets_depth4(0);
      double mult3GeV2ns_Jets_depth1HB(0), mult3GeV2ns_Jets_depth2HB(0), mult3GeV2ns_Jets_depth3HB(0), mult3GeV2ns_Jets_depth4HB(0);
      double mult3GeV_Jets[6] = {0};
      double mult3GeVHE_Jets[6]= {0};
      double mult3GeVHB_Jets[6]= {0};
      double mult2GeV_Jets[6]= {0};
      double mult2GeVHE_Jets[6]= {0};
      double mult2GeVHB_Jets[6]= {0};
      double mult1GeV_Jets[6]= {0};
      double mult1GeVHE_Jets[6]= {0};
      double mult1GeVHB_Jets[6]= {0};
      double mult0GeVHB_Jets[6]= {0};

      double mult3GeV3nsHB_Jet0(0), mult3GeV3nsHB_Jet1(0), mult3GeV3nsHB_Jet2(0), mult3GeV3nsHB_Jet3(0);
      double mult0GeV5nsHB_Jet0(0), mult0GeV5nsHB_Jet1(0), mult0GeV5nsHB_Jet2(0), mult0GeV5nsHB_Jet3(0);
      double DepthVariableMult(0);
      double JetEta1(0), JetEta2(0), JetEta3(0), JetEta4(0);
      double JetPhi1(0), JetPhi2(0), JetPhi3(0), JetPhi4(0);
      int genJetRequirementPassed = 0; // keeps track of for this event, were there L1 jets that passed the gen matching requirement of L1 jet to parton?                                                       
      double DelayedHitCounter1GeV[26] = {0};
      double DelayedHitCounter2GeV[26] = {0};
      double DelayedHitCounter3GeV[26] = {0};
      double DelayedHitCounter4GeV[26] = {0};
      double DelayedHitCounter1GeV_jet2[26] = {0};
      double DelayedHitCounter2GeV_jet2[26] = {0};
      double DelayedHitCounter3GeV_jet2[26] = {0};
      double DelayedHitCounter4GeV_jet2[26] = {0};
      double DelayedHitCounter1GeV_jet3[26] = {0};
      double DelayedHitCounter2GeV_jet3[26] = {0};
      double DelayedHitCounter3GeV_jet3[26] = {0};
      double DelayedHitCounter4GeV_jet3[26] = {0};
      double DelayedHitCounter1GeV_jet4[26] = {0};
      double DelayedHitCounter2GeV_jet4[26] = {0};
      double DelayedHitCounter3GeV_jet4[26] = {0};
      double DelayedHitCounter4GeV_jet4[26] = {0};

      double AllHitCounter0GeV(0), AllHitCounter1GeV(0), AllHitCounter2GeV(0), AllHitCounter3GeV(0), AllHitCounter4GeV(0), AllHitCounter0GeV_jet2(0), AllHitCounter1GeV_jet2(0), AllHitCounter2GeV_jet2(0), AllHitCounter3GeV_jet2(0), AllHitCounter4GeV_jet2(0), AllHitCounter0GeV_jet3(0), AllHitCounter1GeV_jet3(0), AllHitCounter2GeV_jet3(0), AllHitCounter3GeV_jet3(0), AllHitCounter4GeV_jet3(0), AllHitCounter0GeV_jet4(0), AllHitCounter1GeV_jet4(0), AllHitCounter2GeV_jet4(0), AllHitCounter3GeV_jet4(0), AllHitCounter4GeV_jet4(0);
      double avgTimeJet1[6] = {0};
      double numHitsJet1[6] = {0};
      double avgTimeJet2[6] = {0};
      double numHitsJet2[6] = {0};
      double avgTimeJet3[6] = {0};
      double numHitsJet3[6] = {0};
      double avgTimeJet4[6] = {0};
      double numHitsJet4[6] = {0};

      // loop over L1 jets, and only do first four (4 highest energy L1 jets from 4 leptons)
      for(uint jetIt=0; jetIt < nJetemu && jetIt < 4; jetIt++){
	hJetEt->Fill(l1emu_->jetEt[jetIt]); // these are already in order of highest E_T
	//	if ((jetIt == 0) && (l1emu_->jetEt[jetIt] < 1000) ) ET_Jet1 = l1emu_->jetEt[jetIt];
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
	    //	    if (generator_->partPt[partonN] < 20) continue; // require pT > 20 GeV -- remove requirement since some jets could come from low energy partons
	    if ( abs(intersect(generator_->partVx[partonN],generator_->partVy[partonN],generator_->partVz[partonN], generator_->partPx[partonN],generator_->partPy[partonN],generator_->partPz[partonN])[0]) < 3 ) { // require parton is in the HE HB region
	      //	      std::cout << generator_->partId[partonN] << std::endl;
	      partonEta = intersect(generator_->partVx[partonN],generator_->partVy[partonN],generator_->partVz[partonN], generator_->partPx[partonN],generator_->partPy[partonN],generator_->partPz[partonN])[0];
	      partonPhi = intersect(generator_->partVx[partonN],generator_->partVy[partonN],generator_->partVz[partonN], generator_->partPx[partonN],generator_->partPy[partonN],generator_->partPz[partonN])[1];
	      // see if LLP decayed in HCAL volume
	      // requiring quark gen vertex to be in the HCAL volume
	      // HE region 3.88 - 5.68 m and out to 2.95 m radius
	      // HB region z below 3.88m, radius 1.79 - 2.95m
	      double radius = sqrt(generator_->partVx[partonN]*generator_->partVx[partonN] + generator_->partVy[partonN]*generator_->partVy[partonN]);
	      if ( abs(generator_->partVz[partonN]) < 568 && (radius < 295) && ( (abs(generator_->partVz[partonN]) > 388) || radius > 179 ) ) {
		partonEtaHCAL = partonEta;
		partonPhiHCAL = partonPhi;
	      }
	    } // end parton Pt cut loop 
	  } // end loop restricting to quarks or gluons 
	  // calculate the delta R between the L1 jet and the generator parton that passed all cuts
	  deltaR_partonGen_L1jet = deltaR(Jet_eta,Jet_phi,partonEta,partonPhi);
	  if (deltaR_partonGen_L1jet < min_deltaR_partonGen_L1jet ) min_deltaR_partonGen_L1jet = deltaR_partonGen_L1jet; // find the min dR between each L1 jet and the gen partons. Use this for gen matching criterion
	  deltaR_partonGenHCAL_L1jet = deltaR(Jet_eta,Jet_phi,partonEtaHCAL,partonPhiHCAL);
          if (deltaR_partonGenHCAL_L1jet < min_deltaR_partonGenHCAL_L1jet ) min_deltaR_partonGenHCAL_L1jet = deltaR_partonGenHCAL_L1jet; // find the min dR between each L1 jet and the gen partons. Use this for gen matching criterion, when LLP decays in HCAL volume
	} // end parton loop
	if (jetIt == 0) DeltaR_partons_L1Jet_1->Fill(min_deltaR_partonGen_L1jet); // minimum deltaR plot of distance between L1 jet 1 and partons passing cuts
	if (jetIt == 1) DeltaR_partons_L1Jet_2->Fill(min_deltaR_partonGen_L1jet);
	if (jetIt == 2) DeltaR_partons_L1Jet_3->Fill(min_deltaR_partonGen_L1jet);
	if (jetIt == 3) DeltaR_partons_L1Jet_4->Fill(min_deltaR_partonGen_L1jet);
	if (jetIt == 0) DeltaR_partonsHCAL_L1Jet_1->Fill(min_deltaR_partonGenHCAL_L1jet); // minimum deltaR plot of distance between L1 jet 1 and partons passing cuts that have a vertex in the HCAL volume
	if (jetIt == 1) DeltaR_partonsHCAL_L1Jet_2->Fill(min_deltaR_partonGenHCAL_L1jet);
	if (jetIt == 2) DeltaR_partonsHCAL_L1Jet_3->Fill(min_deltaR_partonGenHCAL_L1jet);
	if (jetIt == 3) DeltaR_partonsHCAL_L1Jet_4->Fill(min_deltaR_partonGenHCAL_L1jet);

	// DeltaR cut between L1 jet and generator partons -- this is where GEN MATCHING is enforced
	//	if (min_deltaR_partonGen_L1jet > 0.3) continue; 
	genJetRequirementPassed += 1; // keeps track of for this event, were there L1 jets that passed the gen matching requirement of L1 jet to parton?        

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

	  /*
	  // ****************************** GEN MATCHED SECTION FOR MULTIPLICITY, TIMING AND ENERGY PROFILES  *********************************** 
	  // other option preventing double counting later in HCAL TP loop                                                                             
	  // ************************************************************************************************************************************ 
	  // this is saved in a file in HCAL research folder in case needed

	  // *************************** END OF GEN MATCHED SECTION *********************************
	  // other option preventing double counting later in HCAL TP loop
	  // ****************************************************************************************
	  */
	} // closing HCAL TP loop

        // minimum DeltaR between L1 jet and HCAL TP on a per event basis    
        if (jetIt == 0) DeltaR_TP_L1Jet_1->Fill(min_DeltaR);
        if (jetIt == 1) DeltaR_TP_L1Jet_2->Fill(min_DeltaR); 
        if (jetIt == 2) DeltaR_TP_L1Jet_3->Fill(min_DeltaR); 
        if (jetIt == 3) DeltaR_TP_L1Jet_4->Fill(min_DeltaR);
	min_DeltaR += 0;
      } // closing L1 Jets loop
      // fill DeltaR histograms for DR between each L1 Jet to understand spatial distribution
      DeltaR_L1Jets_1_2->Fill(deltaR(JetEta1,JetPhi1,JetEta2,JetPhi2));
      DeltaR_L1Jets_1_3->Fill(deltaR(JetEta1,JetPhi1,JetEta3,JetPhi3));
      DeltaR_L1Jets_1_4->Fill(deltaR(JetEta1,JetPhi1,JetEta4,JetPhi4));
      DeltaR_L1Jets_2_3->Fill(deltaR(JetEta2,JetPhi2,JetEta3,JetPhi3));
      DeltaR_L1Jets_2_4->Fill(deltaR(JetEta2,JetPhi2,JetEta4,JetPhi4));
      DeltaR_L1Jets_3_4->Fill(deltaR(JetEta3,JetPhi3,JetEta4,JetPhi4));

      double min_DR_L1_TP = 100;
      int LLPdecayDetAcc = 0; // count how many LLPs decay before end of detector acceptance 

      // HCAL TP loop
      // Used for plots with all HCAL TPs (not matched to L1 jets), and when each HCAL TP is associated to one L1 Jet + apply DR restrictions (in L1 jet loop inside of HCAL TP loop)
      for (int HcalTPIt = 0; HcalTPIt < nCaloTPemu; HcalTPIt++){
	tpEtaemu = l1CaloTPemu_->hcalTPieta[HcalTPIt]; // use for HB HE restrictions                                 
	tpPhiemu = l1CaloTPemu_->hcalTPCaliphi[HcalTPIt];  // hcalTPiphi[HcalTPIt]; // use for deltaR calculation
	tpEtemu = l1CaloTPemu_->hcalTPet[HcalTPIt]; // used for energy normalization in the energy depth plots    
	nDepth = l1CaloTPemu_->hcalTPnDepths[HcalTPIt]; // how many HCAL depth layers to go over
	
       	if (nDepth == 0) continue; // skipping events where depth = 0, since here timing = -1 and energy = 0 (invalid event)

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
	/*
        double exp_time = 0;
        if (abs(tpEtaemu) <= 4) exp_time = 5.97;
        if (abs(tpEtaemu) <= 8 && abs(tpEtaemu) > 4) exp_time = 7.4;
        if (abs(tpEtaemu) <= 12 && abs(tpEtaemu) > 8) exp_time = 9.4;
        if (abs(tpEtaemu) <= 16 && abs(tpEtaemu) > 12) exp_time = 11.8;
	*/

	// Energy deposited in each depth layer for every HCAL TP (4 in HB, 7 in HE)  
	hcalTPdepth[0] = l1CaloTPemu_->hcalTPDepth1[HcalTPIt];
	hcalTPdepth[1] = l1CaloTPemu_->hcalTPDepth2[HcalTPIt];
	hcalTPdepth[2] = l1CaloTPemu_->hcalTPDepth3[HcalTPIt];
	hcalTPdepth[3] = l1CaloTPemu_->hcalTPDepth4[HcalTPIt];
	hcalTPdepth[4] = l1CaloTPemu_->hcalTPDepth5[HcalTPIt];
	hcalTPdepth[5] = l1CaloTPemu_->hcalTPDepth6[HcalTPIt];
	hcalTPdepth[6] = l1CaloTPemu_->hcalTPDepth7[HcalTPIt];
	// timing info for each layer, in 25 ns with resolution 0.5 ns, with the expected time subtracted off
	hcalTPtiming[0] = l1CaloTPemu_->hcalTPtiming1[HcalTPIt];// - exp_time;
	hcalTPtiming[1] = l1CaloTPemu_->hcalTPtiming2[HcalTPIt];// - exp_time;
	hcalTPtiming[2] = l1CaloTPemu_->hcalTPtiming3[HcalTPIt];// - exp_time;
	hcalTPtiming[3] = l1CaloTPemu_->hcalTPtiming4[HcalTPIt];// - exp_time;
	hcalTPtiming[4] = l1CaloTPemu_->hcalTPtiming5[HcalTPIt];// - exp_time;
	hcalTPtiming[5] = l1CaloTPemu_->hcalTPtiming6[HcalTPIt];// - exp_time;
	hcalTPtiming[6] = l1CaloTPemu_->hcalTPtiming7[HcalTPIt];// - exp_time;
	for (int i = 0; i < 4; i++ ) {
	  if ( (abs(tpEtaemu) == 1) && (hcalTPtiming[i] > -0.5) ) {
	    centralTiming->Fill(hcalTPtiming[i]-5.5); // filling distribution of time of hits (corrected by expected TOF) in the center of the barrel (ieta = 1)
	  }
	}

	// ratio of energy in first two HCAL layers to total energy in HCAL, only for high energy TPs          
        if ( tpEtemu > 10 ) {
          Ratio_Depth->Fill( (hcalTPdepth[2]+hcalTPdepth[3]) / tpEtemu);
          if (abs(tpEtaemu) < 16) {
            Ratio_DepthHB->Fill( (hcalTPdepth[2]+hcalTPdepth[3]) / tpEtemu);
          }
          if (abs(tpEtaemu) >= 16 && abs(tpEtaemu) < 29 ) {
            Ratio_DepthHE->Fill( (hcalTPdepth[2]+hcalTPdepth[3]) / tpEtemu);
          }
        }

	// loop over HCAL depths for every HCAL TP
        for (int depthIt = 0; depthIt < nDepth-1; depthIt++){
	  // filling energy and time plots for each of 7 HCAL depths  
	  Energy_Depth->Fill(depthIt+1,hcalTPdepth[depthIt]/tpEtemu); // normalized by total energy in event so is fractional energy in each layer      
	  if (hcalTPtiming[depthIt] > 0) Timing_Depth->Fill(depthIt+1,hcalTPtiming[depthIt]); // raw timing value in each layer  
          if (tpEtemu > 10 ) {
            Energy_Depth_HighE->Fill(depthIt+1,hcalTPdepth[depthIt]/tpEtemu); // energy depth for high energy HCAL TPs
          }
	  if (abs(tpEtaemu) < 16) {
	    Energy_DepthHB->Fill(depthIt+1,hcalTPdepth[depthIt]/tpEtemu);
	    if (hcalTPtiming[depthIt] > 0) Timing_DepthHB->Fill(depthIt+1,hcalTPtiming[depthIt]);
	    if (tpEtemu > 10 ) {
	      Energy_DepthHB_HighE->Fill(depthIt+1,hcalTPdepth[depthIt]/tpEtemu);
	    }
	  }
	  if (abs(tpEtaemu) >= 16 && abs(tpEtaemu) < 29 ) {
	    Energy_DepthHE->Fill(depthIt+1,hcalTPdepth[depthIt]/tpEtemu);
	    if (hcalTPtiming[depthIt] > 0) Timing_DepthHE->Fill(depthIt+1,hcalTPtiming[depthIt]);
	    if (tpEtemu > 10 ) {
	      Energy_DepthHE_HighE->Fill(depthIt+1,hcalTPdepth[depthIt]/tpEtemu);
	    }
	  }

          // count multiplicity of layers given a timing and energy threshold   
          // 3 GeV energy cut
          if (hcalTPdepth[depthIt] > 3){
	    for (int TDCns = 1; TDCns <= 5; TDCns++) {
	      if (hcalTPtiming[depthIt] > TDCns ) mult3GeV[TDCns] += 1;
              // 3 GeV HB regions
	      if ( (abs(tpEtaemu) < 16) ) {
		if (hcalTPtiming[depthIt] > TDCns ) mult3GeVHB[TDCns] += 1;
	      }
	      // 3 GeV HE regions
	      if ( (abs(tpEtaemu) >= 16) && (abs(tpEtaemu) < 29) ){
		if (hcalTPtiming[depthIt] > TDCns ) mult3GeVHE[TDCns] += 1;
	      }
	    }
	    if (abs(tpEtaemu) < 5 ) mult3GeVcaloT1 += 1;
	    if ((abs(tpEtaemu) > 4) && (abs(tpEtaemu) < 9) ) mult3GeVcaloT2 += 1;
	    if ((abs(tpEtaemu) > 8) && (abs(tpEtaemu) < 13) ) mult3GeVcaloT3 += 1;
	    if ((abs(tpEtaemu) > 12) && (abs(tpEtaemu) < 17) ) mult3GeVcaloT4 += 1;
	  } // closing 3 GeV energy cut loop
	  // 2 GeV energy cut
	  if (hcalTPdepth[depthIt] > 2){
            for (int TDCns = 1; TDCns <= 5; TDCns++) {
	      if (hcalTPtiming[depthIt] > TDCns ) mult2GeV[TDCns] += 1;
	      // 2 GeV HB regions                                
	      if ( (abs(tpEtaemu) < 16) ){
		if (hcalTPtiming[depthIt] > TDCns ) mult2GeVHB[TDCns] += 1;
	      }
	      // 2 GeV HE regions
	      if ( (abs(tpEtaemu) >= 16) && (abs(tpEtaemu) < 29) ){
		if (hcalTPtiming[depthIt] > TDCns ) mult2GeVHE[TDCns] += 1;
	      }
	    }
	    if (abs(tpEtaemu) < 5 ) mult2GeVcaloT1 += 1;
	    if ((abs(tpEtaemu) > 4) && (abs(tpEtaemu) < 9) ) mult2GeVcaloT2 += 1;
	    if ((abs(tpEtaemu) > 8) && (abs(tpEtaemu) < 13) ) mult2GeVcaloT3 += 1;
	    if ((abs(tpEtaemu) > 12) && (abs(tpEtaemu) < 17) ) mult2GeVcaloT4 += 1;
	  } // closing 2 GeV energy cut loop
	  // 1 GeV energy cut
	  if (hcalTPdepth[depthIt] > 1){
	    for (int TDCns = 1; TDCns <= 5; TDCns++) {
	      if (hcalTPtiming[depthIt] > TDCns ) mult1GeV[TDCns] += 1;
	      // 1 GeV HB regions                                                      
	      if ( (abs(tpEtaemu) < 16) ){
                if (hcalTPtiming[depthIt] > TDCns ) mult1GeVHB[TDCns] += 1;
	      }
	      // 1 GeV HE regions
	      if ( (abs(tpEtaemu) >= 16) && (abs(tpEtaemu) < 29) ){
		if (hcalTPtiming[depthIt] > TDCns ) mult1GeVHE[TDCns] += 1;
	      }
	    }
	    if (abs(tpEtaemu) < 5 ) mult1GeVcaloT1 += 1;
            if ((abs(tpEtaemu) > 4) && (abs(tpEtaemu) < 9) ) mult1GeVcaloT2 += 1;
            if ((abs(tpEtaemu) > 8) && (abs(tpEtaemu) < 13) ) mult1GeVcaloT3 += 1;
	    if ((abs(tpEtaemu) > 12) && (abs(tpEtaemu) < 17) ) mult1GeVcaloT4 += 1;
	  } // closing 1 GeV energy cut loop
	} // closing HCAL depths loop
	// here have calculated multiplicity for a single HCAL TP in an event

        // *********************** GEN MATCHED PLOTS FOR MULTIPLICITY, TIMING, AND ENERGY DEPTH ************************ 
        // these plots are made from gen matched L1 jets (have a parton within DR<0.3), and only associating each HCAL TP to a single L1 jet
        // ************************************************************************************************************* 

	// in HCAL TP loop, find which L1 jet is closest to the HCAL TP
	uint nJetemu(0);
	nJetemu = l1emu_->nJets;
	double min_DeltaR = 100;
	double DeltaR = 100;
	int closestJet(-1);
	double closest_Jet_eta = 1000;
	double closest_Jet_phi = 1000;

	// loop over L1 jets to get the distance from them to the HCAL TPs, find closet L1 jet to the TP
	for (uint jetIt=0; (jetIt < nJetemu) && (jetIt < 4); jetIt++){ // loop over L1 jets, and only do first four (4 highest energy L1 jets from 4 leptons)
	  if (l1emu_->jetEt[jetIt] < 20 ) continue; // require jet is greater than 20 GeV to attempt matching to HCAL TP
	  double Jet_eta;
	  double Jet_phi;
	  Jet_eta = l1emu_->jetEta[jetIt];
	  Jet_phi = l1emu_->jetPhi[jetIt];

	  DeltaR = deltaR(Jet_eta,Jet_phi,TP_eta,TP_phi);   // this is DeltaR for the HCAL TPs to L1 Jet 
	  if (DeltaR < min_DeltaR) {
	    min_DeltaR = DeltaR; // find min delta R between L1 Jet and HCAL TPs -- this is reset for each HCAL TP, so is which jet is closest to the TP
	    closestJet = jetIt; // record which L1 jet is the closest
	    closest_Jet_eta = Jet_eta; // save eta and phi of the closest L1 jet for use in the gen matching to a parton
	    closest_Jet_phi = Jet_phi;
	  }
	  // find overall smallest DR between leading L1 jet and HCAL TP on a per event basis
	  if ((jetIt == 0) && (DeltaR < min_DR_L1_TP)) min_DR_L1_TP = DeltaR;
	} // closing the L1 jet loop
      
	// ************* GEN PARTICLE MATCHING CODE ******************                                                         
	// already got the entries from genTree at start of event loop                  
	double min_deltaR_partonGen_L1jet = 1000;
	double min_deltaR_partonGenHCAL_L1jet = 1000; // min dR defined per each L1 jet
	double min_deltaR_partonGen1m_L1jet = 1000;
       	double parton1vertex = 0;
      	double parton1number = 0;
        LLPdecayDetAcc = 0; // count how many LLPs decay before end of detector acceptance
	for (int partonN = 0; partonN < generator_->nPart; partonN ++) {
	  // reset values for eta, phi of partons (w/ and w/o requirement of vertex within HCAL volume), and for dR between parton and L1 jet on a per-parton basis
	  double partonEta = 1000;
	  double partonPhi = 1000;
	  double partonEtaHCAL = 1000;
	  double partonPhiHCAL = 1000;
          double partonEta1m = 1000;
          double partonPhi1m = 1000;
	  double deltaR_partonGen_L1jet = 1000;
	  double deltaR_partonGenHCAL_L1jet = 1000;
          double deltaR_partonGen1m_L1jet = 1000;

	  if (generator_->partHardProcess[partonN] == 0 ) continue; // require hard process - this is what LLP results from 
	  if ( (abs(generator_->partId[partonN]) >= 1 && abs(generator_->partId[partonN]) <=5 ) || (abs(generator_->partId[partonN]) == 21) ) { // only consider generator particles that are partons (quarks, gluons) from the LLP decay. Top quark is unstable 
	    // detector cuts for the HB and HE regions
	    //    if (generator_->partPt[partonN] < 20) continue; // require pT > 20 GeV -- remove requirement since some jets could come from low energy partons
	    if ( abs(intersect(generator_->partVx[partonN],generator_->partVy[partonN],generator_->partVz[partonN], generator_->partPx[partonN],generator_->partPy[partonN],generator_->partPz[partonN])[0]) <= 3 ) { // require parton is in the HE HB region
	      //      std::cout << generator_->partId[partonN] << std::endl;
	      partonEta = intersect(generator_->partVx[partonN],generator_->partVy[partonN],generator_->partVz[partonN], generator_->partPx[partonN],generator_->partPy[partonN],generator_->partPz[partonN])[0];
	      partonPhi = intersect(generator_->partVx[partonN],generator_->partVy[partonN],generator_->partVz[partonN], generator_->partPx[partonN],generator_->partPy[partonN],generator_->partPz[partonN])[1];
	      // see if LLP decayed anywhere in detector acceptance (HCAL, ECAL, tracker) since that is what is relevant for studies
	      double radius = sqrt(generator_->partVx[partonN]*generator_->partVx[partonN] + generator_->partVy[partonN]*generator_->partVy[partonN]);
	      if ( radius < 295 && abs(generator_->partVz[partonN]) < 568 ) LLPdecayDetAcc += 1; // if a LLP decays before end of HCAL, increment this variable                                  
	      // specifically choose partons that decay with a certain distance to consider time values from this as a check
	      if (sqrt(generator_->partVx[partonN]*generator_->partVx[partonN] + generator_->partVy[partonN]*generator_->partVy[partonN] + generator_->partVz[partonN]*generator_->partVz[partonN]) > 100 ) {
		partonEta1m = intersect(generator_->partVx[partonN],generator_->partVy[partonN],generator_->partVz[partonN], generator_->partPx[partonN],generator_->partPy[partonN],generator_->partPz[partonN])[0];
		partonPhi1m = intersect(generator_->partVx[partonN],generator_->partVy[partonN],generator_->partVz[partonN], generator_->partPx[partonN],generator_->partPy[partonN],generator_->partPz[partonN])[1];
	      }
	      if (HcalTPIt == 1) {
		double vertex = sqrt(generator_->partVx[partonN]*generator_->partVx[partonN] + generator_->partVy[partonN]*generator_->partVy[partonN] + generator_->partVz[partonN]*generator_->partVz[partonN]);
		//		  std::cout << vertex << " ----- " << parton1vertex << std::endl;
		if ( parton1vertex == vertex && parton1number+1 == partonN ) { // means found second b quark with same vertex 
		  //		  std::cout << generator_->partParent[partonN] << " and entry " << jentry << std::endl;
		  double genLLPBeta = sqrt((generator_->partPx[partonN-1] + generator_->partPx[partonN])*(generator_->partPx[partonN-1] + generator_->partPx[partonN]) + (generator_->partPy[partonN-1] + generator_->partPy[partonN])*(generator_->partPy[partonN-1] + generator_->partPy[partonN]) + (generator_->partPz[partonN-1] + generator_->partPz[partonN])*(generator_->partPz[partonN-1] + generator_->partPz[partonN])) / (generator_->partE[partonN-1] + generator_->partE[partonN]);
		  double genLLPGamma = 1./TMath::Sqrt(1.-genLLPBeta*genLLPBeta);
		  LLPdecayRadiusDetAcceptance->Fill(radius,abs(generator_->partVz[partonN]),1); // fill the radius for events that have LLP decaying anywhere in HCAL or earlier
		  LLPdecayXyzDetAcceptance->Fill(vertex / (genLLPBeta * genLLPGamma));
		  //		    std::cout << radius << " radius, and entry = " << jentry << " and parton number = " << partonN << std::endl;
		}
		else {
		  parton1vertex = vertex;
		  parton1number = partonN;
		}
	      }
	      // see if LLP decayed in HCAL volume
	      // requiring quark gen vertex to be in the HCAL volume
	      // HE region 3.88 - 5.68 m and out to 2.95 m radius
	      // HB region z below 3.88m, radius 1.79 - 2.95m
	      if ( abs(generator_->partVz[partonN]) < 568 && (radius < 295) && ( (abs(generator_->partVz[partonN]) > 388) || radius > 179 ) ) {
		partonEtaHCAL = intersect(generator_->partVx[partonN],generator_->partVy[partonN],generator_->partVz[partonN], generator_->partPx[partonN],generator_->partPy[partonN],generator_->partPz[partonN])[0];
		partonPhiHCAL = intersect(generator_->partVx[partonN],generator_->partVy[partonN],generator_->partVz[partonN], generator_->partPx[partonN],generator_->partPy[partonN],generator_->partPz[partonN])[1];
	      }
	    } // end parton Pt cut loop 
	  } // end loop restricting to quarks or gluons 
	  // calculate the delta R between the L1 jet and the generator parton that passed all cuts
	  deltaR_partonGen_L1jet = deltaR(closest_Jet_eta,closest_Jet_phi,partonEta,partonPhi);
	  if (deltaR_partonGen_L1jet < min_deltaR_partonGen_L1jet ) min_deltaR_partonGen_L1jet = deltaR_partonGen_L1jet; // find the min dR between each L1 jet and the gen partons. Use this for gen matching criterion
	  deltaR_partonGenHCAL_L1jet = deltaR(closest_Jet_eta,closest_Jet_phi,partonEtaHCAL,partonPhiHCAL);
          if (deltaR_partonGenHCAL_L1jet < min_deltaR_partonGenHCAL_L1jet ) min_deltaR_partonGenHCAL_L1jet = deltaR_partonGenHCAL_L1jet; // find the min dR between each L1 jet and the gen partons. Use this for gen matching criterion, when LLP decays in HCAL volume
	  deltaR_partonGen1m_L1jet = deltaR(closest_Jet_eta,closest_Jet_phi,partonEta1m,partonPhi1m);
	  if (deltaR_partonGen1m_L1jet < min_deltaR_partonGen1m_L1jet ) min_deltaR_partonGen1m_L1jet = deltaR_partonGen1m_L1jet;
	} // end parton loop
	// ********************************************* GEN MATCHING REQUIREMENT *****************************************
	//        if (min_deltaR_partonGen_L1jet > 0.3) continue;  // require that the L1 jet has a generator parton nearby
	//      	if (inputFile.substr(0,5) == "../mh" &&  min_deltaR_partonGen1m_L1jet > 0.3 ) continue; // only look at L1 jets resulting from a LLP that lasted 1m or more

	genJetRequirementPassed += 1; // how many L1 jets this event passed the generator matching requirement? without gen matching this should be same as total number of events
        // loop over L1 jets to find DR between L1 jet and each TP, given that they are gen matched
        for (uint jetIt=0; (jetIt < nJetemu) && (jetIt < 4); jetIt++){ // loop over L1 jets, and only do first four (4 highest energy L1 jets from 4 leptons)    
          if (l1emu_->jetEt[jetIt] < 20 ) continue; // require jet is greater than 20 GeV to attempt matching to HCAL TP                                         
          double Jet_eta;
          double Jet_phi;
          Jet_eta = l1emu_->jetEta[jetIt];
          Jet_phi = l1emu_->jetPhi[jetIt];
	  DeltaR = deltaR(Jet_eta,Jet_phi,TP_eta,TP_phi);
          // fill DR plots between each L1 jet and every HCAL TP    
          if (jetIt == 0) DeltaR_TPall_L1Jet_1->Fill(DeltaR);
          if (jetIt == 1) DeltaR_TPall_L1Jet_2->Fill(DeltaR);
          if (jetIt == 2) DeltaR_TPall_L1Jet_3->Fill(DeltaR);
          if (jetIt == 3) DeltaR_TPall_L1Jet_4->Fill(DeltaR);
          // fill DR plots between each L1 jet and HCAL TPs with at least one layer above energy and timing cuts    
          if ( (hcalTPtiming[0]>3 && hcalTPdepth[0]>3) || (hcalTPtiming[1]>3 && hcalTPdepth[1]>3) || (hcalTPtiming[2]>3 && hcalTPdepth[2]>3) || (hcalTPtiming[3]>3 && hcalTPdepth[3]>3) || (hcalTPtiming[4]>3 && hcalTPdepth[4]>3) || (hcalTPtiming[5]>3 && hcalTPdepth[5]>3) || (hcalTPtiming[6]>3 && hcalTPdepth[6]>3) ) {
            if (jetIt == 0) DeltaR_TPet_L1Jet_1->Fill(DeltaR);
            if (jetIt == 1) DeltaR_TPet_L1Jet_2->Fill(DeltaR);
            if (jetIt == 2) DeltaR_TPet_L1Jet_3->Fill(DeltaR);
            if (jetIt == 3) DeltaR_TPet_L1Jet_4->Fill(DeltaR);
          }
	} // closing L1 jet loop plotting distance from each HCAL TP to the L1 jet

        if ( min_DeltaR > DR_threshold ) continue; // don't fill matched multiplicity for TPs that are DR > 0.5 away from their nearest L1 jet
	// Ratio of energy in first two HCAL layers to all HCAL layers. Only consider for high energy TPs > 10 GeV
	if ( tpEtemu > 10 ) {
	  Ratio_Depth_Jets->Fill( (hcalTPdepth[2]+hcalTPdepth[3]) / tpEtemu);
	  if (abs(tpEtaemu) < 16) {
	    Ratio_DepthHB_Jets->Fill( (hcalTPdepth[2]+hcalTPdepth[3]) / tpEtemu);
	  }
	  if ((abs(tpEtaemu) >= 16) && (abs(tpEtaemu) < 29) ) {
	    Ratio_DepthHE_Jets->Fill( (hcalTPdepth[2]+hcalTPdepth[3]) / tpEtemu);
	  }
	}

	// loop over HCAL depths for the HCAL TP
	// filling fractional energy deposit and avg time plots for each of 7 HCAL depths. Done for all HCAL TPs within DR 0.5 of the L1 Jet
	for (int depthIt = 0; depthIt < nDepth-1; depthIt++){
	  if (inputFile.substr(0,5) == "../mh" ) { // HCAL volume decay requirement for energy depth plots for LLPs
	    if (min_deltaR_partonGenHCAL_L1jet < 1) Energy_Depth_Jets->Fill(depthIt+1,hcalTPdepth[depthIt]/tpEtemu); // normalized by total energy in event so is fractional energy in each layer
	    if (hcalTPtiming[depthIt] > 0) Timing_Depth_Jets->Fill(depthIt+1,hcalTPtiming[depthIt]); // raw timing value in each layer
	    if (min_deltaR_partonGenHCAL_L1jet < 1) {
	      if (tpEtemu > 10 ) Energy_Depth_Jets_HighE->Fill(depthIt+1,hcalTPdepth[depthIt]/tpEtemu); // restricting to high energy HCAL TPs
	      if (hcalTPtiming[depthIt] > 0) VolTiming_Depth_Jets->Fill(depthIt+1,hcalTPtiming[depthIt]);
	    }
	    if (abs(tpEtaemu) < 16) {
	      if (hcalTPtiming[depthIt] > 1) Energy_DepthHB_Jets_1ns->Fill(depthIt+1,hcalTPdepth[depthIt]); // not normalized by total energy of TP tower
              if (hcalTPtiming[depthIt] > 2) Energy_DepthHB_Jets_2ns->Fill(depthIt+1,hcalTPdepth[depthIt]);
	      if (hcalTPtiming[depthIt] > 0) Timing_DepthHB_Jets->Fill(depthIt+1,hcalTPtiming[depthIt]);
	      if (min_deltaR_partonGenHCAL_L1jet < 1) {
		Energy_DepthHB_Jets->Fill(depthIt+1,hcalTPdepth[depthIt]/tpEtemu);
		if (tpEtemu > 10 ) Energy_DepthHB_Jets_HighE->Fill(depthIt+1,hcalTPdepth[depthIt]/tpEtemu);
		if (hcalTPtiming[depthIt] > 0) VolTiming_DepthHB_Jets->Fill(depthIt+1,hcalTPtiming[depthIt]);
	      }
	    }
	    if (abs(tpEtaemu) >= 16 && abs(tpEtaemu) < 29 ) {
	      if (hcalTPtiming[depthIt] > 1) Energy_DepthHE_Jets_1ns->Fill(depthIt+1,hcalTPdepth[depthIt]);
	      if (hcalTPtiming[depthIt] > 2) Energy_DepthHE_Jets_2ns->Fill(depthIt+1,hcalTPdepth[depthIt]);
	      if (hcalTPtiming[depthIt] > 0) Timing_DepthHE_Jets->Fill(depthIt+1,hcalTPtiming[depthIt]);
	      if (min_deltaR_partonGenHCAL_L1jet < 1) {
		Energy_DepthHE_Jets->Fill(depthIt+1,hcalTPdepth[depthIt]/tpEtemu);
		if (tpEtemu > 10 ) Energy_DepthHE_Jets_HighE->Fill(depthIt+1,hcalTPdepth[depthIt]/tpEtemu);
		if (hcalTPtiming[depthIt] > 0) VolTiming_DepthHE_Jets->Fill(depthIt+1,hcalTPtiming[depthIt]);
	      }
	    }
	  }
	  else { // no HCAL volume decay requirements for QCD or neutrino gun
            Energy_Depth_Jets->Fill(depthIt+1,hcalTPdepth[depthIt]/tpEtemu); // normalized by total energy in event so is fractional energy in each layer
            if (hcalTPtiming[depthIt] > 0) Timing_Depth_Jets->Fill(depthIt+1,hcalTPtiming[depthIt]); // raw timing value in each layer  
	    if (hcalTPtiming[depthIt] > 0) VolTiming_Depth_Jets->Fill(depthIt+1,hcalTPtiming[depthIt]);
            if (tpEtemu > 10 ) {
              Energy_Depth_Jets_HighE->Fill(depthIt+1,hcalTPdepth[depthIt]/tpEtemu); // restricting to high energy HCAL TPs                                               
            }
            if (abs(tpEtaemu) < 16) {
	      if (hcalTPtiming[depthIt] > 1) Energy_DepthHB_Jets_1ns->Fill(depthIt+1,hcalTPdepth[depthIt]);
              if (hcalTPtiming[depthIt] > 2) Energy_DepthHB_Jets_2ns->Fill(depthIt+1,hcalTPdepth[depthIt]);
              Energy_DepthHB_Jets->Fill(depthIt+1,hcalTPdepth[depthIt]/tpEtemu);
              if (hcalTPtiming[depthIt] > 0) Timing_DepthHB_Jets->Fill(depthIt+1,hcalTPtiming[depthIt]);
	      if (hcalTPtiming[depthIt] > 0) VolTiming_DepthHB_Jets->Fill(depthIt+1,hcalTPtiming[depthIt]);
              if (tpEtemu > 10 ) {
                Energy_DepthHB_Jets_HighE->Fill(depthIt+1,hcalTPdepth[depthIt]/tpEtemu);
              }
            }
            if (abs(tpEtaemu) >= 16 && abs(tpEtaemu) < 29 ) {
	      if (hcalTPtiming[depthIt] > 1) Energy_DepthHE_Jets_1ns->Fill(depthIt+1,hcalTPdepth[depthIt]);
              if (hcalTPtiming[depthIt] > 2) Energy_DepthHE_Jets_2ns->Fill(depthIt+1,hcalTPdepth[depthIt]);
              Energy_DepthHE_Jets->Fill(depthIt+1,hcalTPdepth[depthIt]/tpEtemu);
              if (hcalTPtiming[depthIt] > 0) Timing_DepthHE_Jets->Fill(depthIt+1,hcalTPtiming[depthIt]);
	      if (hcalTPtiming[depthIt] > 0) VolTiming_DepthHE_Jets->Fill(depthIt+1,hcalTPtiming[depthIt]);
              if (tpEtemu > 10 ) {
                Energy_DepthHE_Jets_HighE->Fill(depthIt+1,hcalTPdepth[depthIt]/tpEtemu);
              }
            }
	  }
     
	  if ( ( hcalTPtiming[depthIt] >= 0 ) && ( abs(tpEtaemu) < 29 ) ) {
	    // delayed hit fraction calculation, over entire HCAL, near most energetic L1 jet
	    if (closestJet == 0 ) {
	      if (hcalTPdepth[depthIt] > 0 ) AllHitCounter0GeV +=1;
	      NumberTPtiming_energy->Fill(hcalTPtiming[depthIt], hcalTPdepth[depthIt],1);
	      for (int GeV=0; GeV<6; GeV++) { // scan over energy 0-5 GeV for the cells near L1 jet
		if (hcalTPdepth[depthIt] > GeV && depthIt > 0 ) { // is cell over set energy value?  exclude first depth layer 
		  avgTimeJet1[GeV] += hcalTPtiming[depthIt];
		  numHitsJet1[GeV] += 1;

		  NumberTPtiming[GeV]->Fill(hcalTPtiming[depthIt],1); // add one to the # of TP cells at the timing value
		  if (depthIt == 0 ) NumberTPtiming_depth1[GeV]->Fill(hcalTPtiming[depthIt],1); // depth 1
		  if (depthIt == 1 ) NumberTPtiming_depth2[GeV]->Fill(hcalTPtiming[depthIt],1); 
		  if (depthIt == 2 ) NumberTPtiming_depth3[GeV]->Fill(hcalTPtiming[depthIt],1);
		  if (depthIt == 3 ) NumberTPtiming_depth4[GeV]->Fill(hcalTPtiming[depthIt],1);
		}
 	      }
	      //	      if (depthIt > 0) {
	      if (hcalTPdepth[depthIt] > 1) {
		for (int i=0; i<26; i++)  if (hcalTPtiming[depthIt] > i ) DelayedHitCounter1GeV[i] += 1; // iterator counts the TDC value
		AllHitCounter1GeV += 1;
	      }
	      if (hcalTPdepth[depthIt] > 2) {
		for (int i=0; i<26; i++) if (hcalTPtiming[depthIt] > i ) DelayedHitCounter2GeV[i] += 1; // iterator counts the TDC value
		//		std::cout << "jet 1 timing, energy = " << hcalTPtiming[depthIt] << ", " << hcalTPdepth[depthIt]<< std::endl;
		AllHitCounter2GeV += 1;
	      }
	      if (hcalTPdepth[depthIt] > 3) {
		for (int i=0; i<26; i++) if (hcalTPtiming[depthIt] > i ) DelayedHitCounter3GeV[i] += 1;
		AllHitCounter3GeV += 1;
	      }
	      if (hcalTPdepth[depthIt] > 4) {
		for (int i=0; i<26; i++)  if (hcalTPtiming[depthIt] > i ) DelayedHitCounter4GeV[i] += 1;
		AllHitCounter4GeV += 1;
	      }
	    }
	    
	    if (closestJet == 1) {
	      if (hcalTPdepth[depthIt] > 0 ) AllHitCounter0GeV_jet2 +=1;
	      for (int GeV=0; GeV<6; GeV++) { // scan over energy 0-5 GeV for the cells near L1 jet                                                                      
                if (hcalTPdepth[depthIt] > GeV && depthIt > 0 ) { // is cell over set energy value?  first depth layer is excluded      
                    avgTimeJet2[GeV] += hcalTPtiming[depthIt];
                    numHitsJet2[GeV] += 1;
                  }
		}
	      //	      if (depthIt > 0 ) {
	      if (hcalTPdepth[depthIt] > 1) {
		for (int i=0; i<26; i++) if (hcalTPtiming[depthIt] > i ) DelayedHitCounter1GeV_jet2[i] += 1; // iterator counts the TDC value
		AllHitCounter1GeV_jet2 += 1;
	      }
	      if (hcalTPdepth[depthIt] > 2) {
		for (int i=0; i<26; i++) if (hcalTPtiming[depthIt] > i ) DelayedHitCounter2GeV_jet2[i] += 1; // iterator counts the TDC value 
		//		std::cout << "jet 2 timing, energy = " << hcalTPtiming[depthIt] << ", " << hcalTPdepth[depthIt] << std::endl;
		AllHitCounter2GeV_jet2 += 1;
	      }
	      if (hcalTPdepth[depthIt] > 3) {
		for (int i=0; i<26; i++) if (hcalTPtiming[depthIt] > i ) DelayedHitCounter3GeV_jet2[i] += 1;
		AllHitCounter3GeV_jet2 += 1;
	      }
	      if (hcalTPdepth[depthIt] > 4) {
		for (int i=0; i<26; i++) if (hcalTPtiming[depthIt] > i ) DelayedHitCounter4GeV_jet2[i] += 1;
		AllHitCounter4GeV_jet2 += 1;
	      }
	    }
	    
	    if (closestJet == 2) {
	      if (hcalTPdepth[depthIt] > 0 ) AllHitCounter0GeV_jet3 +=1;
              for (int GeV=0; GeV<6; GeV++) { // scan over energy 0-5 GeV for the cells near L1 jet 
                if (hcalTPdepth[depthIt] > GeV && depthIt > 0 ) { // is cell over set energy value?  first depth layer is excluded 
		  avgTimeJet3[GeV] += hcalTPtiming[depthIt];
		  numHitsJet3[GeV] += 1;
		}
	      }
	      //	      if (depthIt > 0) {
	      if (hcalTPdepth[depthIt] > 1) {
		for (int i=0; i<26; i++) if (hcalTPtiming[depthIt] > i ) DelayedHitCounter1GeV_jet3[i] += 1; // iterator counts the TDC value
		AllHitCounter1GeV_jet3 += 1;
	      }
	      if (hcalTPdepth[depthIt] > 2) {
		for (int i=0; i<26; i++) if (hcalTPtiming[depthIt] > i ) DelayedHitCounter2GeV_jet3[i] += 1; // iterator counts the TDC value
		AllHitCounter2GeV_jet3 += 1;
	      }
	      if (hcalTPdepth[depthIt] > 3) {
		for (int i=0; i<26; i++) if (hcalTPtiming[depthIt] > i ) DelayedHitCounter3GeV_jet3[i] += 1;
		AllHitCounter3GeV_jet3 += 1;
	      }
	      if (hcalTPdepth[depthIt] > 4) {
		for (int i=0; i<26; i++) if (hcalTPtiming[depthIt] > i ) DelayedHitCounter4GeV_jet3[i] += 1;
		AllHitCounter4GeV_jet3 += 1;
	      }
	    }
	    
	    if (closestJet == 3) {
	      if (hcalTPdepth[depthIt] > 0 ) AllHitCounter0GeV_jet4 +=1;
              for (int GeV=0; GeV<6; GeV++) { // scan over energy 0-5 GeV for the cells near L1 jet 
                if (hcalTPdepth[depthIt] > GeV && depthIt > 0 ) { // is cell over set energy value?  first depth layer is excluded
		  avgTimeJet4[GeV] += hcalTPtiming[depthIt];
		  numHitsJet4[GeV] += 1;
		}
	      }
	      //	      if (depthIt > 0) {
	      if (hcalTPdepth[depthIt] > 1) {
		for (int i=0; i<26; i++) if (hcalTPtiming[depthIt] > i ) DelayedHitCounter1GeV_jet4[i] += 1; // iterator counts the TDC value
		AllHitCounter1GeV_jet4 += 1;
	      }
	      if (hcalTPdepth[depthIt] > 2) {
		for (int i=0; i<26; i++) if (hcalTPtiming[depthIt] > i ) DelayedHitCounter2GeV_jet4[i] += 1; // iterator counts the TDC value
		AllHitCounter2GeV_jet4 += 1;
	      }
	      if (hcalTPdepth[depthIt] > 3) {
		for (int i=0; i<26; i++) if (hcalTPtiming[depthIt] > i ) DelayedHitCounter3GeV_jet4[i] += 1;
		AllHitCounter3GeV_jet4 += 1;
	      }
	      if (hcalTPdepth[depthIt] > 4) {
		for (int i=0; i<26; i++) if (hcalTPtiming[depthIt] > i ) DelayedHitCounter4GeV_jet4[i] += 1;
		AllHitCounter4GeV_jet4 += 1;
	      }
	    }
	  }

	  // count multiplicity based on which jet is closest to the HCAL TP
          // multiplicity counter now 
	  if ( (hcalTPdepth[depthIt] > 3) && (hcalTPtiming[depthIt] > 3) && (abs(tpEtaemu) < 16) ) { // 3 GeV 3ns in HB region
	    if (closestJet == 0) mult3GeV3nsHB_Jet0 += 1;
	    //	    std::cout << hcalTPdepth[depthIt] << " and timing " << hcalTPtiming[depthIt] << std::endl;
	    if (closestJet == 1) mult3GeV3nsHB_Jet1 += 1;
	    if (closestJet == 2) mult3GeV3nsHB_Jet2 += 1;
	    if (closestJet == 3) mult3GeV3nsHB_Jet3 += 1;
	  } // close HB region loop   
          if ((hcalTPtiming[depthIt] > 5) && (abs(tpEtaemu) < 16) ) { // 3 GeV 3ns in HB region
            if (closestJet == 0) mult0GeV5nsHB_Jet0 += 1;
            if (closestJet == 1) mult0GeV5nsHB_Jet1 += 1;
            if (closestJet == 2) mult0GeV5nsHB_Jet2 += 1;
            if (closestJet == 3) mult0GeV5nsHB_Jet3 += 1;
          } // close HB region loop

	  // count multiplicity of layers given a timing and energy threshold   
	  // 3 GeV energy cut
	  if (hcalTPdepth[depthIt] > 3){
	    for (int TDCns = 1; TDCns <= 5; TDCns++) {
	      if (hcalTPtiming[depthIt] > TDCns ) {
		if (TDCns == 2)	{
		  if (depthIt == 0 ) mult3GeV2ns_Jets_depth1 += 1; // used for depth comparisons, this is HB + HE
		  if (depthIt == 1 ) mult3GeV2ns_Jets_depth2 += 1;
		  if (depthIt == 2 ) mult3GeV2ns_Jets_depth3 += 1;
		  if (depthIt == 3 ) mult3GeV2ns_Jets_depth4 += 1;
		  if ( (abs(tpEtaemu) < 16) ) {
		    if (depthIt == 0 ) mult3GeV2ns_Jets_depth1HB += 1;
		    if (depthIt == 1 ) mult3GeV2ns_Jets_depth2HB += 1;
		    if (depthIt == 2 ) mult3GeV2ns_Jets_depth3HB += 1;
		    if (depthIt == 3 ) mult3GeV2ns_Jets_depth4HB += 1;
		  } // closing HB only statement
		} // closing TDC = 2ns statement for depth studies
		mult3GeV_Jets[TDCns] += 1;
		if ( (abs(tpEtaemu) < 16) && depthIt != 0 ) mult3GeVHB_Jets[TDCns] += 1; // 3 GeV HB regions
		if ( (abs(tpEtaemu) >= 16) && (abs(tpEtaemu) < 29) && depthIt != 0 ) mult3GeVHE_Jets[TDCns] += 1; // 3 GeV HE regions
	      } // require timing is greater than TDC ns (1-5ns)
	    } // close TDC ns 1-5 loop
	  } // closing 3 GeV energy cut statement

	  // 0 GeV energy cut
	  if ( (abs(tpEtaemu) < 16) ) {
            for (int TDCns = 1; TDCns <= 5; TDCns++) {
	      if (hcalTPtiming[depthIt] > TDCns ) mult0GeVHB_Jets[TDCns] += 1;
	    }
	  }

	  // 2 GeV energy cut
	  if (hcalTPdepth[depthIt] > 2){
	    for (int TDCns = 1; TDCns <= 5; TDCns++) {
              if (hcalTPtiming[depthIt] > TDCns ) {
		mult2GeV_Jets[TDCns] += 1;
		if ( (abs(tpEtaemu) < 16) && depthIt != 0 ) mult2GeVHB_Jets[TDCns] += 1; // 2 GeV HB regions
		if ( (abs(tpEtaemu) >= 16) && (abs(tpEtaemu) < 29) && depthIt != 0 ) mult2GeVHE_Jets[TDCns] += 1;  // 2 GeV HE regions
              }
            }
	  } // closing 2 GeV energy cut loop

	  // 1 GeV energy cut
	  if (hcalTPdepth[depthIt] > 1){
	    for (int TDCns = 1; TDCns <= 5; TDCns++) {
              if (hcalTPtiming[depthIt] > TDCns ) {
		mult1GeV_Jets[TDCns] += 1;
		if ( (abs(tpEtaemu) < 16) ) mult1GeVHB_Jets[TDCns] += 1; // 1 GeV HB regions
		if ( (abs(tpEtaemu) >= 16) && (abs(tpEtaemu) < 29) ) mult1GeVHE_Jets[TDCns] += 1; // 1 GeV HE regions
              }
            }
	    if (abs(tpEtaemu) < 5 ) mult1GeVcaloT1 += 1;
	    if ((abs(tpEtaemu) > 4) && (abs(tpEtaemu) < 9) ) mult1GeVcaloT2 += 1;
	    if ((abs(tpEtaemu) > 8) && (abs(tpEtaemu) < 13) ) mult1GeVcaloT3 += 1;
	    if ((abs(tpEtaemu) > 12) && (abs(tpEtaemu) < 17) ) mult1GeVcaloT4 += 1;
	  } // closing 1 GeV energy cut loop
	} // closing HCAL depth loop

	if ( (hcalTPdepth[3] / hcalTPdepth[0]) > 1.2) DepthVariableMult += 1;

        // ************************ END OF GEN MATCHING WITHOUT TP DOUBLE COUNTING *************************************
        // *************************************************************************************************************

      } // closing HCAL TP loop
      LLPdecayDetAcceptance->Fill(LLPdecayDetAcc);

      //      std::cout << LLPdecayDetAcc << " LLPs decayed in det volume" << std::endl;
      //      std::cout << "entry = " << jentry << std::endl;

      // htSum histogram to see distribution
      htSumDistribution->Fill(htSum);
      totalGlobal += 1; // how many total events were there? used to compare with global multiplicity sum   
      if (genJetRequirementPassed > 0 ) totalJets += 1; // how many events had at least one jet passing gen matching requirements? used to compare with regional multiplicity sum that is jet matched
      if (totalGlobal != totalJets) std::cout << "still using gen matching requirements!" << std::endl;

      // fill number delayed vs delayed fraction for the jet with highest number of delayed hits
      if ( (DelayedHitCounter1GeV[2] > DelayedHitCounter1GeV_jet2[2]) && (DelayedHitCounter1GeV[2] > DelayedHitCounter1GeV_jet3[2]) && (DelayedHitCounter1GeV[2] > DelayedHitCounter1GeV_jet4[2])) NumberEvents_Fraction_Mult->Fill(DelayedHitCounter1GeV[2] / AllHitCounter1GeV, DelayedHitCounter1GeV[2], 1); // fraction, multiplicity, 1
      else if ( (DelayedHitCounter1GeV_jet2[2] > DelayedHitCounter1GeV_jet3[2]) && (DelayedHitCounter1GeV_jet2[2] > DelayedHitCounter1GeV_jet4[2]) )  NumberEvents_Fraction_Mult->Fill(DelayedHitCounter1GeV_jet2[2] / AllHitCounter1GeV_jet2, DelayedHitCounter1GeV_jet2[2], 1);
      else if ( (DelayedHitCounter1GeV_jet3[2] > DelayedHitCounter1GeV_jet4[2] ) ) NumberEvents_Fraction_Mult->Fill(DelayedHitCounter1GeV_jet3[2] / AllHitCounter1GeV_jet3, DelayedHitCounter1GeV_jet3[2], 1);
      else NumberEvents_Fraction_Mult->Fill(DelayedHitCounter1GeV_jet4[2] / AllHitCounter1GeV_jet4, DelayedHitCounter1GeV_jet4[2], 1);

      for (int i=1; i<8; i++) {
        if (DelayedHitCounter1GeV[i] > min_num_delayed) DelayedHitFraction1GeV[i]->Fill(DelayedHitCounter1GeV[i] / AllHitCounter1GeV);
	if (DelayedHitCounter2GeV[i] > min_num_delayed) DelayedHitFraction2GeV[i]->Fill(DelayedHitCounter2GeV[i] / AllHitCounter2GeV);
	if (DelayedHitCounter3GeV[i] > min_num_delayed) DelayedHitFraction3GeV[i]->Fill(DelayedHitCounter3GeV[i] / AllHitCounter3GeV);
	if (DelayedHitCounter4GeV[i] > min_num_delayed) DelayedHitFraction4GeV[i]->Fill(DelayedHitCounter4GeV[i] / AllHitCounter4GeV);
      }

      if (mult1GeVHB_Jets[2] >= min_num_delayed_jetsum) {
	for (int percent = 0; percent < 11; percent++ ) { // 0.3, 0.4, 0.5, 0.6, 0.7
	  for (int number = 0; number < 11; number++ ) { // 4, 5, 6, 7, 8
	    for (int i=0; i<26; i++) {
	      int jets_passed_1GeV = 0;
	      if (DelayedHitCounter1GeV[i] > min_num_delayed_scan[number] && DelayedHitCounter1GeV[i]/AllHitCounter1GeV > frac_delayed_scan[percent]) jets_passed_1GeV += 1;
	      if (DelayedHitCounter1GeV_jet2[i] > min_num_delayed_scan[number] && DelayedHitCounter1GeV_jet2[i]/AllHitCounter1GeV_jet2 > frac_delayed_scan[percent]) jets_passed_1GeV += 1;
	      if (DelayedHitCounter1GeV_jet3[i] > min_num_delayed_scan[number] && DelayedHitCounter1GeV_jet3[i]/AllHitCounter1GeV_jet3 > frac_delayed_scan[percent]) jets_passed_1GeV += 1;
	      if (DelayedHitCounter1GeV_jet4[i] > min_num_delayed_scan[number] && DelayedHitCounter1GeV_jet4[i]/AllHitCounter1GeV_jet4 > frac_delayed_scan[percent]) jets_passed_1GeV += 1;
	      if ( ((jets_passed_1GeV >= min_num_jets_passed) && (htSum > 120)) || (htSum > 360) ) passedDelayedHitFraction1GeV_ht120[percent][number][i] += 1;
	      if (DelayedHitCounter1GeV[i] > min_num_delayed_scan[number] && DelayedHitCounter1GeV[i]/AllHitCounter1GeV > frac_delayed_scan[percent] ) {
		}
	      else if (DelayedHitCounter1GeV_jet2[i] > min_num_delayed_scan[number] && DelayedHitCounter1GeV_jet2[i]/AllHitCounter1GeV_jet2 > frac_delayed_scan[percent] ) {
	      }
	      else if (DelayedHitCounter1GeV_jet3[i] > min_num_delayed_scan[number] && DelayedHitCounter1GeV_jet3[i]/AllHitCounter1GeV_jet3 > frac_delayed_scan[percent] ) {
	      }
	      else if (DelayedHitCounter1GeV_jet4[i] > min_num_delayed_scan[number] && DelayedHitCounter1GeV_jet4[i]/AllHitCounter1GeV_jet4 > frac_delayed_scan[percent] ) {
	      }
	    
	      int jets_passed_2GeV = 0;
              if (DelayedHitCounter2GeV[i] > min_num_delayed_scan[number] && DelayedHitCounter2GeV[i]/AllHitCounter2GeV > frac_delayed_scan[percent]) jets_passed_2GeV += 1;
              if (DelayedHitCounter2GeV_jet2[i] > min_num_delayed_scan[number] && DelayedHitCounter2GeV_jet2[i]/AllHitCounter2GeV_jet2 > frac_delayed_scan[percent]) jets_passed_2GeV += 1;
              if (DelayedHitCounter2GeV_jet3[i] > min_num_delayed_scan[number] && DelayedHitCounter2GeV_jet3[i]/AllHitCounter2GeV_jet3 > frac_delayed_scan[percent]) jets_passed_2GeV += 1;
              if (DelayedHitCounter2GeV_jet4[i] > min_num_delayed_scan[number] && DelayedHitCounter2GeV_jet4[i]/AllHitCounter2GeV_jet4 > frac_delayed_scan[percent]) jets_passed_2GeV += 1;
              if ( ((jets_passed_2GeV >= min_num_jets_passed) && (htSum > 120)) || (htSum > 360) ) passedDelayedHitFraction2GeV_ht120[percent][number][i] += 1;
	   
	      int jets_passed_3GeV = 0;
              if (DelayedHitCounter3GeV[i] > min_num_delayed_scan[number] && DelayedHitCounter3GeV[i]/AllHitCounter3GeV > frac_delayed_scan[percent]) jets_passed_3GeV += 1;
              if (DelayedHitCounter3GeV_jet2[i] > min_num_delayed_scan[number] && DelayedHitCounter3GeV_jet2[i]/AllHitCounter3GeV_jet2 > frac_delayed_scan[percent]) jets_passed_3GeV += 1;
              if (DelayedHitCounter3GeV_jet3[i] > min_num_delayed_scan[number] && DelayedHitCounter3GeV_jet3[i]/AllHitCounter3GeV_jet3 > frac_delayed_scan[percent]) jets_passed_3GeV += 1;
              if (DelayedHitCounter3GeV_jet4[i] > min_num_delayed_scan[number] && DelayedHitCounter3GeV_jet4[i]/AllHitCounter3GeV_jet4 > frac_delayed_scan[percent]) jets_passed_3GeV += 1;
	      if ( ((jets_passed_3GeV >= min_num_jets_passed) && (htSum > 120)) || (htSum > 360) ) passedDelayedHitFraction3GeV_ht120[percent][number][i] += 1;
	      
	      int jets_passed_4GeV = 0;
              if (DelayedHitCounter4GeV[i] > min_num_delayed_scan[number] && DelayedHitCounter4GeV[i]/AllHitCounter4GeV > frac_delayed_scan[percent]) jets_passed_4GeV += 1;
              if (DelayedHitCounter4GeV_jet2[i] > min_num_delayed_scan[number] && DelayedHitCounter4GeV_jet2[i]/AllHitCounter4GeV_jet2 > frac_delayed_scan[percent]) jets_passed_4GeV += 1;
              if (DelayedHitCounter4GeV_jet3[i] > min_num_delayed_scan[number] && DelayedHitCounter4GeV_jet3[i]/AllHitCounter4GeV_jet3 > frac_delayed_scan[percent]) jets_passed_4GeV += 1;
              if (DelayedHitCounter4GeV_jet4[i] > min_num_delayed_scan[number] && DelayedHitCounter4GeV_jet4[i]/AllHitCounter4GeV_jet4 > frac_delayed_scan[percent]) jets_passed_4GeV += 1;
	      if ( ((jets_passed_4GeV >= min_num_jets_passed) && (htSum > 120)) || (htSum > 360) ) passedDelayedHitFraction4GeV_ht120[percent][number][i] += 1;
	      
	      DelayedHit2D_Fraction1GeV->Fill(i,DelayedHitCounter1GeV[i]/AllHitCounter1GeV,1);
	      DelayedHit2D_Number1GeV->Fill(i,DelayedHitCounter1GeV[i],1);
	      DelayedHit2D_Fraction1GeV->Fill(i,DelayedHitCounter1GeV_jet2[i]/AllHitCounter1GeV_jet2,1);
	      DelayedHit2D_Number1GeV->Fill(i,DelayedHitCounter1GeV_jet2[i],1);
	      DelayedHit2D_Fraction1GeV->Fill(i,DelayedHitCounter1GeV_jet3[i]/AllHitCounter1GeV_jet3,1);
	      DelayedHit2D_Number1GeV->Fill(i,DelayedHitCounter1GeV_jet3[i],1);
	      DelayedHit2D_Fraction1GeV->Fill(i,DelayedHitCounter1GeV_jet4[i]/AllHitCounter1GeV_jet4,1);
	      DelayedHit2D_Number1GeV->Fill(i,DelayedHitCounter1GeV_jet4[i],1);
	      DelayedHit2D_Fraction2GeV->Fill(i,DelayedHitCounter2GeV[i]/AllHitCounter2GeV,1);
	      DelayedHit2D_Number2GeV->Fill(i,DelayedHitCounter2GeV[i],1);
	      DelayedHit2D_Fraction2GeV->Fill(i,DelayedHitCounter2GeV_jet2[i]/AllHitCounter2GeV_jet2,1);
	      DelayedHit2D_Number2GeV->Fill(i,DelayedHitCounter2GeV_jet2[i],1);
	      DelayedHit2D_Fraction2GeV->Fill(i,DelayedHitCounter2GeV_jet3[i]/AllHitCounter2GeV_jet3,1);
	      DelayedHit2D_Number2GeV->Fill(i,DelayedHitCounter2GeV_jet3[i],1);
	      DelayedHit2D_Fraction2GeV->Fill(i,DelayedHitCounter2GeV_jet4[i]/AllHitCounter2GeV_jet4,1);
	      DelayedHit2D_Number2GeV->Fill(i,DelayedHitCounter2GeV_jet4[i],1);
	      DelayedHit2D_Fraction3GeV->Fill(i,DelayedHitCounter3GeV[i]/AllHitCounter3GeV,1);
	      DelayedHit2D_Number3GeV->Fill(i,DelayedHitCounter3GeV[i],1);
	      DelayedHit2D_Fraction3GeV->Fill(i,DelayedHitCounter3GeV_jet2[i]/AllHitCounter3GeV_jet2,1);
	      DelayedHit2D_Number3GeV->Fill(i,DelayedHitCounter3GeV_jet2[i],1);
	      DelayedHit2D_Fraction3GeV->Fill(i,DelayedHitCounter3GeV_jet3[i]/AllHitCounter3GeV_jet3,1);
	      DelayedHit2D_Number3GeV->Fill(i,DelayedHitCounter3GeV_jet3[i],1);
	      DelayedHit2D_Fraction3GeV->Fill(i,DelayedHitCounter3GeV_jet4[i]/AllHitCounter3GeV_jet4,1);
	      DelayedHit2D_Number3GeV->Fill(i,DelayedHitCounter3GeV_jet4[i],1);
	      DelayedHit2D_Fraction4GeV->Fill(i,DelayedHitCounter4GeV[i]/AllHitCounter4GeV,1);
	      DelayedHit2D_Number4GeV->Fill(i,DelayedHitCounter4GeV[i],1);
	      DelayedHit2D_Fraction4GeV->Fill(i,DelayedHitCounter4GeV_jet2[i]/AllHitCounter4GeV_jet2,1);
	      DelayedHit2D_Number4GeV->Fill(i,DelayedHitCounter4GeV_jet2[i],1);
	      DelayedHit2D_Fraction4GeV->Fill(i,DelayedHitCounter4GeV_jet3[i]/AllHitCounter4GeV_jet3,1);
	      DelayedHit2D_Number4GeV->Fill(i,DelayedHitCounter4GeV_jet3[i],1);
	      DelayedHit2D_Fraction4GeV->Fill(i,DelayedHitCounter4GeV_jet4[i]/AllHitCounter4GeV_jet4,1);
	      DelayedHit2D_Number4GeV->Fill(i,DelayedHitCounter4GeV_jet4[i],1);
	    }
	  }
	}
      }
      

      if (genJetRequirementPassed > 0 ) { // fill these histograms when at least one jet has gen matched -- these histograms will sum to give mult3GeV3nsHB_Jets gen matched
	DepthVariable->Fill(DepthVariableMult); // depth variable
	dt3GeV3nsHBJet1Mult_emu->Fill(mult3GeV3nsHB_Jet0);
	dt3GeV3nsHBJet2Mult_emu->Fill(mult3GeV3nsHB_Jet1);
	dt3GeV3nsHBJet3Mult_emu->Fill(mult3GeV3nsHB_Jet2);
	dt3GeV3nsHBJet4Mult_emu->Fill(mult3GeV3nsHB_Jet3);
	dt3GeV3nsHBQuadJetMult_emu->Fill(mult3GeV3nsHB_Jet0+mult3GeV3nsHB_Jet1+mult3GeV3nsHB_Jet2+mult3GeV3nsHB_Jet3);
	dt3GeV3nsHBTripleJetMult_emu->Fill(mult3GeV3nsHB_Jet0+mult3GeV3nsHB_Jet1+mult3GeV3nsHB_Jet2);
	dt3GeV3nsHBDoubleJetMult_emu->Fill(mult3GeV3nsHB_Jet0+mult3GeV3nsHB_Jet1);
	dt3GeV3nsHBSingleJetMult_emu->Fill(mult3GeV3nsHB_Jet0);
        dt0GeV5nsHBJet1Mult_emu->Fill(mult0GeV5nsHB_Jet0);
        dt0GeV5nsHBJet2Mult_emu->Fill(mult0GeV5nsHB_Jet1);
        dt0GeV5nsHBJet3Mult_emu->Fill(mult0GeV5nsHB_Jet2);
        dt0GeV5nsHBJet4Mult_emu->Fill(mult0GeV5nsHB_Jet3);
        dt0GeV5nsHBQuadJetMult_emu->Fill(mult0GeV5nsHB_Jet0+mult0GeV5nsHB_Jet1+mult0GeV5nsHB_Jet2+mult0GeV5nsHB_Jet3);
        dt0GeV5nsHBTripleJetMult_emu->Fill(mult0GeV5nsHB_Jet0+mult0GeV5nsHB_Jet1+mult0GeV5nsHB_Jet2);
        dt0GeV5nsHBDoubleJetMult_emu->Fill(mult0GeV5nsHB_Jet0+mult0GeV5nsHB_Jet1);
        dt0GeV5nsHBSingleJetMult_emu->Fill(mult0GeV5nsHB_Jet0);
	dt3GeV2nsHBQuadJet_depth1_Mult_emu->Fill(mult3GeV2ns_Jets_depth1HB);
	dt3GeV2nsHBQuadJet_depth2_Mult_emu->Fill(mult3GeV2ns_Jets_depth2HB);
	dt3GeV2nsHBQuadJet_depth3_Mult_emu->Fill(mult3GeV2ns_Jets_depth3HB);
	dt3GeV2nsHBQuadJet_depth4_Mult_emu->Fill(mult3GeV2ns_Jets_depth4HB);
        dt3GeV2nsHBHEQuadJet_depth1_Mult_emu->Fill(mult3GeV2ns_Jets_depth1);
        dt3GeV2nsHBHEQuadJet_depth2_Mult_emu->Fill(mult3GeV2ns_Jets_depth2);
        dt3GeV2nsHBHEQuadJet_depth3_Mult_emu->Fill(mult3GeV2ns_Jets_depth3);
        dt3GeV2nsHBHEQuadJet_depth4_Mult_emu->Fill(mult3GeV2ns_Jets_depth4);
      }
      min_DR_L1_TP += 0;    

      if (mult3GeVHB[3] > GeV3ns3Global_threshold ) passedMultGlobal += 1;
      if ( (mult3GeVHB[3] > GeV3ns3Global_threshold) && (htSum > 120) ) passedMultGlobal_120 += 1;
      if ( (mult3GeVHB[3] > GeV3ns3Global_threshold) && (htSum > 350) ) passedMultGlobal_350 += 1;
      //      if (mult0GeV3nsHB_Jets > GeV3ns0Jet_threshold ) passedMultJets += 1;
      if (mult3GeVHB_Jets[3] + mult3GeVHE_Jets[3] > GeV3ns3Jet_threshold ) passedMultJets += 1;
      if (mult3GeVHB_Jets[3] + mult3GeVHE_Jets[3] > 1) passedMultJets3GeV3_1 += 1;
      if (mult3GeVHB_Jets[3] + mult3GeVHE_Jets[3] > 2) passedMultJets3GeV3_2 += 1;
      if (mult3GeVHB_Jets[3] + mult3GeVHE_Jets[3] > 3) passedMultJets3GeV3_3 += 1;
      if (mult3GeVHB_Jets[3] + mult3GeVHE_Jets[3] > 4) passedMultJets3GeV3_4 += 1;
      if (mult3GeVHB_Jets[3] + mult3GeVHE_Jets[3] > 5) passedMultJets3GeV3_5 += 1;
      if (mult3GeVHB_Jets[3] + mult3GeVHE_Jets[3] > 6) passedMultJets3GeV3_6 += 1;
      if (mult3GeVHB_Jets[3] + mult3GeVHE_Jets[3] > 7) passedMultJets3GeV3_7 += 1;
      if ( ((htSum > 120) && mult2GeVHB_Jets[2] + mult2GeVHE_Jets[2] > 1) || (htSum > 360)) passedMultJets3GeV3_ht120_1 += 1;//2] > 1) || (htSum > 360)) passedMultJets3GeV3_ht120_1 += 1; // 2] + mult3GeVHE_Jets[2] > 1+7) || (htSum > 360)) passedMultJets3GeV3_ht120_1 += 1;
      if ( ((htSum > 120) && mult2GeVHB_Jets[2] + mult2GeVHE_Jets[2] > 2) || (htSum > 360)) passedMultJets3GeV3_ht120_2 += 1;//2] > 2) || (htSum > 360)) passedMultJets3GeV3_ht120_2 += 1; // 2] + mult3GeVHE_Jets[2] > 2+7) || (htSum > 360)) passedMultJets3GeV3_ht120_2 += 1;
      if ( ((htSum > 120) && mult2GeVHB_Jets[2] + mult2GeVHE_Jets[2] > 3) || (htSum > 360)) passedMultJets3GeV3_ht120_3 += 1;//2] > 3) || (htSum > 360)) passedMultJets3GeV3_ht120_3 += 1; // 2] + mult3GeVHE_Jets[2] > 3+7) || (htSum > 360)) passedMultJets3GeV3_ht120_3 += 1;
      if ( ((htSum > 120) && mult2GeVHB_Jets[2] + mult2GeVHE_Jets[2] > 4) || (htSum > 360)) passedMultJets3GeV3_ht120_4 += 1;//2] > 4) || (htSum > 360)) passedMultJets3GeV3_ht120_4 += 1; // 2] + mult3GeVHE_Jets[2] > 4+7) || (htSum > 360)) passedMultJets3GeV3_ht120_4 += 1;
      if ( ((htSum > 120) && mult2GeVHB_Jets[2] + mult2GeVHE_Jets[2] > 5) || (htSum > 360)) passedMultJets3GeV3_ht120_5 += 1;//2] > 5) || (htSum > 360)) passedMultJets3GeV3_ht120_5 += 1; // 2] + mult3GeVHE_Jets[2] > 5+7) || (htSum > 360)) passedMultJets3GeV3_ht120_5 += 1;
      if ( ((htSum > 120) && mult2GeVHB_Jets[2] + mult2GeVHE_Jets[2] > 6) || (htSum > 360)) passedMultJets3GeV3_ht120_6 += 1;//2] > 6) || (htSum > 360)) passedMultJets3GeV3_ht120_6 += 1; // 2] + mult3GeVHE_Jets[2] > 6+7) || (htSum > 360)) passedMultJets3GeV3_ht120_6 += 1;
      if ( ((htSum > 120) && mult2GeVHB_Jets[2] + mult2GeVHE_Jets[2] > 7) || (htSum > 360)) passedMultJets3GeV3_ht120_7 += 1;//2] > 7) || (htSum > 360)) passedMultJets3GeV3_ht120_7 += 1; // 2] + mult3GeVHE_Jets[2] > 7+7) || (htSum > 360)) passedMultJets3GeV3_ht120_7 += 1;
      if ( ((mult3GeVHB_Jets[3] + mult3GeVHE_Jets[3] > GeV3ns3Jet_threshold) && (htSum > 120) ) || (htSum > 360)) passedMultJets_120 += 1;
      if ( ((mult3GeVHB_Jets[3] + mult3GeVHE_Jets[3] > GeV3ns3Jet_threshold) && (htSum > 350) ) || (htSum > 360)) passedMultJets_350 += 1;
      if (htSum > 360 ) passedHtSum360 += 1;
      if (htSum > 350 ) passedHtSum350 += 1;
      if (htSum > 120 ) passedHtSum120 += 1;

      //      if (mult3GeVHB_Jets[3] != (mult3GeV3nsHB_Jet0 + mult3GeV3nsHB_Jet1 + mult3GeV3nsHB_Jet2 + mult3GeV3nsHB_Jet3) ) std::cout << "sums are not equal!" << std::endl; // though note, these are only put into histograms when at least one jet is gen matched!

      // after HCAL depth and HCAL TP loops fill the histograms with multiplicity variables. The multiplicity counter is reset on each loop iteration 
      // 3 GeV histograms
      for (int i=1; i<=5; i++) {
	dt3GeVMult_emu[i]->Fill(mult3GeV[i]);
	dt3GeVHEMult_emu[i]->Fill(mult3GeVHE[i]);
	dt3GeVHBMult_emu[i]->Fill(mult3GeVHB[i]);
	dt2GeVMult_emu[i]->Fill(mult2GeV[i]);
        dt2GeVHEMult_emu[i]->Fill(mult2GeVHE[i]);
        dt2GeVHBMult_emu[i]->Fill(mult2GeVHB[i]);
	dt1GeVMult_emu[i]->Fill(mult1GeV[i]);
        dt1GeVHEMult_emu[i]->Fill(mult1GeVHE[i]);
        dt1GeVHBMult_emu[i]->Fill(mult1GeVHB[i]);
      }

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
      if (genJetRequirementPassed > 0 ) {
	for (int i=1; i<=5; i++) {
	  dt3GeVJetMult_emu[i]->Fill(mult3GeV_Jets[i]);
	  dt3GeVHEJetMult_emu[i]->Fill(mult3GeVHE_Jets[i]);
	  dt3GeVHBJetMult_emu[i]->Fill(mult3GeVHB_Jets[i]);
	  dt2GeVJetMult_emu[i]->Fill(mult2GeV_Jets[i]);
          dt2GeVHEJetMult_emu[i]->Fill(mult2GeVHE_Jets[i]);
          dt2GeVHBJetMult_emu[i]->Fill(mult2GeVHB_Jets[i]);
	  dt1GeVJetMult_emu[i]->Fill(mult1GeV_Jets[i]);
          dt1GeVHEJetMult_emu[i]->Fill(mult1GeVHE_Jets[i]);
          dt1GeVHBJetMult_emu[i]->Fill(mult1GeVHB_Jets[i]);
          dt0GeVHBJetMult_emu[i]->Fill(mult0GeVHB_Jets[i]); // 0 GeV
	}
      }

      for (int GeV=0; GeV<6; GeV++) {
	NumberTPavgtiming_jet1[GeV]->Fill(avgTimeJet1[GeV]/numHitsJet1[GeV]);
	NumberTPtotalhits_jet1[GeV]->Fill(numHitsJet1[GeV]);
        NumberTPavgtiming_jet2[GeV]->Fill(avgTimeJet2[GeV]/numHitsJet2[GeV]);
        NumberTPtotalhits_jet2[GeV]->Fill(numHitsJet2[GeV]);
        NumberTPavgtiming_jet3[GeV]->Fill(avgTimeJet3[GeV]/numHitsJet3[GeV]);
        NumberTPtotalhits_jet3[GeV]->Fill(numHitsJet3[GeV]);
        NumberTPavgtiming_jet4[GeV]->Fill(avgTimeJet4[GeV]/numHitsJet4[GeV]);
        NumberTPtotalhits_jet4[GeV]->Fill(numHitsJet4[GeV]);
	if (avgTimeJet1[GeV] > avgTimeJet2[GeV] && avgTimeJet1[GeV] > avgTimeJet3[GeV] && avgTimeJet1[GeV] > avgTimeJet4[GeV] && numHitsJet1[GeV] > 2 ) NumberTPavgtiming_jetMax[GeV]->Fill(avgTimeJet1[GeV]/numHitsJet1[GeV]);
	else if (avgTimeJet2[GeV] > avgTimeJet1[GeV] && avgTimeJet2[GeV] > avgTimeJet3[GeV] && avgTimeJet2[GeV] > avgTimeJet4[GeV] && numHitsJet2[GeV] > 2 ) NumberTPavgtiming_jetMax[GeV]->Fill(avgTimeJet2[GeV]/numHitsJet2[GeV]);
        else if (avgTimeJet3[GeV] > avgTimeJet1[GeV] && avgTimeJet3[GeV] > avgTimeJet2[GeV] && avgTimeJet3[GeV] > avgTimeJet4[GeV] && numHitsJet3[GeV] > 2 ) NumberTPavgtiming_jetMax[GeV]->Fill(avgTimeJet3[GeV]/numHitsJet3[GeV]);
        else if (avgTimeJet4[GeV] > avgTimeJet1[GeV] && avgTimeJet4[GeV] > avgTimeJet2[GeV] && avgTimeJet4[GeV] > avgTimeJet3[GeV] && numHitsJet4[GeV] > 2 ) NumberTPavgtiming_jetMax[GeV]->Fill(avgTimeJet4[GeV]/numHitsJet4[GeV]);
	if (GeV == 3 && avgTimeJet1[GeV]/numHitsJet1[GeV] >=3 ) passedAvgTimeCut += 1;
      }

      // setting the things for the tree used with TMVA, ET_Jet1 has already been set at the top of the event and jet loop
      x0= mult3GeV3nsHB_Jet0;
      x1= mult3GeV3nsHB_Jet1;
      x2= mult3GeV3nsHB_Jet2;
      x3= mult3GeV3nsHB_Jet3;
      x4= htSum;
      x5= AllHitCounter2GeV;
      x6= AllHitCounter2GeV_jet2;
      x7= AllHitCounter2GeV_jet3;
      x8= AllHitCounter2GeV_jet4;
      x9= DelayedHitCounter2GeV[2];
      x10= DelayedHitCounter2GeV_jet2[2];
      x11= DelayedHitCounter2GeV_jet3[2];
      x12= DelayedHitCounter2GeV_jet4[2];
      x13= JetEta1;
      x14= JetEta2;
      x15= JetEta3;
      x16= JetEta4;
      x17= jetEt_1;
      x18= jentry;
      x19= mult3GeVHE_Jets[2];
      x20= mult3GeVHE_Jets[3];
      tree->Fill();
    
      // for each bin fill according to whether our object has a larger corresponding energy
      // Global. nJetBins = 400, jetBinWidht = 1 so this is jetLo + 1, jetLo + 2...etc
      for(int bin=0; bin<nJetBins; bin++){
        if( ((jetEt_1) >= (jetLo + bin*jetBinWidth)) && (mult3GeVHB[3] > GeV3ns3Global_threshold) ) singleJetGlobalRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
        if( ((jetEt_2) >= (jetLo + bin*jetBinWidth)) && (mult3GeVHB[3] > GeV3ns3Global_threshold) ) doubleJetGlobalRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
	if( ((jetEt_3) >= (jetLo + bin*jetBinWidth)) && (mult3GeVHB[3] > GeV3ns3Global_threshold) ) tripleJetGlobalRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
	if( ((jetEt_4) >= (jetLo + bin*jetBinWidth)) && (mult3GeVHB[3] > GeV3ns3Global_threshold) ) quadJetGlobalRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
	// matched to L1 Jets, cut mult3GeV3nsHB_Jets>1 for 1 L1jet matched, mult3GeV3nsHB_Jets>2 for 4 L1 jets matched
    	//      if( ((jetEt_1) >= (jetLo + bin*jetBinWidth)) && (mult0GeV3nsHB_Jets > GeV3ns0Jet_threshold) ) singleJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
	//	if( ((jetEt_2) >= (jetLo + bin*jetBinWidth)) && (mult0GeV3nsHB_Jets > GeV3ns0Jet_threshold) ) doubleJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV  
	//	if( ((jetEt_3) >= (jetLo + bin*jetBinWidth)) && (mult0GeV3nsHB_Jets > GeV3ns0Jet_threshold) ) tripleJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV 
	//	if( ((jetEt_4) >= (jetLo + bin*jetBinWidth)) && (mult0GeV3nsHB_Jets > GeV3ns0Jet_threshold) ) quadJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV

	/*
	if( ((jetEt_1) >= (jetLo + bin*jetBinWidth)) && (mult3GeVHB_Jets[3] > GeV3ns3Jet_threshold) ) singleJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV     
	if( ((jetEt_2) >= (jetLo + bin*jetBinWidth)) && (mult3GeVHB_Jets[3] > GeV3ns3Jet_threshold) ) doubleJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV     
	if( ((jetEt_3) >= (jetLo + bin*jetBinWidth)) && (mult3GeVHB_Jets[3] > GeV3ns3Jet_threshold) ) tripleJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
	if( ((jetEt_4) >= (jetLo + bin*jetBinWidth)) && (mult3GeVHB_Jets[3] > GeV3ns3Jet_threshold) ) quadJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV 
	*/
	if (mult1GeVHB_Jets[2] >= min_num_delayed_jetsum) {
        if( ((jetEt_1) >= (jetLo + bin*jetBinWidth)) && ( (DelayedHitCounter2GeV[3] > min_num_delayed && DelayedHitCounter2GeV[3]/AllHitCounter2GeV > frac_delayed) || (DelayedHitCounter2GeV_jet2[3] > min_num_delayed && DelayedHitCounter2GeV_jet2[3]/AllHitCounter2GeV_jet2 > frac_delayed) || (DelayedHitCounter2GeV_jet3[3] > min_num_delayed && DelayedHitCounter2GeV_jet3[3]/AllHitCounter2GeV_jet3 > frac_delayed) || (DelayedHitCounter2GeV_jet4[3] > min_num_delayed && DelayedHitCounter2GeV_jet4[3]/AllHitCounter2GeV_jet4 > frac_delayed) )  ) singleJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV  
	if( ((jetEt_2) >= (jetLo + bin*jetBinWidth)) && ( (DelayedHitCounter2GeV[3] > min_num_delayed && DelayedHitCounter2GeV[3]/AllHitCounter2GeV > frac_delayed) || (DelayedHitCounter2GeV_jet2[3] > min_num_delayed && DelayedHitCounter2GeV_jet2[3]/AllHitCounter2GeV_jet2 > frac_delayed) || (DelayedHitCounter2GeV_jet3[3] > min_num_delayed && DelayedHitCounter2GeV_jet3[3]/AllHitCounter2GeV_jet3 > frac_delayed) || (DelayedHitCounter2GeV_jet4[3] > min_num_delayed && DelayedHitCounter2GeV_jet4[3]/AllHitCounter2GeV_jet4 > frac_delayed) ) ) doubleJetRates_emu->Fill(jetLo+(bin*jetBinWidth));
	if( ((jetEt_3) >= (jetLo + bin*jetBinWidth)) && ( (DelayedHitCounter2GeV[3] > min_num_delayed && DelayedHitCounter2GeV[3]/AllHitCounter2GeV > frac_delayed) || (DelayedHitCounter2GeV_jet2[3] > min_num_delayed && DelayedHitCounter2GeV_jet2[3]/AllHitCounter2GeV_jet2 > frac_delayed) || (DelayedHitCounter2GeV_jet3[3] > min_num_delayed && DelayedHitCounter2GeV_jet3[3]/AllHitCounter2GeV_jet3 > frac_delayed) || (DelayedHitCounter2GeV_jet4[3] > min_num_delayed && DelayedHitCounter2GeV_jet4[3]/AllHitCounter2GeV_jet4 > frac_delayed) ) ) tripleJetRates_emu->Fill(jetLo+(bin*jetBinWidth));
	if( ((jetEt_4) >= (jetLo + bin*jetBinWidth)) && ( (DelayedHitCounter2GeV[3] > min_num_delayed && DelayedHitCounter2GeV[3]/AllHitCounter2GeV > frac_delayed) || (DelayedHitCounter2GeV_jet2[3] > min_num_delayed && DelayedHitCounter2GeV_jet2[3]/AllHitCounter2GeV_jet2 > frac_delayed) || (DelayedHitCounter2GeV_jet3[3] > min_num_delayed && DelayedHitCounter2GeV_jet3[3]/AllHitCounter2GeV_jet3 > frac_delayed) || (DelayedHitCounter2GeV_jet4[3] > min_num_delayed && DelayedHitCounter2GeV_jet4[3]/AllHitCounter2GeV_jet4 > frac_delayed) ) ) quadJetRates_emu->Fill(jetLo+(bin*jetBinWidth));
	}

	if( ((jetEt_1) >= (jetLo + bin*jetBinWidth)) && (mult3GeVHB_Jets[3] > 1) ) singleJet_th1_Rates_emu->Fill(jetLo+(bin*jetBinWidth));
	if( ((jetEt_1) >= (jetLo + bin*jetBinWidth)) && (mult3GeVHB_Jets[3] > 2) ) singleJet_th2_Rates_emu->Fill(jetLo+(bin*jetBinWidth));
	if( ((jetEt_1) >= (jetLo + bin*jetBinWidth)) && (mult3GeVHB_Jets[3] > 3) ) singleJet_th3_Rates_emu->Fill(jetLo+(bin*jetBinWidth));
	if( ((jetEt_1) >= (jetLo + bin*jetBinWidth)) && (mult3GeVHB_Jets[3] > 4) ) singleJet_th4_Rates_emu->Fill(jetLo+(bin*jetBinWidth));
	if( ((jetEt_1) >= (jetLo + bin*jetBinWidth)) && (mult3GeVHB_Jets[3] > 5) ) singleJet_th5_Rates_emu->Fill(jetLo+(bin*jetBinWidth));
	if( ((jetEt_1) >= (jetLo + bin*jetBinWidth)) && (mult3GeVHB_Jets[3] > 6) ) singleJet_th6_Rates_emu->Fill(jetLo+(bin*jetBinWidth));
	if( ((jetEt_1) >= (jetLo + bin*jetBinWidth)) && (mult3GeVHB_Jets[3] > 7) ) singleJet_th7_Rates_emu->Fill(jetLo+(bin*jetBinWidth));
	if( ((jetEt_4) >= (jetLo + bin*jetBinWidth)) && (mult3GeVHB_Jets[3] > 1) ) quadJet_th1_Rates_emu->Fill(jetLo+(bin*jetBinWidth));
	if( ((jetEt_4) >= (jetLo + bin*jetBinWidth)) && (mult3GeVHB_Jets[3] > 2) ) quadJet_th2_Rates_emu->Fill(jetLo+(bin*jetBinWidth));
	if( ((jetEt_4) >= (jetLo + bin*jetBinWidth)) && (mult3GeVHB_Jets[3] > 3) ) quadJet_th3_Rates_emu->Fill(jetLo+(bin*jetBinWidth));
	if( ((jetEt_4) >= (jetLo + bin*jetBinWidth)) && (mult3GeVHB_Jets[3] > 4) ) quadJet_th4_Rates_emu->Fill(jetLo+(bin*jetBinWidth));
	if( ((jetEt_4) >= (jetLo + bin*jetBinWidth)) && (mult3GeVHB_Jets[3] > 5) ) quadJet_th5_Rates_emu->Fill(jetLo+(bin*jetBinWidth));
	if( ((jetEt_4) >= (jetLo + bin*jetBinWidth)) && (mult3GeVHB_Jets[3] > 6) ) quadJet_th6_Rates_emu->Fill(jetLo+(bin*jetBinWidth));
	if( ((jetEt_4) >= (jetLo + bin*jetBinWidth)) && (mult3GeVHB_Jets[3] > 7) ) quadJet_th7_Rates_emu->Fill(jetLo+(bin*jetBinWidth));

	if( ((jetEt_1) >= (jetLo + bin*jetBinWidth)) ) singleJetOriginalRates_emu->Fill(jetLo+(bin*jetBinWidth));
	if( ((jetEt_4) >= (jetLo + bin*jetBinWidth)) ) quadJetOriginalRates_emu->Fill(jetLo+(bin*jetBinWidth));
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
	//        if( ((htSum) >= htSumLo+(bin*htSumBinWidth)) && (mult0GeV3nsHB_Jets > GeV3ns0Jet_threshold) ) htSumRates_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV
	//	if( ((htSum) >= htSumLo+(bin*htSumBinWidth)) && (mult3GeVHB_Jets[3] > GeV3ns3Jet_threshold) ) htSumRates_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV   
	if ( ((htSum) >= htSumLo+(bin*htSumBinWidth)) && (mult2GeVHB_Jets[2] + mult2GeVHE_Jets[2] > 10 || htSum>360 ) ) htSum4jetRates_emu->Fill(htSumLo+(bin*htSumBinWidth)); // rates for 4 jet sum approach
//2] > GeV3ns3Jet_threshold || htSum>360 ) ) htSum4jetRates_emu->Fill(htSumLo+(bin*htSumBinWidth));
 // 2] + mult3GeVHE_Jets[2] > 10 || htSum>360 ) ) htSum4jetRates_emu->Fill(htSumLo+(bin*htSumBinWidth)); // rates for 4 jet sum approach

	if( ((htSum) >= htSumLo+(bin*htSumBinWidth)) && (mult2GeVHB_Jets[2] + mult2GeVHE_Jets[2] > 1 || htSum>360) ) htSum_th1_Rates_emu->Fill(htSumLo+(bin*htSumBinWidth)); //2] > 1 || htSum>360) ) htSum_th1_Rates_emu->Fill(htSumLo+(bin*htSumBinWidth)); //2] + mult3GeVHE_Jets[2] > 1+7 || htSum>360) ) htSum_th1_Rates_emu->Fill(htSumLo+(bin*htSumBinWidth));
	if( ((htSum) >= htSumLo+(bin*htSumBinWidth)) && (mult2GeVHB_Jets[2] + mult2GeVHE_Jets[2] > 2 || htSum>360) ) htSum_th2_Rates_emu->Fill(htSumLo+(bin*htSumBinWidth)); //2] > 2 || htSum>360) ) htSum_th2_Rates_emu->Fill(htSumLo+(bin*htSumBinWidth));//2] + mult3GeVHE_Jets[2] > 2+7 || htSum>360) ) htSum_th2_Rates_emu->Fill(htSumLo+(bin*htSumBinWidth));
	if( ((htSum) >= htSumLo+(bin*htSumBinWidth)) && (mult2GeVHB_Jets[2] + mult2GeVHE_Jets[2] > 3 || htSum>360) ) htSum_th3_Rates_emu->Fill(htSumLo+(bin*htSumBinWidth));   //2] > 3 || htSum>360) ) htSum_th3_Rates_emu->Fill(htSumLo+(bin*htSumBinWidth));//2] + mult3GeVHE_Jets[2] > 3+7 || htSum>360) ) htSum_th3_Rates_emu->Fill(htSumLo+(bin*htSumBinWidth));
	if( ((htSum) >= htSumLo+(bin*htSumBinWidth)) && (mult2GeVHB_Jets[2] + mult2GeVHE_Jets[2] > 4 || htSum>360) ) htSum_th4_Rates_emu->Fill(htSumLo+(bin*htSumBinWidth));   //2] > 4 || htSum>360) ) htSum_th4_Rates_emu->Fill(htSumLo+(bin*htSumBinWidth));//2] + mult3GeVHE_Jets[2] > 4+7 || htSum>360) ) htSum_th4_Rates_emu->Fill(htSumLo+(bin*htSumBinWidth));
	if( ((htSum) >= htSumLo+(bin*htSumBinWidth)) && (mult2GeVHB_Jets[2] + mult2GeVHE_Jets[2] > 5 || htSum>360) ) htSum_th5_Rates_emu->Fill(htSumLo+(bin*htSumBinWidth));   //2] > 5 || htSum>360) ) htSum_th5_Rates_emu->Fill(htSumLo+(bin*htSumBinWidth));//2] + mult3GeVHE_Jets[2] > 5+7 || htSum>360) ) htSum_th5_Rates_emu->Fill(htSumLo+(bin*htSumBinWidth));
	if( ((htSum) >= htSumLo+(bin*htSumBinWidth)) && (mult2GeVHB_Jets[2] + mult2GeVHE_Jets[2] > 6 || htSum>360) ) htSum_th6_Rates_emu->Fill(htSumLo+(bin*htSumBinWidth));  //2] > 6 || htSum>360) ) htSum_th6_Rates_emu->Fill(htSumLo+(bin*htSumBinWidth));//2] + mult3GeVHE_Jets[2] > 6+7 || htSum>360) ) htSum_th6_Rates_emu->Fill(htSumLo+(bin*htSumBinWidth));
	if( ((htSum) >= htSumLo+(bin*htSumBinWidth)) && (mult2GeVHB_Jets[2] + mult2GeVHE_Jets[2] > 7 || htSum>360) ) htSum_th7_Rates_emu->Fill(htSumLo+(bin*htSumBinWidth)); //2] > 7 || htSum>360) ) htSum_th7_Rates_emu->Fill(htSumLo+(bin*htSumBinWidth));//2] + mult3GeVHE_Jets[2] > 7+7 || htSum>360) ) htSum_th7_Rates_emu->Fill(htSumLo+(bin*htSumBinWidth));

	if (mult1GeVHB_Jets[2] >= min_num_delayed_jetsum) {
	  for (int TDCns=1; TDCns <8; TDCns++) {
	    int jets_passed_1GeV = 0;
	    if (DelayedHitCounter1GeV[TDCns] > min_num_delayed && DelayedHitCounter1GeV[TDCns]/AllHitCounter1GeV > frac_delayed) jets_passed_1GeV += 1;
	    if (DelayedHitCounter1GeV_jet2[TDCns] > min_num_delayed && DelayedHitCounter1GeV_jet2[TDCns]/AllHitCounter1GeV_jet2 > frac_delayed) jets_passed_1GeV += 1;
	    if (DelayedHitCounter1GeV_jet3[TDCns] > min_num_delayed && DelayedHitCounter1GeV_jet3[TDCns]/AllHitCounter1GeV_jet3 > frac_delayed) jets_passed_1GeV += 1;
	    if (DelayedHitCounter1GeV_jet4[TDCns] > min_num_delayed && DelayedHitCounter1GeV_jet4[TDCns]/AllHitCounter1GeV_jet4 > frac_delayed) jets_passed_1GeV += 1;
	    if ((jets_passed_1GeV >= min_num_jets_passed || htSum>360) && ((htSum) >= htSumLo+(bin*htSumBinWidth))) htSum_timing_1GeV_Rates_emu[TDCns]->Fill(htSumLo+(bin*htSumBinWidth)); // rates for scan over energy, timing for fraction delayed approach cut on one jet
            int jets_passed_2GeV = 0;
            if (DelayedHitCounter2GeV[TDCns] > min_num_delayed && DelayedHitCounter2GeV[TDCns]/AllHitCounter2GeV > frac_delayed) jets_passed_2GeV += 1;
            if (DelayedHitCounter2GeV_jet2[TDCns] > min_num_delayed && DelayedHitCounter2GeV_jet2[TDCns]/AllHitCounter2GeV_jet2 > frac_delayed) jets_passed_2GeV += 1;
            if (DelayedHitCounter2GeV_jet3[TDCns] > min_num_delayed && DelayedHitCounter2GeV_jet3[TDCns]/AllHitCounter2GeV_jet3 > frac_delayed) jets_passed_2GeV += 1;
            if (DelayedHitCounter2GeV_jet4[TDCns] > min_num_delayed && DelayedHitCounter2GeV_jet4[TDCns]/AllHitCounter2GeV_jet4 > frac_delayed) jets_passed_2GeV += 1;
            if ((jets_passed_2GeV >= min_num_jets_passed || htSum>360) && ((htSum) >= htSumLo+(bin*htSumBinWidth))) {
	      htSum_timing_2GeV_Rates_emu[TDCns]->Fill(htSumLo+(bin*htSumBinWidth));
	      if (TDCns == 2) htSum1jetRates_emu->Fill(htSumLo+(bin*htSumBinWidth)); // htSum rates plot for single jet tight restriction approach
	    }
            int jets_passed_3GeV = 0;
            if (DelayedHitCounter3GeV[TDCns] > min_num_delayed && DelayedHitCounter3GeV[TDCns]/AllHitCounter3GeV > frac_delayed) jets_passed_3GeV += 1;
            if (DelayedHitCounter3GeV_jet2[TDCns] > min_num_delayed && DelayedHitCounter3GeV_jet2[TDCns]/AllHitCounter3GeV_jet2 > frac_delayed) jets_passed_3GeV += 1;
            if (DelayedHitCounter3GeV_jet3[TDCns] > min_num_delayed && DelayedHitCounter3GeV_jet3[TDCns]/AllHitCounter3GeV_jet3 > frac_delayed) jets_passed_3GeV += 1;
            if (DelayedHitCounter3GeV_jet4[TDCns] > min_num_delayed && DelayedHitCounter3GeV_jet4[TDCns]/AllHitCounter3GeV_jet4 > frac_delayed) jets_passed_3GeV += 1;
            if ((jets_passed_3GeV >= min_num_jets_passed || htSum>360) && ((htSum) >= htSumLo+(bin*htSumBinWidth))) htSum_timing_3GeV_Rates_emu[TDCns]->Fill(htSumLo+(bin*htSumBinWidth));
            int jets_passed_4GeV = 0;
            if (DelayedHitCounter4GeV[TDCns] > min_num_delayed && DelayedHitCounter4GeV[TDCns]/AllHitCounter4GeV > frac_delayed) jets_passed_4GeV += 1;
            if (DelayedHitCounter4GeV_jet2[TDCns] > min_num_delayed && DelayedHitCounter4GeV_jet2[TDCns]/AllHitCounter4GeV_jet2 > frac_delayed) jets_passed_4GeV += 1;
            if (DelayedHitCounter4GeV_jet3[TDCns] > min_num_delayed && DelayedHitCounter4GeV_jet3[TDCns]/AllHitCounter4GeV_jet3 > frac_delayed) jets_passed_4GeV += 1;
            if (DelayedHitCounter4GeV_jet4[TDCns] > min_num_delayed && DelayedHitCounter4GeV_jet4[TDCns]/AllHitCounter4GeV_jet4 > frac_delayed) jets_passed_4GeV += 1;
            if ((jets_passed_4GeV >= min_num_jets_passed || htSum>360) && ((htSum) >= htSumLo+(bin*htSumBinWidth))) htSum_timing_4GeV_Rates_emu[TDCns]->Fill(htSumLo+(bin*htSumBinWidth));
	  }
	}

        if( ((htSum) >= htSumLo+(bin*htSumBinWidth)) && ((mult3GeVHB[2] + mult3GeVHE[2] > GeV3ns3Global_threshold) || (htSum>360)) ) htSumGlobalRates_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV  // htSum rates plot for four jet multiplicity sum approach (GLOBAL)
	if( ((htSum) >= htSumLo+(bin*htSumBinWidth)) ) htSumOriginalRates_emu->Fill(htSumLo+(bin*htSumBinWidth));
      }

      for(int bin=0; bin<nMhtSumBins; bin++){
        if( (mhtSum) >= mhtSumLo+(bin*mhtSumBinWidth) ) mhtSumRates_emu->Fill(mhtSumLo+(bin*mhtSumBinWidth)); //GeV           
      }

      for(int bin=0; bin<nEtSumBins; bin++){
	//        if( ((etSum) >= etSumLo+(bin*etSumBinWidth)) && (mult0GeVHB_Jets[3] > GeV3ns0Jet_threshold) ) etSumRates_emu->Fill(etSumLo+(bin*etSumBinWidth)); //GeV
	//        if( ((etSum) >= etSumLo+(bin*etSumBinWidth)) && (mult3GeVHB_Jets[3] > GeV3ns3Jet_threshold) ) etSumRates_emu->Fill(etSumLo+(bin*etSumBinWidth)); //GeV           
	if (mult1GeVHB_Jets[2] >= min_num_delayed_jetsum) {
	if( ((etSum) >= etSumLo+(bin*etSumBinWidth)) && ( (DelayedHitCounter1GeV[3] > min_num_delayed && DelayedHitCounter1GeV[3]/AllHitCounter1GeV > frac_delayed) || (DelayedHitCounter1GeV_jet2[3] > min_num_delayed && DelayedHitCounter1GeV_jet2[3]/AllHitCounter1GeV_jet2 > frac_delayed) || (DelayedHitCounter1GeV_jet3[3] > min_num_delayed && DelayedHitCounter1GeV_jet3[3]/AllHitCounter1GeV_jet3 > frac_delayed) || (DelayedHitCounter1GeV_jet4[3] > min_num_delayed && DelayedHitCounter1GeV_jet4[3]/AllHitCounter1GeV_jet4 > frac_delayed) ) ) etSumRates_emu->Fill(etSumLo+(bin*etSumBinWidth)); //GeV
	}

        if( ((etSum) >= etSumLo+(bin*etSumBinWidth)) && (mult3GeVHB[3] > GeV3ns3Global_threshold) ) etSumGlobalRates_emu->Fill(etSumLo+(bin*etSumBinWidth)); //GeV  
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
        if( (jetEt_1) >= jetLo + (bin*jetBinWidth) ) singleJetOriginalRates_hw->Fill(jetLo+(bin*jetBinWidth));  //GeV       
      } 

      for(int bin=0; bin<nJetBins; bin++){
        if( (jetEt_2) >= jetLo + (bin*jetBinWidth) ) doubleJetRates_hw->Fill(jetLo+(bin*jetBinWidth));  //GeV
      }  

      for(int bin=0; bin<nJetBins; bin++){
        if( (jetEt_3) >= jetLo + (bin*jetBinWidth) ) tripleJetRates_hw->Fill(jetLo+(bin*jetBinWidth));  //GeV
      }  

      for(int bin=0; bin<nJetBins; bin++){
        if( (jetEt_4) >= jetLo + (bin*jetBinWidth) ) quadJetRates_hw->Fill(jetLo+(bin*jetBinWidth));  //GeV
	if( (jetEt_4) >= jetLo + (bin*jetBinWidth) ) quadJetOriginalRates_hw->Fill(jetLo+(bin*jetBinWidth));  //GeV      
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
        if( (htSum) >= htSumLo+(bin*htSumBinWidth) ) htSum1jetRates_hw->Fill(htSumLo+(bin*htSumBinWidth)); //GeV
        if( (htSum) >= htSumLo+(bin*htSumBinWidth) ) htSum4jetRates_hw->Fill(htSumLo+(bin*htSumBinWidth)); //GeV
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
  std::cout << "num events passed avg time cut 3GeV, >=3ns = " << passedAvgTimeCut << std::endl;

  // saving efficiencies in txt files to be read by plotting macros
  if (inputFile.substr(0,6) == "../QCD" ) {
    std::ofstream DelayedHitFrac_Background;
    DelayedHitFrac_Background.open("DelayedHitFrac_Background.txt", std::ios_base::trunc);
    for (int i=1; i<8; i++) {
      DelayedHitFrac_Background << passedDelayedHitFraction1GeV_ht120[frac_delayed_x10][min_num_delayed][i] / totalJets << std::endl;
    }
    for (int i=1; i<8; i++) {
      DelayedHitFrac_Background << passedDelayedHitFraction2GeV_ht120[frac_delayed_x10][min_num_delayed][i] / totalJets << std::endl;
    }
    for (int i=1; i<8; i++) {
      DelayedHitFrac_Background << passedDelayedHitFraction3GeV_ht120[frac_delayed_x10][min_num_delayed][i] / totalJets << std::endl;
    }
    for (int i=1; i<8; i++) {
      DelayedHitFrac_Background << passedDelayedHitFraction4GeV_ht120[frac_delayed_x10][min_num_delayed][i] / totalJets << std::endl;
    }
    DelayedHitFrac_Background.close();

    std::ofstream DelayedHitFrac_ht120_Background;
    DelayedHitFrac_ht120_Background.open("DelayedHitFrac_ht120_Background_1GeV.txt", std::ios_base::trunc);
    for (int i=1; i<8; i++) DelayedHitFrac_ht120_Background << passedDelayedHitFraction1GeV_ht120[frac_delayed_x10][min_num_delayed][i] / totalJets << std::endl;
    DelayedHitFrac_ht120_Background << passedHtSum360/totalGlobal << std::endl; // efficiency for just a htSum>360 cut
    DelayedHitFrac_ht120_Background.close();
    DelayedHitFrac_ht120_Background.open("DelayedHitFrac_ht120_Background_2GeV.txt", std::ios_base::trunc);
    for (int i=1; i<8; i++) DelayedHitFrac_ht120_Background << passedDelayedHitFraction2GeV_ht120[frac_delayed_x10][min_num_delayed][i] / totalJets << std::endl;
    DelayedHitFrac_ht120_Background << passedHtSum360/totalGlobal << std::endl; // efficiency for just a htSum>360 cut
    DelayedHitFrac_ht120_Background.close();
    DelayedHitFrac_ht120_Background.open("DelayedHitFrac_ht120_Background_3GeV.txt", std::ios_base::trunc);
    for (int i=1; i<8; i++) DelayedHitFrac_ht120_Background << passedDelayedHitFraction3GeV_ht120[frac_delayed_x10][min_num_delayed][i] / totalJets << std::endl;
    DelayedHitFrac_ht120_Background << passedHtSum360/totalGlobal << std::endl; // efficiency for just a htSum>360 cut
    DelayedHitFrac_ht120_Background.close();
    DelayedHitFrac_ht120_Background.open("DelayedHitFrac_ht120_Background_4GeV.txt", std::ios_base::trunc);
    for (int i=1; i<8; i++) DelayedHitFrac_ht120_Background << passedDelayedHitFraction4GeV_ht120[frac_delayed_x10][min_num_delayed][i] / totalJets << std::endl;
    DelayedHitFrac_ht120_Background << passedHtSum360/totalGlobal << std::endl; // efficiency for just a htSum>360 cut
    DelayedHitFrac_ht120_Background.close();
    for (int percent = 0; percent < 11; percent++ ) {
      for (int number = 0; number < 11; number++ ) {
	DelayedHitFrac_ht120_Background.open(Form("DelayedHitFrac_ht120_Background_Fraction%f_Number%u.txt", frac_delayed_scan[percent], min_num_delayed_scan[number]), std::ios_base::trunc);
	for (int i=1; i<8; i++) {
	  DelayedHitFrac_ht120_Background << passedDelayedHitFraction2GeV_ht120[percent][number][i] / totalJets << std::endl;
	}
	DelayedHitFrac_ht120_Background << passedHtSum360/totalGlobal << std::endl; // efficiency for just a htSum>360 cut
	DelayedHitFrac_ht120_Background.close();
      }
    }
    std::ofstream MultiplicityHits3GeV3ns_Background;
    MultiplicityHits3GeV3ns_Background.open("MultiplicityHits3GeV3ns_Background.txt",std::ios_base::trunc);
    MultiplicityHits3GeV3ns_Background << passedMultJets3GeV3_1/totalJets << std::endl;
    MultiplicityHits3GeV3ns_Background << passedMultJets3GeV3_2/totalJets << std::endl;
    MultiplicityHits3GeV3ns_Background << passedMultJets3GeV3_3/totalJets << std::endl;
    MultiplicityHits3GeV3ns_Background << passedMultJets3GeV3_4/totalJets << std::endl;
    MultiplicityHits3GeV3ns_Background << passedMultJets3GeV3_5/totalJets << std::endl;
    MultiplicityHits3GeV3ns_Background << passedMultJets3GeV3_6/totalJets << std::endl;
    MultiplicityHits3GeV3ns_Background << passedMultJets3GeV3_7/totalJets << std::endl;
    MultiplicityHits3GeV3ns_Background.close();
    std::ofstream MultiplicityHits3GeV3ns_ht120_Background;
    MultiplicityHits3GeV3ns_ht120_Background.open("MultiplicityHits3GeV3ns_ht120_Background.txt",std::ios_base::trunc);
    MultiplicityHits3GeV3ns_ht120_Background << passedMultJets3GeV3_ht120_1/totalJets << std::endl;
    MultiplicityHits3GeV3ns_ht120_Background << passedMultJets3GeV3_ht120_2/totalJets << std::endl;
    MultiplicityHits3GeV3ns_ht120_Background << passedMultJets3GeV3_ht120_3/totalJets << std::endl;
    MultiplicityHits3GeV3ns_ht120_Background << passedMultJets3GeV3_ht120_4/totalJets << std::endl;
    MultiplicityHits3GeV3ns_ht120_Background << passedMultJets3GeV3_ht120_5/totalJets << std::endl;
    MultiplicityHits3GeV3ns_ht120_Background << passedMultJets3GeV3_ht120_6/totalJets << std::endl;
    MultiplicityHits3GeV3ns_ht120_Background << passedMultJets3GeV3_ht120_7/totalJets << std::endl;
    MultiplicityHits3GeV3ns_ht120_Background << passedHtSum360/totalGlobal << std::endl; // efficiency for just a htSum>360 cut
    MultiplicityHits3GeV3ns_ht120_Background.close();
  }

  if (inputFile.substr(0,5) == "../mh" && inputFile.substr(0,21) != "../mh1000_pl500_noPU/" ) {
    std::ofstream DelayedHitFrac_Signal;
    DelayedHitFrac_Signal.open(Form("DelayedHitFrac_Signal_%s.txt", inputFile.substr(3,14).c_str()), std::ios_base::trunc);
    for (int i=1; i<8; i++) {
      DelayedHitFrac_Signal << passedDelayedHitFraction1GeV_ht120[frac_delayed_x10][min_num_delayed][i] / totalJets << std::endl;
    }
    for (int i=1; i<8; i++) {
      DelayedHitFrac_Signal << passedDelayedHitFraction2GeV_ht120[frac_delayed_x10][min_num_delayed][i] / totalJets << std::endl;
    }
    for (int i=1; i<8; i++) {
      DelayedHitFrac_Signal << passedDelayedHitFraction3GeV_ht120[frac_delayed_x10][min_num_delayed][i] / totalJets << std::endl;
    }
    for (int i=1; i<8; i++) {
      DelayedHitFrac_Signal << passedDelayedHitFraction4GeV_ht120[frac_delayed_x10][min_num_delayed][i] / totalJets << std::endl;
    }
    DelayedHitFrac_Signal.close();
    std::ofstream DelayedHitFrac_ht120_Signal;
    DelayedHitFrac_ht120_Signal.open(Form("DelayedHitFrac_ht120_Signal_1GeV_%s.txt",inputFile.substr(3,14).c_str()), std::ios_base::trunc);
    for (int i=1; i<8; i++) DelayedHitFrac_ht120_Signal << passedDelayedHitFraction1GeV_ht120[frac_delayed_x10][min_num_delayed][i] / totalJets << std::endl;
    DelayedHitFrac_ht120_Signal << passedHtSum360/totalGlobal << std::endl;
    DelayedHitFrac_ht120_Signal.close();
    DelayedHitFrac_ht120_Signal.open(Form("DelayedHitFrac_ht120_Signal_2GeV_%s.txt",inputFile.substr(3,14).c_str()), std::ios_base::trunc);
    for (int i=1; i<8; i++) DelayedHitFrac_ht120_Signal << passedDelayedHitFraction2GeV_ht120[frac_delayed_x10][min_num_delayed][i] / totalJets << std::endl;
    DelayedHitFrac_ht120_Signal << passedHtSum360/totalGlobal << std::endl;
    DelayedHitFrac_ht120_Signal.close();
    DelayedHitFrac_ht120_Signal.open(Form("DelayedHitFrac_ht120_Signal_3GeV_%s.txt",inputFile.substr(3,14).c_str()), std::ios_base::trunc);
    for (int i=1; i<8; i++) DelayedHitFrac_ht120_Signal << passedDelayedHitFraction3GeV_ht120[frac_delayed_x10][min_num_delayed][i] / totalJets << std::endl;
    DelayedHitFrac_ht120_Signal << passedHtSum360/totalGlobal << std::endl;
    DelayedHitFrac_ht120_Signal.close();
    DelayedHitFrac_ht120_Signal.open(Form("DelayedHitFrac_ht120_Signal_4GeV_%s.txt",inputFile.substr(3,14).c_str()), std::ios_base::trunc);
    for (int i=1; i<8; i++) DelayedHitFrac_ht120_Signal << passedDelayedHitFraction4GeV_ht120[frac_delayed_x10][min_num_delayed][i] / totalJets << std::endl;
    DelayedHitFrac_ht120_Signal << passedHtSum360/totalGlobal << std::endl;
    DelayedHitFrac_ht120_Signal.close();
    // saving efficiences for scanned fraction and number delayed
    for (int percent = 0; percent < 11; percent++ ) {
      for (int number = 0; number < 11; number++ ) {
	DelayedHitFrac_ht120_Signal.open(Form("DelayedHitFrac_ht120_Signal_%s_Fraction%f_Number%u.txt", inputFile.substr(3,14).c_str(), frac_delayed_scan[percent], min_num_delayed_scan[number]), std::ios_base::trunc);
	for (int i=1; i<8; i++) {
	  DelayedHitFrac_ht120_Signal << passedDelayedHitFraction2GeV_ht120[percent][number][i] / totalJets << std::endl;
	}
	DelayedHitFrac_ht120_Signal << passedHtSum360/totalGlobal << std::endl; // efficiency for just a htSum>360 cut
	DelayedHitFrac_ht120_Signal.close();
      }
    }
    std::ofstream MultiplicityHits3GeV3ns_Signal;
    MultiplicityHits3GeV3ns_Signal.open(Form("MultiplicityHits3GeV3ns_Signal_%s.txt", inputFile.substr(3,14).c_str()),std::ios_base::trunc);
    MultiplicityHits3GeV3ns_Signal << passedMultJets3GeV3_1/totalJets << std::endl;
    MultiplicityHits3GeV3ns_Signal << passedMultJets3GeV3_2/totalJets << std::endl;
    MultiplicityHits3GeV3ns_Signal << passedMultJets3GeV3_3/totalJets << std::endl;
    MultiplicityHits3GeV3ns_Signal << passedMultJets3GeV3_4/totalJets << std::endl;
    MultiplicityHits3GeV3ns_Signal << passedMultJets3GeV3_5/totalJets << std::endl;
    MultiplicityHits3GeV3ns_Signal << passedMultJets3GeV3_6/totalJets << std::endl;
    MultiplicityHits3GeV3ns_Signal << passedMultJets3GeV3_7/totalJets << std::endl;
    MultiplicityHits3GeV3ns_Signal.close();
    std::ofstream MultiplicityHits3GeV3ns_ht120_Signal;
    MultiplicityHits3GeV3ns_ht120_Signal.open(Form("MultiplicityHits3GeV3ns_ht120_Signal_%s.txt", inputFile.substr(3,14).c_str()),std::ios_base::trunc);
    MultiplicityHits3GeV3ns_ht120_Signal << passedMultJets3GeV3_ht120_1/totalJets << std::endl;
    MultiplicityHits3GeV3ns_ht120_Signal << passedMultJets3GeV3_ht120_2/totalJets << std::endl;
    MultiplicityHits3GeV3ns_ht120_Signal << passedMultJets3GeV3_ht120_3/totalJets << std::endl;
    MultiplicityHits3GeV3ns_ht120_Signal << passedMultJets3GeV3_ht120_4/totalJets << std::endl;
    MultiplicityHits3GeV3ns_ht120_Signal << passedMultJets3GeV3_ht120_5/totalJets << std::endl;
    MultiplicityHits3GeV3ns_ht120_Signal << passedMultJets3GeV3_ht120_6/totalJets << std::endl;
    MultiplicityHits3GeV3ns_ht120_Signal << passedMultJets3GeV3_ht120_7/totalJets << std::endl;
    MultiplicityHits3GeV3ns_ht120_Signal << passedHtSum360/totalGlobal << std::endl; // efficiency for just a htSum>360 cut 
    MultiplicityHits3GeV3ns_ht120_Signal.close();
  }

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

    singleJet_th1_Rates_emu->Scale(norm);
    singleJet_th2_Rates_emu->Scale(norm);
    singleJet_th3_Rates_emu->Scale(norm);
    singleJet_th4_Rates_emu->Scale(norm);
    singleJet_th5_Rates_emu->Scale(norm);
    singleJet_th6_Rates_emu->Scale(norm);
    singleJet_th7_Rates_emu->Scale(norm);
    quadJet_th1_Rates_emu->Scale(norm);
    quadJet_th2_Rates_emu->Scale(norm);
    quadJet_th3_Rates_emu->Scale(norm);
    quadJet_th4_Rates_emu->Scale(norm);
    quadJet_th5_Rates_emu->Scale(norm);
    quadJet_th6_Rates_emu->Scale(norm);
    quadJet_th7_Rates_emu->Scale(norm);
    htSum_th1_Rates_emu->Scale(norm);
    htSum_th2_Rates_emu->Scale(norm);
    htSum_th3_Rates_emu->Scale(norm);
    htSum_th4_Rates_emu->Scale(norm);
    htSum_th5_Rates_emu->Scale(norm);
    htSum_th6_Rates_emu->Scale(norm);
    htSum_th7_Rates_emu->Scale(norm);
    for (int TDCns = 1; TDCns < 8; TDCns ++) {
      htSum_timing_1GeV_Rates_emu[TDCns]->Scale(norm);
      htSum_timing_2GeV_Rates_emu[TDCns]->Scale(norm);
      htSum_timing_3GeV_Rates_emu[TDCns]->Scale(norm);
      htSum_timing_4GeV_Rates_emu[TDCns]->Scale(norm);
    }

    singleJetGlobalRates_emu->Scale(norm);
    doubleJetGlobalRates_emu->Scale(norm);
    tripleJetGlobalRates_emu->Scale(norm);
    quadJetGlobalRates_emu->Scale(norm);
    singleJetOriginalRates_emu->Scale(norm);
    quadJetOriginalRates_emu->Scale(norm);
    singleEgRates_emu->Scale(norm);
    doubleEgRates_emu->Scale(norm);
    singleTauRates_emu->Scale(norm);
    doubleTauRates_emu->Scale(norm);
    singleISOEgRates_emu->Scale(norm);
    doubleISOEgRates_emu->Scale(norm);
    singleISOTauRates_emu->Scale(norm);
    doubleISOTauRates_emu->Scale(norm);
    htSum1jetRates_emu->Scale(norm);
    htSum4jetRates_emu->Scale(norm);
    htSumGlobalRates_emu->Scale(norm);
    htSumOriginalRates_emu->Scale(norm);
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
    Energy_DepthHE_avg_Jets_1ns = Energy_DepthHE_Jets_1ns->ProfileX();
    Energy_DepthHE_avg_Jets_2ns = Energy_DepthHE_Jets_2ns->ProfileX();
    Energy_DepthHB_avg_Jets = Energy_DepthHB_Jets->ProfileX();
    Energy_DepthHB_avg_Jets_1ns = Energy_DepthHB_Jets_1ns->ProfileX();
    Energy_DepthHB_avg_Jets_2ns = Energy_DepthHB_Jets_2ns->ProfileX();
    Timing_Depth_avg_Jets = Timing_Depth_Jets->ProfileX();
    Timing_DepthHE_avg_Jets = Timing_DepthHE_Jets->ProfileX();
    Timing_DepthHB_avg_Jets = Timing_DepthHB_Jets->ProfileX();

    VolTiming_Depth_avg_Jets = VolTiming_Depth_Jets->ProfileX();
    VolTiming_DepthHB_avg_Jets = VolTiming_DepthHB_Jets->ProfileX();
    VolTiming_DepthHE_avg_Jets = VolTiming_DepthHE_Jets->ProfileX();

    //set the errors for the rates
    //want error -> error * sqrt(norm) ?

    if (inputFile.substr(0,11) == "../Neutrino" ) {
      int SJet60GeV = singleJetGlobalRates_emu->GetBinContent(singleJetGlobalRates_emu->GetXaxis()->FindBin(60)); // get rate value for single jet at 60 GeV ET threshold
      int SJet60GeV_original = singleJetOriginalRates_emu->GetBinContent(singleJetOriginalRates_emu->GetXaxis()->FindBin(60)); // get rate value for single jet at 60 GeV ET threshold  
      int QJet60GeV = quadJetGlobalRates_emu->GetBinContent(quadJetGlobalRates_emu->GetXaxis()->FindBin(60)); // get rate value for quad jet at 60 GeV ET threshold
      int QJet60GeV_original = quadJetOriginalRates_emu->GetBinContent(quadJetOriginalRates_emu->GetXaxis()->FindBin(60)); // get rate value for quad jet at 60 GeV ET threshold 
      int htSum120GeV = htSumGlobalRates_emu->GetBinContent(htSumGlobalRates_emu->GetXaxis()->FindBin(120)); // get rate value for ht sum rate at 120 GeV ET threshold
      int htSum350GeV = htSumGlobalRates_emu->GetBinContent(htSumGlobalRates_emu->GetXaxis()->FindBin(350)); // get rate value for ht sum rate at 350 GeV ET threshold
      int htSum360GeV = htSumGlobalRates_emu->GetBinContent(htSumGlobalRates_emu->GetXaxis()->FindBin(360)); // get rate value for ht sum rate at 360 GeV ET threshold     
      int htSum120GeV_original = htSumOriginalRates_emu->GetBinContent(htSumOriginalRates_emu->GetXaxis()->FindBin(120)); // get rate value for ht sum rate at 120 GeV ET threshold on the original rates plot
      int htSum350GeV_original = htSumOriginalRates_emu->GetBinContent(htSumOriginalRates_emu->GetXaxis()->FindBin(350)); // get rate value for ht sum rate at 350 GeV ET threshold on the original rates plot
      int htSum360GeV_original = htSumOriginalRates_emu->GetBinContent(htSumOriginalRates_emu->GetXaxis()->FindBin(360)); // get rate value for ht sum rate at 360 GeV ET threshold on the original rates plot
      std::cout << "For global multiplicity threshold of >" << GeV3ns3Global_threshold << " the single jet rate at 60 GeV = " << SJet60GeV << " and the quad jet rate = " << QJet60GeV << std::endl;
      std::cout << "For global multiplicity threshold of >" << GeV3ns3Global_threshold << " the htSum=120 rate = " << htSum120GeV <<  " and htSum=350 rate = " << htSum350GeV << " and htSum=360 rate = " << htSum360GeV  << std::endl;
      std::cout << "Original htSum rate at 360 GeV = " << htSum360GeV_original << std::endl;
      std::cout << "Original htSum rate at 350 GeV = " << htSum350GeV_original << std::endl;     
      std::cout << "Original htSum rate at 120 GeV = " << htSum120GeV_original << std::endl;

      int SJet60GeV_l = singleJetRates_emu->GetBinContent(singleJetRates_emu->GetXaxis()->FindBin(60)); // get rate value for single jet at 60 GeV ET threshold
      int QJet60GeV_l = quadJetRates_emu->GetBinContent(quadJetRates_emu->GetXaxis()->FindBin(60)); // get rate value for quad jet at 60 GeV ET threshold  
      int htSum120GeV_l = htSum4jetRates_emu->GetBinContent(htSum4jetRates_emu->GetXaxis()->FindBin(120)); // get rate value for ht sum rate at 120 GeV ET threshold
      int htSum120GeV_l_single = htSum1jetRates_emu->GetBinContent(htSum1jetRates_emu->GetXaxis()->FindBin(120)); // get rate value for ht sum rate at 120 GeV ET threshold  
      std::cout << "new htSum rate 4 jet at 120 GeV = " << htSum120GeV_l << std::endl;
      std::cout << "new htSum rate 1 jet at 120 GeV = " << htSum120GeV_l_single << std::endl;
      int htSum350GeV_l = htSum4jetRates_emu->GetBinContent(htSum4jetRates_emu->GetXaxis()->FindBin(350)); // get rate value for ht sum rate at 350 GeV ET threshold
      int htSum360GeV_l = htSum4jetRates_emu->GetBinContent(htSum4jetRates_emu->GetXaxis()->FindBin(360)); // get rate value for ht sum rate at 360 GeV ET threshold  

      int SJet60GeV_l_th1 = singleJet_th1_Rates_emu->GetBinContent(singleJet_th1_Rates_emu->GetXaxis()->FindBin(60));
      int SJet60GeV_l_th2 = singleJet_th2_Rates_emu->GetBinContent(singleJet_th2_Rates_emu->GetXaxis()->FindBin(60));
      int SJet60GeV_l_th3 = singleJet_th3_Rates_emu->GetBinContent(singleJet_th3_Rates_emu->GetXaxis()->FindBin(60));
      int SJet60GeV_l_th4 = singleJet_th4_Rates_emu->GetBinContent(singleJet_th4_Rates_emu->GetXaxis()->FindBin(60));
      int SJet60GeV_l_th5 = singleJet_th5_Rates_emu->GetBinContent(singleJet_th5_Rates_emu->GetXaxis()->FindBin(60));
      int SJet60GeV_l_th6 = singleJet_th6_Rates_emu->GetBinContent(singleJet_th6_Rates_emu->GetXaxis()->FindBin(60));
      int SJet60GeV_l_th7 = singleJet_th7_Rates_emu->GetBinContent(singleJet_th7_Rates_emu->GetXaxis()->FindBin(60));
      int QJet60GeV_l_th1 = quadJet_th1_Rates_emu->GetBinContent(quadJet_th1_Rates_emu->GetXaxis()->FindBin(60));
      int QJet60GeV_l_th2 = quadJet_th2_Rates_emu->GetBinContent(quadJet_th2_Rates_emu->GetXaxis()->FindBin(60));
      int QJet60GeV_l_th3 = quadJet_th3_Rates_emu->GetBinContent(quadJet_th3_Rates_emu->GetXaxis()->FindBin(60));
      int QJet60GeV_l_th4 = quadJet_th4_Rates_emu->GetBinContent(quadJet_th4_Rates_emu->GetXaxis()->FindBin(60));
      int QJet60GeV_l_th5 = quadJet_th5_Rates_emu->GetBinContent(quadJet_th5_Rates_emu->GetXaxis()->FindBin(60));
      int QJet60GeV_l_th6 = quadJet_th6_Rates_emu->GetBinContent(quadJet_th6_Rates_emu->GetXaxis()->FindBin(60));
      int QJet60GeV_l_th7 = quadJet_th7_Rates_emu->GetBinContent(quadJet_th7_Rates_emu->GetXaxis()->FindBin(60));
      int htSum120GeV_l_th1 = htSum_th1_Rates_emu->GetBinContent(htSum_th1_Rates_emu->GetXaxis()->FindBin(120));
      int htSum120GeV_l_th2 = htSum_th2_Rates_emu->GetBinContent(htSum_th2_Rates_emu->GetXaxis()->FindBin(120));
      int htSum120GeV_l_th3 = htSum_th3_Rates_emu->GetBinContent(htSum_th3_Rates_emu->GetXaxis()->FindBin(120));
      int htSum120GeV_l_th4 = htSum_th4_Rates_emu->GetBinContent(htSum_th4_Rates_emu->GetXaxis()->FindBin(120));
      int htSum120GeV_l_th5 = htSum_th5_Rates_emu->GetBinContent(htSum_th5_Rates_emu->GetXaxis()->FindBin(120));
      int htSum120GeV_l_th6 = htSum_th6_Rates_emu->GetBinContent(htSum_th6_Rates_emu->GetXaxis()->FindBin(120));
      int htSum120GeV_l_th7 = htSum_th7_Rates_emu->GetBinContent(htSum_th7_Rates_emu->GetXaxis()->FindBin(120));
      // frac mult around one jet
      int htSum120GeV_l_timing_1GeV[8] = {0};
      int htSum120GeV_l_timing_2GeV[8] = {0};
      int htSum120GeV_l_timing_3GeV[8] = {0};
      int htSum120GeV_l_timing_4GeV[8] = {0};
      for (int TDCns = 1; TDCns<8; TDCns++) {
	htSum120GeV_l_timing_1GeV[TDCns] = htSum_timing_1GeV_Rates_emu[TDCns]->GetBinContent(htSum_timing_1GeV_Rates_emu[TDCns]->GetXaxis()->FindBin(120));
        htSum120GeV_l_timing_2GeV[TDCns] = htSum_timing_2GeV_Rates_emu[TDCns]->GetBinContent(htSum_timing_2GeV_Rates_emu[TDCns]->GetXaxis()->FindBin(120));
        htSum120GeV_l_timing_3GeV[TDCns] = htSum_timing_3GeV_Rates_emu[TDCns]->GetBinContent(htSum_timing_3GeV_Rates_emu[TDCns]->GetXaxis()->FindBin(120));
        htSum120GeV_l_timing_4GeV[TDCns] = htSum_timing_4GeV_Rates_emu[TDCns]->GetBinContent(htSum_timing_4GeV_Rates_emu[TDCns]->GetXaxis()->FindBin(120));
      }

      std::cout << "For L1 jet matched multiplicity threshold of > " << GeV3ns3Jet_threshold << " the single jet rate at 60 GeV = " << SJet60GeV_l << " and the quad jet rate = " << QJet60GeV_l << std::endl;
      std::cout << "For L1 jet matched multiplicity threshold of > " << GeV3ns3Jet_threshold << " the htSum=120 rate = " << htSum120GeV_l << " and htSum=350 rate = " << htSum350GeV_l << " and htSum=360 rate = " << htSum360GeV_l << std::endl;
      std::cout << inputFile << std::endl;

      std::ofstream single_quad_htSum120_JetRate;
      single_quad_htSum120_JetRate.open("single_quad_htSum120_JetRate.txt", std::ios_base::trunc);
      single_quad_htSum120_JetRate << SJet60GeV_l_th1/1000 << std::endl;
      single_quad_htSum120_JetRate << SJet60GeV_l_th2/1000 << std::endl;
      single_quad_htSum120_JetRate << SJet60GeV_l_th3/1000 << std::endl;
      single_quad_htSum120_JetRate << SJet60GeV_l_th4/1000 << std::endl;
      single_quad_htSum120_JetRate << SJet60GeV_l_th5/1000 << std::endl;
      single_quad_htSum120_JetRate << SJet60GeV_l_th6/1000 << std::endl;
      single_quad_htSum120_JetRate << SJet60GeV_l_th7/1000 << std::endl;
      single_quad_htSum120_JetRate << QJet60GeV_l_th1/1000 << std::endl;
      single_quad_htSum120_JetRate << QJet60GeV_l_th2/1000 << std::endl;
      single_quad_htSum120_JetRate << QJet60GeV_l_th3/1000 << std::endl;
      single_quad_htSum120_JetRate << QJet60GeV_l_th4/1000 << std::endl;
      single_quad_htSum120_JetRate << QJet60GeV_l_th5/1000 << std::endl;
      single_quad_htSum120_JetRate << QJet60GeV_l_th6/1000 << std::endl;
      single_quad_htSum120_JetRate << QJet60GeV_l_th7/1000 << std::endl;
      single_quad_htSum120_JetRate << htSum120GeV_l_th1/1000 << std::endl;
      single_quad_htSum120_JetRate << htSum120GeV_l_th2/1000 <<std::endl;
      single_quad_htSum120_JetRate << htSum120GeV_l_th3/1000 <<std::endl;
      single_quad_htSum120_JetRate << htSum120GeV_l_th4/1000 <<std::endl;
      single_quad_htSum120_JetRate << htSum120GeV_l_th5/1000 <<std::endl;
      single_quad_htSum120_JetRate << htSum120GeV_l_th6/1000 <<std::endl;
      single_quad_htSum120_JetRate << htSum120GeV_l_th7/1000 <<std::endl;
      single_quad_htSum120_JetRate << htSum360GeV_original/1000 << std::endl;
      single_quad_htSum120_JetRate.close();

      std::cout << "saving to htsum timing file" << std::endl;
      std::ofstream htSum120_Rate_1jet;
      htSum120_Rate_1jet.open("htSum120_Rate_1jet_1GeV.txt", std::ios_base::trunc);
      for (int TDCns = 1; TDCns < 8; TDCns ++) htSum120_Rate_1jet << htSum120GeV_l_timing_1GeV[TDCns]/1000 << std::endl;
      htSum120_Rate_1jet << htSum360GeV_original/1000 << std::endl;
      htSum120_Rate_1jet.close();
      htSum120_Rate_1jet.open("htSum120_Rate_1jet_2GeV.txt", std::ios_base::trunc);
      for (int TDCns = 1; TDCns < 8; TDCns ++) htSum120_Rate_1jet << htSum120GeV_l_timing_2GeV[TDCns]/1000 << std::endl;
      htSum120_Rate_1jet << htSum360GeV_original/1000 << std::endl;
      htSum120_Rate_1jet.close();
      htSum120_Rate_1jet.open("htSum120_Rate_1jet_3GeV.txt", std::ios_base::trunc);
      for (int TDCns = 1; TDCns < 8; TDCns ++) htSum120_Rate_1jet << htSum120GeV_l_timing_3GeV[TDCns]/1000 << std::endl;
      htSum120_Rate_1jet << htSum360GeV_original/1000 << std::endl;
      htSum120_Rate_1jet.close();
      htSum120_Rate_1jet.open("htSum120_Rate_1jet_4GeV.txt", std::ios_base::trunc);
      for (int TDCns = 1; TDCns < 8; TDCns ++) htSum120_Rate_1jet << htSum120GeV_l_timing_4GeV[TDCns]/1000 << std::endl;
      htSum120_Rate_1jet << htSum360GeV_original/1000 << std::endl;
      htSum120_Rate_1jet.close();

      std::cout << "saving factor changes" << std::endl;
      EffRate_file << "REGIONAL" << std::endl;
      EffRate_file << "Factor change in htSum Rate at 120 GeV = " << htSum120GeV_original/htSum120GeV_l << std::endl;
      EffRate_file << "Factor change in htSum Rate at 350 GeV = " << htSum350GeV_original/htSum350GeV_l << std::endl;
      EffRate_file << "Factor change in htSum Rate at 360 GeV = " << htSum360GeV_original/htSum360GeV_l << std::endl;
      EffRate_file << "Factor change in single jet rate at 60 GeV = " << SJet60GeV_original/SJet60GeV_l << std::endl;
      EffRate_file << "Factor change in quad jet rate at 60 GeV = " << QJet60GeV_original/QJet60GeV_l << std::endl;
      EffRate_file << "GLOBAL" << std::endl;
      EffRate_file << "Factor change in htSum Rate at 120 GeV = " << htSum120GeV_original/htSum120GeV << std::endl;
      EffRate_file << "Factor change in htSum Rate at 350 GeV = " << htSum350GeV_original/htSum350GeV << std::endl;
      EffRate_file << "Factor change in htSum Rate at 360 GeV = " << htSum360GeV_original/htSum360GeV << std::endl;
      EffRate_file << "Factor change in single jet rate at 60 GeV = " << SJet60GeV_original/SJet60GeV << std::endl;
      EffRate_file << "Factor change in quad jet rate at 60 GeV = " << QJet60GeV_original/QJet60GeV << std::endl;
      EffRate_file << " " << std::endl;
    }
    if (inputFile.substr(0,16) == "../mh1000_pl500_" || inputFile.substr(0,16) == "../mh350__pl500_" ) {
      EffRate_file << "REGIONAL pl=0.5m, Factor reduction in eff with timing and ht=120 cuts as compared with only ht=360 threshold["<< GeV3ns3Jet_threshold-1 << "] = " << (passedMultJets_120/totalJets)/(passedHtSum360/totalGlobal) << ";" << std::endl;
      EffRate_file << "GLOBAL pl=0.5m, Factor reduction in eff with timing and ht=120 cuts as compared with only ht=360 threshold["<< GeV3ns3Global_threshold-1 << "] = " << (passedMultGlobal/totalGlobal)/(passedHtSum360/totalGlobal) << ";" << std::endl;
      EffRate_file << " " << std::endl;
    }
    if (inputFile.substr(0,17) == "../mh1000_pl1000_" || inputFile.substr(0,17) == "../mh350__pl1000_" ) {
      EffRate_file << "REGIONAL pl=1m, Factor reduction in eff with timing and ht=120 cuts as compared with only ht=360 threshold["<< GeV3ns3Jet_threshold-1 << "] = " << (passedMultJets_120/totalJets)/(passedHtSum360/totalGlobal) << ";" << std::endl;
      EffRate_file << "GLOBAL pl=1m, Factor reduction in eff with timing and ht=120 cuts as compared with only ht=360 threshold["<< GeV3ns3Global_threshold-1 << "] = " << (passedMultGlobal/totalGlobal)/(passedHtSum360/totalGlobal) << ";" << std::endl;
      EffRate_file << " " << std::endl;
    }
    if (inputFile.substr(0,6) == "../QCD" ) {
      EffRate_file << "REGIONAL QCD, Factor reduction in eff with timing and ht=120 cuts as compared with only ht=360 threshold["<< GeV3ns3Jet_threshold-1 << "] = " << (passedMultJets_120/totalJets)/(passedHtSum360/totalGlobal) << ";" << std::endl;
      EffRate_file << "GLOBAL QCD, Factor reduction in eff with timing and ht=120 cuts as compared with only ht=360 threshold["<< GeV3ns3Global_threshold-1 << "] = " << (passedMultGlobal/totalGlobal)/(passedHtSum360/totalGlobal) << ";" << std::endl;
      EffRate_file << " " << std::endl;
    }
    std::cout << "writing histograms" <<std::endl;
    htSumDistribution->Write();
    LLPdecayDetAcceptance->Write();
    LLPdecayRadiusDetAcceptance->Write();
    LLPdecayXyzDetAcceptance->Write();
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
    htSum1jetRates_emu->Write();
    htSum4jetRates_emu->Write();
    htSumGlobalRates_emu->Write();
    mhtSumRates_emu->Write();
    etSumRates_emu->Write();
    etSumGlobalRates_emu->Write();
    metSumRates_emu->Write();
    metHFSumRates_emu->Write();
    // 3 GeV
    for (int i=1; i<=5; i++) {
      dt3GeVMult_emu[i]->Write();
      dt3GeVHEMult_emu[i]->Write();
      dt3GeVHBMult_emu[i]->Write();
      dt3GeVJetMult_emu[i]->Write();
      dt3GeVHEJetMult_emu[i]->Write();
      dt3GeVHBJetMult_emu[i]->Write();

      dt2GeVMult_emu[i]->Write();
      dt2GeVHEMult_emu[i]->Write();
      dt2GeVHBMult_emu[i]->Write();
      dt2GeVJetMult_emu[i]->Write();
      dt2GeVHEJetMult_emu[i]->Write();
      dt2GeVHBJetMult_emu[i]->Write();

      dt1GeVMult_emu[i]->Write();
      dt1GeVHEMult_emu[i]->Write();
      dt1GeVHBMult_emu[i]->Write();
      dt1GeVJetMult_emu[i]->Write();
      dt1GeVHEJetMult_emu[i]->Write();
      dt1GeVHBJetMult_emu[i]->Write();

      dt0GeVHBJetMult_emu[i]->Write();
    }

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

    DepthVariable->Write();
    for (int i=1; i<8; i++) {
      DelayedHitFraction1GeV[i]->Write();
      DelayedHitFraction2GeV[i]->Write();
      DelayedHitFraction3GeV[i]->Write();
      DelayedHitFraction4GeV[i]->Write();
    }
    for (int GeV=0; GeV<6; GeV++) {
      NumberTPtiming[GeV]->Write();
      NumberTPavgtiming_jet1[GeV]->Write();
      NumberTPavgtiming_jet2[GeV]->Write();
      NumberTPavgtiming_jet3[GeV]->Write();
      NumberTPavgtiming_jet4[GeV]->Write();
      NumberTPavgtiming_jetMax[GeV]->Write();
      NumberTPtotalhits_jet1[GeV]->Write();
      NumberTPtotalhits_jet2[GeV]->Write();
      NumberTPtotalhits_jet3[GeV]->Write();
      NumberTPtotalhits_jet4[GeV]->Write();
      NumberTPtiming_depth1[GeV]->Write();
      NumberTPtiming_depth2[GeV]->Write();
      NumberTPtiming_depth3[GeV]->Write();
      NumberTPtiming_depth4[GeV]->Write();
    }

    if (inputFile.substr(3,11) != "Neutrino") {
      TCanvas *c1 = new TCanvas("c1","Graph Draw Options",200,10,1200,500);
      c1->Divide(2,1);
      c1->cd(1);
      //      DelayedHit2D_Number1GeV->Draw("COLZ");
      DelayedHit2D_Number1GeV->ProfileX()->GetYaxis()->SetRangeUser(0,15);
      DelayedHit2D_Number1GeV->ProfileX()->GetYaxis()->SetTitle("Number of Delayed Hits");
      DelayedHit2D_Number1GeV->ProfileX()->Draw();
      gPad->SetLogz();
      c1->cd(2);
      //      DelayedHit2D_Fraction1GeV->Draw("COLZ");
      DelayedHit2D_Fraction1GeV->ProfileX()->GetYaxis()->SetRangeUser(0,1);
      DelayedHit2D_Fraction1GeV->ProfileX()->GetYaxis()->SetTitle("Fraction of Delayed Hits");
      DelayedHit2D_Fraction1GeV->ProfileX()->Draw();
      gPad->SetLogz();
      TCanvas *c2 = new TCanvas("c2","Graph Draw Options",200,10,1200,500);
      c2->Divide(2,1);
      c2->cd(1);
      //      DelayedHit2D_Number2GeV->Draw("COLZ");
      DelayedHit2D_Number2GeV->ProfileX()->GetYaxis()->SetRangeUser(0,15);
      DelayedHit2D_Number2GeV->ProfileX()->GetYaxis()->SetTitle("Number of Delayed Hits");
      DelayedHit2D_Number2GeV->ProfileX()->Draw();
      gPad->SetLogz();
      c2->cd(2);
      //      DelayedHit2D_Fraction2GeV->Draw("COLZ");
      DelayedHit2D_Fraction2GeV->ProfileX()->GetYaxis()->SetRangeUser(0,1);
      DelayedHit2D_Fraction2GeV->ProfileX()->GetYaxis()->SetTitle("Fraction of Delayed Hits");
      DelayedHit2D_Fraction2GeV->ProfileX()->Draw();
      gPad->SetLogz();
      TCanvas *c3 = new TCanvas("c3","Graph Draw Options",200,10,1200,500);
      c3->Divide(2,1);
      c3->cd(1);
      //      DelayedHit2D_Number3GeV->Draw("COLZ");
      DelayedHit2D_Number3GeV->ProfileX()->GetYaxis()->SetRangeUser(0,15);
      DelayedHit2D_Number3GeV->ProfileX()->GetYaxis()->SetTitle("Number of Delayed Hits");
      DelayedHit2D_Number3GeV->ProfileX()->Draw();
      gPad->SetLogz();
      c3->cd(2);
      //      DelayedHit2D_Fraction3GeV->Draw("COLZ");
      DelayedHit2D_Fraction3GeV->ProfileX()->GetYaxis()->SetRangeUser(0,1);
      DelayedHit2D_Fraction3GeV->ProfileX()->GetYaxis()->SetTitle("Fraction of Delayed Hits");
      DelayedHit2D_Fraction3GeV->ProfileX()->Draw();
      gPad->SetLogz();
      TCanvas *c4 = new TCanvas("c4","Graph Draw Options",200,10,1200,500);
      c4->Divide(2,1);
      c4->cd(1);
      //      DelayedHit2D_Number4GeV->Draw("COLZ");
      DelayedHit2D_Number4GeV->ProfileX()->GetYaxis()->SetRangeUser(0,15);
      DelayedHit2D_Number4GeV->ProfileX()->GetYaxis()->SetTitle("Number of Delayed Hits");
      DelayedHit2D_Number4GeV->ProfileX()->Draw();
      gPad->SetLogz();
      c4->cd(2);
      //      DelayedHit2D_Fraction4GeV->Draw("COLZ");
      DelayedHit2D_Fraction4GeV->ProfileX()->GetYaxis()->SetRangeUser(0,1);
      DelayedHit2D_Fraction4GeV->ProfileX()->GetYaxis()->SetTitle("Fraction of Delayed Hits");
      DelayedHit2D_Fraction4GeV->ProfileX()->Draw();
      gPad->SetLogz();
      if (inputFile.substr(3,3) == "QCD") {
        c1->Print(Form("DelayedHit2D_Number1GeV_%s.pdf",inputFile.substr(3,3).c_str()));
        c2->Print(Form("DelayedHit2D_Number2GeV_%s.pdf",inputFile.substr(3,3).c_str()));
        c3->Print(Form("DelayedHit2D_Number3GeV_%s.pdf",inputFile.substr(3,3).c_str()));
        c4->Print(Form("DelayedHit2D_Number4GeV_%s.pdf",inputFile.substr(3,3).c_str()));
      }
      if (inputFile.substr(3,17) == "mh1000_pl500_noPU") {
        c1->Print(Form("DelayedHit2D_Number1GeV_%s.pdf",inputFile.substr(3,17).c_str()));
        c2->Print(Form("DelayedHit2D_Number2GeV_%s.pdf",inputFile.substr(3,17).c_str()));
        c3->Print(Form("DelayedHit2D_Number3GeV_%s.pdf",inputFile.substr(3,17).c_str()));
        c4->Print(Form("DelayedHit2D_Number4GeV_%s.pdf",inputFile.substr(3,17).c_str()));
      }
      else {
        c1->Print(Form("DelayedHit2D_Number1GeV_%s.pdf",inputFile.substr(3,14).c_str()));
        c2->Print(Form("DelayedHit2D_Number2GeV_%s.pdf",inputFile.substr(3,14).c_str()));
        c3->Print(Form("DelayedHit2D_Number3GeV_%s.pdf",inputFile.substr(3,14).c_str()));
        c4->Print(Form("DelayedHit2D_Number4GeV_%s.pdf",inputFile.substr(3,14).c_str()));
      }
    }
    DelayedHit2D_Fraction2GeV->Write();
    DelayedHit2D_Number2GeV->Write();

    if (inputFile.substr(3,11) != "Neutrino") {
      TCanvas *c = new TCanvas("c","Graph Draw Options",200,10,600,500);
      c->cd();
      NumberTPtiming_energy->Draw("COLZ");
      gPad->SetLogz();
      if (inputFile.substr(3,3) == "QCD") {
	NumberTPtiming_energy->SetTitle(Form("Number of Cells vs. Energy vs. TDC Values for %s, leading L1 Jet",inputFile.substr(3,3).c_str()));
	c->Print(Form("NumberTPtiming_energy_%s.pdf",inputFile.substr(3,3).c_str()));
      }
      if (inputFile.substr(3,17) == "mh1000_pl500_noPU") {
	NumberTPtiming_energy->SetTitle(Form("Number of Cells vs. Energy vs. TDC Values for %s, leading L1 Jet",inputFile.substr(3,17).c_str()));
	c->Print(Form("NumberTPtiming_energy_%s.pdf",inputFile.substr(3,17).c_str()));
      }
      else {
	NumberTPtiming_energy->SetTitle(Form("Number of Cells vs. Energy vs. TDC Values for %s, leading L1 Jet",inputFile.substr(3,14).c_str()));
	c->Print(Form("NumberTPtiming_energy_%s.pdf",inputFile.substr(3,14).c_str()));
      }
    }
    NumberTPtiming_energy->Write();

    if (inputFile.substr(3,11) != "Neutrino") {
      TCanvas *c = new TCanvas("c","Graph Draw Options",200,10,600,500);
      c->cd();
      NumberEvents_Fraction_Mult->Draw("COLZ");
      //      gPad->SetLogz();
      if (inputFile.substr(3,3) == "QCD") {
        NumberEvents_Fraction_Mult->SetTitle(Form("Number of Events vs. Delayed Fraction vs. Hit Multiplicity for %s, leading L1 Jet",inputFile.substr(3,3).c_str()));
        c->Print(Form("NumberEvents_Fraction_Mult_%s.pdf",inputFile.substr(3,3).c_str()));
      }
      if (inputFile.substr(3,17) == "mh1000_pl500_noPU") {
        NumberEvents_Fraction_Mult->SetTitle(Form("Number of Events vs. Delayed Fraction vs. Hit Multiplicity for %s, leading L1 Jet",inputFile.substr(3,17).c_str()));
        c->Print(Form("NumberEvents_Fraction_Mult_%s.pdf",inputFile.substr(3,17).c_str()));
      }
      else {
        NumberEvents_Fraction_Mult->SetTitle(Form("Number of Events vs. Delayed Fraction vs. Hit Multiplicity for %s, leading L1 Jet",inputFile.substr(3,14).c_str()));
        c->Print(Form("NumberEvents_Fraction_Mult_%s.pdf",inputFile.substr(3,14).c_str()));
      }
    }
    NumberEvents_Fraction_Mult->Write();


    dt3GeV3nsHBJet1Mult_emu->Write();
    dt3GeV3nsHBJet2Mult_emu->Write();
    dt3GeV3nsHBJet3Mult_emu->Write();
    dt3GeV3nsHBJet4Mult_emu->Write();
    dt3GeV3nsHBQuadJetMult_emu->Write();
    dt3GeV3nsHBTripleJetMult_emu->Write();
    dt3GeV3nsHBDoubleJetMult_emu->Write();
    dt3GeV3nsHBSingleJetMult_emu->Write();

    dt3GeV2nsHBQuadJet_depth1_Mult_emu->Write();
    dt3GeV2nsHBQuadJet_depth2_Mult_emu->Write();
    dt3GeV2nsHBQuadJet_depth3_Mult_emu->Write();
    dt3GeV2nsHBQuadJet_depth4_Mult_emu->Write();
    dt3GeV2nsHBHEQuadJet_depth1_Mult_emu->Write();
    dt3GeV2nsHBHEQuadJet_depth2_Mult_emu->Write();
    dt3GeV2nsHBHEQuadJet_depth3_Mult_emu->Write();
    dt3GeV2nsHBHEQuadJet_depth4_Mult_emu->Write();

    dt0GeV5nsHBJet1Mult_emu->Write();
    dt0GeV5nsHBJet2Mult_emu->Write();
    dt0GeV5nsHBJet3Mult_emu->Write();
    dt0GeV5nsHBJet4Mult_emu->Write();
    dt0GeV5nsHBQuadJetMult_emu->Write();
    dt0GeV5nsHBTripleJetMult_emu->Write();
    dt0GeV5nsHBDoubleJetMult_emu->Write();
    dt0GeV5nsHBSingleJetMult_emu->Write();

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
    Energy_DepthHB_Jets_1ns->Write();
    Energy_DepthHB_Jets_2ns->Write();
    Timing_DepthHB_Jets->Write();
    Energy_DepthHE_Jets->Write();
    Energy_DepthHE_Jets_1ns->Write();
    Energy_DepthHE_Jets_2ns->Write();
    Timing_DepthHE_Jets->Write();

    Energy_Depth_Jets_HighE->Write();
    Energy_DepthHB_Jets_HighE->Write();
    Energy_DepthHE_Jets_HighE->Write();
    VolTiming_Depth_Jets->Write();
    VolTiming_DepthHB_Jets->Write();
    VolTiming_DepthHE_Jets->Write();

    Energy_Depth_avg->Write();
    Energy_DepthHE_avg->Write();
    Energy_DepthHB_avg->Write();
    Timing_Depth_avg->Write();
    Timing_DepthHE_avg->Write();
    Timing_DepthHB_avg->Write();

    Energy_Depth_avg_Jets->Write();
    Energy_DepthHE_avg_Jets->Write();
    Energy_DepthHE_avg_Jets_1ns->Write();
    Energy_DepthHE_avg_Jets_2ns->Write();
    Energy_DepthHB_avg_Jets->Write();
    Energy_DepthHB_avg_Jets_1ns->Write();
    Energy_DepthHB_avg_Jets_2ns->Write();
    Timing_Depth_avg_Jets->Write();
    Timing_DepthHE_avg_Jets->Write();
    Timing_DepthHB_avg_Jets->Write();

    VolTiming_Depth_avg_Jets->Write();
    VolTiming_DepthHE_avg_Jets->Write();
    VolTiming_DepthHB_avg_Jets->Write();

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
    htSum1jetRates_hw->Scale(norm);
    htSum4jetRates_hw->Scale(norm);
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
    htSum1jetRates_hw->Write();
    htSum4jetRates_hw->Write();
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
  EffRate_file.close();
}//closes the function 'rates'
