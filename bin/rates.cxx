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
#include "L1Trigger/L1TNtuples/interface/L1AnalysisEventDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoVertexDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisCaloTPDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisGeneratorDataFormat.h"


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

double deltaPhi(double phi1, double phi2) {  // calculate delta phi given two phi values
  double result = phi1 - phi2;
  if(fabs(result) > 9999) return result;
  while (result > TMath::Pi()) result -= 2*TMath::Pi();
  while (result <= -TMath::Pi()) result += 2*TMath::Pi();
  return result;
}

double deltaR(double eta1, double phi1, double eta2, double phi2) { // calculate deltaR given two eta and phi values
  double deta = eta1 - eta2;
  double dphi = deltaPhi(phi1, phi2);
  return sqrt(deta*deta + dphi*dphi);
}

std::vector<double> intersect(double vx, double vy,double vz, double px, double py, double pz) {
  double lightSpeed = 29979245800; // speed of light in cm/s
  double radius = 179; // 130 for calorimeters (ECAL + HCAL) in cm
  double length = 388; // 300 for calorimeters (ECAL + HCAL) in cm
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
  double xPos = tInter*(px/energy)*lightSpeed + vx; // in cm
  double yPos = tInter*(py/energy)*lightSpeed + vy;
  double zPos = tInter*(pz/energy)*lightSpeed + vz;
  // Find eta/phi of intersection                          
  double phi = atan2(yPos,xPos); // return the arc tan in radians                                                                                                                               
  double theta = acos(zPos/sqrt(xPos*xPos + yPos*yPos + zPos*zPos));
  double eta = -log(tan(theta/2.));
  etaphi.push_back(eta);
  etaphi.push_back(phi);
  etaphi.push_back(xPos);
  etaphi.push_back(yPos);
  etaphi.push_back(zPos);
  return etaphi;
}


std::vector<double> closestParton(int L1Jet, L1Analysis::L1AnalysisL1UpgradeDataFormat *l1emu_, L1Analysis::L1AnalysisGeneratorDataFormat *generator_) { 
  // find DR between L1 jet (argument) and a b quark from the LLP 
  double min_dR = 1000;
  double partonNum = -1;
  for (int partonN = 0; partonN < generator_->nPart; partonN ++) {
    double partonEta = 1000;
    double partonPhi = 1000;
    if (generator_->partHardProcess[partonN] == 0 ) continue; 
    if ( (abs(generator_->partId[partonN]) >= 1 && abs(generator_->partId[partonN]) <=5 ) || (abs(generator_->partId[partonN]) == 21) ) { 
      if ( (generator_->partParent[partonN] == 9000006) || (generator_->partParent[partonN] == 9000007) || (generator_->partParent[partonN] == 6000113) ) {
	partonEta = intersect(generator_->partVx[partonN],generator_->partVy[partonN],generator_->partVz[partonN], generator_->partPx[partonN],generator_->partPy[partonN],generator_->partPz[partonN])[0];
	partonPhi = intersect(generator_->partVx[partonN],generator_->partVy[partonN],generator_->partVz[partonN], generator_->partPx[partonN],generator_->partPy[partonN],generator_->partPz[partonN])[1];
	if (l1emu_->jetEt[L1Jet] >= 20 ) {
	  double distance = deltaR(l1emu_->jetEta[L1Jet], l1emu_->jetPhi[L1Jet], partonEta, partonPhi);
	  if (distance < min_dR ) {
	    min_dR = distance;
	    partonNum = (double)partonN;
	  }
	}
      }
    }
  }
  std::vector<double> DR_partonNum;
  DR_partonNum.push_back(min_dR);
  DR_partonNum.push_back(partonNum);
  return DR_partonNum;
}
// need to write another function returning the z and radius position of the LLP associated with the parton closest to this jet
// then make useful distributions based on this gen matching

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

  // histograms based on hit multiplicity from the timing bit
  TH1F * ADC50_3ns_4JetMultHB_emu = new TH1F("ADC50_3ns_4JetMultHB_emu","Number of cells >=50 ADC and >=3ns near 4 L1 Jets (HB, dR(TP,L1Jet)<0.5);Number of cells;Fraction of Entries (normalized)",40,0,40);
  TH1F * ADC50_3ns_4JetMultHE_emu = new TH1F("ADC50_3ns_4JetMultHE_emu","Number of cells >=50 ADC and >=3ns near 4 L1 Jets (HE, dR(TP,L1Jet)<0.5);Number of cells;Fraction of Entries (normalized)",80,0,80);
  TH1F * ADC50_3ns_4JetMultHBHE_emu = new TH1F("ADC50_3ns_4JetMultHBHE_emu","Number of cells >=50 ADC and >=3ns near 4 L1 Jets (HBHE, dR(TP,L1Jet)<0.5);Number of cells;Fraction of Entries (normalized)",80,0,80);

  // HT sum rate distributions to use in rate vs eff plots. Need HT > 360 rate, and HT > 120 + timing cut rate
  TH1F* htSumRates_original_emu = new TH1F("htSumRates_original_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
  TH1F* htSumRates_120timingOR360_emu = new TH1F("htSumRates_120timingOR360_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);

  // GeV / ADC ratios
  std::map<int, TH1F*> GeV_ADC_ratio_HB, GeV_ADC_ratio_HE;
  for (int iEta = 1; iEta <= 16; iEta++) GeV_ADC_ratio_HB[iEta] = new TH1F(Form("GeV_ADC_ratio_HB_%d",iEta),"Energy (GeV) / ADC measurement in HB;Ratio;Number of Entries",50,0,.1);
  for (int iEta = 17; iEta <= 28; iEta++) GeV_ADC_ratio_HE[iEta] = new TH1F(Form("GeV_ADC_ratio_HE_%d",iEta),"Energy (GeV) / ADC measurement in HE;Ratio;Number of Entries",50,0,.1);
  TH2F* ADC_vs_GeV_percell = new TH2F("ADC_vs_GeV_percell","ADC vs. GeV values per cell;ADC;Energy (GeV)",500,0,500,160,0,40);
  TH2F* ADC_vs_GeV_percell_HBHE20 = new TH2F("ADC_vs_GeV_percell_HBHE20","ADC vs. GeV values per cell, abs(iEta)<=20;ADC;Energy (GeV)",500,0,500,160,0,40);
  TH2F* ADC_vs_GeV_percell_HE21 = new TH2F("ADC_vs_GeV_percell_HE21","ADC vs. GeV values per cell in HE, abs(iEta)>=21;ADC;Energy (GeV)",500,0,500,160,0,40);
  TH2F* GeV_ADC_byTPET = new TH2F("GeV_ADC_byTPET","GeV/ADC per cell as a function of TP ET;TP ET (GeV);GeV/ADC per cell",100,0,100,50,0,0.1);
  TH2F* GeV_ADC_byTPET_HBHE20 = new TH2F("GeV_ADC_byTPET_HBHE20","GeV/ADC per cell as a function of TP ET, abs(iEta)<=20;TP ET (GeV);GeV/ADC per cell",100,0,100,50,0,0.1);
  TH2F* GeV_ADC_byTPET_HE21 = new TH2F("GeV_ADC_byTPET_HE21","GeV/ADC per cell as a function of TP ET in HE, abs(iEta)>=21;TP ET (GeV);GeV/ADC per cell",100,0,100,50,0,0.1);

  // gen matching distrbutions
  TH1F* nJet_pt20GeV = new TH1F("nJet_pt20GeV","Number of Jets with a pT > 20 GeV; Number of Jets; Number of Events",15,0,15);
  TH1F* DeltaR_parton_L1jet = new TH1F("DeltaR_parton_L1jet", "Delta R between parton and closest L1 Jet;DeltaR;Number of Entries",50,0,5);
  std::map<int, TH1F*> DeltaR_parton_L1jet_closest;
  for (int closestJet = 0; closestJet < 15; closestJet++) DeltaR_parton_L1jet_closest[closestJet] = new TH1F(Form("DeltaR_parton_L1jet_closest_Jet%d",closestJet), Form("Delta R between a parton and closest L1 Jet (Jet %d);DeltaR;Number of Entries",closestJet),50,0,5);
  TH2F * LLPdecayRadius = new TH2F("LLPdecayRadius", "Decay Radius;Decay Position in cm (z); Decay Radius (cm)",150,-400,400,100,0,250);
  TH1F * LLPdecayXyz = new TH1F("LLPdecayXyz", "LLP ctau;ctau in cm (displacement / beta*gamma);Number of Events",100,0,600);
  TH2F * LLPdecayRadiusTrigAcceptance = new TH2F("LLPdecayRadiusTrigAcceptance", "Decay Radius for LLPs within trigger acceptance;Decay Position in cm (z); Decay Radius (cm)",150,-400,400,100,0,250);
  TH1F * LLPdecayXyzTrigAcceptance = new TH1F("LLPdecayXyzTrigAcceptance", "ctau of LLPs within trigger acceptance;ctau in cm (displacement / beta * gamma);Number of Events",100,0,600);
  TH1F * TOF_LLP_quark = new TH1F("TOF_LLP_quark", "TOF_LLP + TOF_bQuark - TOF_expected; TOF (ns); Number of Events",100,-1,15);
  TH1F * TOF_expected = new TH1F("TOF_expected", "TOF_expected; TOF (ns); Number of Events",100,-5,20);
  TH2F * TOF_vs_TDC = new TH2F("TOF_vs_TDC", "TOF_LLP + TOF_bQuark - TOF_expected vs HCAL TDC within DR<0.2; TOF_LLP + TOF_bQuark - TOF_expected (ns); TDC; Number of Events",100,-5,20,100,-5,20);

  TH1F * htSumDistribution = new TH1F("htSumDistribution","htSum Distribution;HT Sum;Number of Events",100,0,1000);

  // saving rate and efficiencies 
  double passed4JetMult_HBHE_ht120(0); //, passed4JetMult_HB_ht120(0),passed4JetMult_HE_ht120(0);
  double passedHtSum360(0);
  double totalEvents(0);

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

      /*      for(int bin=0; bin<nHtSumBins; bin++){
        if( (htSum) >= htSumLo+(bin*htSumBinWidth)) htSumRates_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV
	}*/

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
      uint nJetemu = l1emu_->nJets; // number of jets per event
      int SumTimingBitJet1_HB = 0;
      int SumTimingBitJet2_HB = 0;
      int SumTimingBitJet3_HB = 0;
      int SumTimingBitJet4_HB = 0;
      int SumTimingBitJet1_HE = 0;
      int SumTimingBitJet2_HE = 0;
      int SumTimingBitJet3_HE = 0;
      int SumTimingBitJet4_HE = 0;
      for (int HcalTPIt = 0; HcalTPIt < nCaloTPemu; HcalTPIt++){
	// eta and phi of the HCAL TP
	double tpEtaemu = l1CaloTPemu_->hcalTPieta[HcalTPIt]; // ieta
	if (abs(tpEtaemu) > 28 ) continue; // don't consider HCAL TPs outside of HB or HE
	double tpPhiemu = l1CaloTPemu_->hcalTPCaliphi[HcalTPIt]; // iphi
	double TP_Eta = etaVal(tpEtaemu); // eta
	double TP_Phi = phiVal(tpPhiemu); // phi
        // for each HCAL TP, find the closest L1 jet   
	double min_DeltaR = 100;
	double DeltaR = 100;
	int closestJet = -1;
	for (uint jetIt = 0; (jetIt < nJetemu) && (jetIt < 4); jetIt++){ // loop over L1 jets, and only do first four (4 highest energy L1 jets from 4 leptons)
      	  if (l1emu_->jetEt[jetIt] < 20 ) continue; // require jet is greater than 20 GeV to attempt matching to HCAL TP
	  double Jet_eta;
	  double Jet_phi;
	  Jet_eta = l1emu_->jetEta[jetIt];
	  Jet_phi = l1emu_->jetPhi[jetIt];
	  DeltaR = deltaR(Jet_eta,Jet_phi,TP_Eta,TP_Phi);   // this is DeltaR for the HCAL TPs to L1 Jet 
	  if (DeltaR < min_DeltaR) {
	    min_DeltaR = DeltaR; // find min delta R between L1 Jet and HCAL TPs -- this is reset for each HCAL TP, so is which jet is closest to the TP
	    closestJet = jetIt; // record which L1 jet is the closest to the HCAL TP under consideration
	  }
	} // closing L1 jet loop
	// HE and HB GEV / ADC ratios, filled by iEta regions
	if (inputFile.substr(27,3) == "QCD" ) { // currently only QCD and nugun have ADC stored in the L1Ntuples
	  if (abs(tpEtaemu) <= 16 ) GeV_ADC_ratio_HB[abs(tpEtaemu)]->Fill( (l1CaloTPemu_->hcalTPDepth1[HcalTPIt] + l1CaloTPemu_->hcalTPDepth2[HcalTPIt] + l1CaloTPemu_->hcalTPDepth3[HcalTPIt] + l1CaloTPemu_->hcalTPDepth4[HcalTPIt]) / (l1CaloTPemu_->hcalTPADC1[HcalTPIt] + l1CaloTPemu_->hcalTPADC2[HcalTPIt] + l1CaloTPemu_->hcalTPADC3[HcalTPIt] + l1CaloTPemu_->hcalTPADC4[HcalTPIt]));
	  if (abs(tpEtaemu) > 16 && abs(tpEtaemu) <=28 ) GeV_ADC_ratio_HE[abs(tpEtaemu)]->Fill( (l1CaloTPemu_->hcalTPDepth1[HcalTPIt] + l1CaloTPemu_->hcalTPDepth2[HcalTPIt] + l1CaloTPemu_->hcalTPDepth3[HcalTPIt] +l1CaloTPemu_->hcalTPDepth4[HcalTPIt] + l1CaloTPemu_->hcalTPDepth5[HcalTPIt] + l1CaloTPemu_->hcalTPDepth6[HcalTPIt] + l1CaloTPemu_->hcalTPDepth7[HcalTPIt]) / (l1CaloTPemu_->hcalTPADC1[HcalTPIt] + l1CaloTPemu_->hcalTPADC2[HcalTPIt] + l1CaloTPemu_->hcalTPADC3[HcalTPIt] + l1CaloTPemu_->hcalTPADC4[HcalTPIt] + l1CaloTPemu_->hcalTPADC5[HcalTPIt] + l1CaloTPemu_->hcalTPADC6[HcalTPIt] + l1CaloTPemu_->hcalTPADC7[HcalTPIt] ));

	  ADC_vs_GeV_percell->Fill(l1CaloTPemu_->hcalTPADC1[HcalTPIt],l1CaloTPemu_->hcalTPDepth1[HcalTPIt]);
	  ADC_vs_GeV_percell->Fill(l1CaloTPemu_->hcalTPADC2[HcalTPIt],l1CaloTPemu_->hcalTPDepth2[HcalTPIt]);
	  ADC_vs_GeV_percell->Fill(l1CaloTPemu_->hcalTPADC3[HcalTPIt],l1CaloTPemu_->hcalTPDepth3[HcalTPIt]);
	  ADC_vs_GeV_percell->Fill(l1CaloTPemu_->hcalTPADC4[HcalTPIt],l1CaloTPemu_->hcalTPDepth4[HcalTPIt]);
	  ADC_vs_GeV_percell->Fill(l1CaloTPemu_->hcalTPADC5[HcalTPIt],l1CaloTPemu_->hcalTPDepth5[HcalTPIt]);
	  ADC_vs_GeV_percell->Fill(l1CaloTPemu_->hcalTPADC6[HcalTPIt],l1CaloTPemu_->hcalTPDepth6[HcalTPIt]);
	  ADC_vs_GeV_percell->Fill(l1CaloTPemu_->hcalTPADC7[HcalTPIt],l1CaloTPemu_->hcalTPDepth7[HcalTPIt]);
	  
	  GeV_ADC_byTPET->Fill(l1CaloTPemu_->hcalTPet[HcalTPIt], (l1CaloTPemu_->hcalTPDepth1[HcalTPIt] / l1CaloTPemu_->hcalTPADC1[HcalTPIt]));
	  GeV_ADC_byTPET->Fill(l1CaloTPemu_->hcalTPet[HcalTPIt], (l1CaloTPemu_->hcalTPDepth2[HcalTPIt] / l1CaloTPemu_->hcalTPADC2[HcalTPIt]));
	  GeV_ADC_byTPET->Fill(l1CaloTPemu_->hcalTPet[HcalTPIt], (l1CaloTPemu_->hcalTPDepth3[HcalTPIt] / l1CaloTPemu_->hcalTPADC3[HcalTPIt]));
	  GeV_ADC_byTPET->Fill(l1CaloTPemu_->hcalTPet[HcalTPIt], (l1CaloTPemu_->hcalTPDepth4[HcalTPIt] / l1CaloTPemu_->hcalTPADC4[HcalTPIt]));
	  GeV_ADC_byTPET->Fill(l1CaloTPemu_->hcalTPet[HcalTPIt], (l1CaloTPemu_->hcalTPDepth5[HcalTPIt] / l1CaloTPemu_->hcalTPADC5[HcalTPIt]));
	  GeV_ADC_byTPET->Fill(l1CaloTPemu_->hcalTPet[HcalTPIt], (l1CaloTPemu_->hcalTPDepth6[HcalTPIt] / l1CaloTPemu_->hcalTPADC6[HcalTPIt]));
	  GeV_ADC_byTPET->Fill(l1CaloTPemu_->hcalTPet[HcalTPIt], (l1CaloTPemu_->hcalTPDepth7[HcalTPIt] / l1CaloTPemu_->hcalTPADC7[HcalTPIt]));

	  if (abs(tpEtaemu) <= 20 ) {
	    ADC_vs_GeV_percell_HBHE20->Fill(l1CaloTPemu_->hcalTPADC1[HcalTPIt],l1CaloTPemu_->hcalTPDepth1[HcalTPIt]);
	    ADC_vs_GeV_percell_HBHE20->Fill(l1CaloTPemu_->hcalTPADC2[HcalTPIt],l1CaloTPemu_->hcalTPDepth2[HcalTPIt]);
	    ADC_vs_GeV_percell_HBHE20->Fill(l1CaloTPemu_->hcalTPADC3[HcalTPIt],l1CaloTPemu_->hcalTPDepth3[HcalTPIt]);
	    ADC_vs_GeV_percell_HBHE20->Fill(l1CaloTPemu_->hcalTPADC4[HcalTPIt],l1CaloTPemu_->hcalTPDepth4[HcalTPIt]);
	    ADC_vs_GeV_percell_HBHE20->Fill(l1CaloTPemu_->hcalTPADC5[HcalTPIt],l1CaloTPemu_->hcalTPDepth5[HcalTPIt]);
	    ADC_vs_GeV_percell_HBHE20->Fill(l1CaloTPemu_->hcalTPADC6[HcalTPIt],l1CaloTPemu_->hcalTPDepth6[HcalTPIt]);
	    ADC_vs_GeV_percell_HBHE20->Fill(l1CaloTPemu_->hcalTPADC7[HcalTPIt],l1CaloTPemu_->hcalTPDepth7[HcalTPIt]);
	    GeV_ADC_byTPET_HBHE20->Fill(l1CaloTPemu_->hcalTPet[HcalTPIt], (l1CaloTPemu_->hcalTPDepth1[HcalTPIt] / l1CaloTPemu_->hcalTPADC1[HcalTPIt]));
	    GeV_ADC_byTPET_HBHE20->Fill(l1CaloTPemu_->hcalTPet[HcalTPIt], (l1CaloTPemu_->hcalTPDepth2[HcalTPIt] / l1CaloTPemu_->hcalTPADC2[HcalTPIt]));
	    GeV_ADC_byTPET_HBHE20->Fill(l1CaloTPemu_->hcalTPet[HcalTPIt], (l1CaloTPemu_->hcalTPDepth3[HcalTPIt] / l1CaloTPemu_->hcalTPADC3[HcalTPIt]));
	    GeV_ADC_byTPET_HBHE20->Fill(l1CaloTPemu_->hcalTPet[HcalTPIt], (l1CaloTPemu_->hcalTPDepth4[HcalTPIt] / l1CaloTPemu_->hcalTPADC4[HcalTPIt]));
	    GeV_ADC_byTPET_HBHE20->Fill(l1CaloTPemu_->hcalTPet[HcalTPIt], (l1CaloTPemu_->hcalTPDepth5[HcalTPIt] / l1CaloTPemu_->hcalTPADC5[HcalTPIt]));
	    GeV_ADC_byTPET_HBHE20->Fill(l1CaloTPemu_->hcalTPet[HcalTPIt], (l1CaloTPemu_->hcalTPDepth6[HcalTPIt] / l1CaloTPemu_->hcalTPADC6[HcalTPIt]));
	    GeV_ADC_byTPET_HBHE20->Fill(l1CaloTPemu_->hcalTPet[HcalTPIt], (l1CaloTPemu_->hcalTPDepth7[HcalTPIt] / l1CaloTPemu_->hcalTPADC7[HcalTPIt]));
	  }
	  if (abs(tpEtaemu) > 20 && abs(tpEtaemu) <=28 ) {
	    ADC_vs_GeV_percell_HE21->Fill(l1CaloTPemu_->hcalTPADC1[HcalTPIt],l1CaloTPemu_->hcalTPDepth1[HcalTPIt]);
	    ADC_vs_GeV_percell_HE21->Fill(l1CaloTPemu_->hcalTPADC2[HcalTPIt],l1CaloTPemu_->hcalTPDepth2[HcalTPIt]);
	    ADC_vs_GeV_percell_HE21->Fill(l1CaloTPemu_->hcalTPADC3[HcalTPIt],l1CaloTPemu_->hcalTPDepth3[HcalTPIt]);
	    ADC_vs_GeV_percell_HE21->Fill(l1CaloTPemu_->hcalTPADC4[HcalTPIt],l1CaloTPemu_->hcalTPDepth4[HcalTPIt]);
	    ADC_vs_GeV_percell_HE21->Fill(l1CaloTPemu_->hcalTPADC5[HcalTPIt],l1CaloTPemu_->hcalTPDepth5[HcalTPIt]);
	    ADC_vs_GeV_percell_HE21->Fill(l1CaloTPemu_->hcalTPADC6[HcalTPIt],l1CaloTPemu_->hcalTPDepth6[HcalTPIt]);
	    ADC_vs_GeV_percell_HE21->Fill(l1CaloTPemu_->hcalTPADC7[HcalTPIt],l1CaloTPemu_->hcalTPDepth7[HcalTPIt]);
	    GeV_ADC_byTPET_HE21->Fill(l1CaloTPemu_->hcalTPet[HcalTPIt], (l1CaloTPemu_->hcalTPDepth1[HcalTPIt] / l1CaloTPemu_->hcalTPADC1[HcalTPIt]));
	    GeV_ADC_byTPET_HE21->Fill(l1CaloTPemu_->hcalTPet[HcalTPIt], (l1CaloTPemu_->hcalTPDepth2[HcalTPIt] / l1CaloTPemu_->hcalTPADC2[HcalTPIt]));
	    GeV_ADC_byTPET_HE21->Fill(l1CaloTPemu_->hcalTPet[HcalTPIt], (l1CaloTPemu_->hcalTPDepth3[HcalTPIt] / l1CaloTPemu_->hcalTPADC3[HcalTPIt]));
	    GeV_ADC_byTPET_HE21->Fill(l1CaloTPemu_->hcalTPet[HcalTPIt], (l1CaloTPemu_->hcalTPDepth4[HcalTPIt] / l1CaloTPemu_->hcalTPADC4[HcalTPIt]));
	    GeV_ADC_byTPET_HE21->Fill(l1CaloTPemu_->hcalTPet[HcalTPIt], (l1CaloTPemu_->hcalTPDepth5[HcalTPIt] / l1CaloTPemu_->hcalTPADC5[HcalTPIt]));
	    GeV_ADC_byTPET_HE21->Fill(l1CaloTPemu_->hcalTPet[HcalTPIt], (l1CaloTPemu_->hcalTPDepth6[HcalTPIt] / l1CaloTPemu_->hcalTPADC6[HcalTPIt]));
	    GeV_ADC_byTPET_HE21->Fill(l1CaloTPemu_->hcalTPet[HcalTPIt], (l1CaloTPemu_->hcalTPDepth7[HcalTPIt] / l1CaloTPemu_->hcalTPADC7[HcalTPIt]));
	  }
	}

	if (min_DeltaR > 0.5 ) continue; // don't consider HCAL TPs that are greater than DR 0.5 from a L1 jet
	int TP_TimingBit = l1CaloTPemu_->hcalTPTimingBit[HcalTPIt];
	if (abs(tpEtaemu) <= 16 ) {
	  if (closestJet == 0 ) SumTimingBitJet1_HB += TP_TimingBit;
          if (closestJet == 1 ) SumTimingBitJet2_HB += TP_TimingBit;
          if (closestJet == 2 ) SumTimingBitJet3_HB += TP_TimingBit;
          if (closestJet == 3 ) SumTimingBitJet4_HB += TP_TimingBit;
	}
	if (abs(tpEtaemu) > 16 && abs(tpEtaemu) <= 28 ) {
          if (closestJet == 0 ) SumTimingBitJet1_HE += TP_TimingBit;
          if (closestJet == 1 ) SumTimingBitJet2_HE += TP_TimingBit;
          if (closestJet == 2 ) SumTimingBitJet3_HE += TP_TimingBit;
          if (closestJet == 3 ) SumTimingBitJet4_HE += TP_TimingBit;
	}
      	// each TP has the Timing Bit set with the multiplicity for that tower
      } // closing HCAL TP loop
      int Sum4Jet_HBHE = SumTimingBitJet1_HB+SumTimingBitJet2_HB+SumTimingBitJet3_HB+SumTimingBitJet4_HB + SumTimingBitJet1_HE+SumTimingBitJet2_HE+SumTimingBitJet3_HE+SumTimingBitJet4_HE;
      int Sum4Jet_HB = SumTimingBitJet1_HB+SumTimingBitJet2_HB+SumTimingBitJet3_HB+SumTimingBitJet4_HB;
      int Sum4Jet_HE = SumTimingBitJet1_HE+SumTimingBitJet2_HE+SumTimingBitJet3_HE+SumTimingBitJet4_HE;
    
      // ************* GEN PARTICLE MATCHING CODE ******************
      //      for (uint jetIt = 0; jetIt < nJetemu; jetIt++)  std::cout << closestParton(jetIt, l1emu_, generator_)[0] << " " << closestParton(jetIt, l1emu_, generator_)[1] << std::endl;
      //      std::cout << " " << std::endl;

      for (int partonN = 0; partonN < generator_->nPart; partonN ++) {
	double partonEta = 1000;
        double partonPhi = 1000;

        if (generator_->partHardProcess[partonN] == 0 ) continue;                  
        if ( (abs(generator_->partId[partonN]) >= 1 && abs(generator_->partId[partonN]) <=5 ) || (abs(generator_->partId[partonN]) == 21) ) { 
          if ( (generator_->partParent[partonN] == 9000006) || (generator_->partParent[partonN] == 9000007) || (generator_->partParent[partonN] == 6000113) ) { // 6000113 possibly, or 9000007, or 9000006
            partonEta = intersect(generator_->partVx[partonN],generator_->partVy[partonN],generator_->partVz[partonN], generator_->partPx[partonN],generator_->partPy[partonN],generator_->partPz[partonN])[0];
            partonPhi = intersect(generator_->partVx[partonN],generator_->partVy[partonN],generator_->partVz[partonN], generator_->partPx[partonN],generator_->partPy[partonN],generator_->partPz[partonN])[1];
	    double partonEta1 = partonEta; // used for TOF and TDC comparisons, partonEta and partonPhi above use in DeltaR from L1 jet
	    double partonPhi1 = partonPhi;
	    double partonEta2 = intersect(generator_->partVx[partonN-1],generator_->partVy[partonN-1],generator_->partVz[partonN-1], generator_->partPx[partonN-1],generator_->partPy[partonN-1],generator_->partPz[partonN-1])[0];
	    double partonPhi2 = intersect(generator_->partVx[partonN-1],generator_->partVy[partonN-1],generator_->partVz[partonN-1], generator_->partPx[partonN-1],generator_->partPy[partonN-1],generator_->partPz[partonN-1])[1];
	    double parton1HCALx = intersect(generator_->partVx[partonN],generator_->partVy[partonN],generator_->partVz[partonN], generator_->partPx[partonN],generator_->partPy[partonN],generator_->partPz[partonN])[2]; // in cm
	    double parton1HCALy = intersect(generator_->partVx[partonN],generator_->partVy[partonN],generator_->partVz[partonN], generator_->partPx[partonN],generator_->partPy[partonN],generator_->partPz[partonN])[3];
	    double parton1HCALz = intersect(generator_->partVx[partonN],generator_->partVy[partonN],generator_->partVz[partonN], generator_->partPx[partonN],generator_->partPy[partonN],generator_->partPz[partonN])[4];
            double parton2HCALx = intersect(generator_->partVx[partonN-1],generator_->partVy[partonN-1],generator_->partVz[partonN-1], generator_->partPx[partonN-1],generator_->partPy[partonN-1],generator_->partPz[partonN-1])[2];
            double parton2HCALy = intersect(generator_->partVx[partonN-1],generator_->partVy[partonN-1],generator_->partVz[partonN-1], generator_->partPx[partonN-1],generator_->partPy[partonN-1],generator_->partPz[partonN-1])[3];
            double parton2HCALz = intersect(generator_->partVx[partonN-1],generator_->partVy[partonN-1],generator_->partVz[partonN-1], generator_->partPx[partonN-1],generator_->partPy[partonN-1],generator_->partPz[partonN-1])[4];
	    double radius = sqrt(generator_->partVx[partonN]*generator_->partVx[partonN] + generator_->partVy[partonN]*generator_->partVy[partonN]);
	    double vertex = sqrt(generator_->partVx[partonN]*generator_->partVx[partonN] + generator_->partVy[partonN]*generator_->partVy[partonN] + generator_->partVz[partonN]*generator_->partVz[partonN]);
	    if ( generator_->partVz[partonN] == generator_->partVz[partonN-1]  ) { 
	      // makes sure found two quarks resulting from LLP decay since they have the same vertex, two b quarks are always separated by 1 in the partonN loop
	      double genLLPBeta = sqrt((generator_->partPx[partonN-1] + generator_->partPx[partonN])*(generator_->partPx[partonN-1] + generator_->partPx[partonN]) + (generator_->partPy[partonN-1] + generator_->partPy[partonN])*(generator_->partPy[partonN-1] + generator_->partPy[partonN]) + (generator_->partPz[partonN-1] + generator_->partPz[partonN])*(generator_->partPz[partonN-1] + generator_->partPz[partonN])) / (generator_->partE[partonN-1] + generator_->partE[partonN]); // beta of LLP considering momentum and energy of two b quarks
	      double genLLPGamma = 1./TMath::Sqrt(1.-genLLPBeta*genLLPBeta);
	      LLPdecayRadius->Fill(generator_->partVz[partonN],radius); // fill the radius and z position for LLP decay
	      LLPdecayXyz->Fill(vertex / (genLLPBeta * genLLPGamma)); // ctau traveled in lab frame

	      if (Sum4Jet_HBHE >= 10 && htSum > 120) {
		LLPdecayRadiusTrigAcceptance->Fill(generator_->partVz[partonN],radius); // fill the radius and z position for LLP decay given that event passes trigger selection
		LLPdecayXyzTrigAcceptance->Fill(vertex / (genLLPBeta * genLLPGamma));
	      }

	      // calculate TOF
	      double lightSpeed = 29979245800; // in cm / s
	      double TOF_LLP = 1000000000*vertex / (genLLPBeta * lightSpeed); // in ns
	      double TOF_bQuark1 = 1000000000*(TMath::Sqrt((generator_->partVx[partonN] - parton1HCALx)*(generator_->partVx[partonN] - parton1HCALx) + (generator_->partVy[partonN] - parton1HCALy)*(generator_->partVy[partonN] -parton1HCALy) + (generator_->partVz[partonN] - parton1HCALz)*(generator_->partVz[partonN] - parton1HCALz))) / lightSpeed; // in ns
              double TOF_bQuark2 = 1000000000*(TMath::Sqrt((generator_->partVx[partonN-1] - parton2HCALx)*(generator_->partVx[partonN-1] - parton2HCALx) + (generator_->partVy[partonN-1] - parton2HCALy)*(generator_->partVy[partonN-1] -parton2HCALy) + (generator_->partVz[partonN-1] - parton2HCALz)*(generator_->partVz[partonN-1] - parton2HCALz))) / lightSpeed; // in ns
              double TOF_expected1 = 1000000000*(TMath::Sqrt(parton1HCALx*parton1HCALx + parton1HCALy*parton1HCALy + parton1HCALz*parton1HCALz)) / lightSpeed; // this is x,y,z of parton intersection with HCAL
	      double TOF_expected2 = 1000000000*(TMath::Sqrt(parton2HCALx*parton2HCALx + parton2HCALy*parton2HCALy + parton2HCALz*parton2HCALz)) / lightSpeed; // this is x,y,z of parton intersection with HCAL
	      TOF_LLP_quark->Fill((TOF_LLP + TOF_bQuark1 - TOF_expected1)); // TOF in ns, all distances have been in cm
	      TOF_LLP_quark->Fill((TOF_LLP + TOF_bQuark2 - TOF_expected2));
	      TOF_expected->Fill(TOF_expected1);
              TOF_expected->Fill(TOF_expected2);

	      // consider HCAL TPs within DR 0.5 of the parton intersection of the HCAL for correlation of time delay
	      for (int HcalTPIt = 0; HcalTPIt < nCaloTPemu; HcalTPIt++){
		// eta and phi of the HCAL TP    
		double tpEtaemu = l1CaloTPemu_->hcalTPieta[HcalTPIt]; // ieta    
		if (abs(tpEtaemu) > 28 ) continue; // don't consider HCAL TPs outside of HB or HE 
		double tpPhiemu = l1CaloTPemu_->hcalTPCaliphi[HcalTPIt]; // iphi
		double TP_Eta = etaVal(tpEtaemu); // eta
		double TP_Phi = phiVal(tpPhiemu); // phi
		if (deltaR(partonEta1, partonPhi1, TP_Eta, TP_Phi) < 0.2 ) {
		  if ( l1CaloTPemu_->hcalTPtiming1[HcalTPIt] >=0 ) TOF_vs_TDC->Fill(TOF_LLP + TOF_bQuark1 - TOF_expected1,l1CaloTPemu_->hcalTPtiming1[HcalTPIt]);
                  if ( l1CaloTPemu_->hcalTPtiming2[HcalTPIt] >=0 )TOF_vs_TDC->Fill(TOF_LLP + TOF_bQuark1 - TOF_expected1,l1CaloTPemu_->hcalTPtiming2[HcalTPIt]);
                  if ( l1CaloTPemu_->hcalTPtiming3[HcalTPIt] >=0 )TOF_vs_TDC->Fill(TOF_LLP + TOF_bQuark1 - TOF_expected1,l1CaloTPemu_->hcalTPtiming3[HcalTPIt]);
                  if ( l1CaloTPemu_->hcalTPtiming4[HcalTPIt] >=0 )TOF_vs_TDC->Fill(TOF_LLP + TOF_bQuark1 - TOF_expected1,l1CaloTPemu_->hcalTPtiming4[HcalTPIt]);
                  if ( l1CaloTPemu_->hcalTPtiming5[HcalTPIt] >=0 )TOF_vs_TDC->Fill(TOF_LLP + TOF_bQuark1 - TOF_expected1,l1CaloTPemu_->hcalTPtiming5[HcalTPIt]);
                  if ( l1CaloTPemu_->hcalTPtiming6[HcalTPIt] >=0 )TOF_vs_TDC->Fill(TOF_LLP + TOF_bQuark1 - TOF_expected1,l1CaloTPemu_->hcalTPtiming6[HcalTPIt]);
                  if ( l1CaloTPemu_->hcalTPtiming7[HcalTPIt] >=0 )TOF_vs_TDC->Fill(TOF_LLP + TOF_bQuark1 - TOF_expected1,l1CaloTPemu_->hcalTPtiming7[HcalTPIt]);
		}
                if (deltaR(partonEta2, partonPhi2, TP_Eta, TP_Phi) < 0.2 ) {
                  if ( l1CaloTPemu_->hcalTPtiming1[HcalTPIt] >=0 )TOF_vs_TDC->Fill(TOF_LLP + TOF_bQuark2 - TOF_expected2,l1CaloTPemu_->hcalTPtiming1[HcalTPIt]);
                  if ( l1CaloTPemu_->hcalTPtiming2[HcalTPIt] >=0 )TOF_vs_TDC->Fill(TOF_LLP + TOF_bQuark2 - TOF_expected2,l1CaloTPemu_->hcalTPtiming2[HcalTPIt]);
                  if ( l1CaloTPemu_->hcalTPtiming3[HcalTPIt] >=0 )TOF_vs_TDC->Fill(TOF_LLP + TOF_bQuark2 - TOF_expected2,l1CaloTPemu_->hcalTPtiming3[HcalTPIt]);
                  if ( l1CaloTPemu_->hcalTPtiming4[HcalTPIt] >=0 )TOF_vs_TDC->Fill(TOF_LLP + TOF_bQuark2 - TOF_expected2,l1CaloTPemu_->hcalTPtiming4[HcalTPIt]);
                  if ( l1CaloTPemu_->hcalTPtiming5[HcalTPIt] >=0 )TOF_vs_TDC->Fill(TOF_LLP + TOF_bQuark2 - TOF_expected2,l1CaloTPemu_->hcalTPtiming5[HcalTPIt]);
                  if ( l1CaloTPemu_->hcalTPtiming6[HcalTPIt] >=0 )TOF_vs_TDC->Fill(TOF_LLP + TOF_bQuark2 - TOF_expected2,l1CaloTPemu_->hcalTPtiming6[HcalTPIt]);
                  if ( l1CaloTPemu_->hcalTPtiming7[HcalTPIt] >=0 )TOF_vs_TDC->Fill(TOF_LLP + TOF_bQuark2 - TOF_expected2,l1CaloTPemu_->hcalTPtiming7[HcalTPIt]);
		}
	      }
	    }

	    double min_DeltaR = 100;    
	    int closestJet = -1;
            for (uint jetIt = 0; jetIt < nJetemu; jetIt++){ // loop over L1 jets, and only do first four (4 highest energy L1 jets from 4 leptons)
              if (l1emu_->jetEt[jetIt] < 20 ) continue; // require jet is greater than 20 GeV to attempt matching to parton
              double Jet_eta = l1emu_->jetEta[jetIt];                 
              double Jet_phi = l1emu_->jetPhi[jetIt];                      
              double DeltaR = deltaR(Jet_eta,Jet_phi,partonEta,partonPhi); // distance between L1 jet and a parton
              if (DeltaR < min_DeltaR) { 
                min_DeltaR = DeltaR; // find min delta R between L1 Jet and parton -- this is reset for each parton
                closestJet = jetIt; // record which L1 jet is the closest to the parton
              }
            } // closing L1 jet loop   
	    if (min_DeltaR < 100) {
	      DeltaR_parton_L1jet->Fill(min_DeltaR); // filled per parton
	      DeltaR_parton_L1jet_closest[closestJet]->Fill(min_DeltaR); // fill just for which jet is closest
	    }
	    //	    std::cout << "energy of closest jet = " << l1emu_->jetEt[closestJet] << std::endl;  
	    //	    std::cout << "the closest jet is = " << closestJet << " with a delta R to the LLP parton of " << min_DeltaR << std::endl;
          } // LLP PDG ID
        } // parton PDG ID
      } // closing parton loop   
      int nJet = 0;
      for (uint jetIt = 0; jetIt < nJetemu; jetIt++) if (l1emu_->jetEt[jetIt] > 20 )  nJet += 1;
      nJet_pt20GeV->Fill(nJet);
      // ********************** end of gen particle matching ********************

      // filling histograms for number of cells above energy and timing threshold set by the timing bit (50 ADC, 3 ns)
      ADC50_3ns_4JetMultHB_emu->Fill(Sum4Jet_HB);
      ADC50_3ns_4JetMultHE_emu->Fill(Sum4Jet_HE);
      ADC50_3ns_4JetMultHBHE_emu->Fill(Sum4Jet_HBHE);

      // saving number of events passed cell multiplicity cuts. These are efficiency values used in rate vs eff plots
      totalEvents += 1; // counting total events for efficiency calculations
      htSumDistribution->Fill(htSum);

      if ( ((htSum > 120) && ( Sum4Jet_HBHE >= 10 )) || (htSum > 360) ) passed4JetMult_HBHE_ht120 += 1;
      if (htSum > 360) passedHtSum360 += 1;

      for(int bin=0; bin<nHtSumBins; bin++){
        if( (htSum) >= htSumLo+(bin*htSumBinWidth) && (Sum4Jet_HBHE)>=10 ) htSumRates_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV, compare to original rates in plot
        if( (htSum) >= htSumLo+(bin*htSumBinWidth) ) htSumRates_original_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV, use for ht > 360 original rates in rate vs. eff plots
        if( (htSum) >= htSumLo+(bin*htSumBinWidth) && ( (Sum4Jet_HBHE)>=10 || htSum > 360 ) ) htSumRates_120timingOR360_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV, use for ht > 120 + timing OR ht > 360 rates in rate vs. eff plots
      }
    }// closes if 'emuOn' is true


    //do routine for L1 hardware quantities
    if (hwOn){

      treeL1TPhw->GetEntry(jentry);
      double tpEt(0.);
      
      for(int i=0; i < l1TPhw_->nHCALTP; i++){
	tpEt = l1TPhw_->hcalTPet[i];
	hcalTP_hw->Fill(tpEt);
      }
      for(int i=0; i < l1TPhw_->nECALTP; i++){
	tpEt = l1TPhw_->ecalTPet[i];
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
      }

      for(int bin=0; bin<nMhtSumBins; bin++){
        if( (mhtSum) >= mhtSumLo+(bin*mhtSumBinWidth) ) mhtSumRates_hw->Fill(mhtSumLo+(bin*mhtSumBinWidth)); //GeV           
      }

      for(int bin=0; bin<nEtSumBins; bin++){
        if( (etSum) >= etSumLo+(bin*etSumBinWidth) ) etSumRates_hw->Fill(etSumLo+(bin*etSumBinWidth)); //GeV           
      }

      for(int bin=0; bin<nMetSumBins; bin++){
        if( (metSum) >= metSumLo+(bin*metSumBinWidth) ) metSumRates_hw->Fill(metSumLo+(bin*metSumBinWidth)); //GeV           
      } 
      for(int bin=0; bin<nMetHFSumBins; bin++){
        if( (metHFSum) >= metHFSumLo+(bin*metHFSumBinWidth) ) metHFSumRates_hw->Fill(metHFSumLo+(bin*metHFSumBinWidth)); //GeV           
      } 

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

    // write histograms based on timing bit
    ADC50_3ns_4JetMultHB_emu->Write();
    ADC50_3ns_4JetMultHE_emu->Write();
    ADC50_3ns_4JetMultHBHE_emu->Write();
    // gen matched plots
    nJet_pt20GeV->Write();
    DeltaR_parton_L1jet->Write();
    for (int closestJet = 0; closestJet < 15; closestJet++) DeltaR_parton_L1jet_closest[closestJet]->Write();
    LLPdecayRadius->Write();
    LLPdecayXyz->Write();
    LLPdecayRadiusTrigAcceptance->Write();
    LLPdecayXyzTrigAcceptance->Write();
    TOF_LLP_quark->Write();
    TOF_expected->Write();
    TOF_vs_TDC->Write();

    htSumDistribution->Write();

    htSumRates_original_emu->Scale(norm);
    htSumRates_120timingOR360_emu->Scale(norm);

    // write histograms for energy / ADC
    for (int iEta = 1; iEta <= 16; iEta ++) GeV_ADC_ratio_HB[iEta]->Write();
    for (int iEta = 17; iEta <= 28; iEta ++) GeV_ADC_ratio_HE[iEta]->Write();
  }
  ADC_vs_GeV_percell->Write();
  ADC_vs_GeV_percell_HBHE20->Write();
  ADC_vs_GeV_percell_HE21->Write();
  GeV_ADC_byTPET->Write();
  GeV_ADC_byTPET_HBHE20->Write();
  GeV_ADC_byTPET_HE21->Write();


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

  std::cout << passed4JetMult_HBHE_ht120 << std::endl;
  std::cout << passedHtSum360 << std::endl;
  std::cout << totalEvents << std::endl;

  // saving efficiencies and rates in txt files to be read by rate vs eff plotting macros
  // signal efficiencies
  if (inputFile.substr(27,2) == "mh" ) {
    std::ofstream MultiplicityHits50ADC3ns_ht120_Signal;
    MultiplicityHits50ADC3ns_ht120_Signal.open(Form("MultiplicityHits50ADC3ns_ht120_Signal_%s.txt", inputFile.substr(27,20).c_str()),std::ios_base::trunc);
    MultiplicityHits50ADC3ns_ht120_Signal << passed4JetMult_HBHE_ht120 / totalEvents << std::endl; // efficiency at HT 120+timing OR HT 360
    MultiplicityHits50ADC3ns_ht120_Signal << passedHtSum360 / totalEvents << std::endl; // efficiency at HT 360
    MultiplicityHits50ADC3ns_ht120_Signal.close();
  }
  // background efficiencies 
  if (inputFile.substr(27,3) == "QCD" ) {
    std::ofstream MultiplicityHits50ADC3ns_ht120_Background;
    MultiplicityHits50ADC3ns_ht120_Background.open("MultiplicityHits50ADC3ns_ht120_Background.txt");
    MultiplicityHits50ADC3ns_ht120_Background << passed4JetMult_HBHE_ht120 / totalEvents << std::endl; // efficiency at HT 120+timing OR HT 360  
    MultiplicityHits50ADC3ns_ht120_Background << passedHtSum360 / totalEvents << std::endl; // efficiency at HT 360  
    MultiplicityHits50ADC3ns_ht120_Background.close();
  }
  // neutrino gun rates
  if (inputFile.substr(27,8) == "Neutrino" ) {
    int htSum_120timingOR360_120 = htSumRates_120timingOR360_emu->GetBinContent(htSumRates_120timingOR360_emu->GetXaxis()->FindBin(120));
    int htSum_original_360 = htSumRates_original_emu->GetBinContent(htSumRates_original_emu->GetXaxis()->FindBin(360));
    std::ofstream NuGunRates;
    NuGunRates.open("NuGunRates.txt");
    NuGunRates << htSum_120timingOR360_120 << std::endl; // rate at HT 120 
    NuGunRates << htSum_original_360 << std::endl; // rate at HT 360 without timing cuts
    NuGunRates.close();
  }


  myfile << "using the following ntuple: " << inputFile << std::endl;
  myfile << "number of colliding bunches = " << numBunch << std::endl;
  myfile << "run luminosity = " << runLum << std::endl;
  myfile << "expected luminosity = " << expectedLum << std::endl;
  myfile << "norm factor used = " << norm << std::endl;
  myfile << "number of good events = " << goodLumiEventCount << std::endl;
  myfile.close(); 
}//closes the function 'rates'
