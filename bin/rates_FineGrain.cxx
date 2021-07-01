// 2020, 2021: edited by Gillian Kopp for HCAL L1 LLP trigger studies, based on using a delayed jet defined by a TP with delayed cells
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
#include <sstream>
#include <algorithm>
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

void rates_LLPflag(bool newConditions, const std::string& inputFileDirectory, double jetEnergy, double timing_depth_OR_flag, double delta_R_req, double triggerability);

int main(int argc, char *argv[])
{
  bool newConditions = true;
  std::string ntuplePath("");
  double jetEt;
  double timing_depth_OR;
  double delta_R;
  double triggerability;

  if (argc != 7) {
    std::cout << "Usage: rates.exe [new/def] [path to ntuples]\n"
	      << "[new/def] indicates new or default (existing) conditions"
	      << "then jet energy (GeV), timing (1) depth (2) OR (3), deltaR, triggerability yes (1) no (0)" << std::endl;
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
    jetEt = atof(argv[3]);
    timing_depth_OR = atof(argv[4]);
    delta_R = atof(argv[5]);
    triggerability = atof(argv[6]);
  }

  rates_LLPflag(newConditions, ntuplePath, jetEt, timing_depth_OR, delta_R, triggerability);

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

std::vector<double> eta_phi_0index(double ieta, double iphi) {
  double TP_eta_0index = 1000;
  if (ieta>0) TP_eta_0index = ieta+27; // 28 to 55
  else if (ieta<0) TP_eta_0index = ieta+28; // 0 to 27                                                                                                              
  double TP_phi_0index = iphi-1; // 0 to 71
  std::vector<double> eta_phi_0index;
  eta_phi_0index.push_back(TP_eta_0index);
  eta_phi_0index.push_back(TP_phi_0index);
  return eta_phi_0index; 
}

std::vector<double> normal_eta_phi(double eta_0index, double phi_0index) {
  double ieta = 1000;
  ieta = eta_0index - 27;
  if (ieta <= 0) ieta -= 1;
  double iphi = phi_0index + 1;
  std::vector<double> normal_eta_phi;
  normal_eta_phi.push_back(ieta);
  normal_eta_phi.push_back(iphi);
  return normal_eta_phi;
}

std::vector<double> intersect(double vx, double vy,double vz, double px, double py, double pz) {
  double lightSpeed = 29979245800; // speed of light in cm/s
  double radius = 295; // outer HCAL radius (cm)
  double length = 568; // end of HCAL endcap (cm)
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
	double distance = deltaR(l1emu_->jetEta[L1Jet], l1emu_->jetPhi[L1Jet], partonEta, partonPhi);
	if (distance < min_dR ) {
	  min_dR = distance;
	  partonNum = (double)partonN;
	}
      }
    }
  }
  std::vector<double> DR_partonNum;
  DR_partonNum.push_back(min_dR);
  DR_partonNum.push_back(partonNum);
  return DR_partonNum; // DeltaR from given L1 jet to parton, which number this parton is
}

// given a parton number for an event (this is parton closest to L1 jet from "closestParton" function, what is the displacement of the LLP this resulted from?                    
std::vector<double> LLPdecayInfo(int partonN, L1Analysis::L1AnalysisGeneratorDataFormat *generator_) {
  double LLPxDecay = generator_->partVx[partonN]; // this is creation vertex of parton = decay vertex of LLP
  double LLPyDecay = generator_->partVy[partonN];
  double LLPzDecay = generator_->partVz[partonN];
  double genLLPBeta = -1000;
  // LLP decay results in two b quarks, calculate beta and gamma of LLP based on the quantities of b quarks
  if ( (generator_->partVz[partonN] == generator_->partVz[partonN-1] ) && (abs(generator_->partId[partonN]) == abs(generator_->partId[partonN-1])) ) {
    genLLPBeta = sqrt((generator_->partPx[partonN-1] + generator_->partPx[partonN])*(generator_->partPx[partonN-1] + generator_->partPx[partonN]) + (generator_->partPy[partonN-1] + generator_->partPy[partonN])*(generator_->partPy[partonN-1] + generator_->partPy[partonN]) + (generator_->partPz[partonN-1] + generator_->partPz[partonN])*(generator_->partPz[partonN-1] + generator_->partPz[partonN])) / (generator_->partE[partonN-1] + generator_->partE[partonN]); // beta of LLP considering momentum and energy of two b quarks
  }
  if ( (generator_->partVz[partonN] == generator_->partVz[partonN+1] ) && (abs(generator_->partId[partonN]) == abs(generator_->partId[partonN+1])) ) {
    genLLPBeta = sqrt((generator_->partPx[partonN+1] + generator_->partPx[partonN])*(generator_->partPx[partonN+1] + generator_->partPx[partonN]) + (generator_->partPy[partonN+1] + generator_->partPy[partonN])*(generator_->partPy[partonN+1] + generator_->partPy[partonN]) + (generator_->partPz[partonN+1] + generator_->partPz[partonN])*(generator_->partPz[partonN+1] + generator_->partPz[partonN])) / (generator_->partE[partonN+1] + generator_->partE[partonN]); // beta of LLP considering momentum and energy of two b quarks
  }
  double genLLPGamma = 1./TMath::Sqrt(1.-genLLPBeta*genLLPBeta);
  double vertex = sqrt(LLPxDecay*LLPxDecay + LLPyDecay*LLPyDecay + LLPzDecay*LLPzDecay);
  // TOF calculation
  double lightSpeed = 29979245800; // in cm / s
  double TOF_LLP = 1000000000*vertex / (genLLPBeta * lightSpeed); // in ns
  double partonHCALx = intersect(LLPxDecay,LLPyDecay,LLPzDecay, generator_->partPx[partonN],generator_->partPy[partonN],generator_->partPz[partonN])[2]; // in cm
  double partonHCALy = intersect(LLPxDecay,LLPyDecay,LLPzDecay, generator_->partPx[partonN],generator_->partPy[partonN],generator_->partPz[partonN])[3];
  double partonHCALz = intersect(LLPxDecay,LLPyDecay,LLPzDecay, generator_->partPx[partonN],generator_->partPy[partonN],generator_->partPz[partonN])[4];
  double TOF_bQuark = 1000000000*(TMath::Sqrt((LLPxDecay - partonHCALx)*(LLPxDecay - partonHCALx) + (LLPyDecay - partonHCALy)*(LLPyDecay -partonHCALy) + (LLPzDecay - partonHCALz)*(LLPzDecay - partonHCALz))) / lightSpeed; // in ns
  double TOF_expected = 1000000000*(TMath::Sqrt(partonHCALx*partonHCALx + partonHCALy*partonHCALy + partonHCALz*partonHCALz)) / lightSpeed; // this is x,y,z of parton intersection with HCAL
  double TOFdelay = TOF_LLP + TOF_bQuark - TOF_expected;
  std::vector<double> LLPinfo;
  LLPinfo.push_back(LLPxDecay);
  LLPinfo.push_back(LLPyDecay);
  LLPinfo.push_back(LLPzDecay);
  LLPinfo.push_back(vertex / (genLLPGamma * genLLPBeta) );
  LLPinfo.push_back(vertex);
  LLPinfo.push_back(TOFdelay);
  LLPinfo.push_back(genLLPBeta);
  return LLPinfo; // LLP x, y, z decay, LLP ctau in cm in LLP rest frame, TOF delay in ns, LLP beta = v/c
}


void rates_LLPflag(bool newConditions, const std::string& inputFileDirectory, double jetEnergy, double timing_depth_OR_flag, double delta_R_req, double triggerability){
  
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

  std::string axR = ";L1 H_{T} (GeV);rate (Hz)";
  std::string axD = ";L1 H_{T} (GeV);events/bin";

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
  TH1F * DeltaR_L1_delayed_hit_emu = new TH1F("DeltaR_L1_delayed_hit_emu","DeltaR between a Delayed Hit and the Closest L1 jet;Delta R;Fraction of Entries (normalized)",20,0,1);
  TH1F * DeltaR_L1_prompt_hit_emu = new TH1F("DeltaR_L1_prompt_hit_emu","DeltaR between a Prompt Hit and the Closest L1 jet;Delta R;Fraction of Entries (normalized)",20,0,1);
  TH1F * DeltaR_L1_delayed_TP_emu = new TH1F("DeltaR_L1_delayed_TP_emu","DeltaR between a Delayed TT and the Closest L1 Jet;Delta R;Fraction of Entries (normalized)",30,0,1.5);

  TH1F * Mult_prompt_hit_emu =new TH1F("Mult_prompt_hit_emu","Number of prompt TPs near a jet with 2 delayed cells;Number of TPs;Fraction of Entries (normalized)",20,0,20);
  TH1F * Mult_delayed_hit_emu = new TH1F("Mult_delayed_hit_emu","Number of TTs with a delayed hit near a L1 jet;Number of TTs;Fraction of Entries (normalized)",10,0,10);
  TH1F * Mult_delayed_hit_jetET_emu = new TH1F("Mult_delayed_hit_jetET_emu","Number of TTs with a delayed hit near a L1 jet, jet ET>40 GeV;Number of TTs;Fraction of Entries (normalized)",10,0,10);
  TH1F * Mult_delayed_hit_promptV_emu = new TH1F("Mult_delayed_hit_promptV_emu","Number of TTs with a delayed hit near L1 jet, prompt veto applied;Number of TTs;Fraction of Entries (normalized)",10,0,10);
  TH1F * Mult_delayed_hit_jetETpromptV_emu = new TH1F("Mult_delayed_hit_jetETpromptV_emu","Number of TTs with a delayed hit near L1 jet, jet ET>40 GeV and prompt veto applied;Number of TTs;Fraction of Entries (normalized)",10,0,10);

  TH1F * HTdistribution_trig_emu = new TH1F("HTdistribution_trig_emu","HT Distribution of Events Passing Calo Cluster Trigger;HT (GeV);Number of Events",35,0,1200);
  TH1F * HTdistribution_emu = new TH1F("HTdistribution_emu","HT Distribution of Events;HT (GeV);Number of Events",35,0,1200);

  // efficiency by jet PT
  TH1F * JetPTdistribution_trig_emu = new TH1F("JetPTdistribution_trig_emu","Jet pT Distribution passing Delayed Jet Trigger;Jet pT (GeV);Number of Jets",15,0,200);
  TH1F * JetPTdistribution_trig120_emu = new TH1F("JetPTdistribution_trig120_emu","Jet pT Distribution passing Delayed Jet Trigger and HT120;Jet pT (GeV);Number of Jets",15,0,200);
  TH1F * JetPTdistribution_emu = new TH1F("JetPTdistribution_emu","Jet pT Distribution;Jet pT (GeV);Number of Jets",15,0,200);
  TH1F * DelayedSeedEnergy = new TH1F("DelayedSeedEnergy","Energy of Delayed Cells;Energy (GeV); Number of Cells",20,0,20);

  // HT sum rate distributions to use in rate vs eff plots. Need HT > 360 rate, and HT > 120 + timing cut rate
  TH1F* htSumRates_original_emu = new TH1F("htSumRates_original_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);

  TH1F * htSumDistribution = new TH1F("htSumDistribution","htSum Distribution;L1_HT (GeV);Number of Events",200,0,2000); // x bins, x low, x up

  TH1F* ctau_allLLP = new TH1F("ctau_allLLP","ctau of LLP in particle rest frame;ctau (m);Number of LLPs",100,0,50);
  TH1F* ctau_LLP = new TH1F("ctau_LLP","ctau of LLP in particle rest frame (triggerable, double counted);ctau (m);Number of LLPs",100,0,50);
  TH1F* ctau_LLP_trigger = new TH1F("ctau_LLP_trigger","ctau of LLP (passed trigger, double counted) in particle rest frame;ctau (m);Number of LLPs",100,0,50);
  TH1F* path_length = new TH1F("path_length","Path length of particle in lab frame (triggerable, double counted);path length (m);Number of LLPs",16,0,8); 
  TH1F* path_length_trigger = new TH1F("path_length_trigger","Path length of particle (passed trigger, double counted) in lab frame;path length (m);Number of LLPs",16,0,8);
  TH1F* path_length_120trigger = new TH1F("path_length_120trigger","Path length of particle (passed trigger, HT>120, double counted) in lab frame;path length (m);Number of LLPs",16,0,8);
  TH1F* path_length2 = new TH1F("path_length2","Path length of particle in lab frame (triggerable, double counted);path length (m);Number of LLPs",16,0,0.5);
  TH1F* path_length2_trigger = new TH1F("path_length2_trigger","Path length of particle (passed trigger, double counted) in lab frame;path length (m);Number of LLPs",16,0,0.5);
  TH1F* path_length2_120trigger = new TH1F("path_length2_120trigger","Path length of particle (passed trigger, HT>120, double counted) in lab frame;path length (m);Number of LLPs",16,0,0.5);
  TH1F* TOFdelay = new TH1F("TOFdelay","TOF delay = LLP TOF + b TOF - TOF expected (ns);TOF delay (ns); Number of Jets",30,0,15);
  TH1F* TOFdelay_trigger = new TH1F("TOFdelay_trigger","TOF delay = LLP TOF + b TOF - TOF expected (ns, triggered);TOF delay (ns); Number of Jets",30,0,15);
  TH1F* TOFdelay_120trigger = new TH1F("TOFdelay_120trigger","TOF delay = LLP TOF + b TOF - TOF expected (ns, triggered, HT>120);TOF delay (ns); Number of Jets",30,0,15);

  TH1F* beta_LLP = new TH1F("beta_LLP","Beta of LLP particle;Beta;Number of LLPs",100,0,1);

  // saving rate and efficiencies 
  double passed4JetMult_HBHE_ht120_1(0), passed4JetMult_HBHE_ht120_2(0);
  double passed_calo_cluster_trig(0);
  double passed_calo_cluster_trig_120(0), passed_calo_cluster_trig_120_2(0);
  double passedHtSum360(0);
  double totalEvents(0);

  /////////////////////////////////
  //////////// eta_depth_tdc //////
  /////////////////////////////////
  int eta_depth_tdc90[30][8] = {{4}};
  //  std::ifstream file("TDCdistribution_Background90_delayed_QCD_2GeV.txt");
  std::ifstream file("TDCdistribution_Background90_QCD_smoothed_2GeV.txt");
  for(int row = 0; row < 30; ++row) {
    std::string line;
    std::getline(file, line);
    if ( !file.good() )
      break;
    std::stringstream iss(line);
    for (int col = 0; col <= 8; ++col) { 
      std::string val;
      std::getline(iss, val, ',');
      if ( !iss.good() )
	break;
      std::stringstream convertor(val);
      convertor >> eta_depth_tdc90[row+1][col]; // now col = 1 corresponds to depth = 1
    }
  }

  /////////////////////////////////
  // loop through all the entries//
  /////////////////////////////////
  std::ofstream DelayedSeed_event_ieta_iphi_depth;
  DelayedSeed_event_ieta_iphi_depth.open(Form("DelayedSeed_event_ieta_iphi_depth_%s.txt", inputFile.substr(0,7).c_str()),std::ios_base::trunc);

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
        if( (jetEt_2) >= jetLo + (bin*jetBinWidth) ) doubleJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
        if( (jetEt_3) >= jetLo + (bin*jetBinWidth) ) tripleJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
        if( (jetEt_4) >= jetLo + (bin*jetBinWidth) ) quadJetRates_emu->Fill(jetLo+(bin*jetBinWidth));  //GeV
      }  
             
      for(int bin=0; bin<nEgBins; bin++){
        if( (egEt_1) >= egLo + (bin*egBinWidth) ) singleEgRates_emu->Fill(egLo+(bin*egBinWidth));  //GeV
        if( (egEt_2) >= egLo + (bin*egBinWidth) ) doubleEgRates_emu->Fill(egLo+(bin*egBinWidth));  //GeV
      }  

      for(int bin=0; bin<nTauBins; bin++){
        if( (tauEt_1) >= tauLo + (bin*tauBinWidth) ) singleTauRates_emu->Fill(tauLo+(bin*tauBinWidth));  //GeV
        if( (tauEt_2) >= tauLo + (bin*tauBinWidth) ) doubleTauRates_emu->Fill(tauLo+(bin*tauBinWidth));  //GeV
      } 

      for(int bin=0; bin<nEgBins; bin++){
        if( (egISOEt_1) >= egLo + (bin*egBinWidth) ) singleISOEgRates_emu->Fill(egLo+(bin*egBinWidth));  //GeV
        if( (egISOEt_2) >= egLo + (bin*egBinWidth) ) doubleISOEgRates_emu->Fill(egLo+(bin*egBinWidth));  //GeV
      }  

      for(int bin=0; bin<nTauBins; bin++){
        if( (tauISOEt_1) >= tauLo + (bin*tauBinWidth) ) singleISOTauRates_emu->Fill(tauLo+(bin*tauBinWidth));  //GeV
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
      /////////// LLP incident on HB? //////
      //////////////////////////////////////
      for (int partonN = 0; partonN < generator_->nPart; partonN ++) {
        if ( (generator_->partVz[partonN] == generator_->partVz[partonN-1] ) && (abs(generator_->partId[partonN]) == abs(generator_->partId[partonN-1])) ) {
          if (generator_->partParent[partonN] == 6000113 ) ctau_allLLP->Fill(LLPdecayInfo(partonN,generator_)[3]/100);
        }
      }

      uint nJetemu = l1emu_->nJets; // number of jets per event
      double nCaloTPemu = l1CaloTPemu_->nHCALTP; // number of TPs varies from 400-1400 per event, approximately Gaussian                                                  

      // triggerability restrictions
      int triggerableJets[nJetemu] = {0};
      int numLLPdecayHB = 0; // count how many LLP decay products are expected to intersect the HB 

      for (uint jetIt = 0; jetIt < nJetemu; jetIt++) { // loop over jets
	if (abs(l1emu_->jetEta[jetIt]) > 2.5) continue; // consider HB+HE jets, HB extends to 1.4. HE extends to 3. Use values of 1, 2.5
	if (inputFile.substr(0,3) == "QCD" || inputFile.substr(0,13) == "TimingBit/QCD" || inputFile.substr(0,18) == "TimingDepthBit/QCD" ) JetPTdistribution_emu->Fill(l1emu_->jetEt[jetIt]);
	double deltaR_parton_jet = closestParton(jetIt, l1emu_, generator_)[0];
	int partonN = closestParton(jetIt, l1emu_, generator_)[1];
	if (deltaR_parton_jet <= 0.5) { // if closest parton is near a HB L1 jet
	  if (triggerability == 1) {
	    numLLPdecayHB += 1; // how many of the partons expected to intersect HB
	    triggerableJets[jetIt] = 1; // 1 if triggerable jet, 0 otherwise
	  }
	  JetPTdistribution_emu->Fill(l1emu_->jetEt[jetIt]);
	  triggerableJets[jetIt] = 1; // 1 if triggerable jet, 0 otherwise
	  ctau_LLP->Fill(LLPdecayInfo(partonN,generator_)[3]/100);
	  path_length->Fill(LLPdecayInfo(partonN,generator_)[4]/100);
	  path_length2->Fill(LLPdecayInfo(partonN,generator_)[4]/100); // first half meter
	  TOFdelay->Fill(LLPdecayInfo(partonN,generator_)[5]); // can be twice per LLP
	}
      } // end jet loop, have determined if a jet is triggerable now
      
      if (triggerability == 0) {
	for (uint jetIt = 0; jetIt < nJetemu; jetIt++) triggerableJets[jetIt] = 1;
	numLLPdecayHB = 1;
      }

      if ( inputFile.substr(0,14) == "Time3Depth1/mh" && numLLPdecayHB == 0 ) continue; // if no LLPs in HB, skip event

      //////////////////////////////////////
      ////////// HCAL TP Loop //////////
      //////////////////////////////////////
      // per TP, count delayed hits (timing bit) and tower energy for prompt hit consideration. Array is 56, 72 (32, 72 if just HB)

      int LLP_flagged_TTs[nJetemu] = {0}; // LLP flagged TTs in each L1 jet (has delayed cell, no prompt cell)

      for (int HcalTPIt = 0; HcalTPIt < nCaloTPemu; HcalTPIt++){ // loop over HCAL TPs
	double TP_ieta = l1CaloTPemu_->hcalTPieta[HcalTPIt]; // ieta
	double TP_iphi = l1CaloTPemu_->hcalTPCaliphi[HcalTPIt]; // iphi

	if (abs(TP_ieta) > 29 ) continue; 

	int TimingFlag = 0;
	int DepthFlag = 0;

        if (inputFile.substr(0,11) == "Time3Depth1") {
          int Timing3_MIP2_Depth1 = l1CaloTPemu_->hcalTPfineGrain[HcalTPIt]; // 10 01 prompt MIP MIP depth
	  if (timing_depth_OR_flag != 2) {
	    int NveryDelayed = (Timing3_MIP2_Depth1 & 0b100000) / 32;
	    int Ndelayed = (Timing3_MIP2_Depth1 & 0b010000) / 16;
	    int PromptVeto = (Timing3_MIP2_Depth1 & 0b001000) / 8;
	    int TotalDelayed = NveryDelayed + Ndelayed;
	    if (TotalDelayed > 0 && PromptVeto == 0) TimingFlag = 1;
	  }
	  if (timing_depth_OR_flag != 1) {
	    if (Timing3_MIP2_Depth1 & 0b000001) DepthFlag = 1;
	  }
	}

	// find closest jet to the TP
	int closestJet = -1;
	double min_DR = 100;
	double TP_eta = etaVal(TP_ieta);
	double TP_phi = phiVal(TP_iphi);
	for (uint jetIt = 0; jetIt < nJetemu; jetIt++) {
	  if (deltaR(l1emu_->jetEta[jetIt], l1emu_->jetPhi[jetIt], TP_eta, TP_phi) < min_DR) {
	    min_DR = deltaR(l1emu_->jetEta[jetIt], l1emu_->jetPhi[jetIt], TP_eta, TP_phi);
	    closestJet = jetIt;
	  }
	}
	if (min_DR<=delta_R_req) {
          if ((timing_depth_OR_flag == 3 && (TimingFlag == 1 || DepthFlag == 1)) || (timing_depth_OR_flag == 1 && TimingFlag == 1) || (timing_depth_OR_flag == 2 && DepthFlag == 1)) { // LLP flagged TT
	    LLP_flagged_TTs[closestJet] += 1;
	    DeltaR_L1_delayed_TP_emu->Fill(min_DR);
	  }
	}
      } // closing HCAL TP loop

      for (uint jetIt = 0; jetIt < nJetemu; jetIt++) {
	if ((inputFile.substr(0,14) == "Time3Depth1/mh") && (triggerableJets[jetIt] == 0)) {
	  LLP_flagged_TTs[jetIt] = 0; // set to 0 if jet is not triggerable
	  continue;
	}
	if (l1emu_->jetEt[jetIt] < jetEnergy ) { // require jet ET > 40 GeV for jet to be delayed
	  LLP_flagged_TTs[jetIt] = 0; // set to 0 if jet is too low energy
	  continue; 
	}
	Mult_delayed_hit_jetETpromptV_emu->Fill(LLP_flagged_TTs[jetIt]);

	// for jet PT efficiencies
	if (LLP_flagged_TTs[jetIt] >= 2) {
	  JetPTdistribution_trig_emu->Fill(l1emu_->jetEt[jetIt]);
	  if (htSum > 120) JetPTdistribution_trig120_emu->Fill(l1emu_->jetEt[jetIt]);

	  // generator quantities for LLP in a jet that was triggered on
	  if (inputFile.substr(0,23) != "Time3Depth1/RelValNuGun") {
	    double partonN = closestParton(jetIt, l1emu_, generator_)[1];
	    ctau_LLP_trigger->Fill(LLPdecayInfo(partonN,generator_)[3]/100);
	    path_length_trigger->Fill(LLPdecayInfo(partonN,generator_)[4]/100);
	    path_length2_trigger->Fill(LLPdecayInfo(partonN,generator_)[4]/100);
	    TOFdelay_trigger->Fill(LLPdecayInfo(partonN,generator_)[5]);
	    //	    if (LLPdecayInfo(partonN,generator_)[5] >= 11 && LLPdecayInfo(partonN,generator_)[5] < 12) std::cout << "TOFdelay_trigger = " << LLPdecayInfo(partonN,generator_)[5] << " for jet = " << jetIt << " at entry = " << jentry << std::endl;
	    if (htSum > 120) path_length_120trigger->Fill(LLPdecayInfo(partonN,generator_)[4]/100);
            if (htSum > 120) path_length2_120trigger->Fill(LLPdecayInfo(partonN,generator_)[4]/100);
	    if (htSum > 120) TOFdelay_120trigger->Fill(LLPdecayInfo(partonN,generator_)[5]);
	    beta_LLP->Fill(LLPdecayInfo(partonN,generator_)[6]);
	  }
	}
      } // closing jet loop

      // how many jets pass the delayed jet trigger
      int num_delayed_jet = -1;
      std::sort(LLP_flagged_TTs, LLP_flagged_TTs+nJetemu, std::greater<int>());
      if (LLP_flagged_TTs[0] == 1) num_delayed_jet = 0;
      if (LLP_flagged_TTs[0] >= 2) num_delayed_jet = 1;
      if (LLP_flagged_TTs[0] >= 2 && LLP_flagged_TTs[1] > 0) num_delayed_jet = 2;
      if (LLP_flagged_TTs[0] >= 2 && LLP_flagged_TTs[1] > 0  && LLP_flagged_TTs[2] > 0) num_delayed_jet = 3;
      
      // saving number of events passed cell multiplicity cuts. These are efficiency values used in rate vs eff plots
      totalEvents += 1; // counting total events for efficiency calculations
      htSumDistribution->Fill(htSum);
      HTdistribution_emu->Fill(htSum);
      if (num_delayed_jet >= 1) HTdistribution_trig_emu->Fill(htSum); // plot HT dist of events passing calo trigger

      if ( num_delayed_jet >= 0 ) passed_calo_cluster_trig += 1;
      if ( num_delayed_jet >= 0 && htSum > 120 ) passed_calo_cluster_trig_120 += 1;
      if ( num_delayed_jet >= 1 && htSum > 120 ) passed_calo_cluster_trig_120_2 += 1;
      if ( ((htSum > 120) && ( num_delayed_jet >= 0 )) || (htSum >= 360) ) passed4JetMult_HBHE_ht120_1 += 1; // +7
      if ( ((htSum > 120) && ( num_delayed_jet >= 1 )) || (htSum >= 360) ) passed4JetMult_HBHE_ht120_2 += 1;

      if (htSum > 360) passedHtSum360 += 1;

      for(int bin=0; bin<nHtSumBins; bin++){
        if( (htSum) >= htSumLo+(bin*htSumBinWidth) && ( num_delayed_jet)>=1 && htSum < 360 ) htSumRates_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV, compare to original rates in plot
        if( (htSum) >= htSumLo+(bin*htSumBinWidth) ) htSumRates_original_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV, use for ht > 360 original rates in rate vs. eff plots
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
        if( (jetEt_2) >= jetLo + (bin*jetBinWidth) ) doubleJetRates_hw->Fill(jetLo+(bin*jetBinWidth));  //GeV
        if( (jetEt_3) >= jetLo + (bin*jetBinWidth) ) tripleJetRates_hw->Fill(jetLo+(bin*jetBinWidth));  //GeV
        if( (jetEt_4) >= jetLo + (bin*jetBinWidth) ) quadJetRates_hw->Fill(jetLo+(bin*jetBinWidth));  //GeV
      }               
      for(int bin=0; bin<nEgBins; bin++){
        if( (egEt_1) >= egLo + (bin*egBinWidth) ) singleEgRates_hw->Fill(egLo+(bin*egBinWidth));  //GeV
        if( (egEt_2) >= egLo + (bin*egBinWidth) ) doubleEgRates_hw->Fill(egLo+(bin*egBinWidth));  //GeV
      }  
      for(int bin=0; bin<nTauBins; bin++){
        if( (tauEt_1) >= tauLo + (bin*tauBinWidth) ) singleTauRates_hw->Fill(tauLo+(bin*tauBinWidth));  //GeV
        if( (tauEt_2) >= tauLo + (bin*tauBinWidth) ) doubleTauRates_hw->Fill(tauLo+(bin*tauBinWidth));  //GeV
      } 
      for(int bin=0; bin<nEgBins; bin++){
        if( (egISOEt_1) >= egLo + (bin*egBinWidth) ) singleISOEgRates_hw->Fill(egLo+(bin*egBinWidth));  //GeV
        if( (egISOEt_2) >= egLo + (bin*egBinWidth) ) doubleISOEgRates_hw->Fill(egLo+(bin*egBinWidth));  //GeV
      }  
      for(int bin=0; bin<nTauBins; bin++){
        if( (tauISOEt_1) >= tauLo + (bin*tauBinWidth) ) singleISOTauRates_hw->Fill(tauLo+(bin*tauBinWidth));  //GeV
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

  DelayedSeed_event_ieta_iphi_depth.close();


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
    HTdistribution_trig_emu->Write();
    HTdistribution_emu->Write();
    htSumDistribution->Write();

    JetPTdistribution_emu->Write();
    JetPTdistribution_trig_emu->Write();
    JetPTdistribution_trig120_emu->Write();
    DelayedSeedEnergy->Write();

    DeltaR_L1_delayed_hit_emu->Write();
    DeltaR_L1_prompt_hit_emu->Write();
    DeltaR_L1_delayed_TP_emu->Write();

    Mult_prompt_hit_emu->Write();
    Mult_delayed_hit_emu->Write();
    Mult_delayed_hit_jetET_emu->Write();
    Mult_delayed_hit_promptV_emu->Write();
    Mult_delayed_hit_jetETpromptV_emu->Write();

    ctau_allLLP->Write();
    ctau_LLP->Write();
    ctau_LLP_trigger->Write();
    path_length->Write();
    path_length_trigger->Write();
    path_length_120trigger->Write();
    path_length2->Write();
    path_length2_trigger->Write();
    path_length2_120trigger->Write();
    TOFdelay->Write();
    TOFdelay_trigger->Write();
    TOFdelay_120trigger->Write();
    beta_LLP->Write();

    htSumRates_original_emu->Scale(norm);
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

  std::cout << inputFile.substr(0,14) << " triggerable events = " << totalEvents << "; Events with jet>40 and HT>120 that: passed delayed TT = " << passed_calo_cluster_trig_120 << "; passed 2 delayed TT = " << passed_calo_cluster_trig_120_2 << std::endl;
  std::cout << "relative efficiency = " << passed_calo_cluster_trig_120_2 / totalEvents << std::endl;
  if (inputFile.substr(0,14) != "Time3Depth1/mh") std::cout << "background rejection = " << totalEvents / passed_calo_cluster_trig_120_2 << std::endl;
  std::cout << passedHtSum360/totalEvents * 100 << " % passed HT360" << std::endl;
  std::cout << passed_calo_cluster_trig / totalEvents * 100 << " % passed calo trig (delayed TT), no HT cut / all events" << std::endl;
  std::cout << passed_calo_cluster_trig_120 / totalEvents * 100 << " % passed calo trig (delayed TT), HT 120 cut / all events" << std::endl;
  std::cout << passed_calo_cluster_trig_120_2 / totalEvents * 100 << " % passed calo trig (2 delayed TT), HT 120 cut / all events" << std::endl;
  std::cout << passed_calo_cluster_trig  << " passed calo trig (delayed TT) no HT cut" << std::endl;
  std::cout << passed4JetMult_HBHE_ht120_1 << " passed calo trig HT120 (delayed TT) OR HT360" << std::endl;
  std::cout << passed4JetMult_HBHE_ht120_2 << " passed calo trig HT120 (2 delayed TT) OR HT360" << std::endl;
  std::cout << (passed4JetMult_HBHE_ht120_1 - passedHtSum360)*100 / totalEvents << " added efficiency HT 360 (delayed TT)" << std::endl;
  std::cout << (passed4JetMult_HBHE_ht120_2 - passedHtSum360)*100 / totalEvents << " added efficiency HT 360 (2 delayed TT)" << std::endl;
  std::cout << (passed4JetMult_HBHE_ht120_2) / passedHtSum360 << " integrated luminosity gain" << std::endl;

  // neutrino gun rates
  if (inputFile.substr(0,23) == "Time3Depth1/RelValNuGun") {
  //  if (inputFile.substr(0,7) == "MinBias" ) {
    std::cout << "htSum_original120 = " << htSumRates_original_emu->GetBinContent(htSumRates_original_emu->GetXaxis()->FindBin(120)) << std::endl;
    std::cout << "htSum_original360 = " << htSumRates_original_emu->GetBinContent(htSumRates_original_emu->GetXaxis()->FindBin(360)) << std::endl;
    std::cout << "htSum_wtiming120 2 hits = " << htSumRates_emu->GetBinContent(htSumRates_emu->GetXaxis()->FindBin(120)) << std::endl;
    std::cout << "ratio change = timing at 120 / original 360 = " << ( htSumRates_emu->GetBinContent(htSumRates_emu->GetXaxis()->FindBin(120))) / (htSumRates_original_emu->GetBinContent(htSumRates_original_emu->GetXaxis()->FindBin(360))) << std::endl;
    std::cout << "ratio change = original 360 / timing at 120 = " << (htSumRates_original_emu->GetBinContent(htSumRates_original_emu->GetXaxis()->FindBin(360))) / ( htSumRates_emu->GetBinContent(htSumRates_emu->GetXaxis()->FindBin(120))) << std::endl;
    std::cout << "ratio change = original 120 / timing at 120 = " << (htSumRates_original_emu->GetBinContent(htSumRates_original_emu->GetXaxis()->FindBin(120))) / ( htSumRates_emu->GetBinContent(htSumRates_emu->GetXaxis()->FindBin(120))) << std::endl;
  }


  myfile << "using the following ntuple: " << inputFile << std::endl;
  myfile << "number of colliding bunches = " << numBunch << std::endl;
  myfile << "run luminosity = " << runLum << std::endl;
  myfile << "expected luminosity = " << expectedLum << std::endl;
  myfile << "norm factor used = " << norm << std::endl;
  myfile << "number of good events = " << goodLumiEventCount << std::endl;
  myfile.close(); 
}//closes the function 'rates'
