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

void rates_delayed_cluster(bool newConditions, const std::string& inputFileDirectory, double TDC_HB_variable, double TDC_HE_variable, double GeV_HB_variable, double GeV_HE_variable, double prompt_TP_energy_variable, double prompt_2x2_energy_variable);

int main(int argc, char *argv[])
{
  bool newConditions = true;
  std::string ntuplePath("");
  double TDC_HB;
  double TDC_HE;
  double GeV_HB;
  double GeV_HE;
  double prompt_TP_energy;
  double prompt_2x2_energy;

  if (argc != 9) {
    std::cout << "Usage: rates.exe [new/def] [path to ntuples]\n"
	      << "[new/def] indicates new or default (existing) conditions"
	      << "then values for TDC, GeV, prompt reject energies" << std::endl;
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
    TDC_HB = atof(argv[3]);
    TDC_HE = atof(argv[4]);
    GeV_HB = atof(argv[5]);
    GeV_HE = atof(argv[6]);
    prompt_TP_energy = atof(argv[7]);
    prompt_2x2_energy = atof(argv[8]);
  }

  rates_delayed_cluster(newConditions, ntuplePath, TDC_HB, TDC_HE, GeV_HB, GeV_HE, prompt_TP_energy, prompt_2x2_energy);

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

std::vector<double> eta_phi_2x2(double ieta, double iphi) {
  double TP_eta_2x2 = 1000;
  if (ieta>0) TP_eta_2x2 = ieta+27; //+15; // 16 to 31                                                                                                                    
  else if (ieta<0) TP_eta_2x2 = ieta+28; //+1+15; // 0 to 15                                                                                                              
  double TP_phi_2x2 = iphi-1; // 0 to 71
  std::vector<double> eta_phi_2x2;
  eta_phi_2x2.push_back(TP_eta_2x2);
  eta_phi_2x2.push_back(TP_phi_2x2);
  return eta_phi_2x2; 
}

std::vector<double> two_two_eta_phi(double eta_2x2, double phi_2x2) {
  double ieta = 1000;
  ieta = eta_2x2 -27;// - 15;
  if (ieta <= 0) ieta -= 1;
  double iphi = phi_2x2 + 1;
  std::vector<double> two_two_eta_phi;
  two_two_eta_phi.push_back(ieta);
  two_two_eta_phi.push_back(iphi);
  return two_two_eta_phi;
}

std::vector<double> intersect(double vx, double vy,double vz, double px, double py, double pz) {
  double lightSpeed = 29979245800; // speed of light in cm/s
  double radius = 295; //179; // 130 for calorimeters (ECAL + HCAL) in cm
  double length = 568; //388; // 300 for calorimeters (ECAL + HCAL) in cm
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
	//	if (l1emu_->jetEt[L1Jet] >= 20 ) {
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

std::vector<double> closestQCDParton(int L1Jet, L1Analysis::L1AnalysisL1UpgradeDataFormat *l1emu_, L1Analysis::L1AnalysisGeneratorDataFormat *generator_) {
  // find DR between L1 jet (argument) and which particle from QCD is nearby 
  double min_dR = 1000;
  double partonNum = -1;
  for (int partonN = 0; partonN < generator_->nPart; partonN ++) {
    double partonEta = 1000;
    double partonPhi = 1000;
    if (generator_->partHardProcess[partonN] == 0 ) continue;
    partonEta = intersect(generator_->partVx[partonN],generator_->partVy[partonN],generator_->partVz[partonN], generator_->partPx[partonN],generator_->partPy[partonN],generator_->partPz[partonN])[0];
    partonPhi = intersect(generator_->partVx[partonN],generator_->partVy[partonN],generator_->partVz[partonN], generator_->partPx[partonN],generator_->partPy[partonN],generator_->partPz[partonN])[1];
    double distance = deltaR(l1emu_->jetEta[L1Jet], l1emu_->jetPhi[L1Jet], partonEta, partonPhi);
    if (distance < min_dR ) {
      min_dR = distance;
      partonNum = (double)partonN;
    }
  }
  std::vector<double> DR_partonNum;
  DR_partonNum.push_back(min_dR);
  DR_partonNum.push_back(partonNum);
  return DR_partonNum; // DeltaR from given L1 jet to parton, which number this parton is  
}


// function returning the x,y,z positions and ctau of the LLP associated with the parton closest to this jet
// then make useful distributions based on this gen matching
std::vector<double> LLPdecayInfo(int partonN, L1Analysis::L1AnalysisL1UpgradeDataFormat *l1emu_, L1Analysis::L1AnalysisGeneratorDataFormat *generator_) {
  // given a parton number for an event (this is parton closest to L1 jet from "closestParton" function, what is the displacement of the LLP this resulted from?
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
  double vertex = sqrt(generator_->partVx[partonN]*generator_->partVx[partonN] + generator_->partVy[partonN]*generator_->partVy[partonN] + generator_->partVz[partonN]*generator_->partVz[partonN]);
  // TOF calculation
  double lightSpeed = 29979245800; // in cm / s
  double TOF_LLP = 1000000000*vertex / (genLLPBeta * lightSpeed); // in ns
  double partonHCALx = intersect(generator_->partVx[partonN],generator_->partVy[partonN],generator_->partVz[partonN], generator_->partPx[partonN],generator_->partPy[partonN],generator_->partPz[partonN])[2]; // in cm
  double partonHCALy = intersect(generator_->partVx[partonN],generator_->partVy[partonN],generator_->partVz[partonN], generator_->partPx[partonN],generator_->partPy[partonN],generator_->partPz[partonN])[3];
  double partonHCALz = intersect(generator_->partVx[partonN],generator_->partVy[partonN],generator_->partVz[partonN], generator_->partPx[partonN],generator_->partPy[partonN],generator_->partPz[partonN])[4];
  double TOF_bQuark = 1000000000*(TMath::Sqrt((generator_->partVx[partonN] - partonHCALx)*(generator_->partVx[partonN] - partonHCALx) + (generator_->partVy[partonN] - partonHCALy)*(generator_->partVy[partonN] -partonHCALy) + (generator_->partVz[partonN] - partonHCALz)*(generator_->partVz[partonN] - partonHCALz))) / lightSpeed; // in ns
  double TOF_expected = 1000000000*(TMath::Sqrt(partonHCALx*partonHCALx + partonHCALy*partonHCALy + partonHCALz*partonHCALz)) / lightSpeed; // this is x,y,z of parton intersection with HCAL
  double TOFdelay = TOF_LLP + TOF_bQuark - TOF_expected;
  std::vector<double> LLPinfo;
  LLPinfo.push_back(LLPxDecay);
  LLPinfo.push_back(LLPyDecay);
  LLPinfo.push_back(LLPzDecay);
  LLPinfo.push_back(vertex / (genLLPGamma * genLLPBeta) );
  LLPinfo.push_back(TOFdelay);
  return LLPinfo; // LLP x, y, z decay, LLP ctau in cm in LLP rest frame, TOF delay in ns
}

std::vector<double> ctau(int partonN, L1Analysis::L1AnalysisGeneratorDataFormat *generator_) {
  double beta = -1000;
  beta = sqrt( (generator_->partPx[partonN] * generator_->partPx[partonN]) + (generator_->partPy[partonN] * generator_->partPy[partonN]) + (generator_->partPz[partonN] * generator_->partPz[partonN]) ) / (generator_->partE[partonN]);
  double gamma = 1./TMath::Sqrt(1.-beta*beta);
  double xDist = generator_->partVx[partonN];// - generator_->dauVx[partonN];
  double yDist = generator_->partVy[partonN];// - generator_->dauVy[partonN];
  double zDist = generator_->partVz[partonN];// - generator_->dauVz[partonN];
  double path_length = sqrt(xDist*xDist + yDist*yDist + zDist*zDist);
  double ctau_rest_frame = path_length / ( gamma * beta);
  double lightSpeed = 29979245800; // in cm / s
  double TOF_LLP = 1000000000*path_length / (beta * lightSpeed); // in ns
  std::vector<double> ctau;
  ctau.push_back(path_length); // cm
  ctau.push_back(ctau_rest_frame); // cm
  ctau.push_back(TOF_LLP); //ns
  return ctau;
}


void rates_delayed_cluster(bool newConditions, const std::string& inputFileDirectory, double TDC_HB_variable, double TDC_HE_variable, double GeV_HB_variable, double GeV_HE_variable, double prompt_TP_energy_variable, double prompt_2x2_energy_variable){
  
  bool hwOn = true;   //are we using data from hardware? (upgrade trigger had to be running!!!)
  bool emuOn = true;  //are we using data from emulator?

  // variables used to scan parameters -- energy, time, delayed and prompt seeds
  //  double prompt_2x2_energy_variable = 4; // energy to define energetic 2x2 region (non-delayed)
  //  double prompt_TP_energy_variable = 3; // energy for a prompt TP
  //  double prompt_2x2_TP_variable = 1; // how many high energy TPs to reject jet based on (if >= this variable)
  double delayed_4x4_variable = 2; // how many cells to count as a delayed seed 4x4 region (require >= this variable)
  //  double TDC_HB_variable = 4; // TDC value to be considered delayed
  //  double TDC_HE_variable = 4;
  //  double GeV_HB_variable = 2; // GeV value to be considered delayed
  //  double GeV_HE_variable = 1;

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

  TH1F* TimingBit_TPenergy = new TH1F("TimingBit_TPenergy", "Threshold Study: Delayed Cells (timing bit) in TP E_{T} Bins;TP E_{T};Delayed Cells (normalized)   ", 20,0,20);

  // histograms based on hit multiplicity from the timing bit
  TH1F * Delayed_2x2_MultHB_emu = new TH1F("Delayed_2x2_MultHB_emu","Number of cells in 2x2 region >=50 ADC and >=3ns;Number of cells;Fraction of Entries (normalized)",10,0,10);
  TH1F * Prompt_2x2_MultHB_emu = new TH1F("Prompt_2x2_MultHB_emu","Number of TPs in 2x2 region that are prompt and high energy;Number of TPs;Fraction of Entries (normalized)",10,0,10);
  TH1F * Prompt_Energy_2x2_MultHB_emu = new TH1F("Prompt_Energy_2x2_MultHB_emu","Energy of TPs in 2x2 region that are prompt;Energy of TPs;Fraction of Entries (normalized)",20,0,10);

  TH1F * DeltaR_L1_delayed_seed_emu = new TH1F("DeltaR_L1_delayed_seed_emu","DeltaR between a Delayed Seed and the Closest L1 jet;Delta R;Fraction of Entries (normalized)",20,0,1);
  TH1F * DeltaR_L1_prompt_seed_emu = new TH1F("DeltaR_L1_prompt_seed_emu","DeltaR between a Prompt Seed and the Closest L1 jet;Delta R;Fraction of Entries (normalized)",20,0,1);
  TH1F * DeltaR_L1_delayed_hit_emu = new TH1F("DeltaR_L1_delayed_hit_emu","DeltaR between a Delayed Hit and the Closest L1 Jet in a Seeded Jet;Delta R;Fraction of Entries (normalized)",30,0,1.5);
  TH1F * DeltaR_L1_prompt_hit_emu = new TH1F("DeltaR_L1_prompt_hit_emu","DeltaR between a Prompt Hit and the Closest L1 Jet in a Seeded Jet;Delta R;Fraction of Entries (normalized)",30,0,1.5);
  TH1F * Mult_delayed_hit_emu = new TH1F("Mult_delayed_hit_emu","Number of delayed hits near a seeded jet;Number of cells;Fraction of Entries (normalized)",20,0,20);
  TH1F * Mult_prompt_hit_emu =new TH1F("Mult_prompt_hit_emu","Number of prompt TPs near a seeded jet;Number of TPs;Fraction of Entries (normalized)",20,0,20);

  TH1F * HTdistribution_trig_emu = new TH1F("HTdistribution_trig_emu","HT Distribution of Events Passing Calo Cluster Trigger;HT (GeV);Number of Events",35,0,1200);
  TH1F * HTdistribution_emu = new TH1F("HTdistribution_emu","HT Distribution of Events;HT (GeV);Number of Events",35,0,1200);

  TH1F * JetPTdistribution_trig_emu = new TH1F("JetPTdistribution_trig_emu","Jet pT Distribution passing Delayed Jet Trigger;Jet pT (GeV);Number of Jets",15,0,200);
  TH1F * JetPTdistribution_trig120_emu = new TH1F("JetPTdistribution_trig120_emu","Jet pT Distribution passing Delayed Jet Trigger and HT120;Jet pT (GeV);Number of Jets",15,0,200);
  TH1F * JetPTdistribution_emu = new TH1F("JetPTdistribution_emu","Jet pT Distribution;Jet pT (GeV);Number of Jets",15,0,200);

  // HT sum rate distributions to use in rate vs eff plots. Need HT > 360 rate, and HT > 120 + timing cut rate
  TH1F* htSumRates_original_emu = new TH1F("htSumRates_original_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
  // OR of HT360 "OR" HT120+timing trigger
  TH1F* htSumRates_120timingOR360_1_emu = new TH1F("htSumRates_120timingOR360_1_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
  TH1F* htSumRates_120timingOR360_2_emu = new TH1F("htSumRates_120timingOR360_2_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
  TH1F* htSumRates_120timingOR360_3_emu = new TH1F("htSumRates_120timingOR360_3_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
  TH1F* htSumRates_120timingOR360_4_emu = new TH1F("htSumRates_120timingOR360_4_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
  TH1F* htSumRates_120timingOR360_5_emu = new TH1F("htSumRates_120timingOR360_5_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
  // just rates of HT120+timing (not OR)
  TH1F* htSumRates_120timing_1_emu = new TH1F("htSumRates_120timing_1_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
  TH1F* htSumRates_120timing_2_emu = new TH1F("htSumRates_120timing_2_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
  TH1F* htSumRates_120timing_3_emu = new TH1F("htSumRates_120timing_3_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
  TH1F* htSumRates_120timing_4_emu = new TH1F("htSumRates_120timing_4_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);
  TH1F* htSumRates_120timing_5_emu = new TH1F("htSumRates_120timing_5_emu",axR.c_str(), nHtSumBins, htSumLo, htSumHi);

  TH1F * htSumDistribution = new TH1F("htSumDistribution","htSum Distribution;L1_HT (GeV);Number of Events",200,0,2000); // x bins, x low, x up

  TH2F* TDC_vs_amp = new TH2F("TDC_vs_amp","TDC value and Amplitude Dependence;TDC value (ns);Amplitude (GeV);Number of Events",50,0,25,50,0,10);
  TH2F* amp_vs_TDC = new TH2F("amp_vs_TDC","Amplitude Dependence of TDC value;Amplitude (GeV);TDC value (ns);Number of Events",50,0,10,50,0,25);

  TH1F* ctau_LLP = new TH1F("ctau_LLP","ctau of LLP in particle rest frame;ctau (m);Number of Events",100,0,50);
  TH1F* ctau_LLP_trigger = new TH1F("ctau_LLP_trigger","ctau of LLP (passed trigger) in particle rest frame;ctau (m);Number of Events",100,0,50);
  TH1F* path_length = new TH1F("path_length","Path length of particle in lab frame;path length (m);Number of Events",100,0,3); 
  TH2F* path_length_energy = new TH2F("path_length_energy","Path length (lab frame) vs energy of parton;path length (m);Energy (GeV)",100,0,3,100,0,1000);

  TH2F* prompt_delayed_seed = new TH2F("prompt_delayed_seed","Prompt TPs and Delayed 2x2 Seeds in a L1 jet;Prompt Seeds;Delayed Seeds",15,0,15,15,0,15);
  TH2F* PDGid_radius = new TH2F("PDGid_radius","PDG ID and creation vertex of particle in QCD sample near delayed seed;PDG ID;Creation Vertex (cm)",250,0,250,90,0,45);
  TH2F* pion_radius_z = new TH2F("pion_radius_z","Pion radius and z position (creation vertex) in QCD sample;Z position;Creation Radius (cm)",50,0,50,150,0,150);
  TH2F* delayed4x4seed_TDC_GeV_HB = new TH2F("delayed4x4seed_TDC_GeV_HB","TDC and Energy in cells forming delayed 4x4 seed in HB;TDC;GeV",25,0,25,20,0,20);
  TH2F* delayed4x4seed_depth_TDC_HB = new TH2F("delayed4x4seed_depth_TDC_HB","TDC and Depth in cells forming delayed 4x4 seed in HB;Depth;TDC",5,0,5,20,0,20);
  TH2F* delayed4x4seed_TDC_GeV_HE = new TH2F("delayed4x4seed_TDC_GeV_HE","TDC and Energy in cells forming delayed 4x4 seed in HE;TDC;GeV",25,0,25,20,0,20);
  TH2F* delayed4x4seed_depth_TDC_HE = new TH2F("delayed4x4seed_depth_TDC_HE","TDC and Depth in cells forming delayed 4x4 seed in HE;Depth;TDC",8,0,8,20,0,20);

  // saving rate and efficiencies 
  double passed4JetMult_HBHE_ht120_1(0), passed4JetMult_HBHE_ht120_2(0), passed4JetMult_HBHE_ht120_3(0), passed4JetMult_HBHE_ht120_4(0), passed4JetMult_HBHE_ht120_5(0); //, passed4JetMult_HB_ht120(0),passed4JetMult_HE_ht120(0);
  double passed_calo_cluster_trig(0);
  double passed_calo_cluster_trig_120(0), passed_calo_cluster_trig_120_2(0), passed_calo_cluster_trig_120_3(0), passed_calo_cluster_trig_120_4(0), passed_calo_cluster_trig_120_5(0);
  double passedHtSum360(0);
  double totalEvents(0);
  double totalEvents_HBdr05(0);
  double totalEvent_ht120 = 0, totalEvent_ht140 = 0, totalEvent_ht160 = 0, totalEvent_ht180 = 0, totalEvent_ht200 = 0, totalEvent_ht220 = 0, totalEvent_ht240 = 0, totalEvent_ht260 = 0,totalEvent_ht280 = 0, totalEvent_ht300 = 0,totalEvent_ht320 = 0, totalEvent_ht340 = 0, totalEvent_ht360 = 0, totalEvent_ht380 = 0;
  double HBHE4Jet_inBins_ht120 = 0, HBHE4Jet_inBins_ht140 = 0, HBHE4Jet_inBins_ht160 = 0, HBHE4Jet_inBins_ht180 = 0, HBHE4Jet_inBins_ht200 = 0, HBHE4Jet_inBins_ht220 = 0, HBHE4Jet_inBins_ht240 = 0, HBHE4Jet_inBins_ht260 = 0,HBHE4Jet_inBins_ht280 = 0, HBHE4Jet_inBins_ht300 = 0,HBHE4Jet_inBins_ht320 = 0, HBHE4Jet_inBins_ht340 = 0, HBHE4Jet_inBins_ht360 = 0, HBHE4Jet_inBins_ht380 = 0;
  double HEHB4Jet_ht120_wTiming = 0;

  /////////////////////////////////
  //////////// eta_depth_tdc //////
  /////////////////////////////////
  int eta_depth_tdc95[30][8] = {{4}};
  std::ifstream file("TDCdistribution_Background90_QCD.txt");
  //  std::ifstream file("TDCdistribution_Background_4ns.txt");
  //  std::ifstream file("TDCdistribution_Background_2ns.txt");
  for(int row = 0; row < 30; ++row) {
    std::string line;
    std::getline(file, line);
    if ( !file.good() )
      break;
    std::stringstream iss(line);
    for (int col = 0; col < 8; ++col) { 
      std::string val;
      std::getline(iss, val, ',');
      if ( !iss.good() )
	break;
      std::stringstream convertor(val);
      convertor >> eta_depth_tdc95[row+1][col]; // now col = 1 corresponds to depth = 1
    }
  }
  //  for (int depth = 1; depth < 8; depth++) std::cout << eta_depth_tdc95[19][depth] << std::endl;


  /////////////////////////////////
  // loop through all the entries//
  /////////////////////////////////
  std::ofstream DelayedSeed_event_ieta_iphi_depth;
  DelayedSeed_event_ieta_iphi_depth.open(Form("DelayedSeed_event_ieta_iphi_depth_%s.txt", inputFile.substr(0,7).c_str()),std::ios_base::trunc);

  for (Long64_t jentry=0; jentry<nentries; jentry++){
    if((jentry%10000)==0) std::cout << "Done " << jentry  << " events of " << nentries << std::endl;
    //    std::cout << jentry << std::endl;
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

      for (int partonN = 0; partonN < generator_->nPart; partonN ++) {
        if (generator_->partParent[partonN] == 6000113) ctau_LLP->Fill(LLPdecayInfo(partonN,l1emu_,generator_)[3]/100);
      }      
      // change lines 688 and 719 to only include HB, or HE+HB!! Jet eta < 1,3; and HCAL TP ieta < 16, 28
      //////////////////////////////////////
      /////////// LLP incident on HB? //////
      //////////////////////////////////////
      uint nJetemu = l1emu_->nJets; // number of jets per event
      int numLLPdecayHB = 0; // count how many LLP decay products are expected to intersect the HB
      double nCaloTPemu = l1CaloTPemu_->nHCALTP; // number of TPs varies from 400-1400 per event, approximately Gaussian                                                  
      // triggerability restrictions

      int triggerableJets[nJetemu] = {0};
      for (uint jetIt = 0; jetIt < nJetemu; jetIt++) { // loop over jets
	if (inputFile.substr(0,3) == "QCD") JetPTdistribution_emu->Fill(l1emu_->jetEt[jetIt]);
	if (abs(l1emu_->jetEta[jetIt]) > 2.5) continue; // consider HB jets, HB extends to 1.4. HE extends to 3. Use values of 1, 2.5
	if (closestParton(jetIt, l1emu_, generator_)[0] <= 0.5) { // if closest parton is near a HB L1 jet
	  numLLPdecayHB += 1; // how many of the partons expected to intersect HB
	  JetPTdistribution_emu->Fill(l1emu_->jetEt[jetIt]);
	  triggerableJets[jetIt] = 1; // 1 if triggerable jet, 0 otherwise
	}
      }
      
      //      numLLPdecayHB += 1; // for testing without triggerability restrictions
      //      std::cout << numLLPdecayHB << " = number of LLP decay products incident on HB" <<std::endl;
      if (inputFile.substr(0,2) == "mh" && numLLPdecayHB > 0) totalEvents_HBdr05 += 1;
      if (inputFile.substr(0,2) == "mh" && numLLPdecayHB == 0 ) continue; // if no LLPs in HB, skip event
   
      //////////////////////////////////////
      ////////// HCAL TP Loop //////////
      //////////////////////////////////////
      // per TP, count delayed hits (timing bit) and tower energy for prompt hit consideration
      double timingbit_eta_phi[56][72] = {{0}}; // 32 72 if just HB
      double prompt_TP_10GeV_eta_phi[56][72] = {{0}};
      double prompt_TP_energy_eta_phi[56][72] = {{0}};

      for (int HcalTPIt = 0; HcalTPIt < nCaloTPemu; HcalTPIt++){ // loop over HCAL TPs
	double tpEtaemu = l1CaloTPemu_->hcalTPieta[HcalTPIt]; // ieta
	double tpPhiemu = l1CaloTPemu_->hcalTPCaliphi[HcalTPIt]; // iphi
	if (abs(tpEtaemu) > 29 ) continue; // 16 just HB. don't consider HCAL TPs outside of HB or HE -- note, 29 should be summed with 28 still presumably 
	// End of HB is eta=1.479m, end of HE is eta=3 from here http://www.hephy.at/user/friedl/diss/html/node8.html.

	int TP_ieta_2x2 = eta_phi_2x2(tpEtaemu,tpPhiemu)[0]; // convert from ieta to 2x2 eta mapping
        int TP_iphi_2x2 = eta_phi_2x2(tpEtaemu,tpPhiemu)[1]; // 2x2 phi mapping

	// each TP has the Timing Bit set with the multiplicity for that tower 
	//	int TP_TimingBit = l1CaloTPemu_->hcalTPTimingBit[HcalTPIt];
	TimingBit_TPenergy->Fill(l1CaloTPemu_->hcalTPet[HcalTPIt], l1CaloTPemu_->hcalTPTimingBit[HcalTPIt]); // fill bin TP energy with value TP timing bit

	// TDC and amplitude dependence investigations, but only for cells where energy is non zero
	if(l1CaloTPemu_->hcalTPDepth1[HcalTPIt] > 0) TDC_vs_amp->Fill(l1CaloTPemu_->hcalTPtiming1[HcalTPIt],l1CaloTPemu_->hcalTPDepth1[HcalTPIt]);
        if(l1CaloTPemu_->hcalTPDepth2[HcalTPIt] > 0) TDC_vs_amp->Fill(l1CaloTPemu_->hcalTPtiming2[HcalTPIt],l1CaloTPemu_->hcalTPDepth2[HcalTPIt]);
	if(l1CaloTPemu_->hcalTPDepth3[HcalTPIt] > 0) TDC_vs_amp->Fill(l1CaloTPemu_->hcalTPtiming3[HcalTPIt],l1CaloTPemu_->hcalTPDepth3[HcalTPIt]);
	if(l1CaloTPemu_->hcalTPDepth4[HcalTPIt] > 0) TDC_vs_amp->Fill(l1CaloTPemu_->hcalTPtiming4[HcalTPIt],l1CaloTPemu_->hcalTPDepth4[HcalTPIt]);
	if(l1CaloTPemu_->hcalTPDepth5[HcalTPIt] > 0) TDC_vs_amp->Fill(l1CaloTPemu_->hcalTPtiming5[HcalTPIt],l1CaloTPemu_->hcalTPDepth5[HcalTPIt]);
        if(l1CaloTPemu_->hcalTPDepth6[HcalTPIt] > 0) TDC_vs_amp->Fill(l1CaloTPemu_->hcalTPtiming6[HcalTPIt],l1CaloTPemu_->hcalTPDepth6[HcalTPIt]);
        if(l1CaloTPemu_->hcalTPDepth7[HcalTPIt] > 0) TDC_vs_amp->Fill(l1CaloTPemu_->hcalTPtiming7[HcalTPIt],l1CaloTPemu_->hcalTPDepth7[HcalTPIt]);
        if(l1CaloTPemu_->hcalTPDepth1[HcalTPIt] > 0) amp_vs_TDC->Fill(l1CaloTPemu_->hcalTPDepth1[HcalTPIt],l1CaloTPemu_->hcalTPtiming1[HcalTPIt]);
        if(l1CaloTPemu_->hcalTPDepth2[HcalTPIt] > 0) amp_vs_TDC->Fill(l1CaloTPemu_->hcalTPDepth2[HcalTPIt],l1CaloTPemu_->hcalTPtiming2[HcalTPIt]);
	if(l1CaloTPemu_->hcalTPDepth3[HcalTPIt] > 0) amp_vs_TDC->Fill(l1CaloTPemu_->hcalTPDepth3[HcalTPIt],l1CaloTPemu_->hcalTPtiming3[HcalTPIt]);
        if(l1CaloTPemu_->hcalTPDepth4[HcalTPIt] > 0) amp_vs_TDC->Fill(l1CaloTPemu_->hcalTPDepth4[HcalTPIt],l1CaloTPemu_->hcalTPtiming4[HcalTPIt]);
        if(l1CaloTPemu_->hcalTPDepth5[HcalTPIt] > 0) amp_vs_TDC->Fill(l1CaloTPemu_->hcalTPDepth5[HcalTPIt],l1CaloTPemu_->hcalTPtiming5[HcalTPIt]);
        if(l1CaloTPemu_->hcalTPDepth6[HcalTPIt] > 0) amp_vs_TDC->Fill(l1CaloTPemu_->hcalTPDepth6[HcalTPIt],l1CaloTPemu_->hcalTPtiming6[HcalTPIt]);
        if(l1CaloTPemu_->hcalTPDepth7[HcalTPIt] > 0) amp_vs_TDC->Fill(l1CaloTPemu_->hcalTPDepth7[HcalTPIt],l1CaloTPemu_->hcalTPtiming7[HcalTPIt]);

	// instead of using pre-set timing bit (3ns 3GeV), set specifically here based on energy and time values
	if (abs(tpEtaemu) <= 16) { 
	  if (l1CaloTPemu_->hcalTPDepth1[HcalTPIt] >= GeV_HB_variable && l1CaloTPemu_->hcalTPtiming1[HcalTPIt]*2 >= eta_depth_tdc95[static_cast<int>(abs(tpEtaemu))][1] ) timingbit_eta_phi[TP_ieta_2x2][TP_iphi_2x2] += 1;//TDC_HB_variable) timingbit_eta_phi[TP_ieta_2x2][TP_iphi_2x2] += 1;
	  if (l1CaloTPemu_->hcalTPDepth2[HcalTPIt] >= GeV_HB_variable && l1CaloTPemu_->hcalTPtiming2[HcalTPIt]*2 >= eta_depth_tdc95[static_cast<int>(abs(tpEtaemu))][2]) timingbit_eta_phi[TP_ieta_2x2][TP_iphi_2x2] += 1;
	  if (l1CaloTPemu_->hcalTPDepth3[HcalTPIt] >= GeV_HB_variable && l1CaloTPemu_->hcalTPtiming3[HcalTPIt]*2 >= eta_depth_tdc95[static_cast<int>(abs(tpEtaemu))][3]) timingbit_eta_phi[TP_ieta_2x2][TP_iphi_2x2] += 1;
	  if (l1CaloTPemu_->hcalTPDepth4[HcalTPIt] >= GeV_HB_variable && l1CaloTPemu_->hcalTPtiming4[HcalTPIt]*2 >= eta_depth_tdc95[static_cast<int>(abs(tpEtaemu))][4]) timingbit_eta_phi[TP_ieta_2x2][TP_iphi_2x2] += 1;
	}

	if (abs(tpEtaemu) > 16) { // lower energy thresholds in HE, since using transverse energy goes as 1/cosh(eta)
	  // investigate which depth layers to include, timing bit excludes depth layer 1

	  //	  if (l1CaloTPemu_->hcalTPDepth1[HcalTPIt] >= 3 && l1CaloTPemu_->hcalTPtiming1[HcalTPIt] >= 3) timingbit_eta_phi[TP_ieta_2x2][TP_iphi_2x2] += 1;
          if (l1CaloTPemu_->hcalTPDepth2[HcalTPIt] >= GeV_HE_variable && l1CaloTPemu_->hcalTPtiming2[HcalTPIt]*2 >= eta_depth_tdc95[static_cast<int>(abs(tpEtaemu))][2] ) timingbit_eta_phi[TP_ieta_2x2][TP_iphi_2x2] += 1; //TDC_HE_variable) timingbit_eta_phi[TP_ieta_2x2][TP_iphi_2x2] += 1;
	  if (l1CaloTPemu_->hcalTPDepth3[HcalTPIt] >= GeV_HE_variable && l1CaloTPemu_->hcalTPtiming3[HcalTPIt]*2 >= eta_depth_tdc95[static_cast<int>(abs(tpEtaemu))][3]) timingbit_eta_phi[TP_ieta_2x2][TP_iphi_2x2] += 1;
	  if (l1CaloTPemu_->hcalTPDepth4[HcalTPIt] >= GeV_HE_variable && l1CaloTPemu_->hcalTPtiming4[HcalTPIt]*2 >= eta_depth_tdc95[static_cast<int>(abs(tpEtaemu))][4]) timingbit_eta_phi[TP_ieta_2x2][TP_iphi_2x2] += 1;
	  if (l1CaloTPemu_->hcalTPDepth5[HcalTPIt] >= GeV_HE_variable && l1CaloTPemu_->hcalTPtiming5[HcalTPIt]*2 >= eta_depth_tdc95[static_cast<int>(abs(tpEtaemu))][5]) timingbit_eta_phi[TP_ieta_2x2][TP_iphi_2x2] += 1;
	  if (l1CaloTPemu_->hcalTPDepth6[HcalTPIt] >= GeV_HE_variable && l1CaloTPemu_->hcalTPtiming6[HcalTPIt]*2 >= eta_depth_tdc95[static_cast<int>(abs(tpEtaemu))][6]) timingbit_eta_phi[TP_ieta_2x2][TP_iphi_2x2] += 1;
	  if (l1CaloTPemu_->hcalTPDepth7[HcalTPIt] >= GeV_HE_variable && l1CaloTPemu_->hcalTPtiming7[HcalTPIt]*2 >= eta_depth_tdc95[static_cast<int>(abs(tpEtaemu))][7]) timingbit_eta_phi[TP_ieta_2x2][TP_iphi_2x2] += 1;
	}

	// prompt veto, counts high energy prompt TPs in region with no delayed hits
	if (l1CaloTPemu_->hcalTPet[HcalTPIt] >= 0) prompt_TP_energy_eta_phi[TP_ieta_2x2][TP_iphi_2x2] = l1CaloTPemu_->hcalTPet[HcalTPIt];
	if (l1CaloTPemu_->hcalTPet[HcalTPIt] >= prompt_TP_energy_variable) prompt_TP_10GeV_eta_phi[TP_ieta_2x2][TP_iphi_2x2] = 1; // sum of depths, transverse energy. Track if TP is above energy threshold (1 if above, 0 if not)
      } // closing HCAL TP loop

      int delayed_calo_objects[nJetemu] = {0}; // delayed calo objects in each L1 jet
      int delayed[nJetemu] = {0}; // delayed hits in each L1 jet
      int prompt_TP[nJetemu] = {0}; // prompt towers in each L1 jet
      int prompt_energy[nJetemu] = {0}; // number of 2x2 non-delayed regions that have high energy
      //      int delayed_2x2_count[32][72] = {{0}}; // per 2x2, number of delayed hits. Odd eta, phi entries will be not be filled -- not corner of a 2x2!
      //      int prompt_2x2_count[32][72] = {{0}}; // per 2x2, number of prompt high energy TPs (when 2x2 has no delayed hits)

      // check 2x2 regions for seeds
      for (int TP_ieta_2x2 = 0; TP_ieta_2x2<56; TP_ieta_2x2+=2) { // 32 for HB only
	for (int TP_iphi_2x2 = 0; TP_iphi_2x2<72; TP_iphi_2x2+=2) {
	  int delayed_2x2 = 0;
	  int delayed_4x4 = 0;
	  int prompt_2x2 = 0;
	  double prompt_energy_2x2 = 0;
	  // grid of 56x72, now break into 2x2 tower sums for a grid fo 28x36 // grid of 32x72, now break into 2x2 tower sums for a new grid of 16x36
	  // even, even = top left corner of 2x2 regions
	  delayed_2x2 += timingbit_eta_phi[TP_ieta_2x2][TP_iphi_2x2] + timingbit_eta_phi[TP_ieta_2x2][TP_iphi_2x2+1] + timingbit_eta_phi[TP_ieta_2x2+1][TP_iphi_2x2] + timingbit_eta_phi[TP_ieta_2x2+1][TP_iphi_2x2+1];
	  if (TP_ieta_2x2 > 0 && TP_iphi_2x2 > 0) delayed_4x4 += delayed_2x2 
					      + timingbit_eta_phi[TP_ieta_2x2][TP_iphi_2x2-2] + timingbit_eta_phi[TP_ieta_2x2][TP_iphi_2x2-1] + timingbit_eta_phi[TP_ieta_2x2+1][TP_iphi_2x2-2] + timingbit_eta_phi[TP_ieta_2x2+1][TP_iphi_2x2-1]
					      + timingbit_eta_phi[TP_ieta_2x2-2][TP_iphi_2x2] + timingbit_eta_phi[TP_ieta_2x2-2][TP_iphi_2x2+1] + timingbit_eta_phi[TP_ieta_2x2-1][TP_iphi_2x2] + timingbit_eta_phi[TP_ieta_2x2-1][TP_iphi_2x2+1]
					      + timingbit_eta_phi[TP_ieta_2x2-2][TP_iphi_2x2-2] + timingbit_eta_phi[TP_ieta_2x2-2][TP_iphi_2x2-1] + timingbit_eta_phi[TP_ieta_2x2-1][TP_iphi_2x2-2] + timingbit_eta_phi[TP_ieta_2x2-1][TP_iphi_2x2-1];
	  if (delayed_2x2 == 0) { // no delayed hits in this 2x2 region
	    prompt_2x2 += prompt_TP_10GeV_eta_phi[TP_ieta_2x2][TP_iphi_2x2] + prompt_TP_10GeV_eta_phi[TP_ieta_2x2][TP_iphi_2x2+1] + prompt_TP_10GeV_eta_phi[TP_ieta_2x2+1][TP_iphi_2x2] + prompt_TP_10GeV_eta_phi[TP_ieta_2x2+1][TP_iphi_2x2+1]; // how many high energy TPs in non-delayed 2x2 region (>3 GeV currently)
	    prompt_energy_2x2 += prompt_TP_energy_eta_phi[TP_ieta_2x2][TP_iphi_2x2] + prompt_TP_energy_eta_phi[TP_ieta_2x2][TP_iphi_2x2+1] + prompt_TP_energy_eta_phi[TP_ieta_2x2+1][TP_iphi_2x2] + prompt_TP_energy_eta_phi[TP_ieta_2x2+1][TP_iphi_2x2+1]; // total energy of non-delayed 2x2
	  }

	  // find standard ieta iphi value from 2x2 ieta and iphi used to define prompt and delayed objects
	  double seed_eta = etaVal(two_two_eta_phi(TP_ieta_2x2, TP_iphi_2x2)[0]);
	  double seed_ieta_1 = two_two_eta_phi(TP_ieta_2x2, TP_iphi_2x2)[0];
          double seed_ieta_2 = two_two_eta_phi(TP_ieta_2x2+1, TP_iphi_2x2)[0];
          double seed_ieta_3 = two_two_eta_phi(TP_ieta_2x2-1, TP_iphi_2x2)[0];
          double seed_ieta_4 = two_two_eta_phi(TP_ieta_2x2-2, TP_iphi_2x2)[0];
	  double seed_phi = phiVal(two_two_eta_phi(TP_ieta_2x2, TP_iphi_2x2)[1]);
	  double seed_iphi_1 = two_two_eta_phi(TP_ieta_2x2, TP_iphi_2x2)[1];
          double seed_iphi_2 = two_two_eta_phi(TP_ieta_2x2, TP_iphi_2x2+1)[1];
          double seed_iphi_3 = two_two_eta_phi(TP_ieta_2x2, TP_iphi_2x2-1)[1];
          double seed_iphi_4 = two_two_eta_phi(TP_ieta_2x2, TP_iphi_2x2-2)[1];

	  int closestJet = -1;
	  double min_DR = 100;
	  // assign corner of 2x2 to closest L1 jet
	  for (uint jetIt = 0; jetIt < nJetemu; jetIt++) {
	    if (deltaR(l1emu_->jetEta[jetIt], l1emu_->jetPhi[jetIt], seed_eta, seed_phi)<min_DR) {
	      min_DR = deltaR(l1emu_->jetEta[jetIt], l1emu_->jetPhi[jetIt], seed_eta, seed_phi);
	      closestJet = jetIt;
	    }
	  }

          if (delayed_2x2 > 0) DeltaR_L1_delayed_hit_emu->Fill(min_DR); // DR from L1 jet center to a delayed hit
          if (prompt_2x2 > 0) DeltaR_L1_prompt_hit_emu->Fill(min_DR); // DR from L1 jet center to a prompt TP
	  if (delayed_4x4 >= delayed_4x4_variable) DeltaR_L1_delayed_seed_emu->Fill(min_DR); // DR from L1 jet center to delayed seed
	  if (prompt_energy_2x2 >= prompt_2x2_energy_variable) DeltaR_L1_prompt_seed_emu->Fill(min_DR); // DR from L1 jet center to prompt seed
	
	  if (min_DR<=0.5) { // check if the jet has a delayed or prompt seed
	    delayed[closestJet] += delayed_2x2; // assign delayed 2x2 hits to the nearest jet
	    prompt_TP[closestJet] += prompt_2x2;
	    if (prompt_energy_2x2 >= prompt_2x2_energy_variable) prompt_energy[closestJet] += 1; // how many 2x2 are non delayed but energetic
	    if (delayed_4x4 >= delayed_4x4_variable) {
	      delayed_calo_objects[closestJet] += 1; // how many 4x4 are delayed near each jet
	      //	      std::cout << seed_eta << " = seed eta, seed phi = " << seed_phi << std::endl;
	      //	      std::cout << seed_ieta_1 << " = seed ieta, seed iphi = " << seed_iphi_1 << std::endl;
	      // find actual time and energy values contributing to the delayed 4x4 seed

	      DelayedSeed_event_ieta_iphi_depth << "Event = " << jentry << " HT = " << htSum << std::endl;

	      for (int HcalTPIt = 0; HcalTPIt < nCaloTPemu; HcalTPIt++){ // loop over HCAL TPs
		double tpEtaemu = l1CaloTPemu_->hcalTPieta[HcalTPIt]; // ieta
		double tpPhiemu = l1CaloTPemu_->hcalTPCaliphi[HcalTPIt]; // iphi
		if ( (seed_ieta_1 == tpEtaemu) || (seed_ieta_2 == tpEtaemu) || (seed_ieta_3 == tpEtaemu) || (seed_ieta_4 == tpEtaemu) ) { // ieta
		  if ( (seed_iphi_1 == tpPhiemu) || (seed_iphi_2 == tpPhiemu) || (seed_iphi_3 == tpPhiemu) || (seed_iphi_4 == tpPhiemu) ) { // iphi
		    if (abs(tpEtaemu) > 16) { // lower energy thresholds in HE, since using transverse energy goes as 1/cosh(eta)
		      if (l1CaloTPemu_->hcalTPDepth2[HcalTPIt] >= GeV_HE_variable && l1CaloTPemu_->hcalTPtiming2[HcalTPIt]*2 >= eta_depth_tdc95[static_cast<int>(abs(tpEtaemu))][2]) {
			delayed4x4seed_TDC_GeV_HE->Fill(l1CaloTPemu_->hcalTPtiming2[HcalTPIt], l1CaloTPemu_->hcalTPDepth2[HcalTPIt]);
                        delayed4x4seed_depth_TDC_HE->Fill(2, l1CaloTPemu_->hcalTPtiming2[HcalTPIt]);
			DelayedSeed_event_ieta_iphi_depth << "ieta, iphi, depth, GeV, ns " << tpEtaemu << ", " << tpPhiemu << ", " << 2 << ", " << l1CaloTPemu_->hcalTPDepth2[HcalTPIt] << ", " << l1CaloTPemu_->hcalTPtiming2[HcalTPIt] << std::endl;
		      }
		      if (l1CaloTPemu_->hcalTPDepth3[HcalTPIt] >= GeV_HE_variable && l1CaloTPemu_->hcalTPtiming3[HcalTPIt]*2 >= eta_depth_tdc95[static_cast<int>(abs(tpEtaemu))][3]) {
			delayed4x4seed_TDC_GeV_HE->Fill(l1CaloTPemu_->hcalTPtiming3[HcalTPIt], l1CaloTPemu_->hcalTPDepth3[HcalTPIt]);
                        delayed4x4seed_depth_TDC_HE->Fill(3, l1CaloTPemu_->hcalTPtiming3[HcalTPIt]);
                        DelayedSeed_event_ieta_iphi_depth << "ieta, iphi, depth, GeV, ns " << tpEtaemu << ", " << tpPhiemu << ", " << 3 << ", " << l1CaloTPemu_->hcalTPDepth3[HcalTPIt] << ", " << l1CaloTPemu_->hcalTPtiming3[HcalTPIt] << std::endl;
		      }
		      if (l1CaloTPemu_->hcalTPDepth4[HcalTPIt] >= GeV_HE_variable && l1CaloTPemu_->hcalTPtiming4[HcalTPIt]*2 >= eta_depth_tdc95[static_cast<int>(abs(tpEtaemu))][4]) {
			delayed4x4seed_TDC_GeV_HE->Fill(l1CaloTPemu_->hcalTPtiming4[HcalTPIt], l1CaloTPemu_->hcalTPDepth4[HcalTPIt]);
                        delayed4x4seed_depth_TDC_HE->Fill(4, l1CaloTPemu_->hcalTPtiming4[HcalTPIt]);
                        DelayedSeed_event_ieta_iphi_depth << "ieta, iphi, depth, GeV, ns " << tpEtaemu << ", " << tpPhiemu << ", " << 4 << ", " << l1CaloTPemu_->hcalTPDepth4[HcalTPIt] << ", " << l1CaloTPemu_->hcalTPtiming4[HcalTPIt] << std::endl;
		      }
		      if (l1CaloTPemu_->hcalTPDepth5[HcalTPIt] >= GeV_HE_variable && l1CaloTPemu_->hcalTPtiming5[HcalTPIt]*2 >= eta_depth_tdc95[static_cast<int>(abs(tpEtaemu))][5]) {
			delayed4x4seed_TDC_GeV_HE->Fill(l1CaloTPemu_->hcalTPtiming5[HcalTPIt], l1CaloTPemu_->hcalTPDepth5[HcalTPIt]);
                        delayed4x4seed_depth_TDC_HE->Fill(5, l1CaloTPemu_->hcalTPtiming5[HcalTPIt]);
                        DelayedSeed_event_ieta_iphi_depth << "ieta, iphi, depth, GeV, ns " << tpEtaemu << ", " << tpPhiemu << ", " << 5 << ", " << l1CaloTPemu_->hcalTPDepth5[HcalTPIt] << ", " << l1CaloTPemu_->hcalTPtiming5[HcalTPIt] << std::endl;
		      }
		      if (l1CaloTPemu_->hcalTPDepth6[HcalTPIt] >= GeV_HE_variable && l1CaloTPemu_->hcalTPtiming6[HcalTPIt]*2 >= eta_depth_tdc95[static_cast<int>(abs(tpEtaemu))][6]) {
			delayed4x4seed_TDC_GeV_HE->Fill(l1CaloTPemu_->hcalTPtiming6[HcalTPIt], l1CaloTPemu_->hcalTPDepth6[HcalTPIt]);
                        delayed4x4seed_depth_TDC_HE->Fill(6, l1CaloTPemu_->hcalTPtiming6[HcalTPIt]);
                        DelayedSeed_event_ieta_iphi_depth << "ieta, iphi, depth, GeV, ns " << tpEtaemu << ", " << tpPhiemu << ", " << 6 << ", " << l1CaloTPemu_->hcalTPDepth6[HcalTPIt] << ", " << l1CaloTPemu_->hcalTPtiming6[HcalTPIt] << std::endl;
		      }
		      if (l1CaloTPemu_->hcalTPDepth7[HcalTPIt] >= GeV_HE_variable && l1CaloTPemu_->hcalTPtiming7[HcalTPIt]*2 >= eta_depth_tdc95[static_cast<int>(abs(tpEtaemu))][7]) {
			delayed4x4seed_TDC_GeV_HE->Fill(l1CaloTPemu_->hcalTPtiming7[HcalTPIt], l1CaloTPemu_->hcalTPDepth7[HcalTPIt]);
                        delayed4x4seed_depth_TDC_HE->Fill(7, l1CaloTPemu_->hcalTPtiming7[HcalTPIt]);
                        DelayedSeed_event_ieta_iphi_depth << "ieta, iphi, depth, GeV, ns " << tpEtaemu << ", " << tpPhiemu << ", " << 7 << ", " << l1CaloTPemu_->hcalTPDepth7[HcalTPIt] << ", " << l1CaloTPemu_->hcalTPtiming7[HcalTPIt] << std::endl;
		      }
		    } // HE
		    if (abs(tpEtaemu) <= 16) {
                      if (l1CaloTPemu_->hcalTPDepth1[HcalTPIt] >= GeV_HB_variable && l1CaloTPemu_->hcalTPtiming1[HcalTPIt]*2 >= eta_depth_tdc95[static_cast<int>(abs(tpEtaemu))][1]) {
			delayed4x4seed_TDC_GeV_HB->Fill(l1CaloTPemu_->hcalTPtiming1[HcalTPIt], l1CaloTPemu_->hcalTPDepth1[HcalTPIt]);
			delayed4x4seed_depth_TDC_HB->Fill(1, l1CaloTPemu_->hcalTPtiming1[HcalTPIt]);
                        DelayedSeed_event_ieta_iphi_depth << "ieta, iphi, depth, GeV, ns " << tpEtaemu << ", " << tpPhiemu << ", " << 1 << ", " << l1CaloTPemu_->hcalTPDepth1[HcalTPIt] << ", " << l1CaloTPemu_->hcalTPtiming1[HcalTPIt] << std::endl;
		      }
		      if (l1CaloTPemu_->hcalTPDepth2[HcalTPIt] >= GeV_HB_variable && l1CaloTPemu_->hcalTPtiming2[HcalTPIt]*2 >= eta_depth_tdc95[static_cast<int>(abs(tpEtaemu))][2]) {
			delayed4x4seed_TDC_GeV_HB->Fill(l1CaloTPemu_->hcalTPtiming2[HcalTPIt], l1CaloTPemu_->hcalTPDepth2[HcalTPIt]);
                        delayed4x4seed_depth_TDC_HB->Fill(2, l1CaloTPemu_->hcalTPtiming2[HcalTPIt]);
                        DelayedSeed_event_ieta_iphi_depth << "ieta, iphi, depth, GeV, ns " << tpEtaemu << ", " << tpPhiemu << ", " << 2 << ", " << l1CaloTPemu_->hcalTPDepth2[HcalTPIt] << ", " << l1CaloTPemu_->hcalTPtiming2[HcalTPIt] << std::endl;
		      }
                      if (l1CaloTPemu_->hcalTPDepth3[HcalTPIt] >= GeV_HB_variable && l1CaloTPemu_->hcalTPtiming3[HcalTPIt]*2 >= eta_depth_tdc95[static_cast<int>(abs(tpEtaemu))][3]) {
			delayed4x4seed_TDC_GeV_HB->Fill(l1CaloTPemu_->hcalTPtiming3[HcalTPIt], l1CaloTPemu_->hcalTPDepth3[HcalTPIt]);
                        delayed4x4seed_depth_TDC_HB->Fill(3, l1CaloTPemu_->hcalTPtiming3[HcalTPIt]);
                        DelayedSeed_event_ieta_iphi_depth << "ieta, iphi, depth, GeV, ns " << tpEtaemu << ", " << tpPhiemu << ", " << 3 << ", " << l1CaloTPemu_->hcalTPDepth3[HcalTPIt] << ", " << l1CaloTPemu_->hcalTPtiming3[HcalTPIt] << std::endl;
		      }
                      if (l1CaloTPemu_->hcalTPDepth4[HcalTPIt] >= GeV_HB_variable && l1CaloTPemu_->hcalTPtiming4[HcalTPIt]*2 >= eta_depth_tdc95[static_cast<int>(abs(tpEtaemu))][4]) {
			delayed4x4seed_TDC_GeV_HB->Fill(l1CaloTPemu_->hcalTPtiming4[HcalTPIt], l1CaloTPemu_->hcalTPDepth4[HcalTPIt]);
                        delayed4x4seed_depth_TDC_HB->Fill(4, l1CaloTPemu_->hcalTPtiming4[HcalTPIt]);
                        DelayedSeed_event_ieta_iphi_depth << "ieta, iphi, depth, GeV, ns " << tpEtaemu << ", " << tpPhiemu << ", " << 4 << ", " << l1CaloTPemu_->hcalTPDepth4[HcalTPIt] << ", " << l1CaloTPemu_->hcalTPtiming4[HcalTPIt] << std::endl;
		      }
		    } // HB
		  } // iphi for 4x4 delayed seed
		} // ieta for 4x4 delayed seed
	      } // closing HCAL TP loop
	    }
	    Delayed_2x2_MultHB_emu->Fill(delayed_4x4); // number of cells in 2x2 that are delayed
	    if (delayed_2x2 == 0 && prompt_2x2 > 0) Prompt_2x2_MultHB_emu->Fill(prompt_2x2);
            if (delayed_2x2 == 0 && prompt_energy_2x2 > 0) Prompt_Energy_2x2_MultHB_emu->Fill(prompt_energy_2x2);
	  }
	}
      }

      // how many jets pass the delayed trigger
      int num_delayed_jet = 0;
      int num_delayed_obj[nJetemu] = {0};

      for (uint jetIt = 0; jetIt < nJetemu; jetIt++) {
	if (l1emu_->jetEt[jetIt] < 40 ) continue; // for jet pt efficiencies
	if ((inputFile.substr(0,2) == "mh") && (triggerableJets[jetIt] == 0)) continue;
	if (delayed_calo_objects[jetIt] >= 1) { // make sure a jet is seeded!	
	  Mult_delayed_hit_emu->Fill(delayed[jetIt]);
	  Mult_prompt_hit_emu->Fill(prompt_TP[jetIt]);
	}
	if (delayed_calo_objects[jetIt] > 0) prompt_delayed_seed->Fill(prompt_TP[jetIt], delayed_calo_objects[jetIt]);

	// num_delayed_jet is 1 if has delayed seed, 2 if delayed seed and passed prompt TP veto, 3 if delayed seed and passed prompt TP veto and prompt 2x2 veto. Used for ROC curves instead of scanning n_delayed_jets. Progressively add on requirements, in this order
	if (delayed_calo_objects[jetIt] >= 1) num_delayed_obj[jetIt] += 1;
	//	if (delayed_calo_objects[jetIt] >= 1 && prompt_TP[jetIt] < prompt_2x2_TP_variable) num_delayed_obj[jetIt] += 2;
	if (delayed_calo_objects[jetIt] >= 1 && prompt_energy[jetIt] == 0) num_delayed_obj[jetIt] += 2;
	//	if (delayed_calo_objects[jetIt] >= 1 && prompt_TP[jetIt] < prompt_2x2_TP_variable && prompt_energy[jetIt] == 0) num_delayed_obj[jetIt] += 1;

	// for jet PT efficiencies
	if (num_delayed_obj[jetIt] >= 3) JetPTdistribution_trig_emu->Fill(l1emu_->jetEt[jetIt]);
        if (num_delayed_obj[jetIt] >= 3 && htSum > 120) JetPTdistribution_trig120_emu->Fill(l1emu_->jetEt[jetIt]);
      }
      std::sort(num_delayed_obj, num_delayed_obj+nJetemu, std::greater<int>());
      if (num_delayed_obj[0] == 1) num_delayed_jet = 1;
      if (num_delayed_obj[0] == 2) num_delayed_jet = 2;
      if (num_delayed_obj[0] == 3) num_delayed_jet = 3;
      if (num_delayed_obj[0] == 3 && num_delayed_obj[1] > 0) num_delayed_jet = 4;
      if (num_delayed_obj[0] == 3 && num_delayed_obj[1] > 0  && num_delayed_obj[2] > 0) num_delayed_jet = 5;
      
      if (num_delayed_jet >= 3) HTdistribution_trig_emu->Fill(htSum); // plot HT dist of events passing calo trigger 
      HTdistribution_emu->Fill(htSum);

      // saving number of events passed cell multiplicity cuts. These are efficiency values used in rate vs eff plots
      totalEvents += 1; // counting total events for efficiency calculations
      htSumDistribution->Fill(htSum);

      if (htSum > 120 && htSum <=140 ) totalEvent_ht120 += 1;
      if (htSum > 140 && htSum <=160 ) totalEvent_ht140 += 1;
      if (htSum > 160 && htSum <=180 ) totalEvent_ht160 += 1;
      if (htSum > 180 && htSum <=200 ) totalEvent_ht180 += 1;
      if (htSum > 200 && htSum <=220 ) totalEvent_ht200 += 1;
      if (htSum > 220 && htSum <=240 ) totalEvent_ht220 += 1;
      if (htSum > 240 && htSum <=260 ) totalEvent_ht240 += 1;
      if (htSum > 260 && htSum <=280 ) totalEvent_ht260 += 1;
      if (htSum > 280 && htSum <=300 ) totalEvent_ht280 += 1;
      if (htSum > 300 && htSum <=320 ) totalEvent_ht300 += 1;
      if (htSum > 320 && htSum <=340 ) totalEvent_ht320 += 1;
      if (htSum > 340 && htSum <=360 ) totalEvent_ht340 += 1;
      if (htSum > 360 && htSum <=380 ) totalEvent_ht360 += 1;
      if (htSum > 380 && htSum <=400 ) totalEvent_ht380 += 1;

      //      if (Sum4Jet_HBHE >= 2) {
      if (num_delayed_jet >= 3) {
	if (htSum > 120 ) HEHB4Jet_ht120_wTiming += 1; // events passing 2 hits over 3ns, 50 ADC // events passing delayed jet and HT 120
	if (htSum > 120 && htSum <=140 ) HBHE4Jet_inBins_ht120 += 1;
	if (htSum > 140 && htSum <=160 ) HBHE4Jet_inBins_ht140 += 1;
	if (htSum > 160 && htSum <=180 ) HBHE4Jet_inBins_ht160 += 1;
	if (htSum > 180 && htSum <=200 ) HBHE4Jet_inBins_ht180 += 1;
	if (htSum > 200 && htSum <=220 ) HBHE4Jet_inBins_ht200 += 1;
	if (htSum > 220 && htSum <=240 ) HBHE4Jet_inBins_ht220 += 1;
	if (htSum > 240 && htSum <=260 ) HBHE4Jet_inBins_ht240 += 1;
	if (htSum > 260 && htSum <=280 ) HBHE4Jet_inBins_ht260 += 1;
	if (htSum > 280 && htSum <=300 ) HBHE4Jet_inBins_ht280 += 1;
	if (htSum > 300 && htSum <=320 ) HBHE4Jet_inBins_ht300 += 1;
	if (htSum > 320 && htSum <=340 ) HBHE4Jet_inBins_ht320 += 1;
	if (htSum > 340 && htSum <=360 ) HBHE4Jet_inBins_ht340 += 1;
	if (htSum > 360 && htSum <=380 ) HBHE4Jet_inBins_ht360 += 1;
	if (htSum > 380 && htSum <=400 ) HBHE4Jet_inBins_ht380 += 1;
      }
      
      if ( num_delayed_jet >= 1 ) passed_calo_cluster_trig += 1;
      if ( num_delayed_jet >= 1 && htSum > 120 ) passed_calo_cluster_trig_120 += 1;
      if ( num_delayed_jet >= 2 && htSum > 120 ) passed_calo_cluster_trig_120_2 += 1;
      if ( num_delayed_jet >= 3 && htSum > 120 ) passed_calo_cluster_trig_120_3 += 1;
      if ( num_delayed_jet >= 4 && htSum > 120 ) passed_calo_cluster_trig_120_4 += 1;
      if ( num_delayed_jet >= 5 && htSum > 120 ) passed_calo_cluster_trig_120_5 += 1;
      if ( ((htSum > 120) && ( num_delayed_jet >= 1 )) || (htSum >= 360) ) passed4JetMult_HBHE_ht120_1 += 1; // +7
      if ( ((htSum > 120) && ( num_delayed_jet >= 2 )) || (htSum >= 360) ) passed4JetMult_HBHE_ht120_2 += 1;
      if ( ((htSum > 120) && ( num_delayed_jet >= 3 )) || (htSum >= 360) ) passed4JetMult_HBHE_ht120_3 += 1;
      if ( ((htSum > 120) && ( num_delayed_jet >= 4 )) || (htSum >= 360) ) passed4JetMult_HBHE_ht120_4 += 1;
      if ( ((htSum > 120) && ( num_delayed_jet >= 5 )) || (htSum >= 360) ) passed4JetMult_HBHE_ht120_5 += 1;

      if (htSum > 360) passedHtSum360 += 1;

      for(int bin=0; bin<nHtSumBins; bin++){
        if( (htSum) >= htSumLo+(bin*htSumBinWidth) && ( num_delayed_jet)>=3 ) htSumRates_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV, compare to original rates in plot
        if( (htSum) >= htSumLo+(bin*htSumBinWidth) ) htSumRates_original_emu->Fill(htSumLo+(bin*htSumBinWidth)); //GeV, use for ht > 360 original rates in rate vs. eff plots

	//GeV, use for ht > 120 + timing OR ht > 360 rates in rate vs. eff plots
        if( (htSum) >= htSumLo+(bin*htSumBinWidth) && ( (num_delayed_jet)>=1 || htSum >= 360 ) ) htSumRates_120timingOR360_1_emu->Fill(htSumLo+(bin*htSumBinWidth));
	if( (htSum) >= htSumLo+(bin*htSumBinWidth) && ( (num_delayed_jet)>=2 || htSum >= 360 ) ) htSumRates_120timingOR360_2_emu->Fill(htSumLo+(bin*htSumBinWidth));
	if( (htSum) >= htSumLo+(bin*htSumBinWidth) && ( (num_delayed_jet)>=3 || htSum >= 360 ) ) htSumRates_120timingOR360_3_emu->Fill(htSumLo+(bin*htSumBinWidth));
	if( (htSum) >= htSumLo+(bin*htSumBinWidth) && ( (num_delayed_jet)>=4 || htSum >= 360 ) ) htSumRates_120timingOR360_4_emu->Fill(htSumLo+(bin*htSumBinWidth));
	if( (htSum) >= htSumLo+(bin*htSumBinWidth) && ( (num_delayed_jet)>=5 || htSum >= 360 ) ) htSumRates_120timingOR360_5_emu->Fill(htSumLo+(bin*htSumBinWidth));
	// use for ht120+timing rates in rate vs. eff plots
        if( (htSum) >= htSumLo+(bin*htSumBinWidth) && ( (num_delayed_jet)>=1 ) ) htSumRates_120timing_1_emu->Fill(htSumLo+(bin*htSumBinWidth));
        if( (htSum) >= htSumLo+(bin*htSumBinWidth) && ( (num_delayed_jet)>=2 ) ) htSumRates_120timing_2_emu->Fill(htSumLo+(bin*htSumBinWidth));
        if( (htSum) >= htSumLo+(bin*htSumBinWidth) && ( (num_delayed_jet)>=3 ) ) htSumRates_120timing_3_emu->Fill(htSumLo+(bin*htSumBinWidth));
        if( (htSum) >= htSumLo+(bin*htSumBinWidth) && ( (num_delayed_jet)>=4 ) ) htSumRates_120timing_4_emu->Fill(htSumLo+(bin*htSumBinWidth));
        if( (htSum) >= htSumLo+(bin*htSumBinWidth) && ( (num_delayed_jet)>=5 ) ) htSumRates_120timing_5_emu->Fill(htSumLo+(bin*htSumBinWidth));
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

    TDC_vs_amp->Write();
    amp_vs_TDC->Write();

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
    TimingBit_TPenergy->Write();

    // write histograms based on timing bit
    Delayed_2x2_MultHB_emu->Write();
    Prompt_2x2_MultHB_emu->Write();
    Prompt_Energy_2x2_MultHB_emu->Write();
    HTdistribution_trig_emu->Write();
    HTdistribution_emu->Write();

    JetPTdistribution_emu->Write();
    JetPTdistribution_trig_emu->Write();
    JetPTdistribution_trig120_emu->Write();

    htSumDistribution->Write();

    prompt_delayed_seed->Write();
    DeltaR_L1_delayed_seed_emu->Write();
    DeltaR_L1_prompt_seed_emu->Write();
    DeltaR_L1_delayed_hit_emu->Write();
    DeltaR_L1_prompt_hit_emu->Write();
    Mult_delayed_hit_emu->Write();
    Mult_prompt_hit_emu->Write();

    ctau_LLP->Write();
    ctau_LLP_trigger->Write();
    path_length->Write();
    path_length_energy->Write();
    PDGid_radius->Write();
    pion_radius_z->Write();
    delayed4x4seed_TDC_GeV_HB->Write();
    delayed4x4seed_TDC_GeV_HE->Write();
    delayed4x4seed_depth_TDC_HB->Write();
    delayed4x4seed_depth_TDC_HE->Write();

    htSumRates_original_emu->Scale(norm);
    // HT120+timing OR HT360
    htSumRates_120timingOR360_1_emu->Scale(norm);
    htSumRates_120timingOR360_2_emu->Scale(norm);
    htSumRates_120timingOR360_3_emu->Scale(norm);
    htSumRates_120timingOR360_4_emu->Scale(norm);
    htSumRates_120timingOR360_5_emu->Scale(norm);
    // HT120+timing
    htSumRates_120timing_1_emu->Scale(norm);
    htSumRates_120timing_2_emu->Scale(norm);
    htSumRates_120timing_3_emu->Scale(norm);
    htSumRates_120timing_4_emu->Scale(norm);
    htSumRates_120timing_5_emu->Scale(norm);

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

  //  std::cout << passed_calo_cluster_trig / totalEvents_HBdr05 * 100 << " passed calo trig, no HT cut / HB events" << std::endl;
  std::cout << passed_calo_cluster_trig_120 / totalEvents * 100 << " passed calo trig (delayed object), HT 120 cut / all events" << std::endl;
  //  std::cout << passed_calo_cluster_trig_120_2 / totalEvents * 100 << " passed calo trig (prompt TP veto), HT 120 cut / all events" << std::endl;
  std::cout << passed_calo_cluster_trig_120_3 / totalEvents * 100 << " passed calo trig (prompt 2x2 veto), HT 120 cut / all events" << std::endl;
  std::cout << passed_calo_cluster_trig / totalEvents * 100 << " passed calo trig (delayed object), no HT cut / all events" << std::endl;
  std::cout << passed_calo_cluster_trig  << " passed calo trig (delayed object), no HT cut" << std::endl;
  std::cout << passed4JetMult_HBHE_ht120_1 << " passed calo trig + HT120 OR HT360 (delayed object)" << std::endl;
  //  std::cout << passed4JetMult_HBHE_ht120_2 << " passed calo trig + HT120 OR HT360 (prompt TP veto)" << std::endl;
  std::cout << passed4JetMult_HBHE_ht120_3 << " passed calo trig + HT120 OR HT360 (prompt 2x2 veto)" << std::endl;
  std::cout << (passed4JetMult_HBHE_ht120_1 - passedHtSum360)*100 / totalEvents << " added efficiency HT 360 (delayed object)" << std::endl;
  //  std::cout << (passed4JetMult_HBHE_ht120_2 - passedHtSum360)*100 / totalEvents << " added efficiency HT 360 (prompt TP veto)" << std::endl;
  std::cout << (passed4JetMult_HBHE_ht120_3 - passedHtSum360)*100 / totalEvents << " added efficiency HT 360 (prompt 2x2 veto)" << std::endl;
  std::cout << passedHtSum360/totalEvents * 100 << " % passed HT360" << std::endl;
  std::cout << (passed4JetMult_HBHE_ht120_3) / passedHtSum360 << " integrated luminosity gain" << std::endl;
  std::cout << totalEvents << " all events" << std::endl;

  // saving efficiencies and rates in txt files to be read by rate vs eff plotting macros
  // signal efficiencies
  if ( (inputFile.substr(0,2) == "mh") ) {
    std::ofstream MultiplicityHits50ADC3ns_ht120_Signal;
    MultiplicityHits50ADC3ns_ht120_Signal.open(Form("MultiplicityHits50ADC3ns_ht120_Signal_%s.txt", inputFile.substr(0,15).c_str()),std::ios_base::trunc);
    std::cout << Form("MultiplicityHits50ADC3ns_ht120_Signal_%s.txt", inputFile.substr(0,15).c_str()) << std::endl;
    MultiplicityHits50ADC3ns_ht120_Signal << passed4JetMult_HBHE_ht120_1 / totalEvents << std::endl; // efficiency at HT 120+timing OR HT 360, delayed seed
    MultiplicityHits50ADC3ns_ht120_Signal << passed4JetMult_HBHE_ht120_2 / totalEvents << std::endl; // prompt TP veto
    MultiplicityHits50ADC3ns_ht120_Signal << passed4JetMult_HBHE_ht120_3 / totalEvents << std::endl; // prompt 2x2 veto
    MultiplicityHits50ADC3ns_ht120_Signal << passed4JetMult_HBHE_ht120_4 / totalEvents << std::endl;
    MultiplicityHits50ADC3ns_ht120_Signal << passed4JetMult_HBHE_ht120_5 / totalEvents << std::endl;
    MultiplicityHits50ADC3ns_ht120_Signal << passedHtSum360 / totalEvents << std::endl; // efficiency at HT 360
    MultiplicityHits50ADC3ns_ht120_Signal << "" << std::endl;
    MultiplicityHits50ADC3ns_ht120_Signal << "" << std::endl;
    //    MultiplicityHits50ADC3ns_ht120_Signal << "Added efficiency at nhit = 1,2,3,4,5 " << std::endl;
    MultiplicityHits50ADC3ns_ht120_Signal << (passed4JetMult_HBHE_ht120_1 - passedHtSum360)*100 / totalEvents << std::endl;
    MultiplicityHits50ADC3ns_ht120_Signal << (passed4JetMult_HBHE_ht120_2 - passedHtSum360)*100 / totalEvents << std::endl;
    MultiplicityHits50ADC3ns_ht120_Signal << (passed4JetMult_HBHE_ht120_3 - passedHtSum360)*100 / totalEvents << std::endl;
    MultiplicityHits50ADC3ns_ht120_Signal << (passed4JetMult_HBHE_ht120_4 - passedHtSum360)*100 / totalEvents << std::endl;
    MultiplicityHits50ADC3ns_ht120_Signal << (passed4JetMult_HBHE_ht120_5 - passedHtSum360)*100 / totalEvents << std::endl;
    MultiplicityHits50ADC3ns_ht120_Signal << "" << std::endl;
    MultiplicityHits50ADC3ns_ht120_Signal << "" << std::endl;
    //    MultiplicityHits50ADC3ns_ht120_Signal << "Efficiency at HT120 + timing, increasing number of delayed jets " << std::endl;
    MultiplicityHits50ADC3ns_ht120_Signal << passed_calo_cluster_trig_120 / totalEvents << std::endl; // delayed seed
    MultiplicityHits50ADC3ns_ht120_Signal << passed_calo_cluster_trig_120_2 / totalEvents << std::endl; // prompt TP veto
    MultiplicityHits50ADC3ns_ht120_Signal << passed_calo_cluster_trig_120_3 / totalEvents << std::endl; // prompt 2x2 veto
    MultiplicityHits50ADC3ns_ht120_Signal << passed_calo_cluster_trig_120_4 / totalEvents << std::endl;
    MultiplicityHits50ADC3ns_ht120_Signal << passed_calo_cluster_trig_120_5 / totalEvents << std::endl;
    MultiplicityHits50ADC3ns_ht120_Signal.close();
    std::ofstream Efficiency_HtBins_Signal;
    Efficiency_HtBins_Signal.open(Form("Efficiency_HtBins_Signal_%s.txt", inputFile.substr(0,15).c_str()),std::ios_base::trunc);
    if (totalEvent_ht120>0) Efficiency_HtBins_Signal << HBHE4Jet_inBins_ht120 / totalEvent_ht120 << std::endl; 
    else Efficiency_HtBins_Signal << 0 << std::endl;
    if (totalEvent_ht140>0) Efficiency_HtBins_Signal << HBHE4Jet_inBins_ht140 / totalEvent_ht140 << std::endl;
    else Efficiency_HtBins_Signal << 0 << std::endl;
    if (totalEvent_ht160>0) Efficiency_HtBins_Signal << HBHE4Jet_inBins_ht160 / totalEvent_ht160 << std::endl;
    else Efficiency_HtBins_Signal << 0 << std::endl;
    if (totalEvent_ht180>0) Efficiency_HtBins_Signal << HBHE4Jet_inBins_ht180 / totalEvent_ht180 << std::endl;
    else Efficiency_HtBins_Signal << 0 << std::endl;
    if (totalEvent_ht200>0) Efficiency_HtBins_Signal << HBHE4Jet_inBins_ht200 / totalEvent_ht200 << std::endl;
    else Efficiency_HtBins_Signal << 0 << std::endl;
    if (totalEvent_ht220>0) Efficiency_HtBins_Signal << HBHE4Jet_inBins_ht220 / totalEvent_ht220 << std::endl;
    else Efficiency_HtBins_Signal << 0 << std::endl;
    if (totalEvent_ht240>0) Efficiency_HtBins_Signal << HBHE4Jet_inBins_ht240 / totalEvent_ht240 << std::endl;
    else Efficiency_HtBins_Signal << 0 << std::endl;
    if (totalEvent_ht260>0) Efficiency_HtBins_Signal << HBHE4Jet_inBins_ht260 / totalEvent_ht260 << std::endl;
    else Efficiency_HtBins_Signal << 0 << std::endl;
    if (totalEvent_ht280>0) Efficiency_HtBins_Signal << HBHE4Jet_inBins_ht280 / totalEvent_ht280 << std::endl;
    else Efficiency_HtBins_Signal << 0 << std::endl;
    if (totalEvent_ht300>0) Efficiency_HtBins_Signal << HBHE4Jet_inBins_ht300 / totalEvent_ht300 << std::endl;
    else Efficiency_HtBins_Signal << 0 << std::endl;
    if (totalEvent_ht320>0) Efficiency_HtBins_Signal << HBHE4Jet_inBins_ht320 / totalEvent_ht320 << std::endl;
    else Efficiency_HtBins_Signal << 0 << std::endl;
    if (totalEvent_ht340>0) Efficiency_HtBins_Signal << HBHE4Jet_inBins_ht340 / totalEvent_ht340 << std::endl;
    else Efficiency_HtBins_Signal << 0 << std::endl;
    if (totalEvent_ht360>0) Efficiency_HtBins_Signal << HBHE4Jet_inBins_ht360 / totalEvent_ht360 << std::endl;
    else Efficiency_HtBins_Signal << 0 << std::endl;
    if (totalEvent_ht380>0) Efficiency_HtBins_Signal << HBHE4Jet_inBins_ht380 / totalEvent_ht380 << std::endl;
    else Efficiency_HtBins_Signal << 0 << std::endl;
    Efficiency_HtBins_Signal << HEHB4Jet_ht120_wTiming / totalEvents << std::endl;
    Efficiency_HtBins_Signal << "Efficiency = events with HT>120 + timing hits over 2 / total events, reported in last line " << std::endl;
    Efficiency_HtBins_Signal.close();
  }
  // background efficiencies 
  if (inputFile.substr(0,3) == "QCD" ) {
    std::ofstream MultiplicityHits50ADC3ns_ht120_Background;
    MultiplicityHits50ADC3ns_ht120_Background.open("MultiplicityHits50ADC3ns_ht120_Background.txt");
    MultiplicityHits50ADC3ns_ht120_Background << passed4JetMult_HBHE_ht120_1 / totalEvents << std::endl; // efficiency at HT 120+timing OR HT 360  
    MultiplicityHits50ADC3ns_ht120_Background << passed4JetMult_HBHE_ht120_2 / totalEvents << std::endl;
    MultiplicityHits50ADC3ns_ht120_Background << passed4JetMult_HBHE_ht120_3 / totalEvents << std::endl;
    MultiplicityHits50ADC3ns_ht120_Background << passed4JetMult_HBHE_ht120_4 / totalEvents << std::endl;
    MultiplicityHits50ADC3ns_ht120_Background << passed4JetMult_HBHE_ht120_5 / totalEvents << std::endl;
    MultiplicityHits50ADC3ns_ht120_Background << passedHtSum360 / totalEvents << std::endl; // efficiency at HT 360  
    MultiplicityHits50ADC3ns_ht120_Background << "" << std::endl;
    MultiplicityHits50ADC3ns_ht120_Background << "" << std::endl;
    //    MultiplicityHits50ADC3ns_ht120_Background << "Added efficiency at nhit = 1,2,3,4,5 " << std::endl;
    MultiplicityHits50ADC3ns_ht120_Background << (passed4JetMult_HBHE_ht120_1 - passedHtSum360)*100 / totalEvents << std::endl;
    MultiplicityHits50ADC3ns_ht120_Background << (passed4JetMult_HBHE_ht120_2 - passedHtSum360)*100 / totalEvents << std::endl;
    MultiplicityHits50ADC3ns_ht120_Background << (passed4JetMult_HBHE_ht120_3 - passedHtSum360)*100 / totalEvents << std::endl;
    MultiplicityHits50ADC3ns_ht120_Background << (passed4JetMult_HBHE_ht120_4 - passedHtSum360)*100 / totalEvents << std::endl;
    MultiplicityHits50ADC3ns_ht120_Background << (passed4JetMult_HBHE_ht120_5 - passedHtSum360)*100 / totalEvents << std::endl;
    MultiplicityHits50ADC3ns_ht120_Background << "" << std::endl;
    MultiplicityHits50ADC3ns_ht120_Background << "" << std::endl;
    //    MultiplicityHits50ADC3ns_ht120_Background << "Efficiency at HT120 + timing, increasing number of delayed jets " << std::endl;
    MultiplicityHits50ADC3ns_ht120_Background << passed_calo_cluster_trig_120 / totalEvents << std::endl;
    MultiplicityHits50ADC3ns_ht120_Background << passed_calo_cluster_trig_120_2 / totalEvents << std::endl;
    MultiplicityHits50ADC3ns_ht120_Background << passed_calo_cluster_trig_120_3 / totalEvents << std::endl;
    MultiplicityHits50ADC3ns_ht120_Background << passed_calo_cluster_trig_120_4 / totalEvents << std::endl;
    MultiplicityHits50ADC3ns_ht120_Background << passed_calo_cluster_trig_120_5 / totalEvents << std::endl;
    MultiplicityHits50ADC3ns_ht120_Background.close();
  }
  // neutrino gun rates
  if (inputFile.substr(0,11) == "RelValNuGun" ) {
  //  if (inputFile.substr(0,7) == "MinBias" ) {
    std::cout << "htSum_original120 = " << htSumRates_original_emu->GetBinContent(htSumRates_original_emu->GetXaxis()->FindBin(120)) << std::endl;
    std::cout << "htSum_original360 = " << htSumRates_original_emu->GetBinContent(htSumRates_original_emu->GetXaxis()->FindBin(360)) << std::endl;
    std::cout << "htSum_wtiming120 2 hits = " << htSumRates_emu->GetBinContent(htSumRates_emu->GetXaxis()->FindBin(120)) << std::endl;
    std::cout << "ratio change = timing at 120 / original 360 = " << ( htSumRates_emu->GetBinContent(htSumRates_emu->GetXaxis()->FindBin(120))) / (htSumRates_original_emu->GetBinContent(htSumRates_original_emu->GetXaxis()->FindBin(360))) << std::endl;
    std::cout << "ratio change = original 360 / timing at 120 = " << (htSumRates_original_emu->GetBinContent(htSumRates_original_emu->GetXaxis()->FindBin(360))) / ( htSumRates_emu->GetBinContent(htSumRates_emu->GetXaxis()->FindBin(120))) << std::endl;
    std::cout << "ratio change = original 120 / timing at 120 = " << (htSumRates_original_emu->GetBinContent(htSumRates_original_emu->GetXaxis()->FindBin(120))) / ( htSumRates_emu->GetBinContent(htSumRates_emu->GetXaxis()->FindBin(120))) << std::endl;
    int htSum_120timingOR360_1_120 = htSumRates_120timingOR360_1_emu->GetBinContent(htSumRates_120timingOR360_1_emu->GetXaxis()->FindBin(120));
    int htSum_120timingOR360_2_120 = htSumRates_120timingOR360_2_emu->GetBinContent(htSumRates_120timingOR360_2_emu->GetXaxis()->FindBin(120));
    int htSum_120timingOR360_3_120 = htSumRates_120timingOR360_3_emu->GetBinContent(htSumRates_120timingOR360_3_emu->GetXaxis()->FindBin(120));
    int htSum_120timingOR360_4_120 = htSumRates_120timingOR360_4_emu->GetBinContent(htSumRates_120timingOR360_4_emu->GetXaxis()->FindBin(120));
    int htSum_120timingOR360_5_120 = htSumRates_120timingOR360_5_emu->GetBinContent(htSumRates_120timingOR360_5_emu->GetXaxis()->FindBin(120));

    int htSum_120timing_1_120 = htSumRates_120timing_1_emu->GetBinContent(htSumRates_120timing_1_emu->GetXaxis()->FindBin(120));
    int htSum_120timing_2_120 = htSumRates_120timing_2_emu->GetBinContent(htSumRates_120timing_2_emu->GetXaxis()->FindBin(120));
    int htSum_120timing_3_120 = htSumRates_120timing_3_emu->GetBinContent(htSumRates_120timing_3_emu->GetXaxis()->FindBin(120));
    int htSum_120timing_4_120 = htSumRates_120timing_4_emu->GetBinContent(htSumRates_120timing_4_emu->GetXaxis()->FindBin(120));
    int htSum_120timing_5_120 = htSumRates_120timing_5_emu->GetBinContent(htSumRates_120timing_5_emu->GetXaxis()->FindBin(120));

    int htSum_original_360 = htSumRates_original_emu->GetBinContent(htSumRates_original_emu->GetXaxis()->FindBin(360));
    std::ofstream NuGunRates;
    NuGunRates.open("NuGunRates_360_OR_120timing.txt");
    NuGunRates << htSum_120timingOR360_1_120 << std::endl; // rate at HT 120 for the OR of the two triggers
    NuGunRates << htSum_120timingOR360_2_120 << std::endl;
    NuGunRates << htSum_120timingOR360_3_120 << std::endl;
    NuGunRates << htSum_120timingOR360_4_120 << std::endl;
    NuGunRates << htSum_120timingOR360_5_120 << std::endl;
    NuGunRates << htSum_original_360 << std::endl; // rate at HT 360 without timing cuts
    NuGunRates << htSumRates_original_emu->GetBinContent(htSumRates_original_emu->GetXaxis()->FindBin(120)) << std::endl; // rate at HT 120 without timing cuts
    NuGunRates.close();
    NuGunRates.open("NuGunRates.txt");
    NuGunRates << htSum_120timing_1_120 << std::endl; // rate at HT 120 for the L1_HT120+timing trigger
    NuGunRates << htSum_120timing_2_120 << std::endl;
    NuGunRates << htSum_120timing_3_120 << std::endl;
    NuGunRates << htSum_120timing_4_120 << std::endl;
    NuGunRates << htSum_120timing_5_120 << std::endl;
    NuGunRates << htSum_original_360 << std::endl; // rate at HT 360 without timing cuts  
    NuGunRates << htSumRates_original_emu->GetBinContent(htSumRates_original_emu->GetXaxis()->FindBin(120)) << std::endl; // rate at HT 120 without timing cuts
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
