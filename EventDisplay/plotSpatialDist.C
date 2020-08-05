#include <vector>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TBenchmark.h"
#include "TRandom.h"
#include "TSystem.h"
#include <TLorentzVector.h>
#include <TStyle.h>
#include "TLegend.h"
#include "TEllipse.h"
#include "TPaveText.h"
#include "TLine.h"
#include "TEllipse.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeDataFormat_eventdisplay.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisCaloTPDataFormat.h"
#include <sstream>
//#include "Math/VectorUtil_Cint.h"

#ifdef __MAKECINT__
#pragma link C++ class vector<float>+;
#endif

#ifdef __MAKECINT__
#pragma link C++ class vector<TLorentzVector>+;
#endif

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

void DrawRegionLines(){
  std::vector<TLine*> regionLines;
  float etaValues[17] = { -3, -2.088, -1.74, -1.392, -1.044, -0.696, -0.348, 0,
			  0.348, 0.696, 1.044, 1.392, 1.74, 2.088, 3 };
  float phiValues[18] = {-2.965, -2.617, -2.268, -1.919, -1.570, -1.221, -0.872, -0.523, -0.174, 
			 0.174, 0.523, 0.872, 1.221, 1.570, 1.919, 2.268, 2.617, 2.965};

  //eta lines
  for(int i = 0; i < 17; i++){
    TLine * line = new TLine(etaValues[i], -3.2, etaValues[i], 3.2); 
    line->SetLineColor(kBlue-7);
    line->SetLineStyle(1);
    regionLines.push_back(line);
  }

  //phi lines
  for(int i = 0; i < 18; i++){
    TLine * line = new TLine(-3, phiValues[i], 3, phiValues[i]); 
    line->SetLineColor(kBlue-7);
    line->SetLineStyle(1);
    regionLines.push_back(line);
  }
  
  for(size_t j = 0; j < regionLines.size(); j++){
    regionLines.at(j)->Draw();
  }
}

void DrawTowerLines(){
  std::vector<TLine*> TowerLines;
  float etaValues[59] = { -2.913, -2.739, -2.565, -2.391, -2.217, -2.088, -2.001, -1.914, -1.827, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957, -0.87, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.87, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.74, 1.827, 1.914, 2.001, 2.088, 2.217, 2.391, 2.565, 2.739, 2.913};

  float phiValues[73] = {-3.142, -3.054, -2.967, -2.880, -2.793, -2.705, -2.618, -2.531, -2.443, -2.356, -2.269, -2.182, -2.094, -2.007, -1.920, -1.833, -1.745, -1.658, -1.571, -1.484, -1.396, -1.309, -1.222, -1.134, -1.047, -0.960, -0.873, -0.785, -0.698, -0.611, -0.524, -0.436, -0.349, -0.262, -0.175, -0.087,
			 0.000, 0.087, 0.175, 0.262, 0.349, 0.436, 0.524, 0.611, 0.698, 0.785, 0.873, 0.960, 1.047, 1.134, 1.222, 1.309, 1.396, 1.484, 1.571, 1.658, 1.745, 1.833, 1.920, 2.007, 2.094, 2.182, 2.269, 2.356, 2.443, 2.531, 2.618, 2.705, 2.793, 2.880, 2.967, 3.054, 3.142};
  
  //eta lines
  for(int i = 0; i < 59; i++){
    TLine * line = new TLine(etaValues[i], -3.2, etaValues[i], 3.2); 
    line->SetLineColor(kGray);
    line->SetLineStyle(3);
    TowerLines.push_back(line);
  }
  
  //phi lines
  for(int i = 0; i < 73; i++){
    TLine * line = new TLine(-3, phiValues[i], 3, phiValues[i]); 
    line->SetLineColor(kGray);
    line->SetLineStyle(3);
    TowerLines.push_back(line);
  }
  
  for(size_t j = 0; j < TowerLines.size(); j++){
    TowerLines.at(j)->Draw();
  }
}

//void plotSpatialDist(int iEvent, const char* file){
void plotSpatialDist(int iEvent){
  gStyle->SetOptStat(0);
  TFile *f = TFile::Open("L1Ntuple_mh125_mx50_pl1000.root","READ");
  if (!f) { return; }
  
  TTree *t = (TTree*) f->Get("l1EventTree/L1EventTree"); // saves event info, in event branch
  TTree *t2 = (TTree*) f->Get("l1UpgradeEmuTree/L1UpgradeTree"); // has info about L1 jets, in L1Upgrade branch   
  TTree *t3 = (TTree*) f->Get("l1CaloTowerEmuTree/L1CaloTowerTree"); // info on HCAL TDC and energy, in CaloTP branch

  L1Analysis::L1AnalysisL1UpgradeDataFormat *vL1Jets = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  L1Analysis::L1AnalysisCaloTPDataFormat *vHcalTPdepth = new L1Analysis::L1AnalysisCaloTPDataFormat();
  ULong64_t event =0; // declare as a unsigned long 64 bit integer
  
  // Create a new canvas.
  TCanvas *c1 = new TCanvas("c1","eta vs phi",200,10,700,700);
  c1->SetFillColor(0);
  c1->GetFrame()->SetFillColor(0);
  c1->GetFrame()->SetBorderSize(6);
  c1->GetFrame()->SetBorderMode(-1);

  const Int_t kUPDATE = 1000;
  TBranch *bEvent = 0;
  TBranch *bL1Jets = 0;
  TBranch *bHcalTPdepth = 0;

  t->SetBranchAddress("event",&event,&bEvent);
  t2->SetBranchAddress("L1Upgrade",&vL1Jets,&bL1Jets);
  t3->SetBranchAddress("CaloTP",&vHcalTPdepth,&bHcalTPdepth);

  // Create one histograms
  TH1F *h                = new TH1F("h","This is the eta distribution",100,-4,4);
  TH2F *h2L1Jets        = new TH2F("h2L1Jets","h2 title",68,-3,3,72,-3.142,3.142);
  TH2F *h2HcalTPdepth        = new TH2F("h2HcalTPdepth","h2 title",68,-3,3,72,-3.142,3.142);

  std::vector<TPaveText*> hcalTPdepthText;

  h->SetFillColor(48);
  int i = iEvent;
  Long64_t tentry = t->LoadTree(i);
  std::cout<<"i "<<i<< " tentry "<< tentry << std::endl;
  bEvent->GetEntry(tentry);
  bL1Jets->GetEntry(tentry);
  bHcalTPdepth->GetEntry(tentry);

  //get the event number
  char name[30];
  sprintf(name,"Event %llu",event);
  std::cout<<event<<std::endl;
  std::cout<<name<<std::endl;
  TH2F *h2 = new TH2F("h2",name,68,-3,3,72,-3.142,3.142);

  int k = 0;

  for (int j = 0; j < vHcalTPdepth->nHCALTP; ++j) {
    double eta = etaVal(vHcalTPdepth->hcalTPieta[j]);
    double phi = phiVal(vHcalTPdepth->hcalTPiphi[j]);
    double pt  = vHcalTPdepth->hcalTPet[j];
    h2HcalTPdepth->Fill(eta, phi, pt);

    if(pt>10){
      std::cout<<"vHcalTPdepth->at(j).Pt() "<<vHcalTPdepth->hcalTPet[j]
	       <<" eta "<<etaVal(vHcalTPdepth->hcalTPieta[j])
	       <<" phi "<<phiVal(vHcalTPdepth->hcalTPiphi[j])<<std::endl;
      std::ostringstream strs;
      strs << pt;
      std::string text = strs.str();
      eta += 0.01;
      phi += 0.01;
      TPaveText *tempText = new TPaveText( eta, phi, eta+0.1, phi+0.1 );
      tempText->AddText(text.c_str());
      tempText->SetFillColor(0);
      tempText->SetLineColor(0);
      tempText->SetShadowColor(0);
      tempText->SetTextColor(kBlue);
      hcalTPdepthText.push_back(tempText);
    }
  }  
  
  for (UInt_t j = 0; j < vL1Jets->nJets; ++j) {
    double eta = vL1Jets->jetEta[j];
    double phi = vL1Jets->jetPhi[j];
    double pt  = vL1Jets->jetEt[j];
    h2L1Jets->Fill(eta, phi, pt);
  }  
  
  h2->GetXaxis()->SetAxisColor(17);
  h2->GetYaxis()->SetAxisColor(17);
  
  h->Draw(); 
  h2->Draw("BOX");
  
  DrawRegionLines();
  DrawTowerLines();
  h2HcalTPdepth->SetFillColor(kMagenta);
  h2HcalTPdepth->Draw("SAME BOX");
  h2L1Jets->SetFillStyle(3001);
  h2L1Jets->SetFillColorAlpha(kSpring, 0.75);
  h2L1Jets->Draw("SAME BOX");
  
  double x1=-99., x2=-99., x3=-99., x4=-99., y1=-99., y2=-99., y3=-99., y4=-99.;
  for (UInt_t j = 0; j < vL1Jets->nJets; ++j) {
    double eta = vL1Jets->jetEta[j];
    double phi = vL1Jets->jetPhi[j];
    if(j == 0) { x1 = eta-0.8; x2 = eta+0.8; y1 = phi-0.8; y2 = phi+0.8; }
    if(j == 1) { x3 = eta-0.8; x4 = eta+0.8; y3 = phi-0.8; y4 = phi+0.8; }
    TEllipse *circ = new TEllipse(eta,phi,.8,.8);
    circ->SetFillStyle(0);
    circ->SetLineStyle(2);
    circ->SetLineColor(kViolet+2);
    circ->Draw("SAME");
  }
  
  float xR=0.8;
  TLegend *l = new TLegend(xR,0.8,xR+0.2,1.0);
  l->AddEntry(h2,"Regions","F");  
  l->AddEntry(h2HcalTPdepth,"HCAL TPGs","F");
  l->AddEntry(h2L1Jets,"L1 jets","F");
  l->Draw();
  h2->GetXaxis()->SetTitle("eta");
  h2->GetYaxis()->SetTitle("phi");
  
  for (UInt_t j = 0; j < hcalTPdepthText.size(); ++j) {
    hcalTPdepthText.at(j)->Draw("SAME");
  }
  
  char saveFile[100];
  sprintf(saveFile,"/afs/cern.ch/work/g/gkopp/HCAL_Trigger/L1Ntuples/HCAL_TP_TimingBitEmulator/CMSSW_10_6_0/src/HcalTrigger/Validation/EventDisplay/Event-%llu-test_trigger.png",event);
  c1->SaveAs(saveFile);
  
  if(x1 > -99.){
    TCanvas *c2 = new TCanvas("c2","leading jet",200,10,700,700);
    c2->SetFillColor(0);
    c2->GetFrame()->SetFillColor(0);
    c2->GetFrame()->SetBorderSize(6);
    c2->GetFrame()->SetBorderMode(-1);
    c2->cd();
    
    h->Draw();
    h2->GetXaxis()->SetRangeUser(x1,x2);
    h2->GetYaxis()->SetRangeUser(y1,y2);
    h2->Draw("BOX");
    DrawRegionLines();
    DrawTowerLines();
    h2HcalTPdepth->Draw("SAME BOX");
    h2L1Jets->Draw("SAME BOX");
    l->Draw();
    char saveFile1[100];
    sprintf(saveFile1,"/afs/cern.ch/work/g/gkopp/HCAL_Trigger/L1Ntuples/HCAL_TP_TimingBitEmulator/CMSSW_10_6_0/src/HcalTrigger/Validation/EventDisplay/Event-%llu-test-jet1_trigger.png",event);
    c2->SaveAs(saveFile1);
  }
  
  if(x3 > -99.){
    TCanvas *c3 = new TCanvas("c3","sub-leading jet",200,10,700,700);
    c3->SetFillColor(0);
    c3->GetFrame()->SetFillColor(0);
    c3->GetFrame()->SetBorderSize(6);
    c3->GetFrame()->SetBorderMode(-1);
    c3->cd();
    
    h->Draw();
    h2->GetXaxis()->SetRangeUser(x3,x4);
    h2->GetYaxis()->SetRangeUser(y3,y4);
    h2->Draw("BOX");
    DrawRegionLines();
    DrawTowerLines();
    h2HcalTPdepth->Draw("SAME BOX");
    h2L1Jets->Draw("SAME BOX");
    l->Draw();
    char saveFile1[100];
    sprintf(saveFile1,"/afs/cern.ch/work/g/gkopp/HCAL_Trigger/L1Ntuples/HCAL_TP_TimingBitEmulator/CMSSW_10_6_0/src/HcalTrigger/Validation/EventDisplay/Event-%llu-test-jet2_trigger.png",event); 
    c3->SaveAs(saveFile1);
  }
  
}
