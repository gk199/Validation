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
#include <sstream>
//#include "Math/VectorUtil_Cint.h"

#ifdef __MAKECINT__
#pragma link C++ class vector<float>+;
#endif

#ifdef __MAKECINT__
#pragma link C++ class vector<TLorentzVector>+;
#endif

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

//void plotVector(int iEvent, const char* file){
void plotVector(int iEvent){
  gStyle->SetOptStat(0);
  //const char* tfile = file;
  //TFile *f = TFile::Open(tfile,"READ");
  TFile *f = TFile::Open("l1TNtuple-ggHBB.root","READ");
  if (!f) { return; }
  
  TTree *t = (TTree*) f->Get("l1NtupleProducer/efficiencyTree");
  std::vector<TLorentzVector> *vpx             = 0;
  std::vector<TLorentzVector> *vEcalTpgs       = 0;
  std::vector<TLorentzVector> *vHcalTpgs       = 0;
  std::vector<TLorentzVector> *vCaloClusters   = 0;
  std::vector<TLorentzVector> *vL1Jets         = 0;
  std::vector<TLorentzVector> *vAK8Jets        = 0;
  std::vector<TLorentzVector> *vSubJets        = 0;
  int event =0;
  
  // Create a new canvas.
  TCanvas *c1 = new TCanvas("c1","eta vs phi",200,10,700,700);
  c1->SetFillColor(0);
  c1->GetFrame()->SetFillColor(0);
  c1->GetFrame()->SetBorderSize(6);
  c1->GetFrame()->SetBorderMode(-1);

  const Int_t kUPDATE = 1000;
  TBranch *bEvent = 0;
  TBranch *bvpx = 0;
  TBranch *bEcalTpgs = 0;
  TBranch *bHcalTpgs = 0;
  TBranch *bCaloClusters = 0;
  TBranch *bL1Jets = 0;
  TBranch *bAK8Jets = 0;
  TBranch *bSubJets = 0;

  t->SetBranchAddress("event",&event,&bEvent);
  t->SetBranchAddress("allRegions",&vpx,&bvpx);
  t->SetBranchAddress("ecalTPGs",&vEcalTpgs,&bEcalTpgs);
  t->SetBranchAddress("hcalTPGs",&vHcalTpgs,&bHcalTpgs);
  t->SetBranchAddress("caloClusters",&vCaloClusters,&bCaloClusters);
  t->SetBranchAddress("l1Jets",&vL1Jets,&bL1Jets);
  t->SetBranchAddress("ak8Jets",&vAK8Jets,&bAK8Jets);
  t->SetBranchAddress("subJets",&vSubJets,&bSubJets);

  // Create one histograms
  TH1F *h                = new TH1F("h","This is the eta distribution",100,-4,4);
  //TH2F *h2               = new TH2F("h2","Event 2988846758",68,-3,3,72,-3.142,3.142);
  TH2F *h2EcalTpgs       = new TH2F("h2EcalTpgs","h2 title",68,-3,3,72,-3.142,3.142);
  TH2F *h2HcalTpgs       = new TH2F("h2HcalTpgs","h2 title",68,-3,3,72,-3.142,3.142);
  TH2F *h2CaloClusters   = new TH2F("h2CaloClusters","h2 title",68,-3,3,72,-3.142,3.142);
  TH2F *h2L1Jets         = new TH2F("h2L1Jets","h2 title",68,-3,3,72,-3.142,3.142);
  TH2F *h2AK8Jets        = new TH2F("h2AK8Jets","h2 title",68,-3,3,72,-3.142,3.142);
  TH2F *h2SubJets        = new TH2F("h2SubJets","h2 title",68,-3,3,72,-3.142,3.142);

  std::vector<TPaveText*> ecalTpgText;
  std::vector<TPaveText*> hcalTpgText;
  std::vector<TPaveText*> caloClusterText;

  h->SetFillColor(48);
  int i = iEvent;
  Long64_t tentry = t->LoadTree(i);
  std::cout<<"i "<<i<< " tentry "<< tentry << std::endl;
  bEvent->GetEntry(tentry);
  bvpx->GetEntry(tentry);
  bEcalTpgs->GetEntry(tentry);
  bHcalTpgs->GetEntry(tentry);
  bCaloClusters->GetEntry(tentry);
  bL1Jets->GetEntry(tentry);
  bAK8Jets->GetEntry(tentry);
  bSubJets->GetEntry(tentry);

  //get the event number
  char name[30];
  sprintf(name,"Event %u",event);
  std::cout<<event<<std::endl;
  std::cout<<name<<std::endl;
  TH2F *h2 = new TH2F("h2",name,68,-3,3,72,-3.142,3.142);

  int k = 0;
  for (UInt_t j = 0; j < vpx->size(); ++j) {
    if(vpx->at(j).Pt()>0)
      h->Fill(vpx->at(j).Eta());
    
    if(vpx->at(j).Pt()>1){
      h2->Fill(vpx->at(j).Eta(),vpx->at(j).Phi(),vpx->at(j).Pt());
      k++; 
      if(vpx->at(j).Pt()>10)
	std::cout<<"region eta "<<vpx->at(j).Eta()<<" region phi "<<vpx->at(j).Phi()<< " pt "<<vpx->at(j).Pt()/2 <<" GeV"<<std::endl;
    }
  }
  std::cout<<"n Regions above pt 0 :"<<k<<std::endl;

  for (UInt_t j = 0; j < vEcalTpgs->size(); ++j) {
    double eta = vEcalTpgs->at(j).Eta();
    double phi = vEcalTpgs->at(j).Phi();
    double pt  = vEcalTpgs->at(j).Pt();
    h2EcalTpgs->Fill(eta, phi, pt);
  
    if(pt>10){
      std::cout<<"vEcalTpgs->at(j).Pt() "<<vEcalTpgs->at(j).Pt()
	       <<" eta "<<vEcalTpgs->at(j).Eta()
	       <<" phi "<<vEcalTpgs->at(j).Phi()<<std::endl;  
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
      ecalTpgText.push_back(tempText);
    }
  }
  
  for (UInt_t j = 0; j < vHcalTpgs->size(); ++j) {
    double eta = vHcalTpgs->at(j).Eta();
    double phi = vHcalTpgs->at(j).Phi();
    double pt  = vHcalTpgs->at(j).Pt();
    h2HcalTpgs->Fill(eta, phi, pt);

    if(pt>10){
      std::cout<<"vHcalTpgs->at(j).Pt() "<<vHcalTpgs->at(j).Pt()
	       <<" eta "<<vHcalTpgs->at(j).Eta()
	       <<" phi "<<vHcalTpgs->at(j).Phi()<<std::endl;
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
      hcalTpgText.push_back(tempText);
    }
  }  
  
  for (UInt_t j = 0; j < vCaloClusters->size(); ++j) {
    double eta = vCaloClusters->at(j).Eta();
    double phi = vCaloClusters->at(j).Phi();
    double pt  = vCaloClusters->at(j).Pt();
    h2CaloClusters->Fill(eta, phi, pt);
    
    if(pt>10){
      std::cout<<"vCaloClusters->at(j).Pt() "<<vCaloClusters->at(j).Pt()
	       <<" eta "<<vCaloClusters->at(j).Eta()
	       <<" phi "<<vCaloClusters->at(j).Phi()<<std::endl;
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
      caloClusterText.push_back(tempText);
    }
  }
  
  for (UInt_t j = 0; j < vL1Jets->size(); ++j) {
    double eta = vL1Jets->at(j).Eta();
    double phi = vL1Jets->at(j).Phi();
    double pt  = vL1Jets->at(j).Pt();
    h2L1Jets->Fill(eta, phi, pt);
  }  
  
  for (UInt_t j = 0; j < vAK8Jets->size(); ++j) {
    double eta = vAK8Jets->at(j).Eta();
    double phi = vAK8Jets->at(j).Phi();
    double pt  = vAK8Jets->at(j).Pt();
    h2AK8Jets->Fill(eta, phi, pt);
  }
  
  for (UInt_t j = 0; j < vSubJets->size(); ++j) {
    double eta = vSubJets->at(j).Eta();
    double phi = vSubJets->at(j).Phi();
    double pt  = vSubJets->at(j).Pt();
    h2SubJets->Fill(eta, phi, pt);
  } 
  
  h2->GetXaxis()->SetAxisColor(17);
  h2->GetYaxis()->SetAxisColor(17);
  
  h->Draw(); 
  h2->Draw("BOX");
  
  DrawRegionLines();
  DrawTowerLines();
  h2HcalTpgs->SetFillColor(kMagenta);
  h2HcalTpgs->Draw("SAME BOX");
  h2EcalTpgs->SetFillColorAlpha(kRed, 0.75);
  h2EcalTpgs->Draw("SAME BOX");
  h2CaloClusters->SetFillStyle(3144);
  h2CaloClusters->SetFillColorAlpha(kOrange, 0.75);
  h2CaloClusters->Draw("SAME BOX");
  h2L1Jets->SetFillStyle(3001);
  h2L1Jets->SetFillColorAlpha(kSpring, 0.75);
  h2L1Jets->Draw("SAME BOX");
  h2AK8Jets->SetFillStyle(3001);
  h2AK8Jets->SetFillColorAlpha(kViolet+2, 0.75);
  h2AK8Jets->Draw("SAME BOX");
  h2SubJets->SetFillStyle(3001);
  h2SubJets->SetFillColorAlpha(kAzure+10, 0.75);
  h2SubJets->Draw("SAME BOX");
  
  double x1=-99., x2=-99., x3=-99., x4=-99., y1=-99., y2=-99., y3=-99., y4=-99.;
  for (UInt_t j = 0; j < vAK8Jets->size(); ++j) {
    double eta = vAK8Jets->at(j).Eta();
    double phi = vAK8Jets->at(j).Phi();
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
  l->AddEntry(h2EcalTpgs,"ECAL TPGs","F");
  l->AddEntry(h2HcalTpgs,"HCAL TPGs","F");
  l->AddEntry(h2CaloClusters,"caloClusters","F");
  l->AddEntry(h2L1Jets,"L1 jets","F");
  l->AddEntry(h2AK8Jets,"AK8 jets","F");
  l->AddEntry(h2SubJets,"Sub jets","F");
  l->Draw();
  h2->GetXaxis()->SetTitle("eta");
  h2->GetYaxis()->SetTitle("phi");
  
  for (UInt_t j = 0; j < ecalTpgText.size(); ++j) {
    ecalTpgText.at(j)->Draw("SAME");
  }
  
  for (UInt_t j = 0; j < hcalTpgText.size(); ++j) {
    hcalTpgText.at(j)->Draw("SAME");
  }
  
  for (UInt_t j = 0; j < caloClusterText.size(); ++j) {
    caloClusterText.at(j)->Draw("SAME");
  }
  
  char saveFile[100];
  sprintf(saveFile,"/afs/cern.ch/work/g/gkopp/HCAL_Trigger/L1Ntuples/HCAL_TP_TimingBitEmulator/CMSSW_10_6_0/src/HcalTrigger/Validation/EventDisplay/Event-%u-test.png",event);//"/afs/cern.ch/work/p/pdas/www/Run3Ntuplizer/EventDisplays/Event-%u-test.png",event);
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
    h2HcalTpgs->Draw("SAME BOX");
    h2EcalTpgs->Draw("SAME BOX");
    h2CaloClusters->Draw("SAME BOX");
    h2L1Jets->Draw("SAME BOX");
    h2AK8Jets->Draw("SAME BOX");
    h2SubJets->Draw("SAME BOX");
    l->Draw();
    char saveFile1[100];
    sprintf(saveFile1,"/afs/cern.ch/work/g/gkopp/HCAL_Trigger/L1Ntuples/HCAL_TP_TimingBitEmulator/CMSSW_10_6_0/src/HcalTrigger/Validation/EventDisplay/Event-%u-test-jet1.png",event);//"/afs/cern.ch/work/p/pdas/www/Run3Ntuplizer/EventDisplays/Event-%u-test-jet1.png",event);
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
    h2HcalTpgs->Draw("SAME BOX");
    h2EcalTpgs->Draw("SAME BOX");
    h2CaloClusters->Draw("SAME BOX");
    h2L1Jets->Draw("SAME BOX");
    h2AK8Jets->Draw("SAME BOX");
    h2SubJets->Draw("SAME BOX");
    l->Draw();
    char saveFile1[100];
    sprintf(saveFile1,"/afs/cern.ch/work/g/gkopp/HCAL_Trigger/L1Ntuples/HCAL_TP_TimingBitEmulator/CMSSW_10_6_0/src/HcalTrigger/Validation/EventDisplay/Event-%u-test-jet2.png",event); //"/afs/cern.ch/work/p/pdas/www/Run3Ntuplizer/EventDisplays/Event-%u-test-jet2.png",event);
    c3->SaveAs(saveFile1);
  }
  
}
