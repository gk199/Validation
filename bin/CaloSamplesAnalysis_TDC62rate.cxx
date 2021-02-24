#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TEfficiency.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TROOT.h"
#include "TGraphAsymmErrors.h"

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

#include "DataFormats/Common/interface/EDProductGetter.h" // mustBeNonZero
#include "CalibFormats/CaloObjects/interface/CaloSamples.h"
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

int main() {
  //  TFile *f = new TFile("/afs/cern.ch/work/g/gkopp/MC_GenProduction/PionGun/CMSSW_11_0_2/src/SinglePion211_E10_eta1phi0_PU00_tdc9pt35_step1.root");
  //  TFile *f = new TFile("/afs/cern.ch/work/g/gkopp/MC_GenProduction/PionGun/CMSSW_11_0_2/src/SinglePion211_E10_eta1phi0_PU00_tdc74pt8_timephase6_TDCflat_step1.root");
  //  TFile *f = new TFile("/afs/cern.ch/work/g/gkopp/MC_GenProduction/PionGun/CMSSW_11_0_2/src/SinglePion211_E10_eta1phi0_PU00_tdc149pt6_timephase6_TDCflat_step1.root");
  TFile *f = new TFile("/afs/cern.ch/work/g/gkopp/MC_GenProduction/PionGun/CMSSW_11_0_2/src/Injected_Energy_0pt001_to_1pt3_HB111_tdc74pt8_tof_step1_TDCflat.root");
  //  TFile *f = new TFile("file:/eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/112X_TDC74pt8/QCD_Pt-15to7000_TuneCP5_Flat_13TeV_pythia8-digi_noPU.root");

  TTreeReader myReader("Events",f);
  TTreeReaderValue<std::vector<CaloSamples>> AllCaloSamples(myReader, "CaloSampless_mix_HcalSamples_HLT.obj");
  TTreeReaderValue<HcalDataFrameContainer<QIE11DataFrame>> FullQIE11DataFrame(myReader, "QIE11DataFrameHcalDataFrameContainer_simHcalUnsuppressedDigis_HBHEQIE11DigiCollection_HLT.obj");

  Int_t tdc_of_interest = 62;

  gROOT->SetBatch(true); // hide plots as they are made

  // TDC = 62 rates from PU
  TH1F *HB = new TH1F("HB","Error Code 62 Rates in HB;Depth;Rate of TDC=62", 5,0,5);
  TH1F *HB62 = new TH1F("HB62","Error Code 62 Rates in HB;Depth;Rate of TDC=62", 5,0,5);
  TH1F *HE = new TH1F("HE","Error Code 62 Rates in HE;Depth;Rate of TDC=62", 8,0,8);
  TH1F *HE62 = new TH1F("HE62","Error Code 62 Rates in HE;Depth;Rate of TDC=62", 8,0,8);
  TCanvas *cHB = new TCanvas();
  TCanvas *cHE = new TCanvas();
  TEfficiency *effHB = 0;
  TEfficiency *effHE = 0;

  // TDC = 62 rates from noise
  TH1F *HBnoise = new TH1F("HBnoise","Error Code 62 Rates in HB (noise);Depth;Rate of TDC=62", 5,0,5);
  TH1F *HBnoise62 = new TH1F("HBnoise62","Error Code 62 Rates in HB (noise);Depth;Rate of TDC=62", 5,0,5);
  TH1F *HEnoise = new TH1F("HEnoise","Error Code 62 Rates in HE (noise);Depth;Rate of TDC=62", 8,0,8);
  TH1F *HEnoise62 = new TH1F("HEnoise62","Error Code 62 Rates in HE (noise);Depth;Rate of TDC=62", 8,0,8);
  TCanvas *cHBnoise = new TCanvas();
  TCanvas *cHEnoise = new TCanvas();
  TEfficiency *effHBnoise = 0;
  TEfficiency *effHEnoise = 0;

  int evtCounter = 0;
  int cells = 0;
  int cells62 = 0;
  int cells_noise = 0;
  int cells62_noise = 0;
  while (myReader.Next()){
    for (QIE11DataFrame frame:*FullQIE11DataFrame) { // loop over QIE11 data frame in HcalDataFrameContainer, this goes over FullQIE11DataFrame->size()
      HcalDetId QIEdetectorID = HcalDetId(frame.id());
      for (int i=0; i<frame.samples(); i++) { // loop over samples in QIE11 data frame, frame.samples is 8 in 110X samples
	if ( i == 2) { // TDC = 62 rates from noise
	  cells_noise += 1;
          if (abs(QIEdetectorID.ieta()) < 16 ) HBnoise->Fill(QIEdetectorID.depth());
          if (abs(QIEdetectorID.ieta()) >= 16 && abs(QIEdetectorID.ieta()) < 30 ) HEnoise->Fill(QIEdetectorID.depth());
	  if ( frame[i].tdc() == tdc_of_interest ) {
	    if (abs(QIEdetectorID.ieta()) < 16 ) HBnoise62->Fill(QIEdetectorID.depth());
	    if (abs(QIEdetectorID.ieta()) >= 16 && abs(QIEdetectorID.ieta()) < 30 ) HEnoise62->Fill(QIEdetectorID.depth());
	    cells62_noise += 1;
	  }
	}
	if ( i == 4) { // TDC = 62 rates from PU
	  cells += 1;
	  if (abs(QIEdetectorID.ieta()) < 16 ) HB->Fill(QIEdetectorID.depth());
          if (abs(QIEdetectorID.ieta()) >= 16 && abs(QIEdetectorID.ieta()) < 30 ) HE->Fill(QIEdetectorID.depth());
	  if ( frame[i].tdc() == tdc_of_interest ) {
	    if (abs(QIEdetectorID.ieta()) < 16 ) HB62->Fill(QIEdetectorID.depth());
	    if (abs(QIEdetectorID.ieta()) >= 16 && abs(QIEdetectorID.ieta()) < 30 ) HE62->Fill(QIEdetectorID.depth());
	    cells62 += 1;
	  }
	}
      }
    }
    evtCounter++; // increment to next event
  }

  std::cout << cells62 << " cells over TDC = 62, total cells = " << cells << std::endl;
  std::cout << cells62_noise << " cells over TDC = 62 (noise), total cells (noise) = " << cells_noise << std::endl;

  cHB->cd();
  if (TEfficiency::CheckConsistency(*HB62,*HB)) {
    effHB = new TEfficiency(*HB62,*HB);
    effHB->SetTitle("Error Code TDC=62 Rates in HB, threshold = 149.6");
    effHB->SetLineWidth(3.);
    effHB->SetLineColor(kBlack);
    effHB->Draw();
    gPad->Update();
    effHB->GetPaintedGraph()->SetMaximum(0.035);
    effHB->GetPaintedGraph()->SetMinimum(0.);
    gPad->Update();
  }
  cHB->SaveAs(Form("HB_TDCerror%i_149pt6.png", tdc_of_interest));

  cHE->cd();
  if (TEfficiency::CheckConsistency(*HE62,*HE)) {
    effHE = new TEfficiency(*HE62,*HE);
    effHE->SetTitle("Error Code TDC=62 Rates in HE, threshold = 149.6");
    effHE->SetLineWidth(3.);
    effHE->SetLineColor(kBlack);
    effHE->Draw();
    gPad->Update();
    effHE->GetPaintedGraph()->SetMaximum(0.45);
    effHE->GetPaintedGraph()->SetMinimum(0.);
    gPad->Update();
  }
  cHE->SaveAs(Form("HE_TDCerror%i_149pt6.png",tdc_of_interest));

  cHBnoise->cd();
  if (TEfficiency::CheckConsistency(*HBnoise62,*HBnoise)) {
    effHBnoise = new TEfficiency(*HBnoise62,*HBnoise);
    effHBnoise->SetTitle("Error Code TDC=62 Rates in HB (noise), threshold = 149.6");
    effHBnoise->SetLineWidth(3.);
    effHBnoise->SetLineColor(kBlack);
    effHBnoise->Draw();
    gPad->Update();
    effHBnoise->GetPaintedGraph()->SetMaximum(0.005);
    effHBnoise->GetPaintedGraph()->SetMinimum(0.);
    gPad->Update();
  }
  cHBnoise->SaveAs(Form("HBnoise_TDCerror%i_149pt6.png", tdc_of_interest));

  cHEnoise->cd();
  if (TEfficiency::CheckConsistency(*HEnoise62,*HEnoise)) {
    effHEnoise = new TEfficiency(*HEnoise62,*HEnoise);
    effHEnoise->SetTitle("Error Code TDC=62 Rates in HE (noise), threshold = 149.6");
    effHEnoise->SetLineWidth(3.);
    effHEnoise->SetLineColor(kBlack);
    effHEnoise->Draw();
    gPad->Update();
    effHEnoise->GetPaintedGraph()->SetMaximum(0.005);
    effHEnoise->GetPaintedGraph()->SetMinimum(0.);
    gPad->Update();
  }
  cHEnoise->SaveAs(Form("HEnoise_TDCerror%i_149pt6.png",tdc_of_interest));

  f->Close();
  return 0;
}
