#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TChain.h"
#include "TROOT.h"

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>

#include "DataFormats/Common/interface/EDProductGetter.h" // mustBeNonZero
#include "CalibFormats/CaloObjects/interface/CaloSamples.h"
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

int main(int argc, char *argv[]) {
//Injected_Energy_0pt001_to_1pt3_HB111_tdc2pt34_tof_step1.root");
//Injected_Energy_0pt001_to_1pt3_HB111_tdc4pt675_tof_step1.root");
//Injected_Energy_0pt001_to_1pt3_HB111_tdc9pt35_tof_step1.root");
//Injected_Energy_0pt001_to_1pt3_HB111_tdc9pt35_tof_timephase7_step1.root");
  TFile *f = new TFile("/afs/cern.ch/work/g/gkopp/MC_GenProduction/PionGun/CMSSW_11_0_2/src/Injected_Energy_0pt001_to_1pt3_HB111_tdc18pt7_tof_step1.root");

  TTreeReader myReader("Events",f);
  TTreeReaderValue<std::vector<CaloSamples>> AllCaloSamples(myReader, "CaloSampless_mix_HcalSamples_HLT.obj");
  TTreeReaderValue<HcalDataFrameContainer<QIE11DataFrame>> FullQIE11DataFrame(myReader, "QIE11DataFrameHcalDataFrameContainer_simHcalUnsuppressedDigis_HBHEQIE11DigiCollection_HLT.obj");

  Int_t ieta_of_interest = atoi(argv[1]); //1;
  Int_t iphi_of_interest = atoi(argv[2]); //1;
  Int_t depth_of_interest = atoi(argv[3]); //1;

  gROOT->SetBatch(true); // hide plots as they are made

  TGraph *gr_adc = new TGraph();
  std::vector< TGraph* > gr_e;
  for (int i = 0; i<40; i++) gr_e.push_back(new TGraph());
  TCanvas *c1_adc = new TCanvas();
  int point_index_adc = 0;
  int point_index_e[40] = {0};

  gr_adc->SetTitle(Form("ADC vs TDC, with TDC threshold = 18.7, timephase = 6, ieta%i, iphi%i, depth%i",ieta_of_interest, iphi_of_interest, depth_of_interest));
  gr_adc->GetYaxis()->SetTitle("TDC value, 0.5 ns steps");
  gr_adc->GetXaxis()->SetTitle("ADC value, 50 = 3 GeV");

  int evtCounter = 0;
  while (myReader.Next()){

    for (QIE11DataFrame frame:*FullQIE11DataFrame) { // loop over QIE11 data frame in HcalDataFrameContainer, this goes over FullQIE11DataFrame->size()
      HcalDetId QIEdetectorID = HcalDetId(frame.id());
      for (int i=0; i<frame.samples(); i++) { // loop over samples in QIE11 data frame
	if ( (frame[i].tdc() <= 60) && (frame[i].soi() == true) && (abs(QIEdetectorID.ieta()) == ieta_of_interest) && (QIEdetectorID.iphi() == iphi_of_interest) && (QIEdetectorID.depth() == depth_of_interest) ) {
	  gr_adc->SetPoint(point_index_adc,frame[i].adc(),frame[i].tdc());
	  for (int j=0; j<40; j++) if (evtCounter < 100 * (j+1) && evtCounter >= 100 * j) {
	      gr_e[j]->SetPoint(point_index_e[j],frame[i].adc(),frame[i].tdc());
	      point_index_e[j] += 1;
	    }
	  point_index_adc += 1;
	}
      }
    }
    evtCounter++; // increment to next event
  }
  gr_adc->Draw("AP");
  gr_adc->SetMarkerStyle(20);
  gr_adc->SetMarkerSize(0.5);
  c1_adc->SaveAs(Form("TDC_vs_ADC_full_18pt7_timephase6_ieta%i_iphi%i_depth%i.pdf",ieta_of_interest, iphi_of_interest, depth_of_interest));
  gr_adc->SetMaximum(30);
  gr_adc->SetMinimum(0);
  gr_adc->Draw("P");
  c1_adc->SaveAs(Form("TDC_vs_ADC_18pt7_timephase6_ieta%i_iphi%i_depth%i.pdf",ieta_of_interest, iphi_of_interest, depth_of_interest));
  for (int i = 0; i < 40; i++) {
    gr_e[i]->SetMarkerColor(kRainBow+1.1*i);
    gr_e[i]->SetMarkerStyle(20);
    gr_e[i]->SetMarkerSize(0.5);
    gr_e[i]->Draw("P");
  }
  c1_adc->SaveAs(Form("TDC_vs_ADC_colors_18pt7_timephase6_ieta%i_iphi%i_depth%i.pdf",ieta_of_interest, iphi_of_interest, depth_of_interest));

  f->Close();
  return 0;
}
