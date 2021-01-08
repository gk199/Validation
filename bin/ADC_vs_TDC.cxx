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

int main() {
  TFile *f = new TFile("/afs/cern.ch/work/g/gkopp/MC_GenProduction/PionGun/CMSSW_11_0_2/src/Injected_Energy_0pt01_to_1pt2_HB111_tdc9pt35_tof_step1.root");

  TTreeReader myReader("Events",f);
  TTreeReaderValue<std::vector<CaloSamples>> AllCaloSamples(myReader, "CaloSampless_mix_HcalSamples_HLT.obj");
  TTreeReaderValue<HcalDataFrameContainer<QIE11DataFrame>> FullQIE11DataFrame(myReader, "QIE11DataFrameHcalDataFrameContainer_simHcalUnsuppressedDigis_HBHEQIE11DigiCollection_HLT.obj");

  Int_t ieta_of_interest = 1;
  Int_t iphi_of_interest = 1;

  gROOT->SetBatch(true); // hide plots as they are made

  TGraph *gr_adc = new TGraph();
  std::vector< TGraph* > gr_e;
  for (int i = 0; i<30; i++) gr_e.push_back(new TGraph());
  TCanvas *c1_adc = new TCanvas();
  int point_index_adc = 0;
  int point_index_e[30] = {0};

  gr_adc->SetTitle("ADC vs TDC, with TDC threshold = 9.35");
  gr_adc->GetYaxis()->SetTitle("TDC value, 0.5 ns steps");
  gr_adc->GetXaxis()->SetTitle("ADC value, 50 = 3 GeV");

  int evtCounter = 0;
  while (myReader.Next()){

    for (QIE11DataFrame frame:*FullQIE11DataFrame) { // loop over QIE11 data frame in HcalDataFrameContainer, this goes over FullQIE11DataFrame->size()
      HcalDetId QIEdetectorID = HcalDetId(frame.id());
      for (int i=0; i<frame.samples(); i++) { // loop over samples in QIE11 data frame
	if ( (frame[i].tdc() <= 60) && (frame[i].soi() == true) && (abs(QIEdetectorID.ieta()) == ieta_of_interest) && (QIEdetectorID.iphi() == iphi_of_interest) ) {
	  gr_adc->SetPoint(point_index_adc,frame[i].adc(),frame[i].tdc());
	  for (int j=0; j<30; j++) if (evtCounter < 100 * (j+1) && evtCounter >= 100 * j) {
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
  c1_adc->SaveAs("TDC_vs_ADC_full.pdf");
  gr_adc->SetMaximum(16);
  gr_adc->SetMinimum(0);
  gr_adc->Draw("P");
  c1_adc->SaveAs("TDC_vs_ADC.pdf");
  for (int i = 0; i < 30; i++) {
    gr_e[i]->SetMarkerColor(kRainBow+1.4*i);
    gr_e[i]->SetMarkerStyle(20);
    gr_e[i]->SetMarkerSize(0.5);
    gr_e[i]->Draw("P");
  }
  c1_adc->SaveAs("TDC_vs_ADC_colors.pdf");

  f->Close();
  return 0;
}
