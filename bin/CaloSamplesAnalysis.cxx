#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TROOT.h"

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

#include "DataFormats/Common/interface/EDProductGetter.h" // mustBeNonZero
#include "CalibFormats/CaloObjects/interface/CaloSamples.h"
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

int main() {
  //  TH1F *h_preciseData = new TH1F("h_preciseData","preciseData;preciseData;Entries",100,0,10);
  //  TFile *f = new TFile("/eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/PionGun/SinglePion211_E10_PU_00_eta1phi0_timeslew-false_step1_CaloSamples_10events.root"); // looked at event 29
  //  TFile *f = new TFile("/eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/PionGun/SinglePion211_E10_PU_00_eta1phi0_step1_threshold1x_CaloSamples_100events.root");
  //  TFile *f = new TFile("/eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/TDC_threshold_18pt7x3/MH-125_MFF-50_CTau-10000mm_step1_CaloSamples.root"); // looked at event 39
  TFile *f = new TFile("/afs/cern.ch/work/g/gkopp/MC_GenProduction/PionGun/CMSSW_11_0_2/src/SinglePion211_E10_eta1phi0_PU_00_step1.root");

  TTreeReader myReader("Events",f);
  TTreeReaderValue<std::vector<CaloSamples>> AllCaloSamples(myReader, "CaloSampless_mix_HcalSamples_HLT.obj");
  TTreeReaderValue<HcalDataFrameContainer<QIE11DataFrame>> FullQIE11DataFrame(myReader, "QIE11DataFrameHcalDataFrameContainer_simHcalUnsuppressedDigis_HBHEQIE11DigiCollection_HLT.obj");

  int events_of_interest[20] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
  Int_t tdc_of_interest = 63;
  std::vector<int> event_withTDC;
  std::vector<int> ieta_withTDC;
  std::vector<int> iphi_withTDC;
  std::vector<int> depth_withTDC;

  gROOT->SetBatch(true); // hide plots as they are made

  TGraph *gr = new TGraph();
  TCanvas *c1 = new TCanvas();
  int point_index = 0;

  gr->SetTitle(Form("Precise Data for TDC=%i", tdc_of_interest));
  gr->GetYaxis()->SetTitle("PreciseData pulse shape");
  gr->GetXaxis()->SetTitle("Time, 0.5 ns steps");

  int evtCounter = 0;
  while (myReader.Next()){
    if (std::find(std::begin(events_of_interest), std::end(events_of_interest), evtCounter) == std::end(events_of_interest)) { // if event is not in event of interest list, skip
      evtCounter++;
      continue;
    }

    std::cout << FullQIE11DataFrame->size() << " QIE11 size for event " << evtCounter << std::endl; // 15840 per event
    for (QIE11DataFrame frame:*FullQIE11DataFrame) { // loop over QIE11 data frame in HcalDataFrameContainer, this goes over FullQIE11DataFrame->size()
      HcalDetId QIEdetectorID = HcalDetId(frame.id());
      for (int i=0; i<frame.samples(); i++) { // loop over samples in QIE11 data frame
	if ( (frame[i].soi() == true) && (frame[i].tdc() == tdc_of_interest) ) {
	  //	  std::cout << frame[i].tdc() << " = TDC value in cell of interest" << std::endl;
	  event_withTDC.push_back(evtCounter);
	  ieta_withTDC.push_back(QIEdetectorID.ieta());
	  iphi_withTDC.push_back(QIEdetectorID.iphi());
	  depth_withTDC.push_back(QIEdetectorID.depth());
	}
      }
    }
    std::cout << ieta_withTDC.size() << " length of ieta, iphi, depth_withTDC vectors" << std::endl;

    //    TGraph *gr = new TGraph();
    //    TCanvas *c1 = new TCanvas();
    //    int point_index = 0; // so don't overwrite TGraph each time with the last pulse shape

    /*
    gr->SetTitle(Form("Precise Data for Event=%i with TDC=%i",evtCounter, tdc_of_interest));
    gr->GetYaxis()->SetTitle("PreciseData pulse shape");
    gr->GetXaxis()->SetTitle("Time, 0.5 ns steps");
    */

    //    std::cout << AllCaloSamples->size() << std::endl; // how many caloSamples do we have? ~6k per event
    for (CaloSamples CaloSample:*AllCaloSamples) { // loop over everything in AllCaloSamples, call it CaloSample
      HcalDetId detectorID = HcalDetId(CaloSample.id()); // convert raw ID to detector ID 
      for (unsigned int of_interest = 0; of_interest < ieta_withTDC.size(); of_interest++) {
	if ( (evtCounter == event_withTDC[of_interest]) && (detectorID.ieta() == ieta_withTDC[of_interest]) && (detectorID.iphi() == iphi_withTDC[of_interest]) && (detectorID.depth() == depth_withTDC[of_interest]) ) { // for looking at the pulse shape of a certain TDC value, check looking at right event, ieta, iphi, depth cell
	  for (int i = 0; i < CaloSample.preciseSize(); i++) { // assume length of precise data is given by precise size
	    gr->SetPoint(point_index,i,CaloSample.preciseAt(i)); // point index, x coordinate (time)
	    point_index += 1;
	    if (i+1 == CaloSample.preciseSize()) {
	      gr->SetPoint(point_index,399,0);
	      point_index += 1;
	    }
	  }
	}
      }
    }
    /*
    gr->Draw();
    c1->SaveAs(Form("CaloSamplesAnalysis_event%i_tdc%i.png",evtCounter, tdc_of_interest));
    gr_adc->Draw();
    c1_adc->SaveAs(Form("CaloSamplesAnalysis_event%i_adc%i.png",evtCounter, adc_of_interest));
    */

    evtCounter++; // increment to next event
  }
  gr->Draw();
  c1->SaveAs(Form("CaloSamplesAnalysis_tdc%i.png", tdc_of_interest));

  f->Close();
  return 0;
}
