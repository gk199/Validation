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
  //  TFile *f = new TFile("/afs/cern.ch/work/g/gkopp/MC_GenProduction/PionGun/CMSSW_11_0_2/src/SinglePion211_E10_eta1phi0_PU_00_step1.root");
  //  TFile *f = new TFile("/afs/cern.ch/work/g/gkopp/MC_GenProduction/PionGun/CMSSW_11_0_2/src/SinglePion211_E0pt05_HB111_tdc9pt35_tof_injected_step1.root");
  TFile *f = new TFile("/afs/cern.ch/work/g/gkopp/MC_GenProduction/PionGun/CMSSW_11_0_2/src/Injected_Energy_0pt01_to_1pt2_HB111_tdc9pt35_tof_step1.root");

  TTreeReader myReader("Events",f);
  TTreeReaderValue<std::vector<CaloSamples>> AllCaloSamples(myReader, "CaloSampless_mix_HcalSamples_HLT.obj");
  TTreeReaderValue<HcalDataFrameContainer<QIE11DataFrame>> FullQIE11DataFrame(myReader, "QIE11DataFrameHcalDataFrameContainer_simHcalUnsuppressedDigis_HBHEQIE11DigiCollection_HLT.obj");

  //  int events_of_interest[20] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
  Int_t adc_of_interest = 50;
  Int_t adc_of_interest2 = 50;
  Int_t ieta_of_interest = 1;
  Int_t iphi_of_interest = 1;
  std::vector<int> event_withADC;
  std::vector<int> ieta_withADC;
  std::vector<int> iphi_withADC;
  std::vector<int> depth_withADC;
  std::vector<int> event_withADC2;
  std::vector<int> ieta_withADC2;
  std::vector<int> iphi_withADC2;
  std::vector<int> depth_withADC2;

  gROOT->SetBatch(true); // hide plots as they are made

  TGraph *gr_adc = new TGraph();
  TGraph *gr_adc2 = new TGraph();
  TCanvas *c1_adc = new TCanvas("c","c",1200,600);
  c1_adc->Divide(2,1);
  int point_index_adc = 0;
  int point_index_adc2 = 0;

  gr_adc->SetTitle(Form("Precise Data for ADC<=%i, ieta=%i", adc_of_interest, ieta_of_interest));
  gr_adc->GetYaxis()->SetTitle("PreciseData pulse shape");
  gr_adc->GetXaxis()->SetTitle("Time, 0.5 ns steps");
  gr_adc2->SetTitle(Form("Precise Data for ADC>=%i, ieta=%i", adc_of_interest2, ieta_of_interest));
  gr_adc2->GetYaxis()->SetTitle("PreciseData pulse shape");
  gr_adc2->GetXaxis()->SetTitle("Time, 0.5 ns steps");

  int evtCounter = 0;
  while (myReader.Next()){
    /*
    if (std::find(std::begin(events_of_interest), std::end(events_of_interest), evtCounter) == std::end(events_of_interest)) { // if event is not in event of interest list, skip
      evtCounter++;
      continue;
    }
    */

    //    std::cout << FullQIE11DataFrame->size() << " QIE11 size for event " << evtCounter << std::endl; // 15840 per event
    for (QIE11DataFrame frame:*FullQIE11DataFrame) { // loop over QIE11 data frame in HcalDataFrameContainer, this goes over FullQIE11DataFrame->size()
      HcalDetId QIEdetectorID = HcalDetId(frame.id());
      for (int i=0; i<frame.samples(); i++) { // loop over samples in QIE11 data frame
	if ( (frame[i].tdc() <= 60) && (frame[i].soi() == true) && (abs(QIEdetectorID.ieta()) == ieta_of_interest) && (QIEdetectorID.iphi() == iphi_of_interest) ) {
	  if ( (frame[i].adc() <= adc_of_interest) ) {
	    if (frame[i].adc() < 20 ) std::cout << "event, ADC, TDC = " << evtCounter << ", " << frame[i].adc() << ", " << frame[i].tdc() << std::endl;
	    event_withADC.push_back(evtCounter);
	    ieta_withADC.push_back(QIEdetectorID.ieta());
	    iphi_withADC.push_back(QIEdetectorID.iphi());
	    depth_withADC.push_back(QIEdetectorID.depth());
	  }
	  if ( (frame[i].adc() >= adc_of_interest2) ) {
	    //	    std::cout << "ADC value = " << frame[i].adc() << std::endl;
	    //	    std::cout << "TDC value = " << frame[i].tdc() << std::endl;
	    event_withADC2.push_back(evtCounter);
	    ieta_withADC2.push_back(QIEdetectorID.ieta());
	    iphi_withADC2.push_back(QIEdetectorID.iphi());
	    depth_withADC2.push_back(QIEdetectorID.depth());
	  }
	}
      }
    }

    //    TGraph *gr_adc = new TGraph();
    //    TCanvas *c1_adc = new TCanvas();
    //    int point_index_adc = 0;

    /*
    gr_adc->SetTitle(Form("Precise Data for Event=%i with ADC=%i",evtCounter, adc_of_interest));
    gr_adc->GetYaxis()->SetTitle("PreciseData pulse shape");
    gr_adc->GetXaxis()->SetTitle("Time, 0.5 ns steps");
    */

    //    std::cout << AllCaloSamples->size() << std::endl; // how many caloSamples do we have? ~6k per event
    for (CaloSamples CaloSample:*AllCaloSamples) { // loop over everything in AllCaloSamples, call it CaloSample
      HcalDetId detectorID = HcalDetId(CaloSample.id()); // convert raw ID to detector ID 
      for (unsigned int of_interest = 0; of_interest < ieta_withADC.size(); of_interest++) {
	if ( (evtCounter == event_withADC[of_interest]) && (detectorID.ieta() == ieta_withADC[of_interest]) && (detectorID.iphi() == iphi_withADC[of_interest]) && (detectorID.depth() == depth_withADC[of_interest]) ) { // for looking at the pulse shape of a certain ADC value, check looking at right event, ieta, iphi, depth cell
          for (int i = 0; i < CaloSample.preciseSize(); i++) { // assume length of precise data is given by precise size
            gr_adc->SetPoint(point_index_adc,i,CaloSample.preciseAt(i)); // point index, x coordinate (time)
            point_index_adc += 1;
            if (i+1 == CaloSample.preciseSize()) {
              gr_adc->SetPoint(point_index_adc,399,0);
              point_index_adc += 1;
            }
          }
        }
      }
      for (unsigned int of_interest = 0; of_interest < ieta_withADC2.size(); of_interest++) {
	if ( (evtCounter == event_withADC2[of_interest]) && (detectorID.ieta() == ieta_withADC2[of_interest]) && (detectorID.iphi() == iphi_withADC2[of_interest]) && (detectorID.depth() == depth_withADC2[of_interest]) ) { // for looking at the pulse shape of a certain ADC value, check looking at right event, ieta, iphi, depth cell
	  for (int i = 0; i < CaloSample.preciseSize(); i++) { // assume length of precise data is given by precise size
	    gr_adc2->SetPoint(point_index_adc2,i,CaloSample.preciseAt(i)); // point index, x coordinate (time)
	    point_index_adc2 += 1;
	    if (i+1 == CaloSample.preciseSize()) {
	      gr_adc2->SetPoint(point_index_adc2,399,0);
	      point_index_adc2 += 1;
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
  c1_adc->cd(1);
  gr_adc->Draw();
  c1_adc->cd(2);
  gr_adc2->Draw();
  c1_adc->SaveAs(Form("CaloSamplesAnalysis_adc%i_adc%i.png", adc_of_interest, adc_of_interest2));

  f->Close();
  return 0;
}
