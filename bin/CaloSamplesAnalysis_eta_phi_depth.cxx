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
//#include "DataFormats/Common/interface/Ref.h"
#include "CalibFormats/CaloObjects/interface/CaloSamples.h"
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
//#include "DataFormats/HcalDetId/interface/HcalDetId.h"
//#include "DataFormats/HcalDigi/interface/QIE11DataFrame.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

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

int main() {
  TFile *f = new TFile("/afs/cern.ch/work/g/gkopp/MC_GenProduction/PionGun/CMSSW_11_0_2/src/SinglePion211_E10_eta1phi0_PU_00_step1.root");
  //  TFile *f = new TFile("/afs/cern.ch/work/g/gkopp/MC_GenProduction/PionGun/CMSSW_11_0_2/src/SinglePion211_E10_eta1phi0_PU_step1.root");
  //  TFile *f = new TFile("/eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/PionGun/SinglePion211_E10_PU_00_eta1phi0_step1_threshold1x_CaloSamples_100events.root");
  //  TFile *f = new TFile("/eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/TDC_threshold_18pt7x3/MH-125_MFF-50_CTau-10000mm_step1_CaloSamples.root");

  TTreeReader myReader("Events",f);
  TTreeReaderValue<std::vector<CaloSamples>> AllCaloSamples(myReader, "CaloSampless_mix_HcalSamples_HLT.obj");
  TTreeReaderValue<HcalDataFrameContainer<QIE11DataFrame>> FullQIE11DataFrame(myReader, "QIE11DataFrameHcalDataFrameContainer_simHcalUnsuppressedDigis_HBHEQIE11DigiCollection_HLT.obj");

  //  int events_of_interest[40] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40};
  Int_t event_of_interest = 37;
  Int_t ieta_of_interest = -17;
  Int_t iphi_of_interest = 24;
  Int_t depth_of_interest = 2;
  double rDepth[4] = {183.60, 190.20, 214.20, 244.80};
  double tofadj[4] = {-3.5, -3.5, -4.5, -5.5};

  gROOT->SetBatch(true); // hide plots as they are made

  int evtCounter = 0;
  while (myReader.Next()){
    //    if (std::find(std::begin(events_of_interest), std::end(events_of_interest), evtCounter) == std::end(events_of_interest)) { // if event is not in list of event of interest, skip
    if (evtCounter != event_of_interest) { // if event is not in event of interest list, skip
      evtCounter++;
      continue;
    }

    std::cout << FullQIE11DataFrame->size() << " QIE11 size for event " << evtCounter << std::endl; // 15840 per event
    for (QIE11DataFrame frame:*FullQIE11DataFrame) { // loop over QIE11 data frame in HcalDataFrameContainer, this goes over FullQIE11DataFrame->size()
      HcalDetId QIEdetectorID = HcalDetId(frame.id());
      if ( (QIEdetectorID.ieta() == ieta_of_interest) && (QIEdetectorID.iphi() == iphi_of_interest) && (QIEdetectorID.depth() == depth_of_interest) ) {
	for (int i=0; i<frame.samples(); i++) { // loop over samples in QIE11 data frame
	  if ( (frame[i].soi() == true) ) {
	    std::cout << frame[i].tdc() << " = TDC (in QIE11, 0-50) value in cell of interest" << std::endl;
	  }
	}
      }
    }
    
    TGraph *gr = new TGraph();
    TCanvas *c1 = new TCanvas();

    //    std::cout << AllCaloSamples->size() << std::endl; // how many caloSamples do we have? ~6k per event
    for (CaloSamples CaloSample:*AllCaloSamples) { // loop over everything in AllCaloSamples, call it CaloSample
      HcalDetId detectorID = HcalDetId(CaloSample.id()); // convert raw ID to detector ID 
      if ( (detectorID.ieta() == ieta_of_interest) && (detectorID.iphi() == iphi_of_interest) && (detectorID.depth() == depth_of_interest) ) { // for a particular ieta, iphi, depth pulse shape
	gr->SetTitle(Form("Precise Data for ieta=%i, iphi=%i, depth=%i, Event=%i",detectorID.ieta(),detectorID.iphi(),detectorID.depth(),evtCounter));
	gr->GetYaxis()->SetTitle("PreciseData pulse shape");
	gr->GetXaxis()->SetTitle("Time, 0.5 ns steps");	  
	for (int i = 0; i < CaloSample.preciseSize(); i++) { // assume length of precise data is given by precise size
	  gr->SetPoint(i,i,CaloSample.preciseAt(i)); // point index, x coordinate (time)
	}
	gr->Draw();
	c1->SaveAs(Form("CaloSamplesAnalysis_event%i_ieta%i_iphi%i_depth%i.png",evtCounter,detectorID.ieta(),detectorID.iphi(),detectorID.depth()));
      }
    }

    // calculate TOF. The reported TDC should be QIE11.tdc() - tof = QIE11.tdc() - (distance + tofadjustment)                             
    // taken from hcalqie11tdc.py code                                                                                                    
    double eta = etaVal(ieta_of_interest);
    double theta = 2 * atan(exp(-eta));
    double distance = rDepth[depth_of_interest - 1] / sin(theta); // distance to collision vertex                                         
    double tof = 1e9*distance/2.99793e10;  // cm/sec,  *10^9 to get ns                                                                    
    tof += tofadj[depth_of_interest-1];
    std::cout << "time of flight = distance + tof adjustment = " << tof << std::endl;
    
    evtCounter++; // increment to next event
  }
  f->Close();
  return 0;
}
