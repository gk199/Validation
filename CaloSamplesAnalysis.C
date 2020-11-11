#include "CalibFormats/CaloObjects/interface/CaloSamples.h"
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"

void CaloSamplesAnalysis() { 
  TH1F *h_preciseData = new TH1F("h_preciseData","preciseData;preciseData;Entries",100,0,10);

  TFile *f = new TFile("/eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/PionGun/SinglePion211_E10_PU_00_eta1phi0_timeslew-false_step1_CaloSamples_10events.root");

  TTreeReader myReader("Events",f);
  TTreeReaderValue<std::vector<CaloSamples>> AllCaloSamples(myReader, "CaloSampless_mix_HcalSamples_HLT.obj");

  Int_t event_of_interest = 1;
  Int_t ieta_of_interest = 1;
  Int_t iphi_of_interest = 1;
  Int_t depth_of_interest = 1;
  std::vector<float> preciseData;

  unsigned int evtCounter = 0;
  while (myReader.Next()){
    //    if (evtCounter != event_of_interest) continue;
    //    std::cout << AllCaloSamples->size() << std::endl; // how many caloSamples do we have? ~6k per event
    for (CaloSamples CaloSample:*AllCaloSamples) { // loop over everything in AllCaloSamples, call it CaloSample
      HcalDetId detectorID = HcalDetId(CaloSample.id());
      if (detectorID.ieta() == ieta_of_interest){ // convert raw ID to detector ID
	if (detectorID.iphi() == iphi_of_interest){
	  if (detectorID.depth() == depth_of_interest){
      	    std::cout<< "yay" << std::endl; // about 10 in this file total

	    TGraph *gr = new TGraph();
      	    for (int i = 0; i < CaloSample.preciseSize(); i++) { // assume length of precise data is given by precise size
	      gr->SetPoint(i,i,CaloSample.preciseAt(i));
	    }

	    TCanvas *c1 = new TCanvas();
	    gr->Draw();
	    c1->SaveAs(Form("CaloSamplesAnalysis_%u.png",evtCounter));
	  }
	}
      }
    }
    evtCounter++; // increment to next event
  }

  f->Close();
}
