#include "CalibFormats/CaloObjects/interface/CaloSamples.h"
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
//#include "DataFormats/HcalDetId/interface/HcalDetId.h"
//#include "DataFormats/HcalDigi/interface/QIE11DataFrame.h"

void CaloSamplesAnalysis() { 
  TH1F *h_preciseData = new TH1F("h_preciseData","preciseData;preciseData;Entries",100,0,10);

  TFile *f = new TFile("/eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/PionGun/SinglePion211_E10_PU_00_eta1phi0_timeslew-false_step1_CaloSamples_10events.root");

  TTreeReader myReader("Events",f);
  TTreeReaderValue<std::vector<CaloSamples>> AllCaloSamples(myReader, "CaloSampless_mix_HcalSamples_HLT.obj");
  //  TTreeReaderValue<HcalDataFrameContainer<QIE11DataFrame>> FullQIE11DataFrame(myReader, "QIE11DataFrameHcalDataFrameContainer_simHcalUnsuppressedDigis_HBHEQIE11DigiCollection_HLT.obj");

  Int_t event_of_interest = 29;
  Int_t ieta_of_interest = 18;
  Int_t iphi_of_interest = 24;
  Int_t depth_of_interest = 2;

  gROOT->SetBatch(true); // hide plots as they are made

  int evtCounter = 0;
  while (myReader.Next()){
    //    for (QIE11DataFrame QIE11Digis:*FullQIE11DataFrame) std::cout << FullQIE11DataFrame->size() << " QIE11 size for event " << evtCounter << std::endl;

    if (evtCounter != event_of_interest) {
      evtCounter++;
      continue;
    }
    //    std::cout << AllCaloSamples->size() << std::endl; // how many caloSamples do we have? ~6k per event
    for (CaloSamples CaloSample:*AllCaloSamples) { // loop over everything in AllCaloSamples, call it CaloSample
      HcalDetId detectorID = HcalDetId(CaloSample.id());
      //      std::cout << detectorID.depth() << "  " << detectorID.ieta() << "  " << detectorID.iphi() << std::endl;

      if (detectorID.ieta() == ieta_of_interest){ // convert raw ID to detector ID
	if (detectorID.iphi() == iphi_of_interest){
	  if (detectorID.depth() == depth_of_interest){
      	    std::cout<< "yay" << std::endl; 

	    TGraph *gr = new TGraph();
      	    for (int i = 0; i < CaloSample.preciseSize(); i++) { // assume length of precise data is given by precise size
	      gr->SetPoint(i,i,CaloSample.preciseAt(i));
	    }

	    TCanvas *c1 = new TCanvas();
	    gr->SetTitle(Form("Precise Data for ieta=%i, iphi=%i, depth=%i, Event=%i",detectorID.ieta(),detectorID.iphi(),detectorID.depth(),evtCounter));
	    gr->GetYaxis()->SetTitle("PreciseData pulse shape");
            gr->GetXaxis()->SetTitle("Time");
	    gr->Draw();
	    c1->SaveAs(Form("CaloSamplesAnalysis_event%i_ieta%i_iphi%i_depth%i.png",evtCounter,detectorID.ieta(),detectorID.iphi(),detectorID.depth()));
	  }
	}
      }
    }
    evtCounter++; // increment to next event
  }

  f->Close();
}
