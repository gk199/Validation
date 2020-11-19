#include "CalibFormats/CaloObjects/interface/CaloSamples.h"
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
//#include "DataFormats/HcalDetId/interface/HcalDetId.h"
//#include "DataFormats/HcalDigi/interface/QIE11DataFrame.h"

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

void CaloSamplesAnalysis() { 
  TH1F *h_preciseData = new TH1F("h_preciseData","preciseData;preciseData;Entries",100,0,10);

  //  TFile *f = new TFile("/eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/PionGun/SinglePion211_E10_PU_00_eta1phi0_timeslew-false_step1_CaloSamples_10events.root");
  TFile *f = new TFile("/eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/PionGun/SinglePion211_E10_PU_00_eta1phi0_step1_threshold1x_CaloSamples_100events.root");

  TTreeReader myReader("Events",f);
  TTreeReaderValue<std::vector<CaloSamples>> AllCaloSamples(myReader, "CaloSampless_mix_HcalSamples_HLT.obj");
  //  TTreeReaderValue<HcalDataFrameContainer<QIE11DataFrame>> FullQIE11DataFrame(myReader, "QIE11DataFrameHcalDataFrameContainer_simHcalUnsuppressedDigis_HBHEQIE11DigiCollection_HLT.obj");

  Int_t event_of_interest = 73;
  Int_t ieta_of_interest = 19;
  Int_t iphi_of_interest = 26;
  Int_t depth_of_interest = 2;
  double rDepth[4] = {183.60, 190.20, 214.20, 244.80};
  double tofadj[4] = {-3.5, -3.5, -4.5, -5.5};

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

	    // calculate TOF. The reported TDC should be QIE11.tdc() - tof = QIE11.tdc() - (distance + tofadjustment)
	    // taken from hcalqie11tdc.py code
	    double eta = etaVal(detectorID.ieta());
	    double theta = 2 * atan(exp(-eta));
	    double distance = rDepth[depth_of_interest - 1] / sin(theta); // distance to collision vertex
	    double tof = 1e9*distance/2.99793e10;  // cm/sec,  *10^9 to get ns
	    tof += tofadj[depth_of_interest-1];
	    std::cout << "time of flight = distance + tof adjustment = " << tof << std::endl;
	  }
	}
      }
    }
    evtCounter++; // increment to next event
  }

  f->Close();
}
