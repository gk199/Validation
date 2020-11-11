#include <vector>
#include "CalibFormats/CaloObjects/interface/CaloSamples.h"
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"

void CaloSamplesAnalysis() { 
  TH1F *h_preciseData = new TH1F("h_preciseData","preciseData;preciseData;Entries",100,0,100);

  TFile *f = new TFile("/eos/cms/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/PionGun/SinglePion211_E10_PU_00_eta1phi0_timeslew-false_step1_CaloSamples_10events.root");
  TTree *t1 = (TTree*)f->Get("Events");

  std::vector<CaloSamples> *AllCaloSamples = new std::vector<CaloSamples>();

  t1->SetBranchAddress("CaloSampless_mix_HcalSamples_HLT.obj",&AllCaloSamples);

  Int_t nentries = (Int_t)t1->GetEntries();
  std::cout << nentries << std::endl;

  Int_t event_of_interest = 1;
  Int_t ieta_of_interest = 1;
  Int_t iphi_of_interest = 1;
  Int_t depth_of_interest = 1;

  for (Int_t i=0; i<nentries; i++) {
    t1->GetEntry(i);
    if (i==1) { 
      std::cout << "hi" <<std::endl;
      std::cout << AllCaloSamples->size() << std::endl;
    }
    for (CaloSamples CaloSample:*AllCaloSamples) { // loop over everything in AllCaloSamples, call it CaloSample
      std::cout << HcalTrigTowerDetId(CaloSample.id()).ieta() << std::endl;
      if (HcalTrigTowerDetId(CaloSample.id()).ieta() == ieta_of_interest) std::cout<< "yay" << std::endl;
    }
  }

  //  TCanvas *c1 = new TCanvas();
  //  h_preciseData->Draw();
  //  c1->SaveAs("CaloSamplesAnalysis.png");

  f->Close();
}
