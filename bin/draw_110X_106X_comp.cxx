#include "PhysicsTools/Utilities/macros/setTDRStyle.C"

#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TPaveStats.h"
#include <map>
#include <string>
#include <vector>

#include <iostream>

int main()
{
  setTDRStyle();
  gROOT->ForceStyle();

  std::vector<std::string> filenames = {"rates_new_cond_LLP_mh1000_pl10000.root", "/afs/cern.ch/work/g/gkopp/HCAL_Trigger/L1Ntuples/HCAL_TP_TimingBitEmulator/CMSSW_11_0_2/src/HcalTrigger/Validation/rates_new_cond_LLP_mh1000_pl10000.root"};
  //  std::vector<std::string> filenames = {"rates_new_cond_106X_nugun.root", "/afs/cern.ch/work/g/gkopp/HCAL_Trigger/L1Ntuples/HCAL_TP_TimingBitEmulator/CMSSW_11_0_2/src/HcalTrigger/Validation/rates_new_cond_110X_nugun.root"};
  std::vector<std::string> multTypes = {"htSumDistribution","L1_Jet1_ET","L1_Jet2_ET","L1_Jet3_ET","L1_Jet4_ET"};//ADC50_3ns_4JetMultHB", "ADC50_3ns_4JetMultHE", "ADC50_3ns_4JetMultHBHE"};

  std::map<std::string, TH1F*> multHists_106X;
  std::map<std::string, TH1F*> multHists_110X;

  std::vector<TFile*> files;
  for(auto file : filenames) {
    files.push_back(TFile::Open(file.c_str()));
  }

  for(auto multType : multTypes) {
    std::string histName(multType);
    //    histName += "_emu";

    multHists_106X[multType] = dynamic_cast<TH1F*>(files.at(0)->Get(histName.c_str()));
    multHists_110X[multType] = dynamic_cast<TH1F*>(files.at(1)->Get(histName.c_str()));
  }

  for(auto pair : multHists_106X) pair.second->SetLineWidth(2);
  for(auto pair : multHists_110X) pair.second->SetLineWidth(2);

  std::vector<TCanvas*> canvases;
  std::vector<TPad*> pad1;
  std::vector<TPad*> pad2;

  for(auto hist : multTypes) {
    canvases.push_back(new TCanvas);
    canvases.back()->SetWindowSize(1.3*canvases.back()->GetWw(), 0.9*canvases.back()->GetWh());
    pad1.push_back(new TPad("pad1", "pad1", 0, 0, 1, 0.9));
    //    pad1.back()->SetGrid();
    pad1.back()->Draw();
    pad1.back()->cd();

    multHists_110X[hist]->SetLineColor(kBlue); 
    multHists_110X[hist]->Scale(1./multHists_110X[hist]->Integral());
    TString name(multHists_110X[hist]->GetName());
    //    multHists_110X[hist]->GetXaxis()->SetRangeUser(0,900);
    //    if ( strcmp(hist.substr(0,2).c_str(), "L1") == 0 ) multHists_110X[hist]->GetXaxis()->SetRangeUser(0,120);

    multHists_110X[hist]->Draw("hist");
    gStyle->SetOptStat(1);
    gPad->Update();
    TPaveStats *st1 = (TPaveStats*)multHists_110X[hist]->FindObject("stats");
    st1->SetName("110X");
    st1->SetY1NDC(0.5);
    st1->SetY2NDC(0.7);
    st1->SetTextColor(kBlue); 

    multHists_106X[hist]->SetLineColor(kRed);
    multHists_106X[hist]->Scale(1./multHists_106X[hist]->Integral());
    multHists_106X[hist]->Draw("hist sames");
    gStyle->SetOptStat(1);
    gPad->Update();
    TPaveStats *st = (TPaveStats*)multHists_106X[hist]->FindObject("stats");
    st->SetName("106X");
    st->SetY1NDC(0.2);
    st->SetY2NDC(0.4);
    st->SetTextColor(kRed); 

    TLegend *leg = new TLegend(0.6, 0.75, 0.9, 0.9);
    leg->AddEntry(multHists_106X[hist],"106X, mh=1TeV, pl=10m","L"); //neutrino gun", "L");
    leg->AddEntry(multHists_110X[hist], "110X, mh=1TeV, pl=10m","L"); //neutrino gun", "L");
    multHists_110X[hist]->SetTitle(Form("Comparison between 110X and 106X MC distributions -- %s", hist.substr(0).c_str())); 
    multHists_110X[hist]->GetXaxis()->SetLabelSize(0.03);
    multHists_110X[hist]->GetYaxis()->SetLabelSize(0.03);
    multHists_110X[hist]->GetXaxis()->SetTitleSize(0.04);
    multHists_110X[hist]->GetYaxis()->CenterTitle(true);
    multHists_110X[hist]->GetYaxis()->SetTitleSize(0.04);
    multHists_110X[hist]->GetXaxis()->SetTitleOffset(1.2);
    multHists_110X[hist]->GetYaxis()->SetTitleOffset(1.5);
    leg->SetBorderSize(0);
    leg->Draw();
    canvases.back()->Print(Form("plots/%sLLPOverlay.pdf", hist.substr(0).c_str()));
  }
   
  return 0;
}
