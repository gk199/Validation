#include "PhysicsTools/Utilities/macros/setTDRStyle.C"

#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include "TLegend.h"
#include "TROOT.h"

#include <map>
#include <string>
#include <vector>

#include <iostream>

int main()
{
  // include comparisons between HW and data TPs
  //  bool includeHW = false;
  //  int rebinFactor = 1;

  setTDRStyle();
  gROOT->ForceStyle();

  // default, then new conditions
  std::vector<std::string> filenames = {"rates_new_cond_110X_nugun.root"};
  std::vector<std::string> multTypes = {"HTdistribution","HTdistribution_trig","DeltaR_L1_delayed_seed","DeltaR_L1_prompt_seed","DeltaR_L1_delayed_hit","DeltaR_L1_prompt_hit","Mult_delayed_hit","Mult_delayed_hit_jetET","Mult_delayed_hit_promptV","Mult_delayed_hit_jetETpromptV","Mult_prompt_hit"};//,"mhit1","mhit2","mhit3"}; // list the names of the TH1Fs that will be overlayed

  std::map<std::string, TH1F*> multHists_nugun;

  std::vector<TFile*> files;
  for(auto file : filenames) {
    files.push_back(TFile::Open(file.c_str()));
  }

  for(auto multType : multTypes) {
    std::string histName(multType);
    histName += "_emu";

    multHists_nugun[multType] = dynamic_cast<TH1F*>(files.at(0)->Get(histName.c_str()));
  }

  for(auto pair : multHists_nugun) pair.second->SetLineWidth(2);

  std::vector<TCanvas*> canvases;
  std::vector<TPad*> pad1;
  std::vector<TPad*> pad2;

  for(auto hist : multTypes) {
    canvases.push_back(new TCanvas);
    canvases.back()->SetWindowSize(canvases.back()->GetWw(), 1.3*canvases.back()->GetWh());
    pad1.push_back(new TPad("pad1", "pad1", 0, 0, 1, 1));
    pad1.back()->SetGrid();
    pad1.back()->Draw();
    pad1.back()->cd();

    multHists_nugun[hist]->SetLineColor(kBlack); 
    TString name(multHists_nugun[hist]->GetName());
    //    multHists_nugun[hist]->SetFillStyle(3005);
    multHists_nugun[hist]->SetLineColor(kBlack);
    multHists_nugun[hist]->Scale(1./multHists_nugun[hist]->Integral());
    multHists_nugun[hist]->Draw("hist pfc");
    multHists_nugun[hist]->SetAxisRange(0,1,"Y");
    TLegend *leg = new TLegend(0.6, 0.75, 0.95, 0.9);
    leg->AddEntry(multHists_nugun[hist],"NuGun", "F");
    multHists_nugun[hist]->GetYaxis()->CenterTitle(true);
    multHists_nugun[hist]->GetXaxis()->SetLabelSize(0.03);
    multHists_nugun[hist]->GetYaxis()->SetLabelSize(0.03);
    multHists_nugun[hist]->GetXaxis()->SetTitleSize(0.04);
    multHists_nugun[hist]->GetYaxis()->SetTitleSize(0.04);
    multHists_nugun[hist]->GetXaxis()->SetTitleOffset(1.2);
    multHists_nugun[hist]->GetYaxis()->SetTitleOffset(1.5);
    leg->SetBorderSize(0);
    leg->Draw();
    canvases.back()->Print(Form("plots/%sOverlay_nugun.pdf", hist.substr(0).c_str()));
  }

  for(auto hist : multTypes) {
    canvases.push_back(new TCanvas);
    canvases.back()->SetWindowSize(canvases.back()->GetWw(), 1.3*canvases.back()->GetWh());
    pad1.push_back(new TPad("pad1", "pad1", 0, 0, 1, 1));
    pad1.back()->SetGrid();
    pad1.back()->Draw();
    pad1.back()->cd();

    multHists_nugun[hist]->SetLineColor(kBlack);
    TString name(multHists_nugun[hist]->GetName());
    //    multHists_nugun[hist]->SetFillStyle(3005);
    multHists_nugun[hist]->SetLineColor(kBlack);
    multHists_nugun[hist]->Scale(1./multHists_nugun[hist]->Integral());
    multHists_nugun[hist]->Draw("hist pfc");
    gPad->SetLogy();
    TLegend *leg = new TLegend(0.6, 0.75, 0.95, 0.9);
    leg->AddEntry(multHists_nugun[hist],"NuGun", "F");
    multHists_nugun[hist]->GetYaxis()->CenterTitle(true);
    multHists_nugun[hist]->GetXaxis()->SetLabelSize(0.03);
    multHists_nugun[hist]->GetYaxis()->SetLabelSize(0.03);
    multHists_nugun[hist]->GetXaxis()->SetTitleSize(0.04);
    multHists_nugun[hist]->GetYaxis()->SetTitleSize(0.04);
    multHists_nugun[hist]->GetXaxis()->SetTitleOffset(1.2);
    multHists_nugun[hist]->GetYaxis()->SetTitleOffset(1.5);
    leg->SetBorderSize(0);
    leg->Draw();
    canvases.back()->Print(Form("plots/%sOverlay_Log_nugun.pdf", hist.substr(0).c_str()));
  }
  return 0;
}
