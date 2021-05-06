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
  std::vector<std::string> filenames = {"rates_new_cond_QCD.root", "rates_new_cond_LLP_pl500.root", "rates_new_cond_LLP_pl1000.root", "rates_new_cond_LLP_pl10000.root"};
  std::vector<std::string> multTypes = {"HTdistribution","HTdistribution_trig","DeltaR_L1_delayed_seed","DeltaR_L1_prompt_seed","DeltaR_L1_delayed_hit","DeltaR_L1_prompt_hit","Mult_delayed_hit","Mult_delayed_hit_jetET","Mult_delayed_hit_promptV","Mult_delayed_hit_jetETpromptV","Mult_prompt_hit"};

  std::map<std::string, TH1F*> multHists_QCD;
  std::map<std::string, TH1F*> multHists_LLPpl500;
  std::map<std::string, TH1F*> multHists_LLPpl1000;
  std::map<std::string, TH1F*> multHists_LLPpl10000;

  std::vector<TFile*> files;
  for(auto file : filenames) {
    files.push_back(TFile::Open(file.c_str()));
  }

  for(auto multType : multTypes) {
    std::string histName(multType);
    histName += "_emu";

    multHists_QCD[multType] = dynamic_cast<TH1F*>(files.at(0)->Get(histName.c_str()));
    multHists_LLPpl500[multType] = dynamic_cast<TH1F*>(files.at(1)->Get(histName.c_str()));
    multHists_LLPpl1000[multType] = dynamic_cast<TH1F*>(files.at(2)->Get(histName.c_str()));
    multHists_LLPpl10000[multType] = dynamic_cast<TH1F*>(files.at(3)->Get(histName.c_str()));
  }

  for(auto pair : multHists_QCD) pair.second->SetLineWidth(2);
  for(auto pair : multHists_LLPpl500) pair.second->SetLineWidth(2);
  for(auto pair : multHists_LLPpl1000) pair.second->SetLineWidth(2);
  for(auto pair : multHists_LLPpl10000) pair.second->SetLineWidth(2);

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

    multHists_QCD[hist]->SetLineColor(kBlack); 
    TString name(multHists_QCD[hist]->GetName());
    multHists_QCD[hist]->SetFillStyle(3005);
    multHists_QCD[hist]->Scale(1./multHists_QCD[hist]->Integral());
    multHists_QCD[hist]->Draw("hist pfc");
    multHists_QCD[hist]->SetAxisRange(0,1,"Y");
    TLegend *leg = new TLegend(0.6, 0.65, 0.95, 0.9);
    multHists_LLPpl500[hist]->SetLineColor(kBlue);
    multHists_LLPpl500[hist]->Scale(1./multHists_LLPpl500[hist]->Integral());
    multHists_LLPpl500[hist]->Draw("hist same");
    multHists_LLPpl1000[hist]->SetLineColor(kGreen+1);
    multHists_LLPpl1000[hist]->Scale(1./multHists_LLPpl1000[hist]->Integral());
    multHists_LLPpl1000[hist]->Draw("hist same");
    multHists_LLPpl10000[hist]->SetLineColor(kRed);
    multHists_LLPpl10000[hist]->Scale(1./multHists_LLPpl10000[hist]->Integral());
    multHists_LLPpl10000[hist]->Draw("hist same");
    leg->AddEntry(multHists_QCD[hist],"QCD", "F");
    leg->AddEntry(multHists_LLPpl500[hist], "LLP, mh=250, c#scale[1.2]{#tau}=0.5m", "L");
    leg->AddEntry(multHists_LLPpl1000[hist], "LLP, mh=250, c#scale[1.2]{#tau}=1m", "L");
    leg->AddEntry(multHists_LLPpl10000[hist], "LLP, mh=250, c#scale[1.2]{#tau}=10m", "L");
    multHists_QCD[hist]->GetYaxis()->CenterTitle(true);
    multHists_QCD[hist]->GetXaxis()->SetLabelSize(0.03);
    multHists_QCD[hist]->GetYaxis()->SetLabelSize(0.03);
    multHists_QCD[hist]->GetXaxis()->SetTitleSize(0.04);
    multHists_QCD[hist]->GetYaxis()->SetTitleSize(0.04);
    multHists_QCD[hist]->GetXaxis()->SetTitleOffset(1.2);
    multHists_QCD[hist]->GetYaxis()->SetTitleOffset(1.5);
    leg->SetBorderSize(0);
    leg->Draw();
    canvases.back()->Print(Form("plots/%sOverlay.pdf", hist.substr(0).c_str()));
  }

  for(auto hist : multTypes) {
    canvases.push_back(new TCanvas);
    canvases.back()->SetWindowSize(canvases.back()->GetWw(), 1.3*canvases.back()->GetWh());
    pad1.push_back(new TPad("pad1", "pad1", 0, 0, 1, 1));
    pad1.back()->SetGrid();
    pad1.back()->Draw();
    pad1.back()->cd();
    multHists_QCD[hist]->SetLineColor(kBlack);
    TString name(multHists_QCD[hist]->GetName());
    multHists_QCD[hist]->SetFillStyle(3005);
    multHists_QCD[hist]->Scale(1./multHists_QCD[hist]->Integral());
    multHists_QCD[hist]->Draw("hist pfc");
    gPad->SetLogy();                                                                                                                                                              
    TLegend *leg = new TLegend(0.6, 0.65, 0.95, 0.9);
    multHists_LLPpl500[hist]->SetLineColor(kBlue);
    multHists_LLPpl500[hist]->Scale(1./multHists_LLPpl500[hist]->Integral());
    multHists_LLPpl500[hist]->Draw("hist same");
    multHists_LLPpl1000[hist]->SetLineColor(kGreen+1);
    multHists_LLPpl1000[hist]->Scale(1./multHists_LLPpl1000[hist]->Integral());
    multHists_LLPpl1000[hist]->Draw("hist same");
    multHists_LLPpl10000[hist]->SetLineColor(kRed);
    multHists_LLPpl10000[hist]->Scale(1./multHists_LLPpl10000[hist]->Integral());
    multHists_LLPpl10000[hist]->Draw("hist same");
    leg->AddEntry(multHists_QCD[hist],"QCD", "F");
    leg->AddEntry(multHists_LLPpl500[hist], "LLP, mh=250, c#scale[1.2]{#tau}=0.5m", "L");
    leg->AddEntry(multHists_LLPpl1000[hist], "LLP, mh=250, c#scale[1.2]{#tau}=1m", "L");
    leg->AddEntry(multHists_LLPpl10000[hist], "LLP, mh=250, c#scale[1.2]{#tau}=10m", "L");
    multHists_QCD[hist]->GetYaxis()->CenterTitle(true);
    multHists_QCD[hist]->GetXaxis()->SetLabelSize(0.03);
    multHists_QCD[hist]->GetYaxis()->SetLabelSize(0.03);
    multHists_QCD[hist]->GetXaxis()->SetTitleSize(0.04);
    multHists_QCD[hist]->GetYaxis()->SetTitleSize(0.04);
    multHists_QCD[hist]->GetXaxis()->SetTitleOffset(1.2);
    multHists_QCD[hist]->GetYaxis()->SetTitleOffset(1.5);
    leg->SetBorderSize(0);
    leg->Draw();
    canvases.back()->Print(Form("plots/%sOverlay_Log.pdf", hist.substr(0).c_str()));
  }
   
  return 0;
}
