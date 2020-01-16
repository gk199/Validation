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
  bool includeHW = false;
  int rebinFactor = 1;

  setTDRStyle();
  gROOT->ForceStyle();

  // default, then new conditions
  // files for L1 rates
  //  std::vector<std::string> filenames = {"rates_def.root", "rates_new_cond.root"};
  std::vector<std::string> filenames = {"rates_new_cond_QCD.root", "rates_new_cond_pl1000.root"};
  std::vector<std::string> rateTypes = {"singleJet", "doubleJet", "tripleJet", "quadJet",
					"singleEg", "singleISOEg", "doubleEg", "doubleISOEg",
					"singleTau", "singleISOTau", "doubleTau", "doubleISOTau",
					"htSum", "etSum", "metSum", "metHFSum"};
  // files for multiplicity overlay plots
  std::vector<std::string> mult_filenames = {"rates_new_cond_pl10000.root", "rates_new_cond_pl1000.root", "rates_new_cond_pl500.root", "rates_new_cond_QCD.root"};
  std::vector<std::string> multTypes = {"dt3GeV1ns","dt3GeV2ns","dt3GeV3ns","dt3GeV4ns","dt3GeV5ns","dt2GeV1ns","dt2GeV2ns","dt2GeV3ns","dt2GeV4ns","dt2GeV5ns","dt1GeV1ns","dt1GeV2ns","dt1GeV3ns","dt1GeV4ns","dt1GeV5ns"};

  std::map<std::string, int> histColor;
  histColor["singleJet"] = histColor["singleEg"] = histColor["singleTau"] = histColor["etSum"] = histColor["metSum"] = histColor["dt3GeV1ns"] = histColor["dt2GeV1ns"] = histColor["dt1GeV1ns"] = kRed;
  histColor["doubleJet"] = histColor["singleISOEg"] = histColor["singleISOTau"] = histColor["htSum"] = histColor["metHFSum"] = histColor["dt3GeV2ns"] = histColor["dt2GeV2ns"] = histColor["dt1GeV2ns"] = kBlue;
  histColor["tripleJet"] = histColor["doubleEg"] = histColor["doubleTau"] = histColor["dt3GeV3ns"] = histColor["dt2GeV3ns"] = histColor["dt1GeV3ns"] = kGreen;
  histColor["quadJet"] = histColor["doubleISOEg"] = histColor["doubleISOTau"] = histColor["dt3GeV4ns"] = histColor["dt2GeV4ns"] = histColor["dt1GeV4ns"] = kBlack;
  histColor["dt3GeV5ns"] = histColor["dt2GeV5ns"] = histColor["dt1GeV5ns"] = kCyan;

  std::map<std::string, TH1F*> rateHists_def;
  std::map<std::string, TH1F*> rateHists_new_cond;
  std::map<std::string, TH1F*> rateHists_hw;
  std::map<std::string, TH1F*> rateHistsRatio;
  
  std::map<std::string, TH1F*> multHists_QCD;
  std::map<std::string, TH1F*> multHists_LLP10000;
  std::map<std::string, TH1F*> multHists_LLP1000;
  std::map<std::string, TH1F*> multHists_LLP500;
  std::map<std::string, TH1F*> multHistsRatio;
  std::map<std::string, TH1F*> multHists_hw;

  std::vector<TFile*> files;
  for(auto file : filenames) {
    files.push_back(TFile::Open(file.c_str()));
  }
  // making rate plots for current and new conditions
  for(auto rateType : rateTypes) {
    std::string histName(rateType);
    std::string histNameHw(histName);
    histName += "Rates_emu";
    histNameHw += "Rates_hw";
    rateHists_def[rateType] = dynamic_cast<TH1F*>(files.at(0)->Get(histName.c_str()));
    rateHists_hw[rateType] = dynamic_cast<TH1F*>(files.at(0)->Get(histNameHw.c_str()));
    rateHists_new_cond[rateType] = dynamic_cast<TH1F*>(files.at(1)->Get(histName.c_str())); 
    rateHists_def[rateType]->Rebin(rebinFactor);
    rateHists_hw[rateType]->Rebin(rebinFactor);
    rateHists_new_cond[rateType]->Rebin(rebinFactor);

    rateHists_def[rateType]->SetLineColor(histColor[rateType]);
    rateHists_hw[rateType]->SetLineColor(histColor[rateType]);
    rateHists_new_cond[rateType]->SetLineColor(histColor[rateType]);
    TString name(rateHists_new_cond[rateType]->GetName());
    name += "_ratio";
    if(includeHW) {
      rateHistsRatio[rateType] = dynamic_cast<TH1F*>(rateHists_def[rateType]->Clone(name));
      rateHistsRatio[rateType]->Divide(rateHists_hw[rateType]);
    }
    else {
      rateHistsRatio[rateType] = dynamic_cast<TH1F*>(rateHists_new_cond[rateType]->Clone(name));
      rateHistsRatio[rateType]->Divide(rateHists_def[rateType]);
    }
    rateHistsRatio[rateType]->SetMinimum(0.6);    
    rateHistsRatio[rateType]->SetMaximum(1.4);    
    rateHistsRatio[rateType]->SetLineWidth(2);    
  }

  // opening the files for multiplicity plots     
  std::vector<TFile*> mult_files;
  for(auto file : mult_filenames) {
    mult_files.push_back(TFile::Open(file.c_str()));
  }
  for(auto multType : multTypes) {
    std::string histName(multType);
    std::string histNameHw(histName);
    histName += "Mult_emu";
    //    histNameHw += "Mult_hw";
    multHists_QCD[multType]  = dynamic_cast<TH1F*>(mult_files.at(3)->Get(histName.c_str()));
    //    multHists_hw[multType]  = dynamic_cast<TH1F*>(mult_files.at(3)->Get(histNameHw.c_str()));
    multHists_LLP10000[multType] = dynamic_cast<TH1F*>(mult_files.at(0)->Get(histName.c_str()));
    multHists_LLP1000[multType] = dynamic_cast<TH1F*>(mult_files.at(1)->Get(histName.c_str()));
    multHists_LLP500[multType] = dynamic_cast<TH1F*>(mult_files.at(2)->Get(histName.c_str()));
    
    multHists_QCD[multType]->Rebin(rebinFactor);
    multHists_LLP10000[multType]->Rebin(rebinFactor);
    multHists_LLP10000[multType]->Rebin(rebinFactor);
    multHists_LLP10000[multType]->Rebin(rebinFactor);

    multHists_QCD[multType]->SetLineColor(histColor[multType]);
    //    multHists_hw[multType]->SetLineColor(histColor[multType]);
    multHists_LLP10000[multType]->SetLineColor(histColor[multType]);
    multHists_LLP1000[multType]->SetLineColor(histColor[multType]);
    multHists_LLP500[multType]->SetLineColor(histColor[multType]);
  }

  for(auto pair : rateHists_new_cond) pair.second->SetLineWidth(2);
  for(auto pair : rateHists_hw) pair.second->SetLineStyle(kDashed);
  for(auto pair : rateHists_def) pair.second->SetLineStyle(kDotted);

  for(auto pair : multHists_LLP10000) pair.second->SetLineWidth(2);
  for(auto pair : multHists_LLP1000) pair.second->SetLineWidth(2);
  for(auto pair : multHists_LLP500) pair.second->SetLineWidth(2);
  //  for(auto pair : multHists_hw) pair.second->SetLineStyle(kDashed);
  //  for(auto pair : multHists_QCD) pair.second->SetLineStyle(kDotted); // analogous to rateHist_def
  for(auto pair : multHists_QCD) pair.second->SetLineWidth(2);

  std::vector<std::string> jetPlots = {"singleJet", "doubleJet", "tripleJet", "quadJet"};
  std::vector<std::string> egPlots = {"singleEg", "singleISOEg", "doubleEg", "doubleISOEg"};
  std::vector<std::string> tauPlots = {"singleTau", "singleISOTau", "doubleTau", "doubleISOTau"};
  std::vector<std::string> scalarSumPlots = {"etSum", "htSum"};
  std::vector<std::string> vectorSumPlots = {"metSum", "metHFSum"};
  // multiplicity plot types 
  std::vector<std::string> multPlots3GeV = {"dt3GeV1ns","dt3GeV2ns","dt3GeV3ns","dt3GeV4ns","dt3GeV5ns"};
  std::vector<std::string> multPlots2GeV = {"dt2GeV1ns","dt2GeV2ns","dt2GeV3ns","dt2GeV4ns","dt2GeV5ns"};
  std::vector<std::string> multPlots1GeV = {"dt1GeV1ns","dt1GeV2ns","dt1GeV3ns","dt1GeV4ns","dt1GeV5ns"}; 
  // used for overlays
  std::vector<std::string> o31 = {"dt3GeV1ns"};
  std::vector<std::string> o32 = {"dt3GeV2ns"};
  std::vector<std::string> o33 = {"dt3GeV3ns"};
  std::vector<std::string> o34 = {"dt3GeV4ns"};
  std::vector<std::string> o35 = {"dt3GeV5ns"};
  std::vector<std::string> o21 = {"dt2GeV1ns"};
  std::vector<std::string> o22 = {"dt2GeV2ns"};
  std::vector<std::string> o23 = {"dt2GeV3ns"};
  std::vector<std::string> o24 = {"dt2GeV4ns"};
  std::vector<std::string> o25 = {"dt2GeV5ns"};
  std::vector<std::string> o11 = {"dt1GeV1ns"};
  std::vector<std::string> o12 = {"dt1GeV2ns"};
  std::vector<std::string> o13 = {"dt1GeV3ns"};
  std::vector<std::string> o14 = {"dt1GeV4ns"};
  std::vector<std::string> o15 = {"dt1GeV5ns"};

  std::vector<TCanvas*> canvases;
  std::vector<TPad*> pad1;
  std::vector<TPad*> pad2;
  std::map<std::string, std::vector<std::string> > plots;
  plots["jet"] = jetPlots;
  plots["eg"] = egPlots;
  plots["tau"] = tauPlots;
  plots["scalarSum"] = scalarSumPlots;
  plots["vectorSum"] = vectorSumPlots;

  std::map<std::string, std::vector<std::string> > mult_plots;
  mult_plots["mult3"] = multPlots3GeV;
  mult_plots["mult2"] = multPlots2GeV;
  mult_plots["mult1"] = multPlots1GeV;

  std::map<std::string, std::vector<std::string> > overlays;
  overlays["3GeV1ns"] = o31;
  overlays["3GeV2ns"] = o32;
  overlays["3GeV3ns"] = o33;
  overlays["3GeV4ns"] = o34;
  overlays["3GeV5ns"] = o35;
  overlays["2GeV1ns"] = o21;
  overlays["2GeV2ns"] = o22;
  overlays["2GeV3ns"] = o23;
  overlays["2GeV4ns"] = o24;
  overlays["2GeV5ns"] = o25;
  overlays["1GeV1ns"] = o11;
  overlays["1GeV2ns"] = o12;
  overlays["1GeV3ns"] = o13;
  overlays["1GeV4ns"] = o14;
  overlays["1GeV5ns"] = o15;

  // looping through all plot collections (jets, eg, tau, scalar, vector)
  for(auto iplot : plots) {
    canvases.push_back(new TCanvas);
    canvases.back()->SetWindowSize(canvases.back()->GetWw(), 1.3*canvases.back()->GetWh());
    pad1.push_back(new TPad("pad1", "pad1", 0, 0.3, 1, 1));
    pad1.back()->SetLogy();
    pad1.back()->SetGrid();
    pad1.back()->Draw();
    pad2.push_back(new TPad("pad2", "pad2", 0, 0, 1, 0.3));
    pad2.back()->SetGrid();
    pad2.back()->Draw();
    
    pad1.back()->cd();
    
    rateHists_def[iplot.second.front()]->Draw("hist");
    TLegend *leg = new TLegend(0.55, 0.9 - 0.1*iplot.second.size(), 0.95, 0.93);
    for(auto hist : iplot.second) {
      rateHists_def[hist]->Draw("hist same");
      rateHists_def[hist]->GetYaxis()->SetRangeUser(10.01, 100000000); // setting the range of the Y axis to show low rates
      if(includeHW) rateHists_hw[hist]->Draw("hist same");
      rateHists_new_cond[hist]->Draw("hist same");
      rateHists_def[hist]->GetYaxis()->SetRangeUser(10.01, 100000000); // setting the range of the Y axis to show low rates
      TString name(rateHists_def[hist]->GetName());
      TString nameHw(rateHists_hw[hist]->GetName());
      std::cout << name << std::endl;
      leg->AddEntry(rateHists_def[hist], name + " (current)", "L");
      if(includeHW) leg->AddEntry(rateHists_hw[hist], name + " (hw)", "L");
      leg->AddEntry(rateHists_new_cond[hist], name + " (new)", "L"); 
    }
    leg->SetBorderSize(0);
    leg->Draw();
    
    pad2.back()->cd();
    rateHistsRatio[iplot.second.front()]->Draw("hist");
    if(includeHW) rateHistsRatio[iplot.second.front()]->GetYaxis()->SetTitle("Current/HW");
    else rateHistsRatio[iplot.second.front()]->GetYaxis()->SetTitle("New/Current");
    for(auto hist : iplot.second) {
      rateHistsRatio[hist]->Draw("hist same");
    }

    if(includeHW) canvases.back()->Print(Form("plots/%sRates_hw.pdf", iplot.first.c_str()));
    //    else canvases.back()->Print(Form("plots/%sRates_emu.pdf", iplot.first.c_str()));
  }

  // multiplicity plot loop
  // QCD
  for (auto iplot : mult_plots){
    canvases.push_back(new TCanvas);
    canvases.back()->SetWindowSize(canvases.back()->GetWw(), 1.3*canvases.back()->GetWh());
    pad1.push_back(new TPad("pad1", "pad1", 0, 0.15, 1, 1));
    //    pad1.back()->SetLogy();
    pad1.back()->SetGrid();
    pad1.back()->Draw();
    pad1.back()->cd();
    multHists_QCD[iplot.second.front()]->Draw("hist");
    TLegend *leg = new TLegend(0.65, 1.1 - 0.1*iplot.second.size(), 0.95, 0.93);
    for (auto hist : iplot.second) {
      multHists_QCD[hist]->GetYaxis()->SetRangeUser(0,2500000);
      // rebin histograms for 2 GeV energy cut, as the x-axis extends further as compared to 3 GeV
      if ((hist == "dt2GeV1ns") || (hist == "dt2GeV2ns") || (hist == "dt2GeV3ns") || (hist == "dt2GeV4ns") || (hist == "dt2GeV5ns")){
	multHists_QCD[hist]->Rebin(rebinFactor*2);
      }
      // rebin histograms for 1 GeV energy cut, as the x-axis extends further here
      if ((hist == "dt1GeV1ns") || (hist == "dt1GeV2ns") || (hist == "dt1GeV3ns") || (hist == "dt1GeV4ns") || (hist == "dt1GeV5ns")){
        multHists_QCD[hist]->Rebin(rebinFactor*4);
      }
      multHists_QCD[hist]->Draw("hist same");
      TString name(multHists_QCD[hist]->GetName());
      leg->AddEntry(multHists_QCD[hist], name(2,7) + " (QCD)", "L");
      multHists_QCD[hist]->SetTitle("Multiplicity for QCD, timing scan at " + name(2,4));
    }
    leg->SetBorderSize(0);
    leg->Draw();
    canvases.back()->Print(Form("plots/%sMult_emu_QCD.pdf", iplot.first.c_str()));
  }

  // LLP 10000
  for (auto iplot : mult_plots){
    canvases.push_back(new TCanvas);
    canvases.back()->SetWindowSize(canvases.back()->GetWw(), 1.3*canvases.back()->GetWh());
    pad1.push_back(new TPad("pad1", "pad1", 0, 0.15, 1, 1));
    //    pad1.back()->SetLogy();
    pad1.back()->SetGrid();
    pad1.back()->Draw();
    pad1.back()->cd();
    multHists_LLP10000[iplot.second.front()]->Draw("hist");
    TLegend *leg = new TLegend(0.65, 1.1 - 0.1*iplot.second.size(), 0.95, 0.93);
    for (auto hist : iplot.second) {
      multHists_LLP10000[hist]->GetYaxis()->SetRangeUser(0,2500000);
      // rebin histograms for 2 GeV energy cut, as the x-axis extends further as compared to 3 GeV                                                                                                                            
      if ((hist == "dt2GeV1ns") || (hist == "dt2GeV2ns") || (hist == "dt2GeV3ns") || (hist == "dt2GeV4ns") || (hist == "dt2GeV5ns")){
        multHists_LLP10000[hist]->Rebin(rebinFactor*2);
      }
      // rebin histograms for 1 GeV energy cut, as the x-axis extends further here                                                                                                           
      if ((hist == "dt1GeV1ns") || (hist == "dt1GeV2ns") || (hist == "dt1GeV3ns") || (hist == "dt1GeV4ns") || (hist == "dt1GeV5ns")){
        multHists_LLP10000[hist]->Rebin(rebinFactor*4);
      }
      multHists_LLP10000[hist]->Draw("hist same");
      TString name(multHists_QCD[hist]->GetName());
      leg->AddEntry(multHists_LLP10000[hist], name(2,7) + " (LLP, pl=10m)", "L");
      multHists_LLP10000[hist]->SetTitle("Multiplicity for LLP pl=10m, timing scan at " + name(2,4));
    }
    leg->SetBorderSize(0);
    leg->Draw();
    canvases.back()->Print(Form("plots/%sMult_emu_LLP10000.pdf", iplot.first.c_str()));
  }

  // LLP 1000
  for (auto iplot : mult_plots){
    canvases.push_back(new TCanvas);
    canvases.back()->SetWindowSize(canvases.back()->GetWw(), 1.3*canvases.back()->GetWh());
    pad1.push_back(new TPad("pad1", "pad1", 0, 0.15, 1, 1));
    //    pad1.back()->SetLogy();
    pad1.back()->SetGrid();
    pad1.back()->Draw();
    pad1.back()->cd();
    multHists_LLP1000[iplot.second.front()]->Draw("hist");
    TLegend *leg = new TLegend(0.65, 1.1 - 0.1*iplot.second.size(), 0.95, 0.93);
    for (auto hist : iplot.second) {
      multHists_LLP1000[hist]->GetYaxis()->SetRangeUser(0,2500000);
      // rebin histograms for 2 GeV energy cut, as the x-axis extends further as compared to 3 GeV
      if ((hist == "dt2GeV1ns") || (hist == "dt2GeV2ns") || (hist == "dt2GeV3ns") || (hist == "dt2GeV4ns") || (hist == "dt2GeV5ns")){
        multHists_LLP1000[hist]->Rebin(rebinFactor*2);
      }
      // rebin histograms for 1 GeV energy cut, as the x-axis extends further here
      if ((hist == "dt1GeV1ns") || (hist == "dt1GeV2ns") || (hist == "dt1GeV3ns") || (hist == "dt1GeV4ns") || (hist == "dt1GeV5ns")){
        multHists_LLP1000[hist]->Rebin(rebinFactor*4);
      }
      multHists_LLP1000[hist]->Draw("hist same");
      TString name(multHists_QCD[hist]->GetName());
      leg->AddEntry(multHists_LLP1000[hist], name(2,7) + " (LLP, pl=1m)", "L");
      multHists_LLP1000[hist]->SetTitle("Multiplicity for LLP pl=1m, timing scan at " + name(2,4));
    }
    leg->SetBorderSize(0);
    leg->Draw();
    canvases.back()->Print(Form("plots/%sMult_emu_LLP1000.pdf", iplot.first.c_str()));
  }

  // LLP 500                                                                
  for (auto iplot : mult_plots){
    canvases.push_back(new TCanvas);
    canvases.back()->SetWindowSize(canvases.back()->GetWw(), 1.3*canvases.back()->GetWh());
    pad1.push_back(new TPad("pad1", "pad1", 0, 0.15, 1, 1));
    //    pad1.back()->SetLogy();
    pad1.back()->SetGrid();
    pad1.back()->Draw();
    pad1.back()->cd();
    multHists_LLP500[iplot.second.front()]->Draw("hist");
    TLegend *leg = new TLegend(0.65, 1.1 - 0.1*iplot.second.size(), 0.95, 0.93);
    for (auto hist : iplot.second) {
      multHists_LLP500[hist]->GetYaxis()->SetRangeUser(0,2500000);
      // rebin histograms for 2 GeV energy cut, as the x-axis extends further as compared to 3 GeV
      if ((hist == "dt2GeV1ns") || (hist == "dt2GeV2ns") || (hist == "dt2GeV3ns") || (hist == "dt2GeV4ns") || (hist == "dt2GeV5ns")){
        multHists_LLP500[hist]->Rebin(rebinFactor*2);
      }
      // rebin histograms for 1 GeV energy cut, as the x-axis extends further here
      if ((hist == "dt1GeV1ns") || (hist == "dt1GeV2ns") || (hist == "dt1GeV3ns") || (hist == "dt1GeV4ns") || (hist == "dt1GeV5ns")){
        multHists_LLP500[hist]->Rebin(rebinFactor*4);
      }
      multHists_LLP500[hist]->Draw("hist same");
      TString name(multHists_QCD[hist]->GetName());
      leg->AddEntry(multHists_LLP500[hist], name(2,7) + " (LLP, pl=0.5m)", "L");
      multHists_LLP500[hist]->SetTitle("Multiplicity for LLP pl=0.5m, timing scan at " + name(2,4));
    }
    leg->SetBorderSize(0);
    leg->Draw();
    canvases.back()->Print(Form("plots/%sMult_emu_LLP500.pdf", iplot.first.c_str()));
  }

  // overlay LLP and QCD at a single energy / timing cut value
  for (auto iplot : overlays ) {
    canvases.push_back(new TCanvas);
    canvases.back()->SetWindowSize(canvases.back()->GetWw(), 1.3*canvases.back()->GetWh());
    pad1.push_back(new TPad("pad1", "pad1", 0, 0.3, 1, 1));
    pad1.back()->SetGrid();
    pad1.back()->Draw();
    pad1.back()->cd();
    multHists_QCD[iplot.second.front()]->Draw("hist");
    TLegend *leg = new TLegend(0.55, 0.9 - 0.1*iplot.second.size(), 0.95, 0.93);
    for (auto hist : iplot.second){
      multHists_QCD[hist]->GetYaxis()->SetRangeUser(0,1000000);
      multHists_QCD[hist]->SetLineColor(kBlack);
      multHists_QCD[hist]->Draw("hist same");
      multHists_LLP500[hist]->SetLineColor(kBlue);
      multHists_LLP500[hist]->Draw("hist same");
      multHists_LLP1000[hist]->SetLineColor(kGreen);
      multHists_LLP1000[hist]->Draw("hist same");
      multHists_LLP10000[hist]->SetLineColor(kRed);
      multHists_LLP10000[hist]->Draw("hist same");
      TString name (multHists_QCD[hist]->GetName());
      leg->AddEntry(multHists_QCD[hist],"QCD", "L");
      leg->AddEntry(multHists_LLP500[hist],"LLP, pl=0.5m", "L");
      leg->AddEntry(multHists_LLP1000[hist], "LLP, pl=1m", "L");
      leg->AddEntry(multHists_LLP10000[hist], "LLP, pl=10m", "L");
      multHists_QCD[hist]->SetTitle("Multiplicity Overlay of QCD and LLPs at " + name(2,4) + " and " + name(6,3));
    }
    leg->SetBorderSize(0);
    leg->Draw();
    canvases.back()->Print(Form("plots/%sOverlay.pdf", iplot.first.c_str()));
  }
  return 0;
}
