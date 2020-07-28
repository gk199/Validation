#include<iostream>
#include<fstream>
#include<vector>
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMarker.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TFile.h"


using namespace std;

int main() {
  // read in signal and background files 
  // first line is efficiency at HT>120 and mult HEHB >= 10
  // second line is efficiency at HT>360 with no timing restrictions
  // Signals
  // mh=125 pl=0.5m
  double signal_mh125_pl500[15];
  ifstream mh125_pl500;
  mh125_pl500.open("Efficiency_HtBins_Signal_mh125__mx50__pl500__.txt");
  int n=0;
  while (mh125_pl500 >> signal_mh125_pl500[n]) n++;
  mh125_pl500.close();
  // mh=125 pl=1m
  double signal_mh125_pl1000[15];
  ifstream mh125_pl1000;
  mh125_pl1000.open("Efficiency_HtBins_Signal_mh125__mx50__pl1000_.txt");
  n=0;
  while (mh125_pl1000 >> signal_mh125_pl1000[n]) n++;
  mh125_pl1000.close();
  // mh=125 pl=10m
  double signal_mh125_pl10000[15];
  ifstream mh125_pl10000;
  mh125_pl10000.open("Efficiency_HtBins_Signal_mh125__mx50__pl10000.txt");
  n=0;
  while (mh125_pl10000 >> signal_mh125_pl10000[n]) n++;
  mh125_pl10000.close();

  // mh=1000 pl=0.5m
  double signal_mh1000_pl500[15];
  ifstream mh1000_pl500;
  mh1000_pl500.open("Efficiency_HtBins_Signal_mh1000_mx450_pl500__.txt");
  n=0;
  while (mh1000_pl500 >> signal_mh1000_pl500[n]) n++;
  mh1000_pl500.close();
  // mh=1000 pl=1m 
  double signal_mh1000_pl1000[15];
  ifstream mh1000_pl1000;
  mh1000_pl1000.open("Efficiency_HtBins_Signal_mh1000_mx450_pl1000_.txt");
  n=0;
  while (mh1000_pl1000 >> signal_mh1000_pl1000[n]) n++;
  mh1000_pl1000.close();
  // mh=1000 pl=10m
  double signal_mh1000_pl10000[15];
  ifstream mh1000_pl10000;
  mh1000_pl10000.open("Efficiency_HtBins_Signal_mh1000_mx450_pl10000.txt");
  n=0;
  while (mh1000_pl10000 >> signal_mh1000_pl10000[n]) n++;
  mh1000_pl10000.close();

  // mh=250 pl=1m 
  double signal_mh250_pl1000[15];
  ifstream mh250_pl1000;
  mh250_pl1000.open("Efficiency_HtBins_Signal_mh250__mx120_pl1000_.txt");
  n=0;
  while (mh250_pl1000 >> signal_mh250_pl1000[n]) n++;
  mh250_pl1000.close();

  // mh=350 pl=0.5m 
  double signal_mh350_pl500[15];
  ifstream mh350_pl500;
  mh350_pl500.open("Efficiency_HtBins_Signal_mh350__mx160_pl500__.txt");
  n=0;
  while (mh350_pl500 >> signal_mh350_pl500[n]) n++;
  mh350_pl500.close();
  // mh=350 pl=1m
  double signal_mh350_pl1000[15];
  ifstream mh350_pl1000;
  mh350_pl1000.open("Efficiency_HtBins_Signal_mh350__mx160_pl1000_.txt");
  n=0;
  while (mh350_pl1000 >> signal_mh350_pl1000[n]) n++;
  mh350_pl1000.close();
  // mh=350 pl=10m
  double signal_mh350_pl10000[15];
  ifstream mh350_pl10000;
  mh350_pl10000.open("Efficiency_HtBins_Signal_mh350__mx160_pl10000.txt");
  n=0;
  while (mh350_pl10000 >> signal_mh350_pl10000[n]) n++;
  mh350_pl10000.close();
  /*
  // Background
  double background[15];
  ifstream QCD;
  QCD.open("Efficiency_HtBins_Background.txt");
  n=0;
  while (QCD >> background[n]) n++;
  QCD.close();
  */

  double mh125_eff[3], mh250_eff[1], mh350_eff[3], mh1000_eff[3];
  mh125_eff[0] = signal_mh125_pl500[14];
  mh125_eff[1] = signal_mh125_pl1000[14];
  mh125_eff[2] = signal_mh125_pl10000[14];
  mh250_eff[0] = signal_mh250_pl1000[14];
  mh350_eff[0] = signal_mh350_pl500[14];
  mh350_eff[1] = signal_mh350_pl1000[14];
  mh350_eff[2] = signal_mh350_pl10000[14];
  mh1000_eff[0] = signal_mh1000_pl500[14];
  mh1000_eff[1] = signal_mh1000_pl1000[14];
  mh1000_eff[2] = signal_mh1000_pl10000[14];
  double ctau[3] = {500, 1000, 10000};
  double ctau_10[1] = {1000};

  // mh = 125 GeV
  TGraph *gr_LLP_mh125 = new TGraph(3, ctau, mh125_eff);
  TGraph *gr_LLP_mh250 = new TGraph(1, ctau_10, mh250_eff);
  TGraph *gr_LLP_mh350 = new TGraph(3, ctau, mh350_eff);
  TGraph *gr_LLP_mh1000 = new TGraph(3, ctau, mh1000_eff);

  TCanvas *c1_LLP_mh125 = new TCanvas("c1_LLP_mh125","Graph Draw Options",200,10,600,400);
  gr_LLP_mh125->SetMarkerColor(kRed);
  gr_LLP_mh125->SetMarkerStyle(21);
  gr_LLP_mh125->SetLineColor(kRed);
  gr_LLP_mh125->Draw("AC*");
  gr_LLP_mh125->SetTitle("Signal Efficiency vs. c#scale[1.2]{#tau} of LLP; c#scale[1.2]{#tau} of LLP (mm); Trigger Efficiency (events passed HT120+timing / all events)");
  gr_LLP_mh125->GetHistogram()->SetMinimum(0.);
  gr_LLP_mh125->GetHistogram()->SetMaximum(1.);
  c1_LLP_mh125->SetGrid();
  gr_LLP_mh250->SetMarkerColor(kGreen);
  gr_LLP_mh250->SetMarkerStyle(21);
  gr_LLP_mh250->SetLineColor(kGreen);
  gr_LLP_mh250->Draw("C*");
  gr_LLP_mh350->SetMarkerColor(kBlue);
  gr_LLP_mh350->SetMarkerStyle(21);
  gr_LLP_mh350->SetLineColor(kBlue);
  gr_LLP_mh350->Draw("C*");
  gr_LLP_mh1000->SetMarkerColor(kMagenta);
  gr_LLP_mh1000->SetMarkerStyle(21);
  gr_LLP_mh1000->SetLineColor(kMagenta);
  gr_LLP_mh1000->Draw("C*");
  auto legend_htSum = new TLegend(0.6,0.15,0.75,0.3);
  legend_htSum->AddEntry(gr_LLP_mh125,"m_{H}=125 GeV");
  legend_htSum->AddEntry(gr_LLP_mh250,"m_{H}=250 GeV");
  legend_htSum->AddEntry(gr_LLP_mh350,"m_{H}=350 GeV");
  legend_htSum->AddEntry(gr_LLP_mh1000,"m_{H}=1000 GeV");
  legend_htSum->Draw();
  c1_LLP_mh125->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/TimingBit/Eff_ctau.pdf");

}
