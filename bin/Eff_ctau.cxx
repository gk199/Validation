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
  // mh=125 pl=3m
  double signal_mh125_pl3000[15];
  ifstream mh125_pl3000;
  mh125_pl3000.open("Efficiency_HtBins_Signal_mh125__pl3000__.txt");
  int n=0;
  while (mh125_pl3000 >> signal_mh125_pl3000[n]) n++;
  mh125_pl3000.close();
  // mh=125 pl=30m
  double signal_mh125_pl30000[15];
  ifstream mh125_pl30000;
  mh125_pl3000.open("Efficiency_HtBins_Signal_mh125__pl30000_.txt");
  n=0;
  while (mh125_pl30000 >> signal_mh125_pl30000[n]) n++;
  mh125_pl30000.close();

  // mh=250 pl=0.5m
  double signal_mh250_pl500[15];
  ifstream mh250_pl500;
  mh250_pl500.open("Efficiency_HtBins_Signal_mh250__pl500___.txt");
  n=0;
  while (mh250_pl500 >> signal_mh250_pl500[n]) n++;
  mh250_pl500.close();
  // mh=250 pl=1m
  double signal_mh250_pl1000[15];
  ifstream mh250_pl1000;
  mh250_pl1000.open("Efficiency_HtBins_Signal_mh250__pl1000__.txt");
  n=0;
  while (mh250_pl1000 >> signal_mh250_pl1000[n]) n++;
  mh250_pl1000.close();
  // mh=250 pl=10m
  double signal_mh250_pl10000[15];
  ifstream mh250_pl10000;
  mh250_pl10000.open("Efficiency_HtBins_Signal_mh250__pl10000_.txt");
  n=0;
  while (mh250_pl10000 >> signal_mh250_pl10000[n]) n++;
  mh250_pl10000.close();

  // mh=350 pl=0.5m
  double signal_mh350_pl500[15];
  ifstream mh350_pl500;
  mh350_pl500.open("Efficiency_HtBins_Signal_mh350__pl500___.txt");
  n=0;
  while (mh350_pl500 >> signal_mh350_pl500[n]) n++;
  mh350_pl500.close();
  // mh=350 pl=1m 
  double signal_mh350_pl1000[15];
  ifstream mh350_pl1000;
  mh350_pl1000.open("Efficiency_HtBins_Signal_mh350__pl1000__.txt");
  n=0;
  while (mh350_pl1000 >> signal_mh350_pl1000[n]) n++;
  mh350_pl1000.close();
  // mh=350 pl=10m      
  double signal_mh350_pl10000[15];
  ifstream mh350_pl10000;
  mh350_pl10000.open("Efficiency_HtBins_Signal_mh350__pl10000_.txt");
  n=0;
  while (mh350_pl10000 >> signal_mh350_pl10000[n]) n++;
  mh350_pl10000.close();

  // mh=1000 pl=10m
  double signal_mh1000_pl10000[15];
  ifstream mh1000_pl10000;
  mh1000_pl10000.open("Efficiency_HtBins_Signal_mh1000_pl10000_.txt");
  n=0;
  while (mh1000_pl10000 >> signal_mh1000_pl10000[n]) n++;
  mh1000_pl10000.close();
  // mh=1000 pl=100m
  double signal_mh1000_pl100000[15];
  ifstream mh1000_pl100000;
  mh1000_pl100000.open("Efficiency_HtBins_Signal_mh1000_pl100000.txt");
  n=0;
  while (mh1000_pl100000 >> signal_mh1000_pl100000[n]) n++;
  mh1000_pl100000.close();

  double mh125_eff[2], mh250_eff[3], mh350_eff[3], mh1000_eff[2];
  mh125_eff[0] = signal_mh125_pl3000[14];
  mh125_eff[1] = signal_mh125_pl30000[14];
  mh250_eff[0] = signal_mh250_pl500[14];
  mh250_eff[1] = signal_mh250_pl1000[14];
  mh250_eff[2] = signal_mh250_pl10000[14];
  mh350_eff[0] = signal_mh350_pl500[14];
  mh350_eff[1] = signal_mh350_pl1000[14];
  mh350_eff[2] = signal_mh350_pl10000[14];
  mh1000_eff[0] = signal_mh1000_pl10000[14];
  mh1000_eff[1] = signal_mh1000_pl100000[14];

  double ctau_05_1_10[3] = {500, 1000, 10000};
  double ctau_3_30[2] = {3000, 30000};
  double ctau_10_100[2] = {10000, 100000};

  // mh = 125 GeV
  TGraph *gr_LLP_mh125 = new TGraph(2, ctau_3_30, mh125_eff);
  TGraph *gr_LLP_mh250 = new TGraph(3, ctau_05_1_10, mh250_eff);
  TGraph *gr_LLP_mh350 = new TGraph(3, ctau_05_1_10, mh350_eff);
  TGraph *gr_LLP_mh1000 = new TGraph(2, ctau_10_100, mh1000_eff);

  TCanvas *c1_LLP_mh125 = new TCanvas("c1_LLP_mh125","Graph Draw Options",200,10,600,400);
  gr_LLP_mh250->SetMarkerColor(kGreen);
  gr_LLP_mh250->SetMarkerStyle(21);
  gr_LLP_mh250->SetLineColor(kGreen);
  gr_LLP_mh250->Draw("AC*");
  gr_LLP_mh250->SetTitle("Signal Efficiency vs. c#scale[1.2]{#tau} of LLP; c#scale[1.2]{#tau} of LLP (mm); Trigger Efficiency (events passed HT120+timing / all events)");
  c1_LLP_mh125->SetLogx();
  gr_LLP_mh250->GetHistogram()->SetMinimum(0.);
  gr_LLP_mh250->GetHistogram()->SetMaximum(1.);
  gr_LLP_mh250->GetXaxis()->SetRangeUser(1.0,1002000.0);
  gr_LLP_mh250->Draw("AC*");
  c1_LLP_mh125->Update();
  c1_LLP_mh125->SetGrid();
  gr_LLP_mh125->SetMarkerColor(kRed);
  gr_LLP_mh125->SetMarkerStyle(21);
  gr_LLP_mh125->SetLineColor(kRed);
  gr_LLP_mh125->Draw("C*");
  gr_LLP_mh350->SetMarkerColor(kBlue);
  gr_LLP_mh350->SetMarkerStyle(21);
  gr_LLP_mh350->SetLineColor(kBlue);
  gr_LLP_mh350->Draw("C*");
  gr_LLP_mh1000->SetMarkerColor(kMagenta);
  gr_LLP_mh1000->SetMarkerStyle(21);
  gr_LLP_mh1000->SetLineColor(kMagenta);
  gr_LLP_mh1000->Draw("C*");
  gr_LLP_mh1000->GetXaxis()->SetRangeUser(1.0,1002000.0);
  c1_LLP_mh125->Update();
  auto legend_htSum = new TLegend(0.15,0.75,0.3,0.85);
  legend_htSum->AddEntry(gr_LLP_mh125,"m_{H}=125 GeV");
  legend_htSum->AddEntry(gr_LLP_mh250,"m_{H}=250 GeV");
  legend_htSum->AddEntry(gr_LLP_mh350,"m_{H}=350 GeV");
  legend_htSum->AddEntry(gr_LLP_mh1000,"m_{H}=1000 GeV");
  legend_htSum->Draw();
  gr_LLP_mh1000->GetXaxis()->SetRangeUser(1.0,1002000.0);
  //  c1_LLP_mh125->SetLogx();
  c1_LLP_mh125->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/Eff_ctau.pdf");

}
