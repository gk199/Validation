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
  double signal_mh125_pl500[14];
  ifstream mh125_pl500;
  mh125_pl500.open("Efficiency_HtBins_Signal_mh125__mx50__pl500__.txt");
  int n=0;
  while (mh125_pl500 >> signal_mh125_pl500[n]) n++;
  mh125_pl500.close();
  // mh=125 pl=1m
  double signal_mh125_pl1000[14];
  ifstream mh125_pl1000;
  mh125_pl1000.open("Efficiency_HtBins_Signal_mh125__mx50__pl1000_.txt");
  n=0;
  while (mh125_pl1000 >> signal_mh125_pl1000[n]) n++;
  mh125_pl1000.close();
  // mh=125 pl=10m
  double signal_mh125_pl10000[14];
  ifstream mh125_pl10000;
  mh125_pl10000.open("Efficiency_HtBins_Signal_mh125__mx50__pl10000.txt");
  n=0;
  while (mh125_pl10000 >> signal_mh125_pl10000[n]) n++;
  mh125_pl10000.close();
  /*
  // mh=1000 pl=0.5m
  double signal_mh1000_pl500[14];
  ifstream mh1000_pl500;
  mh1000_pl500.open("Efficiency_HtBins_Signal_mh1000_mx450_pl500__.txt");
  n=0;
  while (mh1000_pl500 >> signal_mh1000_pl500[n]) n++;
  mh1000_pl500.close();
  // mh=1000 pl=1m 
  double signal_mh1000_pl1000[14];
  ifstream mh1000_pl1000;
  mh1000_pl1000.open("Efficiency_HtBins_Signal_mh1000_mx450_pl1000_.txt");
  n=0;
  while (mh1000_pl1000 >> signal_mh1000_pl1000[n]) n++;
  mh1000_pl1000.close();
  // mh=1000 pl=10m
  double signal_mh1000_pl10000[14];
  ifstream mh1000_pl10000;
  mh1000_pl10000.open("Efficiency_HtBins_Signal_mh1000_mx450_pl10000.txt");
  n=0;
  while (mh1000_pl10000 >> signal_mh1000_pl10000[n]) n++;
  mh1000_pl10000.close();

  // mh=250 pl=1m 
  double signal_mh250_pl1000[14];
  ifstream mh250_pl1000;
  mh250_pl1000.open("Efficiency_HtBins_Signal_mh250__mx120_pl1000_.txt");
  n=0;
  while (mh250_pl1000 >> signal_mh250_pl1000[n]) n++;
  mh250_pl1000.close();

  // mh=350 pl=0.5m 
  double signal_mh350_pl500[14];
  ifstream mh350_pl500;
  mh350_pl500.open("Efficiency_HtBins_Signal_mh350__mx160_pl500__.txt");
  n=0;
  while (mh350_pl500 >> signal_mh350_pl500[n]) n++;
  mh350_pl500.close();
  // mh=350 pl=1m
  double signal_mh350_pl1000[14];
  ifstream mh350_pl1000;
  mh350_pl1000.open("Efficiency_HtBins_Signal_mh350__mx160_pl1000_.txt");
  n=0;
  while (mh350_pl1000 >> signal_mh350_pl1000[n]) n++;
  mh350_pl1000.close();
  // mh=350 pl=10m
  double signal_mh350_pl10000[14];
  ifstream mh350_pl10000;
  mh350_pl10000.open("Efficiency_HtBins_Signal_mh350__mx160_pl10000.txt");
  n=0;
  while (mh350_pl10000 >> signal_mh350_pl10000[n]) n++;
  mh350_pl10000.close();
  // Background
  double background[14];
  ifstream QCD;
  QCD.open("Efficiency_HtBins_Background.txt");
  n=0;
  while (QCD >> background[n]) n++;
  QCD.close();
  */

  // ht360
  double L1HT360[14] = {0,0,0,0,0,0,0,0,0,0,0,0,1,1};
  // Middle of HT bin used for determining efficiencies
  double htBin[14] = {130, 150, 170, 190, 210, 230, 250, 270, 290, 310, 330, 350, 370, 390};
  double htBin_er[14] = {10,10,10,10,10,10,10,10,10,10,10,10,10,10};
  double err_y[14]={0};
  
  TGraphErrors *gr_L1HT360 = new TGraphErrors(14, htBin, L1HT360, htBin_er,err_y);
  // 1m
  // mh = 125 GeV
  TGraphErrors *gr_LLP_mh125_pl1000 = new TGraphErrors(14, htBin, signal_mh125_pl1000,htBin_er,err_y);
  // mh = 1000 GeV
  //  TGraphErrors *gr_LLP_mh1000_pl1000 = new TGraphErrors(14, htBin, signal_mh1000_pl1000,htBin_er,err_y);
  // mh = 250 GeV
  //  TGraphErrors *gr_LLP_mh250_pl1000 = new TGraphErrors(14, htBin, signal_mh250_pl1000,htBin_er,err_y);
  // mh = 350 GeV
  //  TGraphErrors *gr_LLP_mh350_pl1000 = new TGraphErrors(14, htBin, signal_mh350_pl1000,htBin_er,err_y);

  // 0.5m            
  // mh = 125 GeV  
  TGraphErrors *gr_LLP_mh125_pl500 = new TGraphErrors(14, htBin, signal_mh125_pl500,htBin_er,err_y);
  // mh = 1000 GeV 
  //  TGraphErrors *gr_LLP_mh1000_pl500 = new TGraphErrors(14, htBin, signal_mh1000_pl500,htBin_er,err_y);
  // mh = 350 GeV  
  //  TGraphErrors *gr_LLP_mh350_pl500 = new TGraphErrors(14, htBin, signal_mh350_pl500,htBin_er,err_y);

  // 10m            
  // mh = 125 GeV
  TGraphErrors *gr_LLP_mh125_pl10000 = new TGraphErrors(14, htBin, signal_mh125_pl10000,htBin_er,err_y);
  // mh = 1000 GeV 
  //  TGraphErrors *gr_LLP_mh1000_pl10000 = new TGraphErrors(14, htBin, signal_mh1000_pl10000,htBin_er,err_y);
  // mh = 350 GeV
  //  TGraphErrors *gr_LLP_mh350_pl10000 = new TGraphErrors(14, htBin, signal_mh350_pl10000,htBin_er,err_y);

  // QCD

  // mh = 125 GeV
  // 0.5 m
  TCanvas *c1_LLP_pl500 = new TCanvas("c1_LLP_pl500","Graph Draw Options",200,10,600,400);
  gr_LLP_mh125_pl500->SetMarkerColor(kBlue);
  gr_LLP_mh125_pl500->SetMarkerStyle(21);
  gr_LLP_mh125_pl500->Draw("AP");
  gr_L1HT360->SetMarkerColor(kRed);
  gr_L1HT360->SetMarkerStyle(21);
  //  gr_L1HT360->Draw("P");
  gr_LLP_mh125_pl500->SetTitle("HT vs. Signal Efficiency for >=2 Cells >=50ADC,3ns near 4 L1 Jets;HT (GeV);LLP Efficiency, mh=125 GeV, c#scale[1.2]{#tau}=0.5m   ");
  gr_LLP_mh125_pl500->GetHistogram()->SetMinimum(0.);
  gr_LLP_mh125_pl500->GetHistogram()->SetMaximum(1.);
  c1_LLP_pl500->SetGrid();
  auto legend_htSum = new TLegend(0.15,0.65,0.45,0.8);
  legend_htSum->AddEntry(gr_LLP_mh125_pl500,"HT>120 with timing cuts");
  //  legend_htSum->AddEntry(gr_L1HT360,"HT>360 at L1");
  legend_htSum->Draw();
  c1_LLP_pl500->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/TimingBit/EffHT_LLP_mh125_pl500.pdf");

  TCanvas *c1_LLP_pl1000 = new TCanvas("c1_LLP_pl1000","Graph Draw Options",200,10,600,400);
  gr_LLP_mh125_pl1000->SetMarkerColor(kBlue);
  gr_LLP_mh125_pl1000->SetMarkerStyle(21);
  gr_LLP_mh125_pl1000->Draw("AP");
  gr_L1HT360->SetMarkerColor(kRed);
  gr_L1HT360->SetMarkerStyle(21);
  //  gr_L1HT360->Draw("P");
  gr_LLP_mh125_pl1000->SetTitle("HT vs. Signal Efficiency for >=2 Cells >=50ADC,3ns near 4 L1 Jets;HT (GeV);LLP Efficiency, mh=125 GeV, c#scale[1.2]{#tau}=1m   ");
  gr_LLP_mh125_pl1000->GetHistogram()->SetMinimum(0.);
  gr_LLP_mh125_pl1000->GetHistogram()->SetMaximum(1.);
  c1_LLP_pl1000->SetGrid();
  auto legend2_htSum = new TLegend(0.15,0.65,0.45,0.8);
  legend2_htSum->AddEntry(gr_LLP_mh125_pl1000,"HT>120 with timing cuts");
  //  legend2_htSum->AddEntry(gr_L1HT360,"HT>360 at L1");
  legend2_htSum->Draw();
  c1_LLP_pl1000->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/TimingBit/EffHT_LLP_mh125_pl1000.pdf");

  TCanvas *c1_LLP_pl10000 = new TCanvas("c1_LLP_pl10000","Graph Draw Options",200,10,600,400);
  gr_LLP_mh125_pl10000->SetMarkerColor(kBlue);
  gr_LLP_mh125_pl10000->SetMarkerStyle(21);
  gr_LLP_mh125_pl10000->Draw("AP");
  gr_L1HT360->SetMarkerColor(kRed);
  gr_L1HT360->SetMarkerStyle(21);
  //  gr_L1HT360->Draw("P");
  gr_LLP_mh125_pl10000->SetTitle("HT vs. Signal Efficiency for >=2 Cells >=50ADC,3ns near 4 L1 Jets;HT (GeV);LLP Efficiency, mh=125 GeV, c#scale[1.2]{#tau}=10m   ");
  gr_LLP_mh125_pl10000->GetHistogram()->SetMinimum(0.);
  gr_LLP_mh125_pl10000->GetHistogram()->SetMaximum(1.);
  c1_LLP_pl10000->SetGrid();
  auto legend3_htSum = new TLegend(0.15,0.65,0.45,0.8);
  legend3_htSum->AddEntry(gr_LLP_mh125_pl10000,"HT>120 with timing cuts");
  //  legend3_htSum->AddEntry(gr_L1HT360,"HT>360 at L1");
  legend3_htSum->Draw();
  c1_LLP_pl10000->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/TimingBit/EffHT_LLP_mh125_pl10000.pdf");
}
