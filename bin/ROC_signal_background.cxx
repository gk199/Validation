#include<iostream>
#include<fstream>
#include<vector>
#include "TGraph.h"
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
  // read in signal and background files as arrays -- these are for the trigger based on prompt and delayed seeds in any L1 jet (HB, HE restrictions)
  // mh = 1000 GeV
  double arr_signal_mh1000_pl10000[20];
  ifstream Signal_mh1000_pl10000;
  Signal_mh1000_pl10000.open("MultiplicityHits50ADC3ns_ht120_Signal_mh1000_mx450_pl10000.txt");
  int n = 0;
  while (Signal_mh1000_pl10000 >> arr_signal_mh1000_pl10000[n]) n++; // signal efficiency
  Signal_mh1000_pl10000.close();
  double arr_signal_mh1000_pl1000[20];
  ifstream Signal_mh1000_pl1000;
  Signal_mh1000_pl1000.open("MultiplicityHits50ADC3ns_ht120_Signal_mh1000_mx450_pl1000_.txt");
  n = 0;
  while (Signal_mh1000_pl1000 >> arr_signal_mh1000_pl1000[n]) n++; // signal efficiency
  Signal_mh1000_pl1000.close();
  double arr_signal_mh1000_pl500[20];
  ifstream Signal_mh1000_pl500;
  Signal_mh1000_pl500.open("MultiplicityHits50ADC3ns_ht120_Signal_mh1000_mx450_pl500__.txt");
  n = 0;
  while (Signal_mh1000_pl500 >> arr_signal_mh1000_pl500[n]) n++; // signal efficiency  
  Signal_mh1000_pl500.close();
  // mh = 350 GeV
  double arr_signal_mh350_pl10000[20];
  ifstream Signal_mh350_pl10000;
  Signal_mh350_pl10000.open("MultiplicityHits50ADC3ns_ht120_Signal_mh350__mx160_pl10000.txt");
  n = 0;
  while (Signal_mh350_pl10000 >> arr_signal_mh350_pl10000[n]) n++; // signal efficiency
  Signal_mh350_pl10000.close();
  double arr_signal_mh350_pl1000[20];
  ifstream Signal_mh350_pl1000;
  Signal_mh350_pl1000.open("MultiplicityHits50ADC3ns_ht120_Signal_mh350__mx160_pl1000_.txt");
  n = 0;
  while (Signal_mh350_pl1000 >> arr_signal_mh350_pl1000[n]) n++; // signal efficiency
  Signal_mh350_pl1000.close();
  double arr_signal_mh350_pl500[20];
  ifstream Signal_mh350_pl500;
  Signal_mh350_pl500.open("MultiplicityHits50ADC3ns_ht120_Signal_mh350__mx160_pl500__.txt");
  n = 0;
  while (Signal_mh350_pl500 >> arr_signal_mh350_pl500[n]) n++; // signal efficiency
  Signal_mh350_pl500.close();
  // mh = 250 GeV     
  double arr_signal_mh250_pl1000[20];
  ifstream Signal_mh250_pl1000;
  Signal_mh250_pl1000.open("MultiplicityHits50ADC3ns_ht120_Signal_mh250__mx120_pl1000_.txt");
  n = 0;
  while (Signal_mh250_pl1000 >> arr_signal_mh250_pl1000[n]) n++; // signal efficiency  
  Signal_mh250_pl1000.close();
  // mh = 125 GeV            
  double arr_signal_mh125_pl10000[20];
  ifstream Signal_mh125_pl10000;
  Signal_mh125_pl10000.open("MultiplicityHits50ADC3ns_ht120_Signal_mh125__mx50__pl10000.txt");
  n = 0;
  while (Signal_mh125_pl10000 >> arr_signal_mh125_pl10000[n]) n++; // signal efficiency
  Signal_mh125_pl10000.close();
  double arr_signal_mh125_pl1000[20];
  ifstream Signal_mh125_pl1000;
  Signal_mh125_pl1000.open("MultiplicityHits50ADC3ns_ht120_Signal_mh125__mx50__pl1000_.txt");
  n = 0;
  while (Signal_mh125_pl1000 >> arr_signal_mh125_pl1000[n]) n++; // signal efficiency  
  Signal_mh125_pl1000.close();
  double arr_signal_mh125_pl500[20];
  ifstream Signal_mh125_pl500;
  Signal_mh125_pl500.open("MultiplicityHits50ADC3ns_ht120_Signal_mh125__mx50__pl500__.txt");
  n = 0;
  while (Signal_mh125_pl500 >> arr_signal_mh125_pl500[n]) n++; // signal efficiency    
  Signal_mh125_pl500.close();

  double arr_background[20];
  ifstream Background;
  Background.open("MultiplicityHits50ADC3ns_ht120_Background.txt");
  n = 0;
  while (Background >> arr_background[n]) n++; // background efficiency
  Background.close();

  // split up these arrays in chunks based on > 120 + timing, >120 + timing OR > 360
  double Signal_120timing_mh1000_pl10000[5],Signal_120timingOR360_mh1000_pl10000[5];
  double Signal_120timing_mh1000_pl1000[5],Signal_120timingOR360_mh1000_pl1000[5];
  double Signal_120timing_mh1000_pl500[5],Signal_120timingOR360_mh1000_pl500[5];
  double Signal_120timing_mh350_pl10000[5],Signal_120timingOR360_mh350_pl10000[5];
  double Signal_120timing_mh350_pl1000[5],Signal_120timingOR360_mh350_pl1000[5];
  double Signal_120timing_mh350_pl500[5],Signal_120timingOR360_mh350_pl500[5];
  double Signal_120timing_mh250_pl1000[5],Signal_120timingOR360_mh250_pl1000[5];
  double Signal_120timing_mh125_pl10000[5],Signal_120timingOR360_mh125_pl10000[5];
  double Signal_120timing_mh125_pl1000[5],Signal_120timingOR360_mh125_pl1000[5];
  double Signal_120timing_mh125_pl500[5],Signal_120timingOR360_mh125_pl500[5];
  double Background_120timing[5],Background_120timingOR360[5];

  for (int i=0; i<5; i++) {
    Signal_120timingOR360_mh1000_pl10000[i] = arr_signal_mh1000_pl10000[i]; // eff at 120+timing OR 360
    Signal_120timing_mh1000_pl10000[i] = arr_signal_mh1000_pl10000[i+11]; // eff at 120+timing, this is last grouping of what is written to txt file
    Signal_120timingOR360_mh1000_pl1000[i] = arr_signal_mh1000_pl1000[i];
    Signal_120timing_mh1000_pl1000[i] = arr_signal_mh1000_pl1000[i+11];
    Signal_120timingOR360_mh1000_pl500[i] = arr_signal_mh1000_pl500[i];
    Signal_120timing_mh1000_pl500[i] = arr_signal_mh1000_pl500[i+11];
    Signal_120timingOR360_mh350_pl10000[i] = arr_signal_mh350_pl10000[i];
    Signal_120timing_mh350_pl10000[i] = arr_signal_mh350_pl10000[i+11];
    Signal_120timingOR360_mh350_pl1000[i] = arr_signal_mh350_pl1000[i];
    Signal_120timing_mh350_pl1000[i] = arr_signal_mh350_pl1000[i+11];
    Signal_120timingOR360_mh350_pl500[i] = arr_signal_mh350_pl500[i];
    Signal_120timing_mh350_pl500[i] = arr_signal_mh350_pl500[i+11];
    Signal_120timingOR360_mh250_pl1000[i] = arr_signal_mh250_pl1000[i];
    Signal_120timing_mh250_pl1000[i] = arr_signal_mh250_pl1000[i+11];
    Signal_120timingOR360_mh125_pl10000[i] = arr_signal_mh125_pl10000[i];
    Signal_120timing_mh125_pl10000[i] = arr_signal_mh125_pl10000[i+11];
    Signal_120timingOR360_mh125_pl1000[i] = arr_signal_mh125_pl1000[i];
    Signal_120timing_mh125_pl1000[i] = arr_signal_mh125_pl1000[i+11];
    Signal_120timingOR360_mh125_pl500[i] = arr_signal_mh125_pl500[i];
    Signal_120timing_mh125_pl500[i] = arr_signal_mh125_pl500[i+11];

    Background_120timingOR360[i] = 1-arr_background[i]; // background rejection
    Background_120timing[i] = 1-arr_background[i+11]; // background rejection                                                                                                                        
  }
  for (int i=0; i<5; i++) {
    std::cout << Background_120timingOR360[i] << " = background OR "
	      << Background_120timing[i] << " = background 120 " 
	      << Signal_120timingOR360_mh1000_pl10000[i] << " = signal OR "
	      << Signal_120timing_mh1000_pl10000[i] << " = signal 120 " << std::endl;
  }



  // make TGraphs with various TDC scans overlayed to compare signal acceptance and background rejection performance
  TGraph *gr_120timing_mh1000_pl10000_timescan = new TGraph (5, Signal_120timing_mh1000_pl10000,Background_120timing);
  TGraph *gr_120timingOR360_mh1000_pl10000_timescan = new TGraph (5, Signal_120timingOR360_mh1000_pl10000,Background_120timingOR360);
  TGraph *gr_120timing_mh1000_pl1000_timescan = new TGraph (5, Signal_120timing_mh1000_pl1000,Background_120timing);
  TGraph *gr_120timingOR360_mh1000_pl1000_timescan = new TGraph (5, Signal_120timingOR360_mh1000_pl1000,Background_120timingOR360);
  TGraph *gr_120timing_mh1000_pl500_timescan = new TGraph (5, Signal_120timing_mh1000_pl500,Background_120timing);
  TGraph *gr_120timingOR360_mh1000_pl500_timescan = new TGraph (5, Signal_120timingOR360_mh1000_pl500,Background_120timingOR360);
  TGraph *gr_120timing_mh350_pl10000_timescan = new TGraph (5, Signal_120timing_mh350_pl10000,Background_120timing);
  TGraph *gr_120timingOR360_mh350_pl10000_timescan = new TGraph (5, Signal_120timingOR360_mh350_pl10000,Background_120timingOR360);
  TGraph *gr_120timing_mh350_pl1000_timescan = new TGraph (5, Signal_120timing_mh350_pl1000,Background_120timing);
  TGraph *gr_120timingOR360_mh350_pl1000_timescan = new TGraph (5, Signal_120timingOR360_mh350_pl1000,Background_120timingOR360);
  TGraph *gr_120timing_mh350_pl500_timescan = new TGraph (5, Signal_120timing_mh350_pl500,Background_120timing);
  TGraph *gr_120timingOR360_mh350_pl500_timescan = new TGraph (5, Signal_120timingOR360_mh350_pl500,Background_120timingOR360);
  TGraph *gr_120timing_mh250_pl1000_timescan = new TGraph (5, Signal_120timing_mh250_pl1000,Background_120timing);
  TGraph *gr_120timingOR360_mh250_pl1000_timescan = new TGraph (5, Signal_120timingOR360_mh250_pl1000,Background_120timingOR360);
  TGraph *gr_120timing_mh125_pl10000_timescan = new TGraph (5, Signal_120timing_mh125_pl10000,Background_120timing);
  TGraph *gr_120timingOR360_mh125_pl10000_timescan = new TGraph (5, Signal_120timingOR360_mh125_pl10000,Background_120timingOR360);
  TGraph *gr_120timing_mh125_pl1000_timescan = new TGraph (5, Signal_120timing_mh125_pl1000,Background_120timing);
  TGraph *gr_120timingOR360_mh125_pl1000_timescan = new TGraph (5, Signal_120timingOR360_mh125_pl1000,Background_120timingOR360);
  TGraph *gr_120timing_mh125_pl500_timescan = new TGraph (5, Signal_120timing_mh125_pl500,Background_120timing);
  TGraph *gr_120timingOR360_mh125_pl500_timescan = new TGraph (5, Signal_120timingOR360_mh125_pl500,Background_120timingOR360);

  TCanvas *ROC_120timingOR360_pl10000 = new TCanvas("ROC_120timingOR360_pl10000","Graph Draw Options",200,10,600,600);
  ROC_120timingOR360_pl10000->SetGrid();
  gr_120timingOR360_mh1000_pl10000_timescan->SetLineColor(4); // blue 
  gr_120timingOR360_mh1000_pl10000_timescan->Draw("AC*");
  gr_120timingOR360_mh1000_pl10000_timescan->GetXaxis()->SetLimits(0.,1.1);
  gr_120timingOR360_mh1000_pl10000_timescan->GetHistogram()->SetMaximum(1.1);
  gr_120timingOR360_mh1000_pl10000_timescan->GetHistogram()->SetMinimum(0.);
  gr_120timingOR360_mh1000_pl10000_timescan->SetTitle("Signal Eff. vs Background Rejection with Prompt Seed Approach;Signal Eff (LLP, mh=1TeV, ctau=10m);QCD Rejection");
  gr_120timingOR360_mh350_pl10000_timescan->SetLineColor(2);
  gr_120timingOR360_mh350_pl10000_timescan->Draw("C*");
  gr_120timingOR360_mh125_pl10000_timescan->SetLineColor(3);
  gr_120timingOR360_mh125_pl10000_timescan->Draw("C*");
  auto legend_pl10000 = new TLegend(0.45,0.65,0.9,0.8);
  legend_pl10000->AddEntry(gr_120timingOR360_mh125_pl10000_timescan,"mh = 125, 120+timing OR 360, delayed objects scanned");
  legend_pl10000->AddEntry(gr_120timingOR360_mh350_pl10000_timescan,"mh = 350, 120+timing OR 360, delayed objects scanned");
  legend_pl10000->AddEntry(gr_120timingOR360_mh1000_pl10000_timescan,"mh = 1000, 120+timing OR 360, delayed objects scanned");
  legend_pl10000->Draw();
  ROC_120timingOR360_pl10000->SaveAs("plots/ROC_120timingOR360_pl10000.pdf");

  TCanvas *ROC_120timingOR360_pl1000 = new TCanvas("ROC_120timingOR360_pl1000","Graph Draw Options",200,10,600,600);
  ROC_120timingOR360_pl1000->SetGrid();
  gr_120timingOR360_mh1000_pl1000_timescan->SetLineColor(4); // blue        
  gr_120timingOR360_mh1000_pl1000_timescan->Draw("AC*");
  gr_120timingOR360_mh1000_pl1000_timescan->GetXaxis()->SetLimits(0.,1.1);
  gr_120timingOR360_mh1000_pl1000_timescan->GetHistogram()->SetMaximum(1.1);
  gr_120timingOR360_mh1000_pl1000_timescan->GetHistogram()->SetMinimum(0.);
  gr_120timingOR360_mh1000_pl1000_timescan->SetTitle("Signal Eff. vs Background Rejection with Prompt Seed Approach;Signal Eff (LLP, ctau=1m);QCD Rejection");
  gr_120timingOR360_mh350_pl1000_timescan->SetLineColor(2);
  gr_120timingOR360_mh350_pl1000_timescan->Draw("C*");
  gr_120timingOR360_mh125_pl1000_timescan->SetLineColor(3);
  gr_120timingOR360_mh125_pl1000_timescan->Draw("C*");
  gr_120timingOR360_mh250_pl1000_timescan->SetLineColor(1);
  gr_120timingOR360_mh250_pl1000_timescan->Draw("C*");
  auto legend_pl1000 = new TLegend(0.45,0.65,0.9,0.8);
  legend_pl1000->AddEntry(gr_120timingOR360_mh125_pl1000_timescan,"mh = 125, 120+timing OR 360, delayed objects scanned");
  legend_pl1000->AddEntry(gr_120timingOR360_mh250_pl1000_timescan,"mh = 250, 120+timing OR 360, delayed objects scanned");
  legend_pl1000->AddEntry(gr_120timingOR360_mh350_pl1000_timescan,"mh = 350, 120+timing OR 360, delayed objects scanned");
  legend_pl1000->AddEntry(gr_120timingOR360_mh1000_pl1000_timescan,"mh = 1000, 120+timing OR 360, delayed objects scanned");
  legend_pl1000->Draw();
  ROC_120timingOR360_pl1000->SaveAs("plots/ROC_120timingOR360_pl1000.pdf");

  TCanvas *ROC_120timingOR360_pl500 = new TCanvas("ROC_120timingOR360_pl500","Graph Draw Options",200,10,600,600);
  ROC_120timingOR360_pl500->SetGrid();
  gr_120timingOR360_mh1000_pl500_timescan->SetLineColor(4); // blue                                                                                                                                  
  gr_120timingOR360_mh1000_pl500_timescan->Draw("AC*");
  gr_120timingOR360_mh1000_pl500_timescan->GetXaxis()->SetLimits(0.,1.1);
  gr_120timingOR360_mh1000_pl500_timescan->GetHistogram()->SetMaximum(1.1);
  gr_120timingOR360_mh1000_pl500_timescan->GetHistogram()->SetMinimum(0.);
  gr_120timingOR360_mh1000_pl500_timescan->SetTitle("Signal Eff. vs Background Rejection with Prompt Seed Approach;Signal Eff (LLP, ctau=0.5m);QCD Rejection");
  gr_120timingOR360_mh350_pl500_timescan->SetLineColor(2);
  gr_120timingOR360_mh350_pl500_timescan->Draw("C*");
  gr_120timingOR360_mh125_pl500_timescan->SetLineColor(3);
  gr_120timingOR360_mh125_pl500_timescan->Draw("C*");
  auto legend_pl500 = new TLegend(0.45,0.65,0.9,0.8);
  legend_pl500->AddEntry(gr_120timingOR360_mh125_pl500_timescan,"mh = 125, 120+timing OR 360, delayed objects scanned");
  legend_pl500->AddEntry(gr_120timingOR360_mh350_pl500_timescan,"mh = 350, 120+timing OR 360, delayed objects scanned");
  legend_pl500->AddEntry(gr_120timingOR360_mh1000_pl500_timescan,"mh = 1000, 120+timing OR 360, delayed objects scanned");
  legend_pl500->Draw();
  ROC_120timingOR360_pl500->SaveAs("plots/ROC_120timingOR360_pl500.pdf");

  // 120 + timing
  TCanvas *ROC_120timing_pl10000 = new TCanvas("ROC_120timing_pl10000","Graph Draw Options",200,10,600,600);
  ROC_120timing_pl10000->SetGrid();
  gr_120timing_mh1000_pl10000_timescan->SetLineColor(4); // blue                                                                                                                                  
  gr_120timing_mh1000_pl10000_timescan->Draw("AC*");
  gr_120timing_mh1000_pl10000_timescan->GetXaxis()->SetLimits(0.,1.1);
  gr_120timing_mh1000_pl10000_timescan->GetHistogram()->SetMaximum(1.1);
  gr_120timing_mh1000_pl10000_timescan->GetHistogram()->SetMinimum(0.);
  gr_120timing_mh1000_pl10000_timescan->SetTitle("Signal Eff. vs Background Rejection with Prompt Seed Approach;Signal Eff (LLP, mh=1TeV, ctau=10m);QCD Rejection");
  gr_120timing_mh350_pl10000_timescan->SetLineColor(2);
  gr_120timing_mh350_pl10000_timescan->Draw("C*");
  gr_120timing_mh125_pl10000_timescan->SetLineColor(3);
  gr_120timing_mh125_pl10000_timescan->Draw("C*");
  auto legend_120timing_pl10000 = new TLegend(0.15,0.15,0.55,0.4);
  legend_120timing_pl10000->AddEntry(gr_120timing_mh125_pl10000_timescan,"mh = 125, 120+timing, delayed objects scanned");
  legend_120timing_pl10000->AddEntry(gr_120timing_mh350_pl10000_timescan,"mh = 350, 120+timing, delayed objects scanned");
  legend_120timing_pl10000->AddEntry(gr_120timing_mh1000_pl10000_timescan,"mh = 1000, 120+timing, delayed objects scanned");
  legend_120timing_pl10000->Draw();
  ROC_120timing_pl10000->SaveAs("plots/ROC_120timing_pl10000.pdf");

  TCanvas *ROC_120timing_pl1000 = new TCanvas("ROC_120timing_pl1000","Graph Draw Options",200,10,600,600);
  ROC_120timing_pl1000->SetGrid();
  gr_120timing_mh1000_pl1000_timescan->SetLineColor(4); // blue                                                                                                                                   
  gr_120timing_mh1000_pl1000_timescan->Draw("AC*");
  gr_120timing_mh1000_pl1000_timescan->GetXaxis()->SetLimits(0.,1.1);
  gr_120timing_mh1000_pl1000_timescan->GetHistogram()->SetMaximum(1.1);
  gr_120timing_mh1000_pl1000_timescan->GetHistogram()->SetMinimum(0.);
  gr_120timing_mh1000_pl1000_timescan->SetTitle("Signal Eff. vs Background Rejection with Prompt Seed Approach;Signal Eff (LLP, ctau=1m);QCD Rejection");
  gr_120timing_mh350_pl1000_timescan->SetLineColor(2);
  gr_120timing_mh350_pl1000_timescan->Draw("C*");
  gr_120timing_mh125_pl1000_timescan->SetLineColor(3);
  gr_120timing_mh125_pl1000_timescan->Draw("C*");
  gr_120timing_mh250_pl1000_timescan->SetLineColor(1);
  gr_120timing_mh250_pl1000_timescan->Draw("C*");
  auto legend_120timing_pl1000 = new TLegend(0.15,0.15,0.55,0.4);
  legend_120timing_pl1000->AddEntry(gr_120timing_mh125_pl1000_timescan,"mh = 125, 120+timing, delayed objects scanned");
  legend_120timing_pl1000->AddEntry(gr_120timing_mh250_pl1000_timescan,"mh = 250, 120+timing, delayed objects scanned");
  legend_120timing_pl1000->AddEntry(gr_120timing_mh350_pl1000_timescan,"mh = 350, 120+timing, delayed objects scanned");
  legend_120timing_pl1000->AddEntry(gr_120timing_mh1000_pl1000_timescan,"mh = 1000, 120+timing, delayed objects scanned");
  legend_120timing_pl1000->Draw();
  ROC_120timing_pl1000->SaveAs("plots/ROC_120timing_pl1000.pdf");

  TCanvas *ROC_120timing_pl500 = new TCanvas("ROC_120timing_pl500","Graph Draw Options",200,10,600,600);
  ROC_120timing_pl500->SetGrid();
  gr_120timing_mh1000_pl500_timescan->SetLineColor(4); // blue                                                                                                                                    
  gr_120timing_mh1000_pl500_timescan->Draw("AC*");
  gr_120timing_mh1000_pl500_timescan->GetXaxis()->SetLimits(0.,1.1);
  gr_120timing_mh1000_pl500_timescan->GetHistogram()->SetMaximum(1.1);
  gr_120timing_mh1000_pl500_timescan->GetHistogram()->SetMinimum(0.);
  gr_120timing_mh1000_pl500_timescan->SetTitle("Signal Eff. vs Background Rejection with Prompt Seed Approach;Signal Eff (LLP, ctau=0.5m);QCD Rejection");
  gr_120timing_mh350_pl500_timescan->SetLineColor(2);
  gr_120timing_mh350_pl500_timescan->Draw("C*");
  gr_120timing_mh125_pl500_timescan->SetLineColor(3);
  gr_120timing_mh125_pl500_timescan->Draw("C*");
  auto legend_120timing_pl500 = new TLegend(0.15,0.15,0.55,0.4);
  legend_120timing_pl500->AddEntry(gr_120timing_mh125_pl500_timescan,"mh = 125, 120+timing, delayed objects scanned");
  legend_120timing_pl500->AddEntry(gr_120timing_mh350_pl500_timescan,"mh = 350, 120+timing, delayed objects scanned");
  legend_120timing_pl500->AddEntry(gr_120timing_mh1000_pl500_timescan,"mh = 1000, 120+timing, delayed objects scanned");
  legend_120timing_pl500->Draw();
  ROC_120timing_pl500->SaveAs("plots/ROC_120timing_pl500.pdf");

  return 0;
}

