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
  // read in signal and background files as arrays -- these are for the trigger based on delayed hit fraction within the first L1 jet
  double arr_signal_mh1000_pl1000[21];
  ifstream Signal_mh1000_pl1000;
  Signal_mh1000_pl1000.open("DelayedHitFrac_Signal_mh1000_pl10.txt");
  int n = 0;
  while (Signal_mh1000_pl1000 >> arr_signal_mh1000_pl1000[n]) n++; // signal efficiency
  Signal_mh1000_pl1000.close();
  double arr_signal_mh1000_pl500[21];
  ifstream Signal_mh1000_pl500;
  Signal_mh1000_pl500.open("DelayedHitFrac_Signal_mh1000_pl50.txt");
  n = 0;
  while (Signal_mh1000_pl500 >> arr_signal_mh1000_pl500[n]) n++; // signal efficiency  
  Signal_mh1000_pl500.close();

  double arr_background[21];
  ifstream Background;
  Background.open("DelayedHitFrac_Background.txt");
  n = 0;
  while (Background >> arr_background[n]) n++; // background efficiency
  Background.close();

  // read in signal and background files as arrays -- these are for the trigger based on quad jet hit multiplicity at 3 GeV, 3ns
  double arr_signal_mult_mh1000_pl1000[7];
  ifstream Signal_mult_mh1000_pl1000;
  Signal_mult_mh1000_pl1000.open("MultiplicityHits3GeV3ns_Signal_mh1000_pl10.txt");
  n = 0;
  while (Signal_mult_mh1000_pl1000 >> arr_signal_mult_mh1000_pl1000[n]) n++; // signal efficiency
  Signal_mult_mh1000_pl1000.close();
  double arr_signal_mult_mh1000_pl500[7];
  ifstream Signal_mult_mh1000_pl500;
  Signal_mult_mh1000_pl500.open("MultiplicityHits3GeV3ns_Signal_mh1000_pl50.txt");
  n = 0;
  while (Signal_mult_mh1000_pl500 >> arr_signal_mult_mh1000_pl500[n]) n++; // signal efficiency
  Signal_mult_mh1000_pl500.close();

  double arr_background_mult[7];
  ifstream Background_mult;
  Background_mult.open("MultiplicityHits3GeV3ns_Background.txt");
  n = 0;
  while (Background_mult >>arr_background_mult[n]) n++; // background efficiency
  Background_mult.close();

  // split up these arrays in chunks based on what energy the TDC scan was taken at
  double Signal_2GeV_mh1000_pl1000[7],Signal_3GeV_mh1000_pl1000[7], Signal_4GeV_mh1000_pl1000[7];
  double Signal_2GeV_mh1000_pl500[7],Signal_3GeV_mh1000_pl500[7], Signal_4GeV_mh1000_pl500[7];
  double Background_2GeV[7], Background_3GeV[7], Background_4GeV[7];
  double Background_3GeV_mult[7];
  for (int i=0; i<7; i++) {
    Signal_2GeV_mh1000_pl1000[i] = arr_signal_mh1000_pl1000[i]; // signal efficiency
    Signal_3GeV_mh1000_pl1000[i] = arr_signal_mh1000_pl1000[i+7];
    Signal_4GeV_mh1000_pl1000[i] = arr_signal_mh1000_pl1000[i+14];
    Signal_2GeV_mh1000_pl500[i] = arr_signal_mh1000_pl500[i]; // signal efficiency
    Signal_3GeV_mh1000_pl500[i] = arr_signal_mh1000_pl500[i+7];
    Signal_4GeV_mh1000_pl500[i] = arr_signal_mh1000_pl500[i+14];
    Background_2GeV[i] = 1-arr_background[i]; // background rejection
    Background_3GeV[i] = 1-arr_background[i+7];
    Background_4GeV[i] = 1-arr_background[i+14];
    Background_3GeV_mult[i] = 1-arr_background_mult[i];
  }

  // make TGraphs with various TDC scans overlayed to compare signal acceptance and background rejection performance
  TGraph *gr_2GeV_mh1000_pl1000_timescan = new TGraph (7, Signal_2GeV_mh1000_pl1000,Background_2GeV);
  TGraph *gr_3GeV_mh1000_pl1000_timescan = new TGraph (7, Signal_3GeV_mh1000_pl1000,Background_3GeV);
  TGraph *gr_4GeV_mh1000_pl1000_timescan = new TGraph (7, Signal_4GeV_mh1000_pl1000,Background_4GeV);
  TGraph *gr_2GeV_mh1000_pl500_timescan = new TGraph (7, Signal_2GeV_mh1000_pl500,Background_2GeV);
  TGraph *gr_3GeV_mh1000_pl500_timescan = new TGraph (7, Signal_3GeV_mh1000_pl500,Background_3GeV);
  TGraph *gr_4GeV_mh1000_pl500_timescan = new TGraph (7, Signal_4GeV_mh1000_pl500,Background_4GeV);
  TGraph *gr_3GeV_timescan_mh1000_pl1000_mult = new TGraph (7, arr_signal_mult_mh1000_pl1000, Background_3GeV_mult);
  TGraph *gr_3GeV_timescan_mh1000_pl500_mult = new TGraph (7, arr_signal_mult_mh1000_pl500, Background_3GeV_mult);

  TCanvas *ROC_DelayedHitFrac_mh1000_pl1000 = new TCanvas("ROC_DelayedHitFrac_mh1000_pl1000","Graph Draw Options",200,10,600,600);
  ROC_DelayedHitFrac_mh1000_pl1000->SetGrid();
  gr_2GeV_mh1000_pl1000_timescan->SetLineColor(4); // blue
  gr_2GeV_mh1000_pl1000_timescan->Draw("AC*");
  gr_2GeV_mh1000_pl1000_timescan->GetXaxis()->SetLimits(0.,1.1);
  gr_2GeV_mh1000_pl1000_timescan->GetHistogram()->SetMaximum(1.1);
  gr_2GeV_mh1000_pl1000_timescan->GetHistogram()->SetMinimum(0.);
  gr_2GeV_mh1000_pl1000_timescan->SetTitle("Signal Eff. vs Background Rejection with Delayed Hit Fraction > 0.5;Signal Eff (LLP, mh=1TeV, ctau=1m);QCD Rejection");
  gr_3GeV_mh1000_pl1000_timescan->SetLineColor(3); // green
  gr_3GeV_mh1000_pl1000_timescan->Draw("C*");
  gr_4GeV_mh1000_pl1000_timescan->SetLineColor(2); // red
  gr_4GeV_mh1000_pl1000_timescan->Draw("C*");
  auto legend_mh1000_pl1000 = new TLegend(0.15,0.15,0.55,0.4);
  legend_mh1000_pl1000->AddEntry(gr_2GeV_mh1000_pl1000_timescan,"2 GeV, TDC scanned");
  legend_mh1000_pl1000->AddEntry(gr_3GeV_mh1000_pl1000_timescan,"3 GeV, TDC scanned");
  legend_mh1000_pl1000->AddEntry(gr_4GeV_mh1000_pl1000_timescan,"4 GeV, TDC scanned");
  legend_mh1000_pl1000->Draw();
  ROC_DelayedHitFrac_mh1000_pl1000->SaveAs("plots/ROC_DelayedHitFrac_mh1000_pl1000.pdf");

  TCanvas *ROC_DelayedHitFrac_mh1000_pl500 = new TCanvas("ROC_DelayedHitFrac_mh1000_pl500","Graph Draw Options",200,10,600,600);
  ROC_DelayedHitFrac_mh1000_pl500->SetGrid();
  gr_2GeV_mh1000_pl500_timescan->SetLineColor(4); // blue
  gr_2GeV_mh1000_pl500_timescan->Draw("AC*");
  gr_2GeV_mh1000_pl500_timescan->GetXaxis()->SetLimits(0.,1.1);
  gr_2GeV_mh1000_pl500_timescan->GetHistogram()->SetMaximum(1.1);
  gr_2GeV_mh1000_pl500_timescan->GetHistogram()->SetMinimum(0.);
  gr_2GeV_mh1000_pl500_timescan->SetTitle("Signal Eff. vs Background Rejection with Delayed Hit Fraction > 0.5;Signal Eff (LLP, mh=1TeV, ctau=0.5m);QCD Rejection");
  gr_3GeV_mh1000_pl500_timescan->SetLineColor(3); // green
  gr_3GeV_mh1000_pl500_timescan->Draw("C*");
  gr_4GeV_mh1000_pl500_timescan->SetLineColor(2); // red 
  gr_4GeV_mh1000_pl500_timescan->Draw("C*");
  auto legend_mh1000_pl500 = new TLegend(0.15,0.15,0.55,0.4);
  legend_mh1000_pl500->AddEntry(gr_2GeV_mh1000_pl1000_timescan,"2 GeV, TDC scanned");
  legend_mh1000_pl500->AddEntry(gr_3GeV_mh1000_pl1000_timescan,"3 GeV, TDC scanned");
  legend_mh1000_pl500->AddEntry(gr_4GeV_mh1000_pl1000_timescan,"4 GeV, TDC scanned");
  legend_mh1000_pl500->Draw();
  ROC_DelayedHitFrac_mh1000_pl500->SaveAs("plots/ROC_DelayedHitFrac_mh1000_pl500.pdf");

  TCanvas *ROC_Frac_Mult_mh1000_pl1000 = new TCanvas("ROC_Frac_Mult_mh1000_pl1000","Graph Draw Options",200,10,600,600);
  ROC_Frac_Mult_mh1000_pl1000->SetGrid();
  gr_3GeV_mh1000_pl1000_timescan->SetLineColor(3); // green 
  gr_3GeV_mh1000_pl1000_timescan->Draw("AC*");
  gr_3GeV_mh1000_pl1000_timescan->GetXaxis()->SetLimits(0.,1.1); 
  gr_3GeV_mh1000_pl1000_timescan->GetHistogram()->SetMaximum(1.1);
  gr_3GeV_mh1000_pl1000_timescan->GetHistogram()->SetMinimum(0.);
  gr_3GeV_mh1000_pl1000_timescan->SetTitle("Signal Eff. vs Background Rejection;Signal Eff (LLP, mh=1TeV, ctau=1m);QCD Rejection");
  gr_3GeV_timescan_mh1000_pl1000_mult->SetLineColor(1);
  gr_3GeV_timescan_mh1000_pl1000_mult->Draw("C*");
  auto legend2_mh1000_pl1000 = new TLegend(0.15,0.15,0.9,0.35);
  legend2_mh1000_pl1000->AddEntry(gr_3GeV_mh1000_pl1000_timescan,"Delayed Hit Frac>0.5 and above 3GeV, TDC scanned");
  legend2_mh1000_pl1000->AddEntry(gr_3GeV_timescan_mh1000_pl1000_mult,"Quad Jet Multiplicity above 3GeV, 3ns, mult scanned");
  legend2_mh1000_pl1000->Draw();
  ROC_Frac_Mult_mh1000_pl1000->SaveAs("plots/ROC_Frac_Mult_3GeV_mh1000_pl1000.pdf");

  TCanvas *ROC_Frac_Mult_mh1000_pl500 = new TCanvas("ROC_Frac_Mult_mh1000_pl500","Graph Draw Options",200,10,600,600);
  ROC_Frac_Mult_mh1000_pl500->SetGrid();
  gr_3GeV_mh1000_pl500_timescan->SetLineColor(3); // green
  gr_3GeV_mh1000_pl500_timescan->Draw("AC*");
  gr_3GeV_mh1000_pl500_timescan->GetXaxis()->SetLimits(0.,1.1);
  gr_3GeV_mh1000_pl500_timescan->GetHistogram()->SetMaximum(1.1);
  gr_3GeV_mh1000_pl500_timescan->GetHistogram()->SetMinimum(0.);
  gr_3GeV_mh1000_pl500_timescan->SetTitle("Signal Eff. vs Background Rejection;Signal Eff (LLP, mh=1TeV, ctau=0.5m);QCD Rejection");
  gr_3GeV_timescan_mh1000_pl500_mult->SetLineColor(1);
  gr_3GeV_timescan_mh1000_pl500_mult->Draw("C*");
  auto legend2_mh1000_pl500 = new TLegend(0.15,0.15,0.9,0.35);
  legend2_mh1000_pl500->AddEntry(gr_3GeV_mh1000_pl500_timescan,"Delayed Hit Frac>0.5 and above 3GeV, TDC scanned");
  legend2_mh1000_pl500->AddEntry(gr_3GeV_timescan_mh1000_pl500_mult,"Quad Jet Multiplicity above 3GeV, 3ns, mult scanned");
  legend2_mh1000_pl500->Draw();
  ROC_Frac_Mult_mh1000_pl500->SaveAs("plots/ROC_Frac_Mult_3GeV_mh1000_pl500.pdf");

  return 0;
}
