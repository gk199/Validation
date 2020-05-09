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
  // read in signal and background efficiency files as arrays -- these are for the trigger based on quad jet hit multiplicity at 3 GeV, 3ns 
  double arr_signal_mult_mh1000_pl1000[7];
  ifstream Signal_mult_mh1000_pl1000;
  Signal_mult_mh1000_pl1000.open("MultiplicityHits3GeV3ns_Signal_mh1000_pl10.txt");
  int n = 0;
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

  // read in neutrino gun rates at the various multiplicity cuts
  double arr_neutrino[21];
  ifstream Neutrino_single_quad_htSum120;
  Neutrino_single_quad_htSum120.open("single_quad_htSum120_JetRate.txt");
  n = 0;
  while (Neutrino_single_quad_htSum120 >> arr_neutrino[n]) n++;
  Neutrino_single_quad_htSum120.close();

  // split up neutrino gun into single jet rate, quad jet rate, and htSum120 GeV rate
  double singleJetRate[7], quadJetRate[7], htSum120Rate[7];
  for (int i=0; i<7; i++) {
    singleJetRate[i] = arr_neutrino[i];
    quadJetRate[i] = arr_neutrino[i+7];
    htSum120Rate[i] = arr_neutrino[i+14];
  }

  // make TGraphs with the efficiency and neutrino gun rates from the multiplicity based trigger
  TGraph *gr_sing_LLP = new TGraph (6, arr_signal_mult_mh1000_pl500, singleJetRate);
  TGraph *gr_quad_LLP = new TGraph (6, arr_signal_mult_mh1000_pl500, quadJetRate);
  TGraph *gr_120_LLP = new TGraph (6, arr_signal_mult_mh1000_pl500, htSum120Rate);
  TGraph *gr_sing_QCD = new TGraph (6, arr_background_mult, singleJetRate);
  TGraph *gr_quad_QCD = new TGraph (6, arr_background_mult, quadJetRate);
  TGraph *gr_120_QCD = new TGraph (6, arr_background_mult, htSum120Rate);

  // single, quad jet rate vs LLP mh=1TeV, ct=0.5m efficiency
  TCanvas *c1_sing_quad = new TCanvas("c1_sing_quad","Graph Draw Options",200,10,600,400);
  gr_sing_LLP->SetLineColor(4); // blue
  gr_sing_LLP->GetHistogram()->SetMinimum(-5.);
  gr_quad_LLP->GetHistogram()->SetMinimum(-5.);
  gr_sing_LLP->Draw("AC*");
  gr_sing_LLP->SetTitle("Jet Rate vs. Signal Efficiency for Regional Multiplicity Cut;LLP, mh=1TeV, c#scale[1.2]{#tau}=0.5m Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  gr_quad_LLP->SetLineColor(2); // red
  gr_quad_LLP->Draw("C*");
  auto legend = new TLegend(0.15,0.7,0.5,0.9);
  legend->AddEntry(gr_sing_LLP,"Single Jet Rate at 60 GeV");
  legend->AddEntry(gr_quad_LLP,"Quad Jet Rate at 60 GeV");
  legend->Draw();
  gr_sing_LLP->GetXaxis()->SetLimits(0.,1.);
  c1_sing_quad->SaveAs("plots/NuGun_JetRates_vs_SignalEff_Pl500Mh1000.pdf");

  // single, quad jet rate vs QCD efficency
  TCanvas *c2_sing_quad = new TCanvas("c2_sing_quad","Graph Draw Options",200,10,600,400);
  gr_sing_QCD->SetLineColor(4); // blue  
  gr_sing_QCD->GetHistogram()->SetMinimum(-5.);
  gr_quad_QCD->GetHistogram()->SetMinimum(-5.);
  gr_sing_QCD->Draw("AC*");
  gr_sing_QCD->SetTitle("Jet Rate vs. Background Efficiency for Regional Multiplicity Cut;QCD Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  gr_quad_QCD->SetLineColor(2); // red  
  gr_quad_QCD->Draw("C*");
  auto legend2 = new TLegend(0.55,0.7,0.9,0.9);
  legend2->AddEntry(gr_sing_QCD,"Single Jet Rate at 60 GeV");
  legend2->AddEntry(gr_quad_QCD,"Quad Jet Rate at 60 GeV");
  legend2->Draw();
  gr_sing_QCD->GetXaxis()->SetLimits(0.,1.);
  c2_sing_quad->SaveAs("plots/NuGun_JetRates_vs_BackgroundEff.pdf");

  // htSum rate, LLP mh=1TeV, ct=0.5m efficiency, regional
  TCanvas *c1_350_120 = new TCanvas("c1_350_120","Graph Draw Options",200,10,600,400);
  gr_120_LLP->SetLineColor(4); // blue 
  gr_120_LLP->GetHistogram()->SetMinimum(-5.);
  gr_120_LLP->Draw("AL*");
  gr_120_LLP->SetTitle("htSum Rate vs. Signal Efficiency for Regional Multiplicity Cut;LLP, mh=1TeV, c#scale[1.2]{#tau}=0.5m Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  //  m_LLP_htSum430->SetMarkerStyle(21);
  //  m_LLP_htSum430->SetMarkerColor(3);
  auto legend_350_120 = new TLegend(0.15,0.7,0.5,0.9);
  legend_350_120->AddEntry(gr_120_LLP,"H_{T}>120 GeV, with timing cuts");
  legend_350_120->Draw();
  gr_120_LLP->GetXaxis()->SetLimits(0.,1.);
  gr_120_LLP->GetHistogram()->SetMinimum(1.);
  gr_120_LLP->GetHistogram()->SetMaximum(10000.);
  c1_350_120->SetLogy(); 
  c1_350_120->SetGrid();
  c1_350_120->SaveAs("plots/NuGun_htSumRates_vs_SignalEff_Pl500Mh1000.pdf");

  // htSum rate, QCD, regional
  TCanvas *c2_350_120 = new TCanvas("c2_350_120","Graph Draw Options",200,10,600,400);
  gr_120_QCD->SetLineColor(4); // blue  
  gr_120_QCD->GetHistogram()->SetMinimum(-5.);
  gr_120_QCD->Draw("AL*");
  gr_120_QCD->SetTitle("htSum Rate vs. Background Efficiency for Regional Multiplicity Cut;QCD Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  //  m_QCD_htSum430->SetMarkerStyle(21);
  //  m_QCD_htSum430->SetMarkerColor(3);
  //  m_QCD_htSum430->Draw();
  gr_120_QCD->GetXaxis()->SetLimits(0.,1.);
  c2_350_120->SetGrid();
  auto legend2_350_120 = new TLegend(0.15,0.7,0.5,0.9);
  legend2_350_120->AddEntry(gr_120_QCD,"H_{T}>120 GeV, with timing cuts");
  legend2_350_120->Draw();
  c2_350_120->SetLogy();
  gr_120_QCD->GetHistogram()->SetMinimum(1.);
  gr_120_QCD->GetHistogram()->SetMaximum(10000.);
  c2_350_120->SaveAs("plots/NuGun_htSumRates_vs_BackgroundEff.pdf");


  return 0;
}
