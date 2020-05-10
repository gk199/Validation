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
  // mh=1TeV, mx=450 GeV
  double arr_signal_mult_mh1000_pl1000[8];
  ifstream Signal_mult_mh1000_pl1000;
  Signal_mult_mh1000_pl1000.open("MultiplicityHits3GeV3ns_ht120_Signal_mh1000_pl10.txt");
  int n = 0;
  while (Signal_mult_mh1000_pl1000 >> arr_signal_mult_mh1000_pl1000[n]) n++; // signal efficiency 
  Signal_mult_mh1000_pl1000.close();
  double arr_signal_mult_mh1000_pl500[8];
  ifstream Signal_mult_mh1000_pl500;
  Signal_mult_mh1000_pl500.open("MultiplicityHits3GeV3ns_ht120_Signal_mh1000_pl50.txt");
  n = 0;
  while (Signal_mult_mh1000_pl500 >> arr_signal_mult_mh1000_pl500[n]) n++; // signal efficiency                
  Signal_mult_mh1000_pl500.close();
  // mh=350 GeV, mx=160 GeV
  double arr_signal_mult_mh350_pl1000[8];
  ifstream Signal_mult_mh350_pl1000;
  Signal_mult_mh350_pl1000.open("MultiplicityHits3GeV3ns_ht120_Signal_mh350_pl100.txt");
  n = 0;
  while (Signal_mult_mh350_pl1000 >> arr_signal_mult_mh350_pl1000[n]) n++; // signal efficiency                                                                                     
  Signal_mult_mh350_pl1000.close();
  double arr_signal_mult_mh350_pl500[8];
  ifstream Signal_mult_mh350_pl500;
  Signal_mult_mh350_pl500.open("MultiplicityHits3GeV3ns_ht120_Signal_mh350_pl500.txt");
  n = 0;
  while (Signal_mult_mh350_pl500 >> arr_signal_mult_mh350_pl500[n]) n++; // signal efficiency                                                                                       
  Signal_mult_mh350_pl500.close();
  // QCD efficiencies
  double arr_background_mult[8];
  ifstream Background_mult;
  Background_mult.open("MultiplicityHits3GeV3ns_ht120_Background.txt");
  n = 0;
  while (Background_mult >>arr_background_mult[n]) n++; // background efficiency                                                             
  Background_mult.close();

  // read in neutrino gun rates at the various multiplicity cuts
  double arr_neutrino[22];
  ifstream Neutrino_single_quad_htSum120;
  Neutrino_single_quad_htSum120.open("single_quad_htSum120_JetRate.txt");
  n = 0;
  while (Neutrino_single_quad_htSum120 >> arr_neutrino[n]) n++;
  Neutrino_single_quad_htSum120.close();

  // split up neutrino gun into single jet rate, quad jet rate, and htSum120 GeV rate
  double singleJetRate[7], quadJetRate[7], htSum120Rate[7];
  double signal_mult_mh1000_pl1000[7], signal_mult_mh1000_pl500[7], signal_mult_mh350_pl1000[7], signal_mult_mh350_pl500[7], background_mult[7];
  double Original_htSumRate430, EffPl1000Mh1000_htSum430_Global, EffPl500Mh1000_htSum430_Global, EffPl1000Mh350_htSum430_Global, EffPl500Mh350_htSum430_Global, EffQCD_htSum430_Global;
  for (int i=0; i<7; i++) {
    singleJetRate[i] = arr_neutrino[i];
    quadJetRate[i] = arr_neutrino[i+7];
    htSum120Rate[i] = arr_neutrino[i+14];
    signal_mult_mh1000_pl1000[i] = arr_signal_mult_mh1000_pl1000[i];
    signal_mult_mh1000_pl500[i] = arr_signal_mult_mh1000_pl500[i];
    signal_mult_mh350_pl1000[i] = arr_signal_mult_mh350_pl1000[i];
    signal_mult_mh350_pl500[i] = arr_signal_mult_mh350_pl500[i];
    background_mult[i] = arr_background_mult[i];
  }
  // comparison points for htsum rates and efficiencies
  Original_htSumRate430 = arr_neutrino[21];
  EffPl1000Mh1000_htSum430_Global = arr_signal_mult_mh1000_pl1000[7];
  EffPl500Mh1000_htSum430_Global = arr_signal_mult_mh1000_pl500[7];
  EffPl1000Mh350_htSum430_Global = arr_signal_mult_mh350_pl1000[7];
  EffPl500Mh350_htSum430_Global = arr_signal_mult_mh350_pl500[7];
  EffQCD_htSum430_Global = arr_background_mult[7];

  // make TGraphs with the efficiency and neutrino gun rates from the multiplicity based trigger
  TGraph *gr_sing_LLP_mh1000_pl500 = new TGraph (6, signal_mult_mh1000_pl500, singleJetRate);
  TGraph *gr_quad_LLP_mh1000_pl500 = new TGraph (6, signal_mult_mh1000_pl500, quadJetRate);
  TGraph *gr_120_LLP_mh1000_pl500 = new TGraph (6, signal_mult_mh1000_pl500, htSum120Rate);
  TGraph *gr_sing_LLP_mh1000_pl1000 = new TGraph (6, signal_mult_mh1000_pl1000, singleJetRate);
  TGraph *gr_quad_LLP_mh1000_pl1000 = new TGraph (6, signal_mult_mh1000_pl1000, quadJetRate);
  TGraph *gr_120_LLP_mh1000_pl1000 = new TGraph (6, signal_mult_mh1000_pl1000, htSum120Rate);
  TGraph *gr_sing_LLP_mh350_pl500 = new TGraph (6, signal_mult_mh350_pl500, singleJetRate);
  TGraph *gr_quad_LLP_mh350_pl500 = new TGraph (6, signal_mult_mh350_pl500, quadJetRate);
  TGraph *gr_120_LLP_mh350_pl500 = new TGraph (6, signal_mult_mh350_pl500, htSum120Rate);
  TGraph *gr_sing_LLP_mh350_pl1000 = new TGraph (6, signal_mult_mh350_pl1000, singleJetRate);
  TGraph *gr_quad_LLP_mh350_pl1000 = new TGraph (6, signal_mult_mh350_pl1000, quadJetRate);
  TGraph *gr_120_LLP_mh350_pl1000 = new TGraph (6, signal_mult_mh350_pl1000, htSum120Rate);
  TGraph *gr_sing_QCD = new TGraph (6, background_mult, singleJetRate);
  TGraph *gr_quad_QCD = new TGraph (6, background_mult, quadJetRate);
  TGraph *gr_120_QCD = new TGraph (6, background_mult, htSum120Rate);
  // ************** Markers for comparison with current benchmark performance ******
  TMarker *m_QCD_htSum430 = new TMarker (EffQCD_htSum430_Global, Original_htSumRate430, 21);
  TMarker *m_LLP_htSum430 = new TMarker (EffPl500Mh1000_htSum430_Global, Original_htSumRate430, 21);
  TMarker *m_LLPpl1000_htSum430 = new TMarker (EffPl1000Mh1000_htSum430_Global, Original_htSumRate430, 21);
  TMarker *m_LLPMh350_htSum430 = new TMarker (EffPl500Mh350_htSum430_Global, Original_htSumRate430, 21);
  TMarker *m_LLPMh350pl1000_htSum430 = new TMarker (EffPl1000Mh350_htSum430_Global, Original_htSumRate430, 21);

  // single, quad jet rate vs LLP mh=1TeV, ct=0.5m efficiency
  TCanvas *c1_sing_quad_mh1000_pl500 = new TCanvas("c1_sing_quad_mh1000_pl500","Graph Draw Options",200,10,600,400);
  gr_sing_LLP_mh1000_pl500->SetLineColor(4); // blue
  gr_sing_LLP_mh1000_pl500->GetHistogram()->SetMinimum(-5.);
  gr_quad_LLP_mh1000_pl500->GetHistogram()->SetMinimum(-5.);
  gr_sing_LLP_mh1000_pl500->Draw("AC*");
  gr_sing_LLP_mh1000_pl500->SetTitle("Jet Rate vs. Signal Efficiency for Regional Multiplicity Cut;LLP, mh=1TeV, c#scale[1.2]{#tau}=0.5m Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  gr_quad_LLP_mh1000_pl500->SetLineColor(2); // red
  gr_quad_LLP_mh1000_pl500->Draw("C*");
  auto legend_sing_quad_mh1000_pl500 = new TLegend(0.15,0.7,0.5,0.9);
  legend_sing_quad_mh1000_pl500->AddEntry(gr_sing_LLP_mh1000_pl500,"Single Jet Rate at 60 GeV");
  legend_sing_quad_mh1000_pl500->AddEntry(gr_quad_LLP_mh1000_pl500,"Quad Jet Rate at 60 GeV");
  legend_sing_quad_mh1000_pl500->Draw();
  gr_sing_LLP_mh1000_pl500->GetXaxis()->SetLimits(0.,1.);
  c1_sing_quad_mh1000_pl500->SaveAs("plots/NuGun_JetRates_vs_SignalEff_Pl500Mh1000.pdf");
  // ct=1m
  TCanvas *c1_sing_quad_mh1000_pl1000 = new TCanvas("c1_sing_quad_mh1000_pl1000","Graph Draw Options",200,10,600,400);
  gr_sing_LLP_mh1000_pl1000->SetLineColor(4); // blue
  gr_sing_LLP_mh1000_pl1000->GetHistogram()->SetMinimum(-5.);
  gr_quad_LLP_mh1000_pl1000->GetHistogram()->SetMinimum(-5.);
  gr_sing_LLP_mh1000_pl1000->Draw("AC*");
  gr_sing_LLP_mh1000_pl1000->SetTitle("Jet Rate vs. Signal Efficiency for Regional Multiplicity Cut;LLP, mh=1TeV, c#scale[1.2]{#tau}=1m Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  gr_quad_LLP_mh1000_pl1000->SetLineColor(2); // red 
  gr_quad_LLP_mh1000_pl1000->Draw("C*");
  auto legend_sing_quad_mh1000_pl1000 = new TLegend(0.15,0.7,0.5,0.9);
  legend_sing_quad_mh1000_pl1000->AddEntry(gr_sing_LLP_mh1000_pl1000,"Single Jet Rate at 60 GeV");
  legend_sing_quad_mh1000_pl1000->AddEntry(gr_quad_LLP_mh1000_pl1000,"Quad Jet Rate at 60 GeV");
  legend_sing_quad_mh1000_pl1000->Draw();
  gr_sing_LLP_mh1000_pl1000->GetXaxis()->SetLimits(0.,1.);
  c1_sing_quad_mh1000_pl1000->SaveAs("plots/NuGun_JetRates_vs_SignalEff_Pl1000Mh1000.pdf");

  // single, quad jet rate vs LLP mh=350GeV, ct=0.5m efficiency 
  TCanvas *c1_sing_quad_mh350_pl500 = new TCanvas("c1_sing_quad_mh350_pl500","Graph Draw Options",200,10,600,400);
  gr_sing_LLP_mh350_pl500->SetLineColor(4); // blue
  gr_sing_LLP_mh350_pl500->GetHistogram()->SetMinimum(-5.);
  gr_quad_LLP_mh350_pl500->GetHistogram()->SetMinimum(-5.);
  gr_sing_LLP_mh350_pl500->Draw("AC*");
  gr_sing_LLP_mh350_pl500->SetTitle("Jet Rate vs. Signal Efficiency for Regional Multiplicity Cut;LLP, mh=350GeV, c#scale[1.2]{#tau}=0.5m Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  gr_quad_LLP_mh350_pl500->SetLineColor(2); // red
  gr_quad_LLP_mh350_pl500->Draw("C*");
  auto legend_sing_quad_mh350_pl500 = new TLegend(0.15,0.7,0.5,0.9);
  legend_sing_quad_mh350_pl500->AddEntry(gr_sing_LLP_mh350_pl500,"Single Jet Rate at 60 GeV");
  legend_sing_quad_mh350_pl500->AddEntry(gr_quad_LLP_mh350_pl500,"Quad Jet Rate at 60 GeV");
  legend_sing_quad_mh350_pl500->Draw();
  gr_sing_LLP_mh350_pl500->GetXaxis()->SetLimits(0.,1.);
  c1_sing_quad_mh350_pl500->SaveAs("plots/NuGun_JetRates_vs_SignalEff_Pl500Mh350.pdf");
  // ct=1m
  TCanvas *c1_sing_quad_mh350_pl1000 = new TCanvas("c1_sing_quad_mh350_pl1000","Graph Draw Options",200,10,600,400);
  gr_sing_LLP_mh350_pl1000->SetLineColor(4); // blue
  gr_sing_LLP_mh350_pl1000->GetHistogram()->SetMinimum(-5.);
  gr_quad_LLP_mh350_pl1000->GetHistogram()->SetMinimum(-5.);
  gr_sing_LLP_mh350_pl1000->Draw("AC*");
  gr_sing_LLP_mh350_pl1000->SetTitle("Jet Rate vs. Signal Efficiency for Regional Multiplicity Cut;LLP, mh=350GeV, c#scale[1.2]{#tau}=1m Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  gr_quad_LLP_mh350_pl1000->SetLineColor(2); // red
  gr_quad_LLP_mh350_pl1000->Draw("C*");
  auto legend_sing_quad_mh350_pl1000 = new TLegend(0.15,0.7,0.5,0.9);
  legend_sing_quad_mh350_pl1000->AddEntry(gr_sing_LLP_mh350_pl1000,"Single Jet Rate at 60 GeV");
  legend_sing_quad_mh350_pl1000->AddEntry(gr_quad_LLP_mh350_pl1000,"Quad Jet Rate at 60 GeV");
  legend_sing_quad_mh350_pl1000->Draw();
  gr_sing_LLP_mh350_pl1000->GetXaxis()->SetLimits(0.,1.);
  c1_sing_quad_mh350_pl1000->SaveAs("plots/NuGun_JetRates_vs_SignalEff_Pl1000Mh350.pdf");

  // single, quad jet rate vs QCD efficency
  TCanvas *c2_sing_quad = new TCanvas("c2_sing_quad","Graph Draw Options",200,10,600,400);
  gr_sing_QCD->SetLineColor(4); // blue  
  gr_sing_QCD->GetHistogram()->SetMinimum(-5.);
  gr_quad_QCD->GetHistogram()->SetMinimum(-5.);
  gr_sing_QCD->Draw("AC*");
  gr_sing_QCD->SetTitle("Jet Rate vs. Background Efficiency for Regional Multiplicity Cut;QCD Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  gr_quad_QCD->SetLineColor(2); // red  
  gr_quad_QCD->Draw("C*");
  auto legend_sing_quad = new TLegend(0.55,0.7,0.9,0.9);
  legend_sing_quad->AddEntry(gr_sing_QCD,"Single Jet Rate at 60 GeV");
  legend_sing_quad->AddEntry(gr_quad_QCD,"Quad Jet Rate at 60 GeV");
  legend_sing_quad->Draw();
  gr_sing_QCD->GetXaxis()->SetLimits(0.,1.);
  c2_sing_quad->SaveAs("plots/NuGun_JetRates_vs_BackgroundEff.pdf");

  // htSum rate, LLP mh=1TeV, ct=0.5m efficiency, regional
  TCanvas *c1_htSum_mh1000_pl500 = new TCanvas("c1_htSum_mh1000_pl500","Graph Draw Options",200,10,600,400);
  gr_120_LLP_mh1000_pl500->SetLineColor(4); // blue 
  gr_120_LLP_mh1000_pl500->GetHistogram()->SetMinimum(-5.);
  gr_120_LLP_mh1000_pl500->Draw("AL*");
  gr_120_LLP_mh1000_pl500->SetTitle("htSum Rate vs. Signal Efficiency for Regional Multiplicity Cut;LLP, mh=1TeV, c#scale[1.2]{#tau}=0.5m Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  m_LLP_htSum430->SetMarkerStyle(21);
  m_LLP_htSum430->SetMarkerColor(3);
  m_LLP_htSum430->Draw();
  auto legend_htSum_mh1000_pl500 = new TLegend(0.15,0.7,0.5,0.9);
  legend_htSum_mh1000_pl500->AddEntry(gr_120_LLP_mh1000_pl500,"H_{T}>120 GeV, with timing cuts");
  legend_htSum_mh1000_pl500->AddEntry(m_LLP_htSum430,"H_{T}>430 GeV, no timing cuts");
  legend_htSum_mh1000_pl500->Draw();
  gr_120_LLP_mh1000_pl500->GetXaxis()->SetLimits(0.,1.);
  gr_120_LLP_mh1000_pl500->GetHistogram()->SetMinimum(1.);
  gr_120_LLP_mh1000_pl500->GetHistogram()->SetMaximum(10000.);
  c1_htSum_mh1000_pl500->SetLogy(); 
  c1_htSum_mh1000_pl500->SetGrid();
  c1_htSum_mh1000_pl500->SaveAs("plots/NuGun_htSumRates_vs_SignalEff_Pl500Mh1000.pdf");
  // ct=1m
  TCanvas *c1_htSum_mh1000_pl1000 = new TCanvas("c1_htSum_mh1000_pl1000","Graph Draw Options",200,10,600,400);
  gr_120_LLP_mh1000_pl1000->SetLineColor(4); // blue 
  gr_120_LLP_mh1000_pl1000->GetHistogram()->SetMinimum(-5.);
  gr_120_LLP_mh1000_pl1000->Draw("AL*");
  gr_120_LLP_mh1000_pl1000->SetTitle("htSum Rate vs. Signal Efficiency for Regional Multiplicity Cut;LLP, mh=1TeV, c#scale[1.2]{#tau}=1m Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  m_LLPpl1000_htSum430->SetMarkerStyle(21);
  m_LLPpl1000_htSum430->SetMarkerColor(3);
  m_LLPpl1000_htSum430->Draw();
  auto legend_htSum_mh1000_pl1000 = new TLegend(0.15,0.7,0.5,0.9);
  legend_htSum_mh1000_pl1000->AddEntry(gr_120_LLP_mh1000_pl1000,"H_{T}>120 GeV, with timing cuts");
  legend_htSum_mh1000_pl1000->AddEntry(m_LLPpl1000_htSum430,"H_{T}>430 GeV, no timing cuts");
  legend_htSum_mh1000_pl1000->Draw();
  gr_120_LLP_mh1000_pl1000->GetXaxis()->SetLimits(0.,1.);
  gr_120_LLP_mh1000_pl1000->GetHistogram()->SetMinimum(1.);
  gr_120_LLP_mh1000_pl1000->GetHistogram()->SetMaximum(10000.);
  c1_htSum_mh1000_pl1000->SetLogy();
  c1_htSum_mh1000_pl1000->SetGrid();
  c1_htSum_mh1000_pl1000->SaveAs("plots/NuGun_htSumRates_vs_SignalEff_Pl1000Mh1000.pdf");

  // htSum rate, LLP mh=350GeV, ct=0.5m efficiency, regional
  TCanvas *c1_htSum_mh350_pl500 = new TCanvas("c1_htSum_mh350_pl500","Graph Draw Options",200,10,600,400);
  gr_120_LLP_mh350_pl500->SetLineColor(4); // blue
  gr_120_LLP_mh350_pl500->GetHistogram()->SetMinimum(-5.);
  gr_120_LLP_mh350_pl500->Draw("AL*");
  gr_120_LLP_mh350_pl500->SetTitle("htSum Rate vs. Signal Efficiency for Regional Multiplicity Cut;LLP, mh=350GeV, c#scale[1.2]{#tau}=0.5m Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  m_LLPMh350_htSum430->SetMarkerStyle(21);
  m_LLPMh350_htSum430->SetMarkerColor(3);
  m_LLPMh350_htSum430->Draw();
  auto legend_htSum_mh350_pl500 = new TLegend(0.15,0.7,0.5,0.9);
  legend_htSum_mh350_pl500->AddEntry(gr_120_LLP_mh350_pl500,"H_{T}>120 GeV, with timing cuts");
  legend_htSum_mh350_pl500->AddEntry(m_LLPMh350_htSum430,"H_{T}>430 GeV, no timing cuts");
  legend_htSum_mh350_pl500->Draw();
  gr_120_LLP_mh350_pl500->GetXaxis()->SetLimits(0.,1.);
  gr_120_LLP_mh350_pl500->GetHistogram()->SetMinimum(1.);
  gr_120_LLP_mh350_pl500->GetHistogram()->SetMaximum(10000.);
  c1_htSum_mh350_pl500->SetLogy();
  c1_htSum_mh350_pl500->SetGrid();
  c1_htSum_mh350_pl500->SaveAs("plots/NuGun_htSumRates_vs_SignalEff_Pl500Mh350.pdf");
  // ct=1m
  TCanvas *c1_htSum_mh350_pl1000 = new TCanvas("c1_htSum_mh350_pl1000","Graph Draw Options",200,10,600,400);
  gr_120_LLP_mh350_pl1000->SetLineColor(4); // blue 
  gr_120_LLP_mh350_pl1000->GetHistogram()->SetMinimum(-5.);
  gr_120_LLP_mh350_pl1000->Draw("AL*");
  gr_120_LLP_mh350_pl1000->SetTitle("htSum Rate vs. Signal Efficiency for Regional Multiplicity Cut;LLP, mh=350GeV, c#scale[1.2]{#tau}=1m Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  m_LLPMh350pl1000_htSum430->SetMarkerStyle(21);
  m_LLPMh350pl1000_htSum430->SetMarkerColor(3);
  m_LLPMh350pl1000_htSum430->Draw();
  auto legend_htSum_mh350_pl1000 = new TLegend(0.15,0.7,0.5,0.9);
  legend_htSum_mh350_pl1000->AddEntry(gr_120_LLP_mh350_pl1000,"H_{T}>120 GeV, with timing cuts");
  legend_htSum_mh350_pl1000->AddEntry(m_LLPMh350pl1000_htSum430,"H_{T}>430 GeV, no timing cuts");
  legend_htSum_mh350_pl1000->Draw();
  gr_120_LLP_mh350_pl1000->GetXaxis()->SetLimits(0.,1.);
  gr_120_LLP_mh350_pl1000->GetHistogram()->SetMinimum(1.);
  gr_120_LLP_mh350_pl1000->GetHistogram()->SetMaximum(10000.);
  c1_htSum_mh350_pl1000->SetLogy();
  c1_htSum_mh350_pl1000->SetGrid();
  c1_htSum_mh350_pl1000->SaveAs("plots/NuGun_htSumRates_vs_SignalEff_Pl1000Mh350.pdf");

  // htSum rate, QCD, regional
  TCanvas *c2_htSum = new TCanvas("c2_htSum","Graph Draw Options",200,10,600,400);
  gr_120_QCD->SetLineColor(4); // blue  
  gr_120_QCD->GetHistogram()->SetMinimum(-5.);
  gr_120_QCD->Draw("AL*");
  gr_120_QCD->SetTitle("htSum Rate vs. Background Efficiency for Regional Multiplicity Cut;QCD Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  m_QCD_htSum430->SetMarkerStyle(21);
  m_QCD_htSum430->SetMarkerColor(3);
  m_QCD_htSum430->Draw();
  gr_120_QCD->GetXaxis()->SetLimits(0.,1.);
  c2_htSum->SetGrid();
  auto legend2_htSum = new TLegend(0.15,0.7,0.5,0.9);
  legend2_htSum->AddEntry(gr_120_QCD,"H_{T}>120 GeV, with timing cuts");
  legend2_htSum->AddEntry(m_QCD_htSum430,"H_{T}>430 GeV, no timing cuts");
  legend2_htSum->Draw();
  c2_htSum->SetLogy();
  gr_120_QCD->GetHistogram()->SetMinimum(1.);
  gr_120_QCD->GetHistogram()->SetMaximum(10000.);
  c2_htSum->SaveAs("plots/NuGun_htSumRates_vs_BackgroundEff.pdf");


  return 0;
}
