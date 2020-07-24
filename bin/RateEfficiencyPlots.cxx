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
  // read in signal and background files 
  // first line is efficiency at HT>120 and mult HEHB >= 10
  // second line is efficiency at HT>360 with no timing restrictions
  // Signals
  // mh=125 pl=3m
  double signal_mh125_pl3000[6];
  ifstream mh125_pl3000;
  mh125_pl3000.open("MultiplicityHits50ADC3ns_ht120_Signal_mh125__mx50__pl3000_.txt");
  int n=0;
  while (mh125_pl3000 >> signal_mh125_pl3000[n]) n++;
  mh125_pl3000.close();

  // mh=1000 pl=10m
  double signal_mh1000_pl10000[6];
  ifstream mh1000_pl10000;
  mh1000_pl10000.open("MultiplicityHits50ADC3ns_ht120_Signal_mh1000_mx450_pl10000.txt");
  n=0;
  while (mh1000_pl10000 >> signal_mh1000_pl10000[n]) n++;
  mh1000_pl10000.close();

  // mh=250 pl=1m 
  double signal_mh250_pl1000[6];
  ifstream mh250_pl1000;
  mh250_pl1000.open("MultiplicityHits50ADC3ns_ht120_Signal_mh250__mx120_pl1000_.txt");
  n=0;
  while (mh250_pl1000 >> signal_mh250_pl1000[n]) n++;
  mh250_pl1000.close();

  // Background
  double background[6];
  ifstream QCD;
  QCD.open("MultiplicityHits50ADC3ns_ht120_Background.txt");
  n=0;
  while (QCD >> background[n]) n++;
  QCD.close();
  // Neutrino gun for rates
  double rates[6];
  ifstream nugun;
  nugun.open("NuGunRates.txt");
  n=0;
  while (nugun >> rates[n]) n++;
  nugun.close();

  // mh = 125 GeV
  TGraph *gr_LLP_mh125_pl3000 = new TGraph(2, signal_mh125_pl3000, rates);
  TMarker *m_mh125_pl3000_ht360 = new TMarker(signal_mh125_pl3000[5], rates[5], 21);
  TMarker *m_mh125_pl3000_ht120 = new TMarker(signal_mh125_pl3000[0], rates[0], 21);
  // mh = 1000 GeV
  TGraph *gr_LLP_mh1000_pl10000 = new TGraph(2, signal_mh1000_pl10000, rates);
  TMarker *m_mh1000_pl10000_ht360 = new TMarker(signal_mh1000_pl10000[5], rates[5], 21);
  TMarker *m_mh1000_pl10000_ht120 = new TMarker(signal_mh1000_pl10000[0], rates[0], 21);
  // mh = 250 GeV
  TGraph *gr_LLP_mh250_pl1000 = new TGraph(2, signal_mh250_pl1000, rates);
  TMarker *m_mh250_pl1000_ht360 = new TMarker(signal_mh250_pl1000[5], rates[5], 21);
  TMarker *m_mh250_pl1000_ht120 = new TMarker(signal_mh250_pl1000[0], rates[0], 21);
  // QCD
  TGraph *gr_background = new TGraph(2, background, rates);
  TMarker *m_background_ht360 = new TMarker(background[5], rates[5], 21);
  TMarker *m_background_ht120 = new TMarker(background[0], rates[0], 21);

  // mh = 125 GeV
  TCanvas *c1_LLP_mh125_pl3000 = new TCanvas("c1_LLP_mh125_pl3000","Graph Draw Options",200,10,600,400);
  gr_LLP_mh125_pl3000->GetHistogram()->SetMinimum(-5.);
  gr_LLP_mh125_pl3000->Draw("A*");
  gr_LLP_mh125_pl3000->SetTitle("htSum Rate vs. Signal Efficiency for >=10 Cells >=50ADC,3ns near 4 L1 Jets;LLP, mh=125GeV, c#scale[1.2]{#tau}=1m Efficiency;Neutrino Gun Rate (Hz, unnorm)   ");
  gr_LLP_mh125_pl3000->GetXaxis()->SetLimits(0.,1.);
  gr_LLP_mh125_pl3000->GetHistogram()->SetMinimum(10000.);
  gr_LLP_mh125_pl3000->GetHistogram()->SetMaximum(100000000.);

  m_mh125_pl3000_ht360->SetMarkerStyle(21);
  m_mh125_pl3000_ht360->SetMarkerColor(3);
  m_mh125_pl3000_ht360->Draw();
  m_mh125_pl3000_ht120->SetMarkerStyle(21);
  m_mh125_pl3000_ht120->SetMarkerColor(2);
  m_mh125_pl3000_ht120->Draw();
  auto legend2_htSum = new TLegend(0.55,0.15,0.9,0.35);
  legend2_htSum->AddEntry(m_mh125_pl3000_ht360,"H_{T}>360 GeV, no timing cuts");
  legend2_htSum->AddEntry(m_mh125_pl3000_ht120,"H_{T}>120 GeV, with timing cuts");
  legend2_htSum->Draw();
  c1_LLP_mh125_pl3000->SetLogy();
  c1_LLP_mh125_pl3000->SetGrid();
  c1_LLP_mh125_pl3000->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/EffRate_mh125_pl3000.pdf");

  // mh = 250 GeV
  TCanvas *c1_LLP_mh250_pl1000 = new TCanvas("c1_LLP_mh250_pl1000","Graph Draw Options",200,10,600,400);
  gr_LLP_mh250_pl1000->GetHistogram()->SetMinimum(-5.);
  gr_LLP_mh250_pl1000->Draw("A*");
  gr_LLP_mh250_pl1000->SetTitle("htSum Rate vs. Signal Efficiency for >=10 Cells >=50ADC,3ns near 4 L1 Jets;LLP, mh=250GeV, c#scale[1.2]{#tau}=1m Efficiency;Neutrino Gun Rate (Hz, unnorm)   ");
  gr_LLP_mh250_pl1000->GetXaxis()->SetLimits(0.,1.);
  gr_LLP_mh250_pl1000->GetHistogram()->SetMinimum(10000.);
  gr_LLP_mh250_pl1000->GetHistogram()->SetMaximum(100000000.);

  m_mh250_pl1000_ht360->SetMarkerStyle(21);
  m_mh250_pl1000_ht360->SetMarkerColor(3);
  m_mh250_pl1000_ht360->Draw();
  m_mh250_pl1000_ht120->SetMarkerStyle(21);
  m_mh250_pl1000_ht120->SetMarkerColor(2);
  m_mh250_pl1000_ht120->Draw();
  auto legend4_htSum = new TLegend(0.55,0.15,0.9,0.35);
  legend4_htSum->AddEntry(m_mh250_pl1000_ht360,"H_{T}>360 GeV, no timing cuts");
  legend4_htSum->AddEntry(m_mh250_pl1000_ht120,"H_{T}>120 GeV, with timing cuts");
  legend4_htSum->Draw();
  c1_LLP_mh250_pl1000->SetLogy();
  c1_LLP_mh250_pl1000->SetGrid();
  c1_LLP_mh250_pl1000->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/EffRate_mh250_pl1000.pdf");

  // mh = 1000 GeV
  TCanvas *c1_LLP_mh1000_pl10000 = new TCanvas("c1_LLP_mh1000_pl10000","Graph Draw Options",200,10,600,400);
  gr_LLP_mh1000_pl10000->GetHistogram()->SetMinimum(-5.);
  gr_LLP_mh1000_pl10000->Draw("A*");
  gr_LLP_mh1000_pl10000->SetTitle("htSum Rate vs. Signal Efficiency for >=10 Cells >=50ADC,3ns near 4 L1 Jets;LLP, mh=1000GeV, c#scale[1.2]{#tau}=10m Efficiency;Neutrino Gun Rate (Hz, unnorm)   ");
  gr_LLP_mh1000_pl10000->GetXaxis()->SetLimits(0.,1.);
  gr_LLP_mh1000_pl10000->GetHistogram()->SetMinimum(10000.);
  gr_LLP_mh1000_pl10000->GetHistogram()->SetMaximum(100000000.);

  m_mh1000_pl10000_ht360->SetMarkerStyle(21);
  m_mh1000_pl10000_ht360->SetMarkerColor(3);
  m_mh1000_pl10000_ht360->Draw();
  m_mh1000_pl10000_ht120->SetMarkerStyle(21);
  m_mh1000_pl10000_ht120->SetMarkerColor(2);
  m_mh1000_pl10000_ht120->Draw();
  auto legend7_htSum = new TLegend(0.55,0.15,0.9,0.35);
  legend7_htSum->AddEntry(m_mh1000_pl10000_ht360,"H_{T}>360 GeV, no timing cuts");
  legend7_htSum->AddEntry(m_mh1000_pl10000_ht120,"H_{T}>120 GeV, with timing cuts");
  legend7_htSum->Draw();
  c1_LLP_mh1000_pl10000->SetLogy();
  c1_LLP_mh1000_pl10000->SetGrid();
  c1_LLP_mh1000_pl10000->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/EffRate_mh1000_pl10000.pdf");

  // background
  TCanvas *c1_background = new TCanvas("c1_background","Graph Draw Options",200,10,600,400);
  gr_background->GetHistogram()->SetMinimum(-5.);
  gr_background->Draw("A*");
  gr_background->SetTitle("htSum Rate vs. Background Efficiency for >=10 Cells >=50ADC,3ns near 4 L1 Jets;QCD Efficiency;Neutrino Gun Rate (Hz, unnorm)   ");
  gr_background->GetXaxis()->SetLimits(0.,1.);
  gr_background->GetHistogram()->SetMinimum(10000.);
  gr_background->GetHistogram()->SetMaximum(100000000.);

  m_background_ht360->SetMarkerStyle(21);
  m_background_ht360->SetMarkerColor(3);
  m_background_ht360->Draw();
  m_background_ht120->SetMarkerStyle(21);
  m_background_ht120->SetMarkerColor(2);
  m_background_ht120->Draw();
  auto legend9_htSum = new TLegend(0.55,0.15,0.9,0.35);
  legend9_htSum->AddEntry(m_background_ht360,"H_{T}>360 GeV, no timing cuts");
  legend9_htSum->AddEntry(m_background_ht120,"H_{T}>120 GeV, with timing cuts");
  legend9_htSum->Draw();
  c1_background->SetLogy();
  c1_background->SetGrid();
  c1_background->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/EffRate_QCDbackground.pdf");
}
