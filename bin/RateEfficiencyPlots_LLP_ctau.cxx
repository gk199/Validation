#include<iostream>
#include<fstream>
#include<vector>
#include "TGraph.h"
#include "TMarker.h"
#include "TLine.h"
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
  mh125_pl3000.open("MultiplicityHits50ADC3ns_ht120_Signal_MH-125_MFF-50_CTau-3000mm_Tun.txt");
  //  mh125_pl3000.open("MultiplicityHits50ADC3ns_ht120_Signal_mh125__mx50__pl3000_.txt");
  int n=0;
  while (mh125_pl3000 >> signal_mh125_pl3000[n]) n++;
  mh125_pl3000.close();
  /*
  // mh=125 pl=30m      
  double signal_mh125_pl30000[6];
  ifstream mh125_pl30000;
  mh125_pl30000.open("MultiplicityHits50ADC3ns_ht120_Signal_MH-125_MFF-50_CTau-30000mm_Tu.txt");  //MultiplicityHits50ADC3ns_ht120_Signal_mh125__mx50__pl30000.txt");
  n=0;
  while (mh125_pl30000 >> signal_mh125_pl30000[n]) n++;
  mh125_pl30000.close();
  */

  // mh=1000 pl=10m
  double signal_mh1000_pl10000[6];
  ifstream mh1000_pl10000;
  mh1000_pl10000.open("MultiplicityHits50ADC3ns_ht120_Signal_MH-1000_MFF-450_CTau-10000mm_.txt"); //MultiplicityHits50ADC3ns_ht120_Signal_mh1000_mx450_pl10000.txt");
  n=0;
  while (mh1000_pl10000 >> signal_mh1000_pl10000[n]) n++;
  mh1000_pl10000.close();
  // mh=1000 pl=100m
  double signal_mh1000_pl100000[6];
  ifstream mh1000_pl100000;
  mh1000_pl100000.open("MultiplicityHits50ADC3ns_ht120_Signal_MH-1000_MFF-450_CTau-100000mm.txt"); //MultiplicityHits50ADC3ns_ht120_Signal_mh1000_mx450_pl100k_.txt");
  n=0;
  while (mh1000_pl100000 >> signal_mh1000_pl100000[n]) n++;
  mh1000_pl100000.close();

  /*
  // mh=350 pl=0.5m 
  double signal_mh350_pl500[6];
  ifstream mh350_pl500;
  mh350_pl500.open("MultiplicityHits50ADC3ns_ht120_Signal_MH-350_MFF-160_CTau-500mm_Tun.txt"); //MultiplicityHits50ADC3ns_ht120_Signal_mh350__mx160_pl500__.txt");
  n=0;
  while (mh350_pl500 >> signal_mh350_pl500[n]) n++;
  mh350_pl500.close();
  */

  // mh=350 pl=1m 
  double signal_mh350_pl1000[6];
  ifstream mh350_pl1000;
  mh350_pl1000.open("MultiplicityHits50ADC3ns_ht120_Signal_MH-350_MFF-160_CTau-1000mm_Tu.txt"); //MultiplicityHits50ADC3ns_ht120_Signal_mh350__mx160_pl1000_.txt");
  n=0;
  while (mh350_pl1000 >> signal_mh350_pl1000[n]) n++;
  mh350_pl1000.close();
  // mh=350 pl=10m      
  double signal_mh350_pl10000[6];
  ifstream mh350_pl10000;
  mh350_pl10000.open("MultiplicityHits50ADC3ns_ht120_Signal_MH-350_MFF-160_CTau-10000mm_T.txt"); //MultiplicityHits50ADC3ns_ht120_Signal_mh350__mx160_pl10000.txt");
  n=0;
  while (mh350_pl10000 >> signal_mh350_pl10000[n]) n++;
  mh350_pl10000.close();

  /*
  // mh=250 pl=0.5m   
  double signal_mh250_pl500[6];
  ifstream mh250_pl500;
  mh250_pl500.open("MultiplicityHits50ADC3ns_ht120_Signal_MH-250_MFF-120_CTau-500mm_Tun.txt"); //MultiplicityHits50ADC3ns_ht120_Signal_mh250__mx160_pl500__.txt");
  n=0;
  while (mh250_pl500 >> signal_mh250_pl500[n]) n++;
  mh250_pl500.close();
  */
  // mh=250 pl=1m 
  double signal_mh250_pl1000[6];
  ifstream mh250_pl1000;
  mh250_pl1000.open("MultiplicityHits50ADC3ns_ht120_Signal_MH-250_MFF-120_CTau-1000mm_Tu.txt"); //MultiplicityHits50ADC3ns_ht120_Signal_mh250__mx120_pl1000_.txt");
  n=0;
  while (mh250_pl1000 >> signal_mh250_pl1000[n]) n++;
  mh250_pl1000.close();
  /*
  // mh=250 pl=10m
  double signal_mh250_pl10000[6];
  ifstream mh250_pl10000;
  mh250_pl10000.open("MultiplicityHits50ADC3ns_ht120_Signal_MH-250_MFF-120_CTau-10000mm_Tun.txt"); //MultiplicityHits50ADC3ns_ht120_Signal_mh250__mx160_pl10000.txt");
  n=0;
  while (mh250_pl10000 >> signal_mh250_pl10000[n]) n++;
  mh250_pl10000.close();
  */

  // Background
  double background[6];
  ifstream QCD;
  QCD.open("MultiplicityHits50ADC3ns_ht120_Background.txt");
  n=0;
  while (QCD >> background[n]) n++;
  QCD.close();

  // Neutrino gun for rates
  double rates[7];
  ifstream nugun;
  nugun.open("NuGunRates.txt");
  n=0;
  while (nugun >> rates[n]) n++;
  nugun.close();

  double nugun_rate[5], cuts_background[5];
  double cuts_mh125_pl3000[5];//, cuts_mh125_pl30000[5];
  double cuts_mh1000_pl10000[5];//, cuts_mh1000_pl100000[5];
  double cuts_mh250_pl1000[5];//cuts_mh250_pl500[5], cuts_mh250_pl1000[5], cuts_mh250_pl10000[5];
  double cuts_mh350_pl10000[5];//cuts_mh350_pl500[5], cuts_mh350_pl1000[5], cuts_mh350_pl10000[5];
  for (int i=0; i<5; i++) {
    nugun_rate[i] = rates[i];
    cuts_mh125_pl3000[i] = signal_mh125_pl3000[i];
    //    cuts_mh125_pl30000[i] = signal_mh125_pl30000[i];
    //    cuts_mh250_pl500[i] = signal_mh250_pl500[i];
    cuts_mh250_pl1000[i] = signal_mh250_pl1000[i];
    //    cuts_mh250_pl10000[i] = signal_mh250_pl10000[i];
    //    cuts_mh350_pl500[i] = signal_mh350_pl500[i];
    //    cuts_mh350_pl1000[i] = signal_mh350_pl1000[i];
    cuts_mh350_pl10000[i] = signal_mh350_pl10000[i];
    cuts_mh1000_pl10000[i] = signal_mh1000_pl10000[i];
    //    cuts_mh1000_pl100000[i] = signal_mh1000_pl100000[i];
    cuts_background[i] = background[i];
  }

  // mh = 125 GeV
  TGraph *gr_LLP_mh125_pl3000 = new TGraph(5, cuts_mh125_pl3000, nugun_rate);
  TMarker *m_mh125_pl3000_ht360 = new TMarker(signal_mh125_pl3000[5], rates[5], 21);
  //  TGraph *gr_LLP_mh125_pl30000 = new TGraph(5, cuts_mh125_pl30000, nugun_rate);
  //  TMarker *m_mh125_pl30000_ht360 = new TMarker(signal_mh125_pl30000[5], rates[5], 21);
  // mh = 250 GeV
  //  TGraph *gr_LLP_mh250_pl500 = new TGraph(5, cuts_mh250_pl500, nugun_rate);
  //  TMarker *m_mh250_pl500_ht360 = new TMarker(signal_mh250_pl500[5], rates[5], 21);
  TGraph *gr_LLP_mh250_pl1000 = new TGraph(5, cuts_mh250_pl1000, nugun_rate);
  TMarker *m_mh250_pl1000_ht360 = new TMarker(signal_mh250_pl1000[5], rates[5], 21);
  //  TGraph *gr_LLP_mh250_pl10000 = new TGraph(5, cuts_mh250_pl10000, nugun_rate);
  //  TMarker *m_mh250_pl10000_ht360 = new TMarker(signal_mh250_pl10000[5], rates[5], 21);
  // mh = 350 GeV
  //  TGraph *gr_LLP_mh350_pl500 = new TGraph(5, cuts_mh350_pl500, nugun_rate);
  //  TMarker *m_mh350_pl500_ht360 = new TMarker(signal_mh350_pl500[5], rates[5], 21);
  //  TGraph *gr_LLP_mh350_pl1000 = new TGraph(5, cuts_mh350_pl1000, nugun_rate);
  //  TMarker *m_mh350_pl1000_ht360 = new TMarker(signal_mh350_pl1000[5], rates[5], 21);
  TGraph *gr_LLP_mh350_pl10000 = new TGraph(5, cuts_mh350_pl10000, nugun_rate);
  TMarker *m_mh350_pl10000_ht360 = new TMarker(signal_mh350_pl10000[5], rates[5], 21);
  // mh = 1000 GeV 
  TGraph *gr_LLP_mh1000_pl10000 = new TGraph(5, cuts_mh1000_pl10000, nugun_rate);
  TMarker *m_mh1000_pl10000_ht360 = new TMarker(signal_mh1000_pl10000[5], rates[5], 21);
  //  TGraph *gr_LLP_mh1000_pl100000 = new TGraph(5, cuts_mh1000_pl100000, nugun_rate);
  //  TMarker *m_mh1000_pl100000_ht360 = new TMarker(signal_mh1000_pl100000[5], rates[5], 21);
  // QCD
  TGraph *gr_background = new TGraph(5, cuts_background, nugun_rate);
  TMarker *m_background_ht360 = new TMarker(background[5], rates[5], 21);
  // comparison to the rate at 120 with no timing cuts  
  TLine *l=new TLine(0.,rates[6],1.,rates[6]);

  // mh = 125 GeV
  TCanvas *c1_LLP_pl3000 = new TCanvas("c1_LLP_pl3000","Graph Draw Options",200,10,600,400);
  gr_LLP_mh125_pl3000->GetHistogram()->SetMinimum(-5.);
  gr_LLP_mh125_pl3000->Draw("AC*");
  gr_LLP_mh125_pl3000->SetTitle("Rate vs. Signal Efficiency for >=1,2,3,4,5 Cells >=50ADC,3ns near 4 L1 Jets;LLP Efficiency;Neutrino Gun Rate (Hz, unnormalized)   ");
  gr_LLP_mh125_pl3000->GetXaxis()->SetLimits(0.,1.);
  gr_LLP_mh125_pl3000->GetHistogram()->SetMinimum(10000.);
  gr_LLP_mh125_pl3000->GetHistogram()->SetMaximum(100000000.);

  m_mh125_pl3000_ht360->SetMarkerStyle(21);
  m_mh125_pl3000_ht360->SetMarkerColor(kRed-9);
  m_mh125_pl3000_ht360->Draw();
  gr_LLP_mh125_pl3000->SetLineColor(kRed);
  gr_LLP_mh125_pl3000->Draw("C*");

  m_mh250_pl1000_ht360->SetMarkerStyle(21);
  m_mh250_pl1000_ht360->SetMarkerColor(kGreen-9);
  m_mh250_pl1000_ht360->Draw();
  gr_LLP_mh250_pl1000->SetLineColor(kGreen);
  gr_LLP_mh250_pl1000->Draw("C*");

  m_mh350_pl10000_ht360->SetMarkerStyle(21);
  m_mh350_pl10000_ht360->SetMarkerColor(kBlue-9);
  m_mh350_pl10000_ht360->Draw();
  gr_LLP_mh350_pl10000->SetLineColor(kBlue);
  gr_LLP_mh350_pl10000->Draw("C*");

  m_mh1000_pl10000_ht360->SetMarkerStyle(21);
  m_mh1000_pl10000_ht360->SetMarkerColor(kMagenta-9);
  m_mh1000_pl10000_ht360->Draw();
  gr_LLP_mh1000_pl10000->SetLineColor(kMagenta);
  gr_LLP_mh1000_pl10000->Draw("C*");
  l->SetLineColor(kBlack);
  l->SetLineStyle(7);
  l->Draw();
  auto legend1_htSum = new TLegend(0.55,0.15,0.9,0.45);
  legend1_htSum->AddEntry(m_mh125_pl3000_ht360,"m_{H}=125, c#scale[1.2]{#tau}=3m; H_{T}>360 GeV, no timing cuts");
  legend1_htSum->AddEntry(gr_LLP_mh125_pl3000,"m_{H}=125, c#scale[1.2]{#tau}=3m; H_{T}>120 GeV, with timing cuts"); 
  legend1_htSum->AddEntry(m_mh250_pl1000_ht360,"m_{H}=250, c#scale[1.2]{#tau}=1m; H_{T}>360 GeV, no timing cuts");
  legend1_htSum->AddEntry(gr_LLP_mh250_pl1000,"m_{H}=250, c#scale[1.2]{#tau}=1m; H_{T}>120 GeV, with timing cuts");
  legend1_htSum->AddEntry(m_mh350_pl10000_ht360,"m_{H}=350, c#scale[1.2]{#tau}=10m; H_{T}>360 GeV, no timing cuts");
  legend1_htSum->AddEntry(gr_LLP_mh350_pl10000,"m_{H}=350, c#scale[1.2]{#tau}=10m; H_{T}>120 GeV, with timing cuts");
  legend1_htSum->AddEntry(m_mh1000_pl10000_ht360,"m_{H}=1000, c#scale[1.2]{#tau}=10m; H_{T}>360 GeV, no timing cuts");
  legend1_htSum->AddEntry(gr_LLP_mh1000_pl10000,"m_{H}=1000, c#scale[1.2]{#tau}=10m; H_{T}>120 GeV, with timing cuts");
  legend1_htSum->AddEntry(l,"Neutrino gun rate at HT=120GeV with no timing cuts");
  legend1_htSum->Draw();
  c1_LLP_pl3000->SetLogy();
  c1_LLP_pl3000->SetGrid();
  c1_LLP_pl3000->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/EffRate_LLP_overlay_mh125_250_350_1000_3_1_10.pdf");

  /*
  // mh = 125 GeV, higher lifetime
  TCanvas *c1_LLP_pl30000 = new TCanvas("c1_LLP_pl30000","Graph Draw Options",200,10,600,400);
  gr_LLP_mh125_pl30000->GetHistogram()->SetMinimum(-5.);
  gr_LLP_mh125_pl30000->Draw("AC*");
  gr_LLP_mh125_pl30000->SetTitle("Rate vs. Signal Efficiency for >=1,2,3,4,5 Cells >=50ADC,3ns near 4 L1 Jets;LLP Efficiency;Neutrino Gun Rate (Hz, unnormalized)   ");
  gr_LLP_mh125_pl30000->GetXaxis()->SetLimits(0.,1.);
  gr_LLP_mh125_pl30000->GetHistogram()->SetMinimum(10000.);
  gr_LLP_mh125_pl30000->GetHistogram()->SetMaximum(100000000.);

  m_mh125_pl30000_ht360->SetMarkerStyle(21);
  m_mh125_pl30000_ht360->SetMarkerColor(kRed-9);
  m_mh125_pl30000_ht360->Draw();
  gr_LLP_mh125_pl30000->SetLineColor(kRed);
  gr_LLP_mh125_pl30000->Draw("C*");
  m_mh1000_pl100000_ht360->SetMarkerStyle(21);
  m_mh1000_pl100000_ht360->SetMarkerColor(kBlue-9);
  m_mh1000_pl100000_ht360->Draw();
  gr_LLP_mh1000_pl100000->SetLineColor(kBlue);
  gr_LLP_mh1000_pl100000->Draw("C*");
  auto legend2_htSum = new TLegend(0.55,0.15,0.9,0.45);
  legend2_htSum->AddEntry(m_mh125_pl30000_ht360,"m_{H}=125, c#scale[1.2]{#tau}=30m; H_{T}>360 GeV, no timing cuts");
  legend2_htSum->AddEntry(gr_LLP_mh125_pl30000,"m_{H}=125, c#scale[1.2]{#tau}=30m; H_{T}>120 GeV, with timing cuts");
  legend2_htSum->AddEntry(m_mh1000_pl100000_ht360,"m_{H}=1000, c#scale[1.2]{#tau}=100m; H_{T}>360 GeV, no timing cuts");
  legend2_htSum->AddEntry(gr_LLP_mh1000_pl100000,"m_{H}=1000, c#scale[1.2]{#tau}=100m; H_{T}>120 GeV, with timing cuts");
  legend2_htSum->Draw();
  c1_LLP_pl30000->SetLogy();
  c1_LLP_pl30000->SetGrid();
  c1_LLP_pl30000->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/EffRate_LLP_overlay_mh125_mh1000_30_100.pdf");

  // mh = 250, 350
  TCanvas *c1_LLP_pl500 = new TCanvas("c1_LLP_pl500","Graph Draw Options",200,10,600,400);
  gr_LLP_mh250_pl500->GetHistogram()->SetMinimum(-5.);
  gr_LLP_mh250_pl500->Draw("AC*");
  gr_LLP_mh250_pl500->SetTitle("Rate vs. Signal Efficiency for >=1,2,3,4,5 Cells >=50ADC,3ns near 4 L1 Jets;LLP Efficiency, c#scale[1.2]{#tau}=0.5m;Neutrino Gun Rate (Hz, unnormalized)   ");
  gr_LLP_mh250_pl500->GetXaxis()->SetLimits(0.,1.);
  gr_LLP_mh250_pl500->GetHistogram()->SetMinimum(10000.);
  gr_LLP_mh250_pl500->GetHistogram()->SetMaximum(100000000.);

  m_mh250_pl500_ht360->SetMarkerStyle(21);
  m_mh250_pl500_ht360->SetMarkerColor(kRed-9);
  m_mh250_pl500_ht360->Draw();
  gr_LLP_mh250_pl500->SetLineColor(kRed);
  gr_LLP_mh250_pl500->Draw("C*");
  m_mh350_pl500_ht360->SetMarkerStyle(21);
  m_mh350_pl500_ht360->SetMarkerColor(kBlue-9);
  m_mh350_pl500_ht360->Draw();
  gr_LLP_mh350_pl500->SetLineColor(kBlue);
  gr_LLP_mh350_pl500->Draw("C*");
  auto legend3_htSum = new TLegend(0.55,0.15,0.9,0.45);
  legend3_htSum->AddEntry(m_mh250_pl500_ht360,"m_{H}=250; H_{T}>360 GeV, no timing cuts");
  legend3_htSum->AddEntry(gr_LLP_mh250_pl500,"m_{H}=250; H_{T}>120 GeV, with timing cuts");
  legend3_htSum->AddEntry(m_mh350_pl500_ht360,"m_{H}=350; H_{T}>360 GeV, no timing cuts");
  legend3_htSum->AddEntry(gr_LLP_mh350_pl500,"m_{H}=350; H_{T}>120 GeV, with timing cuts");
  legend3_htSum->Draw();
  c1_LLP_pl500->SetLogy();
  c1_LLP_pl500->SetGrid();
  c1_LLP_pl500->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/EffRate_LLP_pl500.pdf");

  TCanvas *c1_LLP_pl1000 = new TCanvas("c1_LLP_pl1000","Graph Draw Options",200,10,600,400);
  gr_LLP_mh250_pl1000->GetHistogram()->SetMinimum(-5.);
  gr_LLP_mh250_pl1000->Draw("AC*");
  gr_LLP_mh250_pl1000->SetTitle("Rate vs. Signal Efficiency for >=1,2,3,4,5 Cells >=50ADC,3ns near 4 L1 Jets;LLP Efficiency, c#scale[1.2]{#tau}=1m;Neutrino Gun Rate (Hz, unnormalized)   ");
  gr_LLP_mh250_pl1000->GetXaxis()->SetLimits(0.,1.);
  gr_LLP_mh250_pl1000->GetHistogram()->SetMinimum(10000.);
  gr_LLP_mh250_pl1000->GetHistogram()->SetMaximum(100000000.);

  m_mh250_pl1000_ht360->SetMarkerStyle(21);
  m_mh250_pl1000_ht360->SetMarkerColor(kRed-9);
  m_mh250_pl1000_ht360->Draw();
  gr_LLP_mh250_pl1000->SetLineColor(kRed);
  gr_LLP_mh250_pl1000->Draw("C*");
  m_mh350_pl1000_ht360->SetMarkerStyle(21);
  m_mh350_pl1000_ht360->SetMarkerColor(kBlue-9);
  m_mh350_pl1000_ht360->Draw();
  gr_LLP_mh350_pl1000->SetLineColor(kBlue);
  gr_LLP_mh350_pl1000->Draw("C*");
  auto legend4_htSum = new TLegend(0.55,0.15,0.9,0.45);
  legend4_htSum->AddEntry(m_mh250_pl1000_ht360,"m_{H}=250; H_{T}>360 GeV, no timing cuts");
  legend4_htSum->AddEntry(gr_LLP_mh250_pl1000,"m_{H}=250; H_{T}>120 GeV, with timing cuts");
  legend4_htSum->AddEntry(m_mh350_pl1000_ht360,"m_{H}=350; H_{T}>360 GeV, no timing cuts");
  legend4_htSum->AddEntry(gr_LLP_mh350_pl1000,"m_{H}=350; H_{T}>120 GeV, with timing cuts");
  legend4_htSum->Draw();
  c1_LLP_pl1000->SetLogy();
  c1_LLP_pl1000->SetGrid();
  c1_LLP_pl1000->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/EffRate_LLP_pl1000.pdf");

  TCanvas *c1_LLP_pl10000 = new TCanvas("c1_LLP_pl10000","Graph Draw Options",200,10,600,400);
  gr_LLP_mh250_pl10000->GetHistogram()->SetMinimum(-5.);
  gr_LLP_mh250_pl10000->Draw("AC*");
  gr_LLP_mh250_pl10000->SetTitle("Rate vs. Signal Efficiency for >=1,2,3,4,5 Cells >=50ADC,3ns near 4 L1 Jets;LLP Efficiency, c#scale[1.2]{#tau}=10m;Neutrino Gun Rate (Hz, unnormalized)   ");
  gr_LLP_mh250_pl10000->GetXaxis()->SetLimits(0.,1.);
  gr_LLP_mh250_pl10000->GetHistogram()->SetMinimum(10000.);
  gr_LLP_mh250_pl10000->GetHistogram()->SetMaximum(100000000.);

  m_mh250_pl10000_ht360->SetMarkerStyle(21);
  m_mh250_pl10000_ht360->SetMarkerColor(kRed-9);
  m_mh250_pl10000_ht360->Draw();
  gr_LLP_mh250_pl10000->SetLineColor(kRed);
  gr_LLP_mh250_pl10000->Draw("C*");
  m_mh350_pl10000_ht360->SetMarkerStyle(21);
  m_mh350_pl10000_ht360->SetMarkerColor(kBlue-9);
  m_mh350_pl10000_ht360->Draw();
  gr_LLP_mh350_pl10000->SetLineColor(kBlue);
  gr_LLP_mh350_pl10000->Draw("C*");
  auto legend5_htSum = new TLegend(0.55,0.15,0.9,0.45);
  legend5_htSum->AddEntry(m_mh250_pl10000_ht360,"m_{H}=250; H_{T}>360 GeV, no timing cuts");
  legend5_htSum->AddEntry(gr_LLP_mh250_pl10000,"m_{H}=250; H_{T}>120 GeV, with timing cuts");
  legend5_htSum->AddEntry(m_mh350_pl10000_ht360,"m_{H}=350; H_{T}>360 GeV, no timing cuts");
  legend5_htSum->AddEntry(gr_LLP_mh350_pl10000,"m_{H}=350; H_{T}>120 GeV, with timing cuts");
  legend5_htSum->Draw();
  c1_LLP_pl10000->SetLogy();
  c1_LLP_pl10000->SetGrid();
  c1_LLP_pl10000->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/EffRate_LLP_pl10000.pdf");
  */

  // background
  TCanvas *c1_background = new TCanvas("c1_background","Graph Draw Options",200,10,600,400);
  gr_background->GetHistogram()->SetMinimum(-5.);
  gr_background->Draw("AC*");
  gr_background->SetLineColor(kBlack);
  gr_background->SetTitle("htSum Rate vs. Background Efficiency for >=1,2,3,4,5 Cells >=50ADC,3ns near 4 L1 Jets;QCD Efficiency;Neutrino Gun Rate (Hz, unnormalized)   ");
  gr_background->GetXaxis()->SetLimits(0.,1.);
  gr_background->GetHistogram()->SetMinimum(10000.);
  gr_background->GetHistogram()->SetMaximum(100000000.);

  m_background_ht360->SetMarkerStyle(21);
  m_background_ht360->SetMarkerColor(kGray);
  m_background_ht360->Draw();
  //  m_background_ht120->SetMarkerStyle(21);
  //  m_background_ht120->SetMarkerColor(kBlack);
  //  m_background_ht120->Draw();
  l->SetLineColor(kBlack);
  l->SetLineStyle(7);
  l->Draw();
  auto legend9_htSum = new TLegend(0.55,0.15,0.9,0.35);
  legend9_htSum->AddEntry(m_background_ht360,"H_{T}>360 GeV, no timing cuts");
  legend9_htSum->AddEntry(gr_background,"H_{T}>120 GeV, with timing cuts");        
  legend9_htSum->AddEntry(l,"Neutrino gun rate at HT=120GeV with no timing cuts");
  //  legend9_htSum->AddEntry(m_background_ht120,"H_{T}>120 GeV, with timing cuts");
  legend9_htSum->Draw();
  c1_background->SetLogy();
  c1_background->SetGrid();
  c1_background->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/EffRate_QCDbackground_grey.pdf");
}
