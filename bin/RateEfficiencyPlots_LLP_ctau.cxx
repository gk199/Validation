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
  // mh=125 pl=0.5m
  double signal_mh125_pl500[6];
  ifstream mh125_pl500;
  mh125_pl500.open("MultiplicityHits50ADC3ns_ht120_Signal_mh125__mx50__pl500__.txt");
  int n=0;
  while (mh125_pl500 >> signal_mh125_pl500[n]) n++;
  mh125_pl500.close();
  // mh=125 pl=1m
  double signal_mh125_pl1000[6];
  ifstream mh125_pl1000;
  mh125_pl1000.open("MultiplicityHits50ADC3ns_ht120_Signal_mh125__mx50__pl1000_.txt");
  n=0;
  while (mh125_pl1000 >> signal_mh125_pl1000[n]) n++;
  mh125_pl1000.close();
  // mh=125 pl=10m
  double signal_mh125_pl10000[6];
  ifstream mh125_pl10000;
  mh125_pl10000.open("MultiplicityHits50ADC3ns_ht120_Signal_mh125__mx50__pl10000.txt");
  n=0;
  while (mh125_pl10000 >> signal_mh125_pl10000[n]) n++;
  mh125_pl10000.close();

  // mh=1000 pl=0.5m
  double signal_mh1000_pl500[6];
  ifstream mh1000_pl500;
  mh1000_pl500.open("MultiplicityHits50ADC3ns_ht120_Signal_mh1000_mx450_pl500__.txt");
  n=0;
  while (mh1000_pl500 >> signal_mh1000_pl500[n]) n++;
  mh1000_pl500.close();
  // mh=1000 pl=1m 
  double signal_mh1000_pl1000[6];
  ifstream mh1000_pl1000;
  mh1000_pl1000.open("MultiplicityHits50ADC3ns_ht120_Signal_mh1000_mx450_pl1000_.txt");
  n=0;
  while (mh1000_pl1000 >> signal_mh1000_pl1000[n]) n++;
  mh1000_pl1000.close();
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

  // mh=350 pl=0.5m 
  double signal_mh350_pl500[6];
  ifstream mh350_pl500;
  mh350_pl500.open("MultiplicityHits50ADC3ns_ht120_Signal_mh350__mx160_pl500__.txt");
  n=0;
  while (mh350_pl500 >> signal_mh350_pl500[n]) n++;
  mh350_pl500.close();
  // mh=350 pl=1m
  double signal_mh350_pl1000[6];
  ifstream mh350_pl1000;
  mh350_pl1000.open("MultiplicityHits50ADC3ns_ht120_Signal_mh350__mx160_pl1000_.txt");
  n=0;
  while (mh350_pl1000 >> signal_mh350_pl1000[n]) n++;
  mh350_pl1000.close();
  // mh=350 pl=10m
  double signal_mh350_pl10000[6];
  ifstream mh350_pl10000;
  mh350_pl10000.open("MultiplicityHits50ADC3ns_ht120_Signal_mh350__mx160_pl10000.txt");
  n=0;
  while (mh350_pl10000 >> signal_mh350_pl10000[n]) n++;
  mh350_pl10000.close();

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

  double nugun_rate[4], cuts_mh125_pl1000[4], cuts_mh1000_pl1000[4], cuts_mh250_pl1000[4], cuts_mh350_pl1000[4], cuts_background[4];
  double cuts_mh125_pl500[4], cuts_mh1000_pl500[4], cuts_mh350_pl500[4];
  double cuts_mh125_pl10000[4], cuts_mh1000_pl10000[4], cuts_mh350_pl10000[4];
  for (int i=0; i<4; i++) {
    nugun_rate[i] = rates[i];
    cuts_mh125_pl1000[i] = signal_mh125_pl1000[i];
    cuts_mh1000_pl1000[i] = signal_mh1000_pl1000[i];
    cuts_mh250_pl1000[i] = signal_mh250_pl1000[i];
    cuts_mh350_pl1000[i] = signal_mh350_pl1000[i];
    cuts_mh125_pl500[i] = signal_mh125_pl500[i];
    cuts_mh1000_pl500[i] = signal_mh1000_pl500[i];
    cuts_mh350_pl500[i] = signal_mh350_pl500[i];
    cuts_mh125_pl10000[i] = signal_mh125_pl10000[i];
    cuts_mh1000_pl10000[i] = signal_mh1000_pl10000[i];
    cuts_mh350_pl10000[i] = signal_mh350_pl10000[i];
    cuts_background[i] = background[i];
  }

  // 1m
  // mh = 125 GeV
  TGraph *gr_LLP_mh125_pl1000 = new TGraph(4, cuts_mh125_pl1000, nugun_rate);
  TMarker *m_mh125_pl1000_ht360 = new TMarker(signal_mh125_pl1000[5], rates[5], 21);
  //  TMarker *m_mh125_pl1000_ht120 = new TMarker(signal_mh125_pl1000[0], rates[0], 21);
  // mh = 1000 GeV
  TGraph *gr_LLP_mh1000_pl1000 = new TGraph(4, cuts_mh1000_pl1000, nugun_rate);
  TMarker *m_mh1000_pl1000_ht360 = new TMarker(signal_mh1000_pl1000[5], rates[5], 21);
  // mh = 250 GeV
  TGraph *gr_LLP_mh250_pl1000 = new TGraph(4, cuts_mh250_pl1000, nugun_rate);
  TMarker *m_mh250_pl1000_ht360 = new TMarker(signal_mh250_pl1000[5], rates[5], 21);
  // mh = 350 GeV
  TGraph *gr_LLP_mh350_pl1000 = new TGraph(4, cuts_mh350_pl1000, nugun_rate);
  TMarker *m_mh350_pl1000_ht360 = new TMarker(signal_mh350_pl1000[5], rates[5], 21);

  // 0.5m            
  // mh = 125 GeV  
  TGraph *gr_LLP_mh125_pl500 = new TGraph(4, cuts_mh125_pl500, nugun_rate);
  TMarker *m_mh125_pl500_ht360 = new TMarker(signal_mh125_pl500[5], rates[5], 21);
  // mh = 1000 GeV 
  TGraph *gr_LLP_mh1000_pl500 = new TGraph(4, cuts_mh1000_pl500, nugun_rate);
  TMarker *m_mh1000_pl500_ht360 = new TMarker(signal_mh1000_pl500[5], rates[5], 21);
  // mh = 350 GeV  
  TGraph *gr_LLP_mh350_pl500 = new TGraph(4, cuts_mh350_pl500, nugun_rate);
  TMarker *m_mh350_pl500_ht360 = new TMarker(signal_mh350_pl500[5], rates[5], 21);

  // 10m            
  // mh = 125 GeV  
  TGraph *gr_LLP_mh125_pl10000 = new TGraph(4, cuts_mh125_pl10000, rates);
  TMarker *m_mh125_pl10000_ht360 = new TMarker(signal_mh125_pl10000[5], rates[5], 21);
  // mh = 1000 GeV 
  TGraph *gr_LLP_mh1000_pl10000 = new TGraph(4, cuts_mh1000_pl10000, nugun_rate);
  TMarker *m_mh1000_pl10000_ht360 = new TMarker(signal_mh1000_pl10000[5], rates[5], 21);
  // mh = 350 GeV
  TGraph *gr_LLP_mh350_pl10000 = new TGraph(4, cuts_mh350_pl10000, nugun_rate);
  TMarker *m_mh350_pl10000_ht360 = new TMarker(signal_mh350_pl10000[5], rates[5], 21);

  // QCD
  TGraph *gr_background = new TGraph(4, cuts_background, nugun_rate);
  TMarker *m_background_ht360 = new TMarker(background[5], rates[5], 21);
  // comparison to the rate at 120 with no timing cuts
  TLine *l=new TLine(0.,rates[6],1.,rates[6]);


  // 1m
  TCanvas *c1_LLP_pl1000 = new TCanvas("c1_LLP_pl1000","Graph Draw Options",200,10,600,400);
  gr_LLP_mh125_pl1000->Draw("AL*");
  //  gr_LLP_mh125_pl1000->SetTitle("Rate vs. Signal Efficiency for >=1,2,3,4 Cells >=50ADC,3ns near 4 L1 Jets;Added LLP Efficiency, c#scale[1.2]{#tau}=1m;Neutrino Gun Rate (Hz, unnormalized)   ");
  gr_LLP_mh125_pl1000->SetTitle("Rate vs. Signal Efficiency for >=1,2 Delayed Jets in HB;Added LLP Efficiency, c#scale[1.2]{#tau}=1m;Neutrino Gun Rate (Hz, unnormalized)   ");
  gr_LLP_mh125_pl1000->GetXaxis()->SetLimits(0.,1.);
  gr_LLP_mh125_pl1000->GetHistogram()->SetMinimum(1.);
  gr_LLP_mh125_pl1000->GetHistogram()->SetMaximum(100000000.);

  m_mh125_pl1000_ht360->SetMarkerStyle(21);
  m_mh125_pl1000_ht360->SetMarkerColor(kRed-9);
  m_mh125_pl1000_ht360->Draw();
  gr_LLP_mh125_pl1000->SetLineColor(kRed);
  gr_LLP_mh125_pl1000->Draw("L*");
  m_mh250_pl1000_ht360->SetMarkerStyle(21);
  m_mh250_pl1000_ht360->SetMarkerColor(kGreen-9);
  m_mh250_pl1000_ht360->Draw();
  gr_LLP_mh250_pl1000->SetLineColor(kGreen);
  gr_LLP_mh250_pl1000->Draw("L*");
  m_mh350_pl1000_ht360->SetMarkerStyle(21);
  m_mh350_pl1000_ht360->SetMarkerColor(kBlue-9);
  m_mh350_pl1000_ht360->Draw();
  gr_LLP_mh350_pl1000->SetLineColor(kBlue);
  gr_LLP_mh350_pl1000->Draw("L*");
  m_mh1000_pl1000_ht360->SetMarkerStyle(21);
  m_mh1000_pl1000_ht360->SetMarkerColor(kMagenta-9);
  m_mh1000_pl1000_ht360->Draw();
  gr_LLP_mh1000_pl1000->SetLineColor(kMagenta);
  gr_LLP_mh1000_pl1000->Draw("L*");
  l->SetLineColor(kBlack);
  l->SetLineStyle(7);
  l->Draw();
  auto legend2_htSum = new TLegend(0.35,0.15,0.9,0.45);
  legend2_htSum->AddEntry(m_mh125_pl1000_ht360,"m_{H}=125; H_{T}>360 GeV, no timing cuts");
  legend2_htSum->AddEntry(gr_LLP_mh125_pl1000,"m_{H}=125; H_{T}>120 GeV, with timing cuts"); 
  legend2_htSum->AddEntry(m_mh250_pl1000_ht360,"m_{H}=250; H_{T}>360 GeV, no timing cuts");
  legend2_htSum->AddEntry(gr_LLP_mh250_pl1000,"m_{H}=250; H_{T}>120 GeV, with timing cuts"); 
  legend2_htSum->AddEntry(m_mh350_pl1000_ht360,"m_{H}=350; H_{T}>360 GeV, no timing cuts");
  legend2_htSum->AddEntry(gr_LLP_mh350_pl1000,"m_{H}=350; H_{T}>120 GeV, with timing cuts");  
  legend2_htSum->AddEntry(m_mh1000_pl1000_ht360,"m_{H}=1000; H_{T}>360 GeV, no timing cuts");
  legend2_htSum->AddEntry(gr_LLP_mh1000_pl1000,"m_{H}=1000; H_{T}>120 GeV, with timing cuts");
  legend2_htSum->AddEntry(l,"Neutrino gun rate at HT=120GeV with no timing cuts");
  legend2_htSum->Draw();
  c1_LLP_pl1000->SetGrid();
  //  c1_LLP_pl1000->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/TimingBit/EffRate_LLP_pl1000_linear.pdf");
  c1_LLP_pl1000->SetLogy();
  c1_LLP_pl1000->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/TimingBit/EffRate_LLP_pl1000.pdf");

  // 0.5 m
  TCanvas *c1_LLP_pl500 = new TCanvas("c1_LLP_pl500","Graph Draw Options",200,10,600,400);
  gr_LLP_mh125_pl500->Draw("AL*");
  //  gr_LLP_mh125_pl500->SetTitle("Rate vs. Signal Efficiency for >=1,2,3,4 Cells >=50ADC,3ns near 4 L1 Jets;Added LLP Efficiency, c#scale[1.2]{#tau}=0.5m;Neutrino Gun Rate (Hz, unnormalized)   ");
  gr_LLP_mh125_pl500->SetTitle("Rate vs. Signal Efficiency for >=1,2 Delayed Jets in HB;Added LLP Efficiency, c#scale[1.2]{#tau}=0.5m;Neutrino Gun Rate (Hz, unnormalized)   ");
  gr_LLP_mh125_pl500->GetXaxis()->SetLimits(0.,1.);
  gr_LLP_mh125_pl500->GetHistogram()->SetMinimum(1.);
  gr_LLP_mh125_pl500->GetHistogram()->SetMaximum(100000000.);

  m_mh125_pl500_ht360->SetMarkerStyle(21);
  m_mh125_pl500_ht360->SetMarkerColor(kRed-9);
  m_mh125_pl500_ht360->Draw();
  gr_LLP_mh125_pl500->SetLineColor(kRed);
  gr_LLP_mh125_pl500->Draw("L*");
  m_mh350_pl500_ht360->SetMarkerStyle(21);
  m_mh350_pl500_ht360->SetMarkerColor(kBlue-9);
  m_mh350_pl500_ht360->Draw();
  gr_LLP_mh350_pl500->SetLineColor(kBlue);
  gr_LLP_mh350_pl500->Draw("L*");
  m_mh1000_pl500_ht360->SetMarkerStyle(21);
  m_mh1000_pl500_ht360->SetMarkerColor(kMagenta-9);
  m_mh1000_pl500_ht360->Draw();
  gr_LLP_mh1000_pl500->SetLineColor(kMagenta);
  gr_LLP_mh1000_pl500->Draw("L*");
  l->SetLineColor(kBlack);
  l->SetLineStyle(7);
  l->Draw();
  auto legend3_htSum = new TLegend(0.3,0.15,0.9,0.4);
  legend3_htSum->AddEntry(m_mh125_pl500_ht360,"m_{H}=125; H_{T}>360 GeV, no timing cuts");
  legend3_htSum->AddEntry(gr_LLP_mh125_pl500,"m_{H}=125; H_{T}>120 GeV, with timing cuts");
  legend3_htSum->AddEntry(m_mh350_pl500_ht360,"m_{H}=350; H_{T}>360 GeV, no timing cuts");
  legend3_htSum->AddEntry(gr_LLP_mh350_pl500,"m_{H}=350; H_{T}>120 GeV, with timing cuts");
  legend3_htSum->AddEntry(m_mh1000_pl500_ht360,"m_{H}=1000; H_{T}>360 GeV, no timing cuts");
  legend3_htSum->AddEntry(gr_LLP_mh1000_pl500,"m_{H}=1000; H_{T}>120 GeV, with timing cuts");
  legend3_htSum->AddEntry(l,"Neutrino gun rate at HT=120GeV with no timing cuts");
  legend3_htSum->Draw();
  c1_LLP_pl500->SetGrid();
  //  c1_LLP_pl500->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/TimingBit/EffRate_LLP_pl500_linear.pdf");
  c1_LLP_pl500->SetLogy();
  c1_LLP_pl500->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/TimingBit/EffRate_LLP_pl500.pdf");

  // 10m
  TCanvas *c1_LLP_pl10000 = new TCanvas("c1_LLP_pl10000","Graph Draw Options",200,10,600,400);
  gr_LLP_mh125_pl10000->Draw("AL*");
  //  gr_LLP_mh125_pl10000->SetTitle("Rate vs. Signal Efficiency for >=1,2,3,4 Cells >=50ADC,3ns near 4 L1 Jets;Added LLP Efficiency, c#scale[1.2]{#tau}=10m;Neutrino Gun Rate (Hz, unnormalized)   ");
  gr_LLP_mh125_pl10000->SetTitle("Rate vs. Signal Efficiency for >=1,2 Delayed Jets in HB;Added LLP Efficiency, c#scale[1.2]{#tau}=10m;Neutrino Gun Rate (Hz, unnormalized)   ");
  gr_LLP_mh125_pl10000->GetXaxis()->SetLimits(0.,1.);
  gr_LLP_mh125_pl10000->GetHistogram()->SetMinimum(1.);
  gr_LLP_mh125_pl10000->GetHistogram()->SetMaximum(100000000.);

  m_mh125_pl10000_ht360->SetMarkerStyle(21);
  m_mh125_pl10000_ht360->SetMarkerColor(kRed-9);
  m_mh125_pl10000_ht360->Draw();
  gr_LLP_mh125_pl10000->SetLineColor(kRed);
  gr_LLP_mh125_pl10000->Draw("L*");
  m_mh350_pl10000_ht360->SetMarkerStyle(21);
  m_mh350_pl10000_ht360->SetMarkerColor(kBlue-9);
  m_mh350_pl10000_ht360->Draw();
  gr_LLP_mh350_pl10000->SetLineColor(kBlue);
  gr_LLP_mh350_pl10000->Draw("L*");
  m_mh1000_pl10000_ht360->SetMarkerStyle(21);
  m_mh1000_pl10000_ht360->SetMarkerColor(kMagenta-9);
  m_mh1000_pl10000_ht360->Draw();
  gr_LLP_mh1000_pl10000->SetLineColor(kMagenta);
  gr_LLP_mh1000_pl10000->Draw("L*");
  l->SetLineColor(kBlack);
  l->SetLineStyle(7);
  l->Draw();
  auto legend4_htSum = new TLegend(0.3,0.15,0.9,0.4);
  legend4_htSum->AddEntry(m_mh125_pl10000_ht360,"m_{H}=125; H_{T}>360 GeV, no timing cuts");
  legend4_htSum->AddEntry(gr_LLP_mh125_pl10000,"m_{H}=125; H_{T}>120 GeV, with timing cuts");
  legend4_htSum->AddEntry(m_mh350_pl10000_ht360,"m_{H}=350; H_{T}>360 GeV, no timing cuts");
  legend4_htSum->AddEntry(gr_LLP_mh350_pl10000,"m_{H}=350; H_{T}>120 GeV, with timing cuts");
  legend4_htSum->AddEntry(m_mh1000_pl10000_ht360,"m_{H}=1000; H_{T}>360 GeV, no timing cuts");
  legend4_htSum->AddEntry(gr_LLP_mh1000_pl10000,"m_{H}=1000; H_{T}>120 GeV, with timing cuts");
  legend4_htSum->AddEntry(l,"Neutrino gun rate at HT=120GeV with no timing cuts");
  legend4_htSum->Draw();
  c1_LLP_pl10000->SetGrid();
  //  c1_LLP_pl10000->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/TimingBit/EffRate_LLP_pl10000_linear.pdf");
  c1_LLP_pl10000->SetLogy();
  c1_LLP_pl10000->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/TimingBit/EffRate_LLP_pl10000.pdf");

  // background
  TCanvas *c1_background = new TCanvas("c1_background","Graph Draw Options",200,10,600,400);
  gr_background->Draw("AL*");
  gr_background->SetLineColor(kBlack);
  //  gr_background->SetTitle("Rate vs. Background Efficiency for >=1,2,3,4 Cells >=50ADC,3ns near 4 L1 Jets;Added QCD Efficiency;Neutrino Gun Rate (Hz, unnormalized)   ");
  gr_background->SetTitle("Rate vs. Background Efficiency for >=1,2 Delayed Jets in HB;Added QCD Efficiency, c#scale[1.2]{#tau}=0.5m;Neutrino Gun Rate (Hz, unnormalized)   ");
  gr_background->GetXaxis()->SetLimits(0.,1.);
  gr_background->GetHistogram()->SetMinimum(1.);
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
  auto legend9_htSum = new TLegend(0.35,0.15,0.9,0.35);
  legend9_htSum->AddEntry(m_background_ht360,"H_{T}>360 GeV, no timing cuts");
  legend9_htSum->AddEntry(gr_background,"H_{T}>120 GeV, with timing cuts");                                                                                                                                                        
  legend9_htSum->AddEntry(l,"Neutrino gun rate at HT=120GeV with no timing cuts");
  //  legend9_htSum->AddEntry(m_background_ht120,"H_{T}>120 GeV, with timing cuts");
  legend9_htSum->Draw();
  c1_background->SetGrid();
  //  c1_background->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/TimingBit/EffRate_QCDbackground_grey_linear.pdf");
  c1_background->SetLogy();
  c1_background->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/TimingBit/EffRate_QCDbackground_grey.pdf");
}
