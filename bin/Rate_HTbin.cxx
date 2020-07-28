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
  double signal_mh125_pl3000[14];
  ifstream mh125_pl3000;
  mh125_pl3000.open("Efficiency_HtBins_Signal_MH-125_MFF-50_CTau-3000mm_Tun.txt"); //Efficiency_HtBins_Signal_mh125__mx50__pl500__.txt");
  int n=0;
  while (mh125_pl3000 >> signal_mh125_pl3000[n]) n++;
  mh125_pl3000.close();
  // mh=250 pl=1m
  double signal_mh250_pl1000[14];
  ifstream mh250_pl1000;
  mh250_pl1000.open("Efficiency_HtBins_Signal_MH-250_MFF-120_CTau-1000mm_Tu.txt"); //Efficiency_HtBins_Signal_mh125__mx50__pl1000_.txt");
  n=0;
  while (mh250_pl1000 >> signal_mh250_pl1000[n]) n++;
  mh250_pl1000.close();
  // mh=350 pl=1m
  double signal_mh350_pl1000[14];
  ifstream mh350_pl1000;
  mh350_pl1000.open("Efficiency_HtBins_Signal_MH-350_MFF-160_CTau-1000mm_Tu.txt");
  n=0;
  while (mh350_pl1000 >> signal_mh350_pl1000[n]) n++;
  mh350_pl1000.close();
  // mh=1000 pl=10m
  double signal_mh1000_pl10000[14];
  ifstream mh1000_pl10000;
  mh1000_pl10000.open("Efficiency_HtBins_Signal_MH-1000_MFF-450_CTau-10000mm_.txt");
  n=0;
  while (mh1000_pl10000 >> signal_mh1000_pl10000[n]) n++;
  mh1000_pl10000.close();

  // ht360
  double L1HT360[14] = {0,0,0,0,0,0,0,0,0,0,0,0,1,1};
  // Middle of HT bin used for determining efficiencies
  double htBin[14] = {130, 150, 170, 190, 210, 230, 250, 270, 290, 310, 330, 350, 370, 390};
  double htBin_er[14] = {10,10,10,10,10,10,10,10,10,10,10,10,10,10};
  double err_y[14]={0};
  
  TGraphErrors *gr_L1HT360 = new TGraphErrors(14, htBin, L1HT360, htBin_er,err_y);

  TGraphErrors *gr_LLP_mh250_pl1000 = new TGraphErrors(14, htBin, signal_mh250_pl1000,htBin_er,err_y);
  TGraphErrors *gr_LLP_mh125_pl3000 = new TGraphErrors(14, htBin, signal_mh125_pl3000,htBin_er,err_y);
  TGraphErrors *gr_LLP_mh350_pl1000 = new TGraphErrors(14, htBin, signal_mh350_pl1000,htBin_er,err_y);
  TGraphErrors *gr_LLP_mh1000_pl10000 = new TGraphErrors(14, htBin, signal_mh1000_pl10000,htBin_er,err_y);
  for (int i=0; i<14; i++) std::cout << signal_mh1000_pl10000[i] << " and " << htBin[i] << std::endl;

  TCanvas *c1_LLP_mh125_pl3000 = new TCanvas("c1_LLP_mh125_pl3000","Graph Draw Options",200,10,600,400);
  gr_LLP_mh125_pl3000->SetMarkerColor(kBlue);
  gr_LLP_mh125_pl3000->SetMarkerStyle(21);
  gr_LLP_mh125_pl3000->Draw("AP");
  gr_L1HT360->SetMarkerColor(kRed);
  gr_L1HT360->SetMarkerStyle(21);
  //  gr_L1HT360->Draw("P");
  gr_LLP_mh125_pl3000->SetTitle("HT vs. Signal Efficiency for >=2 Cells >=50ADC,3ns near 4 L1 Jets;HT (GeV);LLP Efficiency, mh=125 GeV, c#scale[1.2]{#tau}=3m   ");
  gr_LLP_mh125_pl3000->GetHistogram()->SetMinimum(0.);
  gr_LLP_mh125_pl3000->GetHistogram()->SetMaximum(1.);
  c1_LLP_mh125_pl3000->SetGrid();
  auto legend_htSum = new TLegend(0.15,0.65,0.45,0.8);
  legend_htSum->AddEntry(gr_LLP_mh125_pl3000,"HT>120 with timing cuts");
  //  legend_htSum->AddEntry(gr_L1HT360,"HT>360 at L1");
  legend_htSum->Draw();
  c1_LLP_mh125_pl3000->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/EffHT_LLP_mh125_pl3000.pdf");

  TCanvas *c1_LLP_mh250_pl1000 = new TCanvas("c1_LLP_mh250_pl1000","Graph Draw Options",200,10,600,400);
  gr_LLP_mh250_pl1000->SetMarkerColor(kBlue);
  gr_LLP_mh250_pl1000->SetMarkerStyle(21);
  gr_LLP_mh250_pl1000->Draw("AP");
  gr_L1HT360->SetMarkerColor(kRed);
  gr_L1HT360->SetMarkerStyle(21);
  gr_LLP_mh250_pl1000->SetTitle("HT vs. Signal Efficiency for >=2 Cells >=50ADC,3ns near 4 L1 Jets;HT (GeV);LLP Efficiency, mh=250 GeV, c#scale[1.2]{#tau}=1m   ");
  gr_LLP_mh250_pl1000->GetHistogram()->SetMinimum(0.);
  gr_LLP_mh250_pl1000->GetHistogram()->SetMaximum(1.);
  c1_LLP_mh250_pl1000->SetGrid();
  auto legend2_htSum = new TLegend(0.15,0.65,0.45,0.8);
  legend2_htSum->AddEntry(gr_LLP_mh250_pl1000,"HT>120 with timing cuts");
  legend2_htSum->Draw();
  c1_LLP_mh250_pl1000->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/EffHT_LLP_mh250_pl1000.pdf");

  TCanvas *c1_LLP_mh350_pl1000 = new TCanvas("c1_LLP_mh350_pl1000","Graph Draw Options",200,10,600,400);
  gr_LLP_mh350_pl1000->SetMarkerColor(kBlue);
  gr_LLP_mh350_pl1000->SetMarkerStyle(21);
  gr_LLP_mh350_pl1000->Draw("AP");
  gr_L1HT360->SetMarkerColor(kRed);
  gr_L1HT360->SetMarkerStyle(21);
  gr_LLP_mh350_pl1000->SetTitle("HT vs. Signal Efficiency for >=2 Cells >=50ADC,3ns near 4 L1 Jets;HT (GeV);LLP Efficiency, mh=350 GeV, c#scale[1.2]{#tau}=1m   ");
  gr_LLP_mh350_pl1000->GetHistogram()->SetMinimum(0.);
  gr_LLP_mh350_pl1000->GetHistogram()->SetMaximum(1.);
  c1_LLP_mh350_pl1000->SetGrid();
  auto legend3_htSum = new TLegend(0.15,0.65,0.45,0.8);
  legend3_htSum->AddEntry(gr_LLP_mh350_pl1000,"HT>120 with timing cuts");
  legend3_htSum->Draw();
  c1_LLP_mh350_pl1000->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/EffHT_LLP_mh350_pl1000.pdf");

  TCanvas *c1_LLP_mh1000_pl10000 = new TCanvas("c1_LLP_mh1000_pl10000","Graph Draw Options",200,10,600,400);
  gr_LLP_mh1000_pl10000->SetMarkerColor(kBlue);
  gr_LLP_mh1000_pl10000->SetMarkerStyle(21);
  gr_LLP_mh1000_pl10000->Draw("AP");
  gr_L1HT360->SetMarkerColor(kRed);
  gr_L1HT360->SetMarkerStyle(21);
  gr_LLP_mh1000_pl10000->SetTitle("HT vs. Signal Efficiency for >=2 Cells >=50ADC,3ns near 4 L1 Jets;HT (GeV);LLP Efficiency, mh=1000 GeV, c#scale[1.2]{#tau}=10m   ");
  gr_LLP_mh1000_pl10000->GetHistogram()->SetMinimum(0.);
  gr_LLP_mh1000_pl10000->GetHistogram()->SetMaximum(1.);
  c1_LLP_mh1000_pl10000->SetGrid();
  auto legend4_htSum = new TLegend(0.15,0.65,0.45,0.8);
  legend4_htSum->AddEntry(gr_LLP_mh1000_pl10000,"HT>120 with timing cuts");
  legend4_htSum->Draw();
  c1_LLP_mh1000_pl10000->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/TimingBit/OfficialProduction/EffHT_LLP_mh1000_pl10000.pdf");

}
