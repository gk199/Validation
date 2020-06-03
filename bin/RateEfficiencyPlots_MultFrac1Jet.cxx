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
  double arr_signal_mult_mh1000_pl10000_1GeV[8];
  ifstream Signal_mult_mh1000_pl10000;
  Signal_mult_mh1000_pl10000.open("DelayedHitFrac_ht120_Signal_1GeV_mh1000_pl10000.txt");
  int n = 0;
  while (Signal_mult_mh1000_pl10000 >> arr_signal_mult_mh1000_pl10000_1GeV[n]) n++; // signal efficiency
  Signal_mult_mh1000_pl10000.close();
  double arr_signal_mult_mh1000_pl10000_2GeV[8];
  Signal_mult_mh1000_pl10000.open("DelayedHitFrac_ht120_Signal_2GeV_mh1000_pl10000.txt");
  n = 0;
  while (Signal_mult_mh1000_pl10000 >> arr_signal_mult_mh1000_pl10000_2GeV[n]) n++; // signal efficiency
  Signal_mult_mh1000_pl10000.close();
  double arr_signal_mult_mh1000_pl10000_3GeV[8];
  Signal_mult_mh1000_pl10000.open("DelayedHitFrac_ht120_Signal_3GeV_mh1000_pl10000.txt");
  n = 0;
  while (Signal_mult_mh1000_pl10000 >> arr_signal_mult_mh1000_pl10000_3GeV[n]) n++; // signal efficiency
  Signal_mult_mh1000_pl10000.close();
  double arr_signal_mult_mh1000_pl10000_4GeV[8];
  Signal_mult_mh1000_pl10000.open("DelayedHitFrac_ht120_Signal_4GeV_mh1000_pl10000.txt");
  n = 0;
  while (Signal_mult_mh1000_pl10000 >> arr_signal_mult_mh1000_pl10000_4GeV[n]) n++; // signal efficiency
  Signal_mult_mh1000_pl10000.close();

  double arr_signal_mult_mh1000_pl1000_1GeV[8];
  ifstream Signal_mult_mh1000_pl1000;
  Signal_mult_mh1000_pl1000.open("DelayedHitFrac_ht120_Signal_1GeV_mh1000_pl1000_.txt");
  n = 0;
  while (Signal_mult_mh1000_pl1000 >> arr_signal_mult_mh1000_pl1000_1GeV[n]) n++; // signal efficiency 
  Signal_mult_mh1000_pl1000.close();
  double arr_signal_mult_mh1000_pl1000_2GeV[8];
  Signal_mult_mh1000_pl1000.open("DelayedHitFrac_ht120_Signal_2GeV_mh1000_pl1000_.txt");
  n = 0;
  while (Signal_mult_mh1000_pl1000 >> arr_signal_mult_mh1000_pl1000_2GeV[n]) n++; // signal efficiency
  Signal_mult_mh1000_pl1000.close();
  double arr_signal_mult_mh1000_pl1000_3GeV[8];
  Signal_mult_mh1000_pl1000.open("DelayedHitFrac_ht120_Signal_3GeV_mh1000_pl1000_.txt");
  n = 0;
  while (Signal_mult_mh1000_pl1000 >> arr_signal_mult_mh1000_pl1000_3GeV[n]) n++; // signal efficiency
  Signal_mult_mh1000_pl1000.close();
  double arr_signal_mult_mh1000_pl1000_4GeV[8];
  Signal_mult_mh1000_pl1000.open("DelayedHitFrac_ht120_Signal_4GeV_mh1000_pl1000_.txt");
  n = 0;
  while (Signal_mult_mh1000_pl1000 >> arr_signal_mult_mh1000_pl1000_4GeV[n]) n++; // signal efficiency
  Signal_mult_mh1000_pl1000.close();

  double arr_signal_mult_mh1000_pl500_1GeV[8];
  ifstream Signal_mult_mh1000_pl500;
  Signal_mult_mh1000_pl500.open("DelayedHitFrac_ht120_Signal_1GeV_mh1000_pl500__.txt");
  n = 0;
  while (Signal_mult_mh1000_pl500 >> arr_signal_mult_mh1000_pl500_1GeV[n]) n++; // signal efficiency
  Signal_mult_mh1000_pl500.close();
  double arr_signal_mult_mh1000_pl500_2GeV[8];
  Signal_mult_mh1000_pl500.open("DelayedHitFrac_ht120_Signal_2GeV_mh1000_pl500__.txt");
  n = 0;
  while (Signal_mult_mh1000_pl500 >> arr_signal_mult_mh1000_pl500_2GeV[n]) n++; // signal efficiency
  Signal_mult_mh1000_pl500.close();
  double arr_signal_mult_mh1000_pl500_3GeV[8];
  Signal_mult_mh1000_pl500.open("DelayedHitFrac_ht120_Signal_3GeV_mh1000_pl500__.txt");
  n = 0;
  while (Signal_mult_mh1000_pl500 >> arr_signal_mult_mh1000_pl500_3GeV[n]) n++; // signal efficiency
  Signal_mult_mh1000_pl500.close();
  double arr_signal_mult_mh1000_pl500_4GeV[8];
  Signal_mult_mh1000_pl500.open("DelayedHitFrac_ht120_Signal_4GeV_mh1000_pl500__.txt");
  n = 0;
  while (Signal_mult_mh1000_pl500 >> arr_signal_mult_mh1000_pl500_4GeV[n]) n++; // signal efficiency
  Signal_mult_mh1000_pl500.close();

  // mh=350 GeV, mx=160 GeV
  double arr_signal_mult_mh350_pl10000_1GeV[8];
  ifstream Signal_mult_mh350_pl10000;
  Signal_mult_mh350_pl10000.open("DelayedHitFrac_ht120_Signal_1GeV_mh350__pl10000.txt");
  n = 0;
  while (Signal_mult_mh350_pl10000 >> arr_signal_mult_mh350_pl10000_1GeV[n]) n++;
  Signal_mult_mh350_pl10000.close();
  double arr_signal_mult_mh350_pl10000_2GeV[8];
  Signal_mult_mh350_pl10000.open("DelayedHitFrac_ht120_Signal_2GeV_mh350__pl10000.txt");
  n = 0;
  while (Signal_mult_mh350_pl10000 >> arr_signal_mult_mh350_pl10000_2GeV[n]) n++;
  Signal_mult_mh350_pl10000.close();
  double arr_signal_mult_mh350_pl10000_3GeV[8];
  Signal_mult_mh350_pl10000.open("DelayedHitFrac_ht120_Signal_3GeV_mh350__pl10000.txt");
  n = 0;
  while (Signal_mult_mh350_pl10000 >> arr_signal_mult_mh350_pl10000_3GeV[n]) n++;
  Signal_mult_mh350_pl10000.close();
  double arr_signal_mult_mh350_pl10000_4GeV[8];
  Signal_mult_mh350_pl10000.open("DelayedHitFrac_ht120_Signal_4GeV_mh350__pl10000.txt");
  n = 0;
  while (Signal_mult_mh350_pl10000 >> arr_signal_mult_mh350_pl10000_4GeV[n]) n++;
  Signal_mult_mh350_pl10000.close();

  double arr_signal_mult_mh350_pl1000_1GeV[8];
  ifstream Signal_mult_mh350_pl1000;
  Signal_mult_mh350_pl1000.open("DelayedHitFrac_ht120_Signal_1GeV_mh350__pl1000_.txt");
  n = 0;
  while (Signal_mult_mh350_pl1000 >> arr_signal_mult_mh350_pl1000_1GeV[n]) n++;
  Signal_mult_mh350_pl1000.close();
  double arr_signal_mult_mh350_pl1000_2GeV[8];
  Signal_mult_mh350_pl1000.open("DelayedHitFrac_ht120_Signal_2GeV_mh350__pl1000_.txt");
  n = 0;
  while (Signal_mult_mh350_pl1000 >> arr_signal_mult_mh350_pl1000_2GeV[n]) n++;
  Signal_mult_mh350_pl1000.close();
  double arr_signal_mult_mh350_pl1000_3GeV[8];
  Signal_mult_mh350_pl1000.open("DelayedHitFrac_ht120_Signal_3GeV_mh350__pl1000_.txt");
  n = 0;
  while (Signal_mult_mh350_pl1000 >> arr_signal_mult_mh350_pl1000_3GeV[n]) n++;
  Signal_mult_mh350_pl1000.close();
  double arr_signal_mult_mh350_pl1000_4GeV[8];
  Signal_mult_mh350_pl1000.open("DelayedHitFrac_ht120_Signal_4GeV_mh350__pl1000_.txt");
  n = 0;
  while (Signal_mult_mh350_pl1000 >> arr_signal_mult_mh350_pl1000_4GeV[n]) n++;
  Signal_mult_mh350_pl1000.close();

  double arr_signal_mult_mh350_pl500_1GeV[8];
  ifstream Signal_mult_mh350_pl500;
  Signal_mult_mh350_pl500.open("DelayedHitFrac_ht120_Signal_1GeV_mh350__pl500__.txt");
  n = 0;
  while (Signal_mult_mh350_pl500 >> arr_signal_mult_mh350_pl500_1GeV[n]) n++;
  Signal_mult_mh350_pl500.close();
  double arr_signal_mult_mh350_pl500_2GeV[8];
  Signal_mult_mh350_pl500.open("DelayedHitFrac_ht120_Signal_2GeV_mh350__pl500__.txt");
  n = 0;
  while (Signal_mult_mh350_pl500 >> arr_signal_mult_mh350_pl500_2GeV[n]) n++;
  Signal_mult_mh350_pl500.close();
  double arr_signal_mult_mh350_pl500_3GeV[8];
  Signal_mult_mh350_pl500.open("DelayedHitFrac_ht120_Signal_3GeV_mh350__pl500__.txt");
  n = 0;
  while (Signal_mult_mh350_pl500 >> arr_signal_mult_mh350_pl500_3GeV[n]) n++;
  Signal_mult_mh350_pl500.close();
  double arr_signal_mult_mh350_pl500_4GeV[8];
  Signal_mult_mh350_pl500.open("DelayedHitFrac_ht120_Signal_4GeV_mh350__pl500__.txt");
  n = 0;
  while (Signal_mult_mh350_pl500 >> arr_signal_mult_mh350_pl500_4GeV[n]) n++;
  Signal_mult_mh350_pl500.close();

  // mh=250 GeV, mx=160 GeV                                              
  double arr_signal_mult_mh250_pl10000_1GeV[8];
  ifstream Signal_mult_mh250_pl10000;
  Signal_mult_mh250_pl10000.open("DelayedHitFrac_ht120_Signal_1GeV_mh250__pl10000.txt");
  n = 0;
  while (Signal_mult_mh250_pl10000 >> arr_signal_mult_mh250_pl10000_1GeV[n]) n++;
  Signal_mult_mh250_pl10000.close();
  double arr_signal_mult_mh250_pl10000_2GeV[8];
  Signal_mult_mh250_pl10000.open("DelayedHitFrac_ht120_Signal_2GeV_mh250__pl10000.txt");
  n = 0;
  while (Signal_mult_mh250_pl10000 >> arr_signal_mult_mh250_pl10000_2GeV[n]) n++;
  Signal_mult_mh250_pl10000.close();
  double arr_signal_mult_mh250_pl10000_3GeV[8];
  Signal_mult_mh250_pl10000.open("DelayedHitFrac_ht120_Signal_3GeV_mh250__pl10000.txt");
  n = 0;
  while (Signal_mult_mh250_pl10000 >> arr_signal_mult_mh250_pl10000_3GeV[n]) n++;
  Signal_mult_mh250_pl10000.close();
  double arr_signal_mult_mh250_pl10000_4GeV[8];
  Signal_mult_mh250_pl10000.open("DelayedHitFrac_ht120_Signal_4GeV_mh250__pl10000.txt");
  n = 0;
  while (Signal_mult_mh250_pl10000 >> arr_signal_mult_mh250_pl10000_4GeV[n]) n++;
  Signal_mult_mh250_pl10000.close();

  double arr_signal_mult_mh250_pl1000_1GeV[8];
  ifstream Signal_mult_mh250_pl1000;
  Signal_mult_mh250_pl1000.open("DelayedHitFrac_ht120_Signal_1GeV_mh250__pl1000_.txt");
  n = 0;
  while (Signal_mult_mh250_pl1000 >> arr_signal_mult_mh250_pl1000_1GeV[n]) n++;
  Signal_mult_mh250_pl1000.close();
  double arr_signal_mult_mh250_pl1000_2GeV[8];
  Signal_mult_mh250_pl1000.open("DelayedHitFrac_ht120_Signal_2GeV_mh250__pl1000_.txt");
  n = 0;
  while (Signal_mult_mh250_pl1000 >> arr_signal_mult_mh250_pl1000_2GeV[n]) n++;
  Signal_mult_mh250_pl1000.close();
  double arr_signal_mult_mh250_pl1000_3GeV[8];
  Signal_mult_mh250_pl1000.open("DelayedHitFrac_ht120_Signal_3GeV_mh250__pl1000_.txt");
  n = 0;
  while (Signal_mult_mh250_pl1000 >> arr_signal_mult_mh250_pl1000_3GeV[n]) n++;
  Signal_mult_mh250_pl1000.close();
  double arr_signal_mult_mh250_pl1000_4GeV[8];
  Signal_mult_mh250_pl1000.open("DelayedHitFrac_ht120_Signal_4GeV_mh250__pl1000_.txt");
  n = 0;
  while (Signal_mult_mh250_pl1000 >> arr_signal_mult_mh250_pl1000_4GeV[n]) n++;
  Signal_mult_mh250_pl1000.close();

  double arr_signal_mult_mh250_pl500_1GeV[8];
  ifstream Signal_mult_mh250_pl500;
  Signal_mult_mh250_pl500.open("DelayedHitFrac_ht120_Signal_1GeV_mh250__pl500__.txt");
  n = 0;
  while (Signal_mult_mh250_pl500 >> arr_signal_mult_mh250_pl500_1GeV[n]) n++;
  Signal_mult_mh250_pl500.close();
  double arr_signal_mult_mh250_pl500_2GeV[8];
  Signal_mult_mh250_pl500.open("DelayedHitFrac_ht120_Signal_2GeV_mh250__pl500__.txt");
  n = 0;
  while (Signal_mult_mh250_pl500 >> arr_signal_mult_mh250_pl500_2GeV[n]) n++;
  Signal_mult_mh250_pl500.close();
  double arr_signal_mult_mh250_pl500_3GeV[8];
  Signal_mult_mh250_pl500.open("DelayedHitFrac_ht120_Signal_3GeV_mh250__pl500__.txt");
  n = 0;
  while (Signal_mult_mh250_pl500 >> arr_signal_mult_mh250_pl500_3GeV[n]) n++;
  Signal_mult_mh250_pl500.close();
  double arr_signal_mult_mh250_pl500_4GeV[8];
  Signal_mult_mh250_pl500.open("DelayedHitFrac_ht120_Signal_4GeV_mh250__pl500__.txt");
  n = 0;
  while (Signal_mult_mh250_pl500 >> arr_signal_mult_mh250_pl500_4GeV[n]) n++;
  Signal_mult_mh250_pl500.close();

  // mh=125 GeV, mx=160 GeV 
  double arr_signal_mult_mh125_pl10000_1GeV[8];
  ifstream Signal_mult_mh125_pl10000;
  Signal_mult_mh125_pl10000.open("DelayedHitFrac_ht120_Signal_1GeV_mh125__pl10000.txt");
  n = 0;
  while (Signal_mult_mh125_pl10000 >> arr_signal_mult_mh125_pl10000_1GeV[n]) n++;
  Signal_mult_mh125_pl10000.close();
  double arr_signal_mult_mh125_pl10000_2GeV[8];
  Signal_mult_mh125_pl10000.open("DelayedHitFrac_ht120_Signal_2GeV_mh125__pl10000.txt");
  n = 0;
  while (Signal_mult_mh125_pl10000 >> arr_signal_mult_mh125_pl10000_2GeV[n]) n++;
  Signal_mult_mh125_pl10000.close();
  double arr_signal_mult_mh125_pl10000_3GeV[8];
  Signal_mult_mh125_pl10000.open("DelayedHitFrac_ht120_Signal_3GeV_mh125__pl10000.txt");
  n = 0;
  while (Signal_mult_mh125_pl10000 >> arr_signal_mult_mh125_pl10000_3GeV[n]) n++;
  Signal_mult_mh125_pl10000.close();
  double arr_signal_mult_mh125_pl10000_4GeV[8];
  Signal_mult_mh125_pl10000.open("DelayedHitFrac_ht120_Signal_4GeV_mh125__pl10000.txt");
  n = 0;
  while (Signal_mult_mh125_pl10000 >> arr_signal_mult_mh125_pl10000_4GeV[n]) n++;
  Signal_mult_mh125_pl10000.close();

  double arr_signal_mult_mh125_pl1000_1GeV[8];
  ifstream Signal_mult_mh125_pl1000;
  Signal_mult_mh125_pl1000.open("DelayedHitFrac_ht120_Signal_1GeV_mh125__pl1000_.txt");
  n = 0;
  while (Signal_mult_mh125_pl1000 >> arr_signal_mult_mh125_pl1000_1GeV[n]) n++;
  Signal_mult_mh125_pl1000.close();
  double arr_signal_mult_mh125_pl1000_2GeV[8];
  Signal_mult_mh125_pl1000.open("DelayedHitFrac_ht120_Signal_2GeV_mh125__pl1000_.txt");
  n = 0;
  while (Signal_mult_mh125_pl1000 >> arr_signal_mult_mh125_pl1000_2GeV[n]) n++;
  Signal_mult_mh125_pl1000.close();
  double arr_signal_mult_mh125_pl1000_3GeV[8];
  Signal_mult_mh125_pl1000.open("DelayedHitFrac_ht120_Signal_3GeV_mh125__pl1000_.txt");
  n = 0;
  while (Signal_mult_mh125_pl1000 >> arr_signal_mult_mh125_pl1000_3GeV[n]) n++;
  Signal_mult_mh125_pl1000.close();
  double arr_signal_mult_mh125_pl1000_4GeV[8];
  Signal_mult_mh125_pl1000.open("DelayedHitFrac_ht120_Signal_4GeV_mh125__pl1000_.txt");
  n = 0;
  while (Signal_mult_mh125_pl1000 >> arr_signal_mult_mh125_pl1000_4GeV[n]) n++;
  Signal_mult_mh125_pl1000.close();

  double arr_signal_mult_mh125_pl500_1GeV[8];
  ifstream Signal_mult_mh125_pl500;
  Signal_mult_mh125_pl500.open("DelayedHitFrac_ht120_Signal_1GeV_mh125__pl500__.txt");
  n = 0;
  while (Signal_mult_mh125_pl500 >> arr_signal_mult_mh125_pl500_1GeV[n]) n++;
  Signal_mult_mh125_pl500.close();
  double arr_signal_mult_mh125_pl500_2GeV[8];
  Signal_mult_mh125_pl500.open("DelayedHitFrac_ht120_Signal_2GeV_mh125__pl500__.txt");
  n = 0;
  while (Signal_mult_mh125_pl500 >> arr_signal_mult_mh125_pl500_2GeV[n]) n++;
  Signal_mult_mh125_pl500.close();
  double arr_signal_mult_mh125_pl500_3GeV[8];
  Signal_mult_mh125_pl500.open("DelayedHitFrac_ht120_Signal_3GeV_mh125__pl500__.txt");
  n = 0;
  while (Signal_mult_mh125_pl500 >> arr_signal_mult_mh125_pl500_3GeV[n]) n++;
  Signal_mult_mh125_pl500.close();
  double arr_signal_mult_mh125_pl500_4GeV[8];
  Signal_mult_mh125_pl500.open("DelayedHitFrac_ht120_Signal_4GeV_mh125__pl500__.txt");
  n = 0;
  while (Signal_mult_mh125_pl500 >> arr_signal_mult_mh125_pl500_4GeV[n]) n++;
  Signal_mult_mh125_pl500.close();

  // QCD efficiencies
  double arr_background_mult_1GeV[8];
  ifstream Background_mult;
  Background_mult.open("DelayedHitFrac_ht120_Background_1GeV.txt");
  n = 0;
  while (Background_mult >>arr_background_mult_1GeV[n]) n++; // background efficiency                                                             
  Background_mult.close();
  double arr_background_mult_2GeV[8];
  Background_mult.open("DelayedHitFrac_ht120_Background_2GeV.txt");
  n = 0;
  while (Background_mult >>arr_background_mult_2GeV[n]) n++; // background efficiency
  Background_mult.close();
  double arr_background_mult_3GeV[8];
  Background_mult.open("DelayedHitFrac_ht120_Background_3GeV.txt");
  n = 0;
  while (Background_mult >>arr_background_mult_3GeV[n]) n++; // background efficiency
  Background_mult.close();
  double arr_background_mult_4GeV[8];
  Background_mult.open("DelayedHitFrac_ht120_Background_4GeV.txt");
  n = 0;
  while (Background_mult >>arr_background_mult_4GeV[n]) n++; // background efficiency
  Background_mult.close();

  // read in neutrino gun rates at the various multiplicity cuts
  double arr_neutrino_1GeV[22];
  ifstream Neutrino_htSum120_rate;
  Neutrino_htSum120_rate.open("htSum120_Rate_1jet_1GeV.txt");
  n = 0;
  while (Neutrino_htSum120_rate >> arr_neutrino_1GeV[n]) n++;
  Neutrino_htSum120_rate.close();
  double arr_neutrino_2GeV[22];
  Neutrino_htSum120_rate.open("htSum120_Rate_1jet_2GeV.txt");
  n = 0;
  while (Neutrino_htSum120_rate >> arr_neutrino_2GeV[n]) n++;
  Neutrino_htSum120_rate.close();
  double arr_neutrino_3GeV[22];
  Neutrino_htSum120_rate.open("htSum120_Rate_1jet_3GeV.txt");
  n = 0;
  while (Neutrino_htSum120_rate >> arr_neutrino_3GeV[n]) n++;
  Neutrino_htSum120_rate.close();
  double arr_neutrino_4GeV[22];
  Neutrino_htSum120_rate.open("htSum120_Rate_1jet_4GeV.txt");
  n = 0;
  while (Neutrino_htSum120_rate >> arr_neutrino_4GeV[n]) n++;
  Neutrino_htSum120_rate.close();

  // split up neutrino gun into single jet rate, quad jet rate, and htSum120 GeV rate
  //  double htSum120Rate[7];
  double htSum120Rate_1GeV[7], htSum120Rate_2GeV[7], htSum120Rate_3GeV[7], htSum120Rate_4GeV[7];
  double signal_mult_mh1000_pl10000_1GeV[7], signal_mult_mh1000_pl10000_2GeV[7], signal_mult_mh1000_pl10000_3GeV[7], signal_mult_mh1000_pl10000_4GeV[7];
  double signal_mult_mh1000_pl1000_1GeV[7], signal_mult_mh1000_pl1000_2GeV[7], signal_mult_mh1000_pl1000_3GeV[7], signal_mult_mh1000_pl1000_4GeV[7];
  double signal_mult_mh1000_pl500_1GeV[7], signal_mult_mh1000_pl500_2GeV[7], signal_mult_mh1000_pl500_3GeV[7], signal_mult_mh1000_pl500_4GeV[7];
  double signal_mult_mh350_pl10000_1GeV[7], signal_mult_mh350_pl10000_2GeV[7], signal_mult_mh350_pl10000_3GeV[7], signal_mult_mh350_pl10000_4GeV[7];
  double signal_mult_mh350_pl1000_1GeV[7], signal_mult_mh350_pl1000_2GeV[7], signal_mult_mh350_pl1000_3GeV[7], signal_mult_mh350_pl1000_4GeV[7];
  double signal_mult_mh350_pl500_1GeV[7], signal_mult_mh350_pl500_2GeV[7], signal_mult_mh350_pl500_3GeV[7], signal_mult_mh350_pl500_4GeV[7];
  double signal_mult_mh250_pl10000_1GeV[7], signal_mult_mh250_pl10000_2GeV[7], signal_mult_mh250_pl10000_3GeV[7], signal_mult_mh250_pl10000_4GeV[7];
  double signal_mult_mh250_pl1000_1GeV[7], signal_mult_mh250_pl1000_2GeV[7], signal_mult_mh250_pl1000_3GeV[7], signal_mult_mh250_pl1000_4GeV[7];
  double signal_mult_mh250_pl500_1GeV[7], signal_mult_mh250_pl500_2GeV[7], signal_mult_mh250_pl500_3GeV[7], signal_mult_mh250_pl500_4GeV[7];
  double signal_mult_mh125_pl10000_1GeV[7], signal_mult_mh125_pl10000_2GeV[7], signal_mult_mh125_pl10000_3GeV[7], signal_mult_mh125_pl10000_4GeV[7];
  double signal_mult_mh125_pl1000_1GeV[7], signal_mult_mh125_pl1000_2GeV[7], signal_mult_mh125_pl1000_3GeV[7], signal_mult_mh125_pl1000_4GeV[7];
  double signal_mult_mh125_pl500_1GeV[7], signal_mult_mh125_pl500_2GeV[7], signal_mult_mh125_pl500_3GeV[7], signal_mult_mh125_pl500_4GeV[7];
  double background_mult_1GeV[7], background_mult_2GeV[7], background_mult_3GeV[7], background_mult_4GeV[7];
  double Original_htSumRate430, EffPl10000Mh1000_htSum430_Global, EffPl1000Mh1000_htSum430_Global, EffPl500Mh1000_htSum430_Global, EffPl10000Mh350_htSum430_Global, EffPl1000Mh350_htSum430_Global, EffPl500Mh350_htSum430_Global, EffPl10000Mh250_htSum430_Global, EffPl1000Mh250_htSum430_Global, EffPl500Mh250_htSum430_Global, EffPl10000Mh125_htSum430_Global, EffPl1000Mh125_htSum430_Global, EffPl500Mh125_htSum430_Global, EffQCD_htSum430_Global;
  for (int i=0; i<7; i++) {
    //    htSum120Rate[i] = arr_neutrino_2GeV[i];
    htSum120Rate_1GeV[i] = arr_neutrino_1GeV[i];
    htSum120Rate_2GeV[i] = arr_neutrino_2GeV[i];
    htSum120Rate_3GeV[i] = arr_neutrino_3GeV[i];
    htSum120Rate_4GeV[i] = arr_neutrino_4GeV[i];
    signal_mult_mh1000_pl10000_1GeV[i] = arr_signal_mult_mh1000_pl10000_1GeV[i];
    signal_mult_mh1000_pl10000_2GeV[i] = arr_signal_mult_mh1000_pl10000_2GeV[i];
    signal_mult_mh1000_pl10000_3GeV[i] = arr_signal_mult_mh1000_pl10000_3GeV[i];
    signal_mult_mh1000_pl10000_4GeV[i] = arr_signal_mult_mh1000_pl10000_4GeV[i];
    signal_mult_mh1000_pl1000_1GeV[i] = arr_signal_mult_mh1000_pl1000_1GeV[i];
    signal_mult_mh1000_pl1000_2GeV[i] = arr_signal_mult_mh1000_pl1000_2GeV[i];
    signal_mult_mh1000_pl1000_3GeV[i] = arr_signal_mult_mh1000_pl1000_3GeV[i];
    signal_mult_mh1000_pl1000_4GeV[i] = arr_signal_mult_mh1000_pl1000_4GeV[i];
    signal_mult_mh1000_pl500_1GeV[i] = arr_signal_mult_mh1000_pl500_1GeV[i];
    signal_mult_mh1000_pl500_2GeV[i] = arr_signal_mult_mh1000_pl500_2GeV[i];
    signal_mult_mh1000_pl500_3GeV[i] = arr_signal_mult_mh1000_pl500_3GeV[i];
    signal_mult_mh1000_pl500_4GeV[i] = arr_signal_mult_mh1000_pl500_4GeV[i];

    signal_mult_mh350_pl10000_1GeV[i] = arr_signal_mult_mh350_pl10000_1GeV[i];
    signal_mult_mh350_pl10000_2GeV[i] = arr_signal_mult_mh350_pl10000_2GeV[i];
    signal_mult_mh350_pl10000_3GeV[i] = arr_signal_mult_mh350_pl10000_3GeV[i];
    signal_mult_mh350_pl10000_4GeV[i] = arr_signal_mult_mh350_pl10000_4GeV[i];
    signal_mult_mh350_pl1000_1GeV[i] = arr_signal_mult_mh350_pl1000_1GeV[i];
    signal_mult_mh350_pl1000_2GeV[i] = arr_signal_mult_mh350_pl1000_2GeV[i];
    signal_mult_mh350_pl1000_3GeV[i] = arr_signal_mult_mh350_pl1000_3GeV[i];
    signal_mult_mh350_pl1000_4GeV[i] = arr_signal_mult_mh350_pl1000_4GeV[i];
    signal_mult_mh350_pl500_1GeV[i] = arr_signal_mult_mh350_pl500_1GeV[i];
    signal_mult_mh350_pl500_2GeV[i] = arr_signal_mult_mh350_pl500_2GeV[i];
    signal_mult_mh350_pl500_3GeV[i] = arr_signal_mult_mh350_pl500_3GeV[i];
    signal_mult_mh350_pl500_4GeV[i] = arr_signal_mult_mh350_pl500_4GeV[i];

    signal_mult_mh250_pl10000_1GeV[i] = arr_signal_mult_mh250_pl10000_1GeV[i];
    signal_mult_mh250_pl10000_2GeV[i] = arr_signal_mult_mh250_pl10000_2GeV[i];
    signal_mult_mh250_pl10000_3GeV[i] = arr_signal_mult_mh250_pl10000_3GeV[i];
    signal_mult_mh250_pl10000_4GeV[i] = arr_signal_mult_mh250_pl10000_4GeV[i];
    signal_mult_mh250_pl1000_1GeV[i] = arr_signal_mult_mh250_pl1000_1GeV[i];
    signal_mult_mh250_pl1000_2GeV[i] = arr_signal_mult_mh250_pl1000_2GeV[i];
    signal_mult_mh250_pl1000_3GeV[i] = arr_signal_mult_mh250_pl1000_3GeV[i];
    signal_mult_mh250_pl1000_4GeV[i] = arr_signal_mult_mh250_pl1000_4GeV[i];
    signal_mult_mh250_pl500_1GeV[i] = arr_signal_mult_mh250_pl500_1GeV[i];
    signal_mult_mh250_pl500_2GeV[i] = arr_signal_mult_mh250_pl500_2GeV[i];
    signal_mult_mh250_pl500_3GeV[i] = arr_signal_mult_mh250_pl500_3GeV[i];
    signal_mult_mh250_pl500_4GeV[i] = arr_signal_mult_mh250_pl500_4GeV[i];

    signal_mult_mh125_pl10000_1GeV[i] = arr_signal_mult_mh125_pl10000_1GeV[i];
    signal_mult_mh125_pl10000_2GeV[i] = arr_signal_mult_mh125_pl10000_2GeV[i];
    signal_mult_mh125_pl10000_3GeV[i] = arr_signal_mult_mh125_pl10000_3GeV[i];
    signal_mult_mh125_pl10000_4GeV[i] = arr_signal_mult_mh125_pl10000_4GeV[i];
    signal_mult_mh125_pl1000_1GeV[i] = arr_signal_mult_mh125_pl1000_1GeV[i];
    signal_mult_mh125_pl1000_2GeV[i] = arr_signal_mult_mh125_pl1000_2GeV[i];
    signal_mult_mh125_pl1000_3GeV[i] = arr_signal_mult_mh125_pl1000_3GeV[i];
    signal_mult_mh125_pl1000_4GeV[i] = arr_signal_mult_mh125_pl1000_4GeV[i];
    signal_mult_mh125_pl500_1GeV[i] = arr_signal_mult_mh125_pl500_1GeV[i];
    signal_mult_mh125_pl500_2GeV[i] = arr_signal_mult_mh125_pl500_2GeV[i];
    signal_mult_mh125_pl500_3GeV[i] = arr_signal_mult_mh125_pl500_3GeV[i];
    signal_mult_mh125_pl500_4GeV[i] = arr_signal_mult_mh125_pl500_4GeV[i];

    background_mult_1GeV[i] = arr_background_mult_1GeV[i];
    background_mult_2GeV[i] = arr_background_mult_2GeV[i];
    background_mult_3GeV[i] = arr_background_mult_3GeV[i];
    background_mult_4GeV[i] = arr_background_mult_4GeV[i];
  }
  // comparison points for htsum rates and efficiencies
  Original_htSumRate430 = arr_neutrino_1GeV[7];
  EffPl10000Mh1000_htSum430_Global = arr_signal_mult_mh1000_pl10000_1GeV[7];
  EffPl1000Mh1000_htSum430_Global = arr_signal_mult_mh1000_pl1000_1GeV[7];
  EffPl500Mh1000_htSum430_Global = arr_signal_mult_mh1000_pl500_1GeV[7];
  EffPl10000Mh350_htSum430_Global = arr_signal_mult_mh350_pl10000_1GeV[7];
  EffPl1000Mh350_htSum430_Global = arr_signal_mult_mh350_pl1000_1GeV[7];
  EffPl500Mh350_htSum430_Global = arr_signal_mult_mh350_pl500_1GeV[7];
  EffPl10000Mh250_htSum430_Global = arr_signal_mult_mh250_pl10000_1GeV[7];
  EffPl1000Mh250_htSum430_Global = arr_signal_mult_mh250_pl1000_1GeV[7];
  EffPl500Mh250_htSum430_Global = arr_signal_mult_mh250_pl500_1GeV[7];
  EffPl10000Mh125_htSum430_Global = arr_signal_mult_mh125_pl10000_1GeV[7];
  EffPl1000Mh125_htSum430_Global = arr_signal_mult_mh125_pl1000_1GeV[7];
  EffPl500Mh125_htSum430_Global = arr_signal_mult_mh125_pl500_1GeV[7];
  EffQCD_htSum430_Global = arr_background_mult_1GeV[7];

  // make TGraphs with the efficiency and neutrino gun rates from the multiplicity based trigger
  TGraph *gr_120_LLP_mh1000_pl500_1GeV = new TGraph (4, signal_mult_mh1000_pl500_1GeV, htSum120Rate_1GeV);
  TGraph *gr_120_LLP_mh1000_pl500_2GeV = new TGraph (4, signal_mult_mh1000_pl500_2GeV, htSum120Rate_2GeV);
  TGraph *gr_120_LLP_mh1000_pl500_3GeV = new TGraph (4, signal_mult_mh1000_pl500_3GeV, htSum120Rate_3GeV);
  TGraph *gr_120_LLP_mh1000_pl500_4GeV = new TGraph (4, signal_mult_mh1000_pl500_4GeV, htSum120Rate_4GeV);
  TGraph *gr_120_LLP_mh1000_pl1000_1GeV = new TGraph (4, signal_mult_mh1000_pl1000_1GeV, htSum120Rate_1GeV);
  TGraph *gr_120_LLP_mh1000_pl1000_2GeV = new TGraph (4, signal_mult_mh1000_pl1000_2GeV, htSum120Rate_2GeV);
  TGraph *gr_120_LLP_mh1000_pl1000_3GeV = new TGraph (4, signal_mult_mh1000_pl1000_3GeV, htSum120Rate_3GeV);
  TGraph *gr_120_LLP_mh1000_pl1000_4GeV = new TGraph (4, signal_mult_mh1000_pl1000_4GeV, htSum120Rate_4GeV);
  TGraph *gr_120_LLP_mh1000_pl10000_1GeV = new TGraph (4, signal_mult_mh1000_pl10000_1GeV, htSum120Rate_1GeV);
  TGraph *gr_120_LLP_mh1000_pl10000_2GeV = new TGraph (4, signal_mult_mh1000_pl10000_2GeV, htSum120Rate_2GeV);
  TGraph *gr_120_LLP_mh1000_pl10000_3GeV = new TGraph (4, signal_mult_mh1000_pl10000_3GeV, htSum120Rate_3GeV);
  TGraph *gr_120_LLP_mh1000_pl10000_4GeV = new TGraph (4, signal_mult_mh1000_pl10000_4GeV, htSum120Rate_4GeV);

  TGraph *gr_120_LLP_mh350_pl500_1GeV = new TGraph (4, signal_mult_mh350_pl500_1GeV, htSum120Rate_1GeV);
  TGraph *gr_120_LLP_mh350_pl500_2GeV = new TGraph (4, signal_mult_mh350_pl500_2GeV, htSum120Rate_2GeV);
  TGraph *gr_120_LLP_mh350_pl500_3GeV = new TGraph (4, signal_mult_mh350_pl500_3GeV, htSum120Rate_3GeV);
  TGraph *gr_120_LLP_mh350_pl500_4GeV = new TGraph (4, signal_mult_mh350_pl500_4GeV, htSum120Rate_4GeV);
  TGraph *gr_120_LLP_mh350_pl1000_1GeV = new TGraph (4, signal_mult_mh350_pl1000_1GeV, htSum120Rate_1GeV);
  TGraph *gr_120_LLP_mh350_pl1000_2GeV = new TGraph (4, signal_mult_mh350_pl1000_2GeV, htSum120Rate_2GeV);
  TGraph *gr_120_LLP_mh350_pl1000_3GeV = new TGraph (4, signal_mult_mh350_pl1000_3GeV, htSum120Rate_3GeV);
  TGraph *gr_120_LLP_mh350_pl1000_4GeV = new TGraph (4, signal_mult_mh350_pl1000_4GeV, htSum120Rate_4GeV);
  TGraph *gr_120_LLP_mh350_pl10000_1GeV = new TGraph (4, signal_mult_mh350_pl10000_1GeV, htSum120Rate_1GeV);
  TGraph *gr_120_LLP_mh350_pl10000_2GeV = new TGraph (4, signal_mult_mh350_pl10000_2GeV, htSum120Rate_2GeV);
  TGraph *gr_120_LLP_mh350_pl10000_3GeV = new TGraph (4, signal_mult_mh350_pl10000_3GeV, htSum120Rate_3GeV);
  TGraph *gr_120_LLP_mh350_pl10000_4GeV = new TGraph (4, signal_mult_mh350_pl10000_4GeV, htSum120Rate_4GeV);

  TGraph *gr_120_LLP_mh250_pl500_1GeV = new TGraph (4, signal_mult_mh250_pl500_1GeV, htSum120Rate_1GeV);
  TGraph *gr_120_LLP_mh250_pl500_2GeV = new TGraph (4, signal_mult_mh250_pl500_2GeV, htSum120Rate_2GeV);
  TGraph *gr_120_LLP_mh250_pl500_3GeV = new TGraph (4, signal_mult_mh250_pl500_3GeV, htSum120Rate_3GeV);
  TGraph *gr_120_LLP_mh250_pl500_4GeV = new TGraph (4, signal_mult_mh250_pl500_4GeV, htSum120Rate_4GeV);
  TGraph *gr_120_LLP_mh250_pl1000_1GeV = new TGraph (4, signal_mult_mh250_pl1000_1GeV, htSum120Rate_1GeV);
  TGraph *gr_120_LLP_mh250_pl1000_2GeV = new TGraph (4, signal_mult_mh250_pl1000_2GeV, htSum120Rate_2GeV);
  TGraph *gr_120_LLP_mh250_pl1000_3GeV = new TGraph (4, signal_mult_mh250_pl1000_3GeV, htSum120Rate_3GeV);
  TGraph *gr_120_LLP_mh250_pl1000_4GeV = new TGraph (4, signal_mult_mh250_pl1000_4GeV, htSum120Rate_4GeV);
  TGraph *gr_120_LLP_mh250_pl10000_1GeV = new TGraph (4, signal_mult_mh250_pl10000_1GeV, htSum120Rate_1GeV);
  TGraph *gr_120_LLP_mh250_pl10000_2GeV = new TGraph (4, signal_mult_mh250_pl10000_2GeV, htSum120Rate_2GeV);
  TGraph *gr_120_LLP_mh250_pl10000_3GeV = new TGraph (4, signal_mult_mh250_pl10000_3GeV, htSum120Rate_3GeV);
  TGraph *gr_120_LLP_mh250_pl10000_4GeV = new TGraph (4, signal_mult_mh250_pl10000_4GeV, htSum120Rate_4GeV);

  TGraph *gr_120_LLP_mh125_pl500_1GeV = new TGraph (4, signal_mult_mh125_pl500_1GeV, htSum120Rate_1GeV);
  TGraph *gr_120_LLP_mh125_pl500_2GeV = new TGraph (4, signal_mult_mh125_pl500_2GeV, htSum120Rate_2GeV);
  TGraph *gr_120_LLP_mh125_pl500_3GeV = new TGraph (4, signal_mult_mh125_pl500_3GeV, htSum120Rate_3GeV);
  TGraph *gr_120_LLP_mh125_pl500_4GeV = new TGraph (4, signal_mult_mh125_pl500_4GeV, htSum120Rate_4GeV);
  TGraph *gr_120_LLP_mh125_pl1000_1GeV = new TGraph (4, signal_mult_mh125_pl1000_1GeV, htSum120Rate_1GeV);
  TGraph *gr_120_LLP_mh125_pl1000_2GeV = new TGraph (4, signal_mult_mh125_pl1000_2GeV, htSum120Rate_2GeV);
  TGraph *gr_120_LLP_mh125_pl1000_3GeV = new TGraph (4, signal_mult_mh125_pl1000_3GeV, htSum120Rate_3GeV);
  TGraph *gr_120_LLP_mh125_pl1000_4GeV = new TGraph (4, signal_mult_mh125_pl1000_4GeV, htSum120Rate_4GeV);
  TGraph *gr_120_LLP_mh125_pl10000_1GeV = new TGraph (4, signal_mult_mh125_pl10000_1GeV, htSum120Rate_1GeV);
  TGraph *gr_120_LLP_mh125_pl10000_2GeV = new TGraph (4, signal_mult_mh125_pl10000_2GeV, htSum120Rate_2GeV);
  TGraph *gr_120_LLP_mh125_pl10000_3GeV = new TGraph (4, signal_mult_mh125_pl10000_3GeV, htSum120Rate_3GeV);
  TGraph *gr_120_LLP_mh125_pl10000_4GeV = new TGraph (4, signal_mult_mh125_pl10000_4GeV, htSum120Rate_4GeV);

  TGraph *gr_120_QCD_1GeV = new TGraph (4, background_mult_1GeV, htSum120Rate_1GeV);
  TGraph *gr_120_QCD_2GeV = new TGraph (4, background_mult_2GeV, htSum120Rate_2GeV);
  TGraph *gr_120_QCD_3GeV = new TGraph (4, background_mult_3GeV, htSum120Rate_3GeV);
  TGraph *gr_120_QCD_4GeV = new TGraph (4, background_mult_4GeV, htSum120Rate_4GeV);

  // ************** Markers for comparison with current benchmark performance ******
  TMarker *m_QCD_htSum430 = new TMarker (EffQCD_htSum430_Global, Original_htSumRate430, 21);
  TMarker *m_LLPMh1000pl500_htSum430 = new TMarker (EffPl500Mh1000_htSum430_Global, Original_htSumRate430, 21);
  TMarker *m_LLPMh1000pl1000_htSum430 = new TMarker (EffPl1000Mh1000_htSum430_Global, Original_htSumRate430, 21);
  TMarker *m_LLPMh1000pl10000_htSum430 = new TMarker (EffPl10000Mh1000_htSum430_Global, Original_htSumRate430, 21);
  TMarker *m_LLPMh350pl500_htSum430 = new TMarker (EffPl500Mh350_htSum430_Global, Original_htSumRate430, 21);
  TMarker *m_LLPMh350pl1000_htSum430 = new TMarker (EffPl1000Mh350_htSum430_Global, Original_htSumRate430, 21);
  TMarker *m_LLPMh350pl10000_htSum430 = new TMarker (EffPl10000Mh350_htSum430_Global, Original_htSumRate430, 21);
  TMarker *m_LLPMh250pl500_htSum430 = new TMarker (EffPl500Mh250_htSum430_Global, Original_htSumRate430, 21);
  TMarker *m_LLPMh250pl1000_htSum430 = new TMarker (EffPl1000Mh250_htSum430_Global, Original_htSumRate430, 21);
  TMarker *m_LLPMh250pl10000_htSum430 = new TMarker (EffPl10000Mh250_htSum430_Global, Original_htSumRate430, 21);
  TMarker *m_LLPMh125pl500_htSum430 = new TMarker (EffPl500Mh125_htSum430_Global, Original_htSumRate430, 21);
  TMarker *m_LLPMh125pl1000_htSum430 = new TMarker (EffPl1000Mh125_htSum430_Global, Original_htSumRate430, 21);
  TMarker *m_LLPMh125pl10000_htSum430 = new TMarker (EffPl10000Mh125_htSum430_Global, Original_htSumRate430, 21);

  // htSum rate, LLP mh=1TeV, ct=0.5m efficiency, regional
  TCanvas *c1_htSum_mh1000_pl500 = new TCanvas("c1_htSum_mh1000_pl500","Graph Draw Options",200,10,600,400);
  gr_120_LLP_mh1000_pl500_1GeV->SetLineColor(4); // blue 
  gr_120_LLP_mh1000_pl500_1GeV->GetHistogram()->SetMinimum(-5.);
  gr_120_LLP_mh1000_pl500_1GeV->Draw("AL*");
  gr_120_LLP_mh1000_pl500_2GeV->SetLineColor(3);
  gr_120_LLP_mh1000_pl500_3GeV->SetLineColor(2);
  gr_120_LLP_mh1000_pl500_4GeV->SetLineColor(1);
  gr_120_LLP_mh1000_pl500_2GeV->Draw("L*");
  gr_120_LLP_mh1000_pl500_3GeV->Draw("L*");
  gr_120_LLP_mh1000_pl500_4GeV->Draw("L*");
  gr_120_LLP_mh1000_pl500_1GeV->SetTitle("htSum Rate vs. Signal Efficiency for One Jet Delayed Hits and Hit Fraction;LLP, mh=1TeV, c#scale[1.2]{#tau}=0.5m Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  m_LLPMh1000pl500_htSum430->SetMarkerStyle(21);
  m_LLPMh1000pl500_htSum430->SetMarkerColor(3);
  m_LLPMh1000pl500_htSum430->Draw();
  auto legend_htSum_mh1000_pl500 = new TLegend(0.45,0.15,0.9,0.35);
  legend_htSum_mh1000_pl500->AddEntry(gr_120_LLP_mh1000_pl500_1GeV,"H_{T}>120 GeV, with timing cuts, 1GeV");
  legend_htSum_mh1000_pl500->AddEntry(gr_120_LLP_mh1000_pl500_2GeV,"H_{T}>120 GeV, with timing cuts, 2GeV");
  legend_htSum_mh1000_pl500->AddEntry(gr_120_LLP_mh1000_pl500_3GeV,"H_{T}>120 GeV, with timing cuts, 3GeV");
  legend_htSum_mh1000_pl500->AddEntry(gr_120_LLP_mh1000_pl500_4GeV,"H_{T}>120 GeV, with timing cuts, 4GeV");
  legend_htSum_mh1000_pl500->AddEntry(m_LLPMh1000pl500_htSum430,"H_{T}>430 GeV, no timing cuts");
  legend_htSum_mh1000_pl500->Draw();
  gr_120_LLP_mh1000_pl500_1GeV->GetXaxis()->SetLimits(0.,1.);
  gr_120_LLP_mh1000_pl500_1GeV->GetHistogram()->SetMinimum(1.);
  gr_120_LLP_mh1000_pl500_1GeV->GetHistogram()->SetMaximum(100000.);
  c1_htSum_mh1000_pl500->SetLogy(); 
  c1_htSum_mh1000_pl500->SetGrid();
  c1_htSum_mh1000_pl500->SaveAs("plots/NuGun_singleJetTrigger_htSumRates_vs_SignalEff_Pl500Mh1000.pdf");
  c1_htSum_mh1000_pl500->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/NuGun_singleJetTrigger_htSumRates_vs_SignalEff_Pl500Mh1000.pdf");
  // ct=1m
  TCanvas *c1_htSum_mh1000_pl1000 = new TCanvas("c1_htSum_mh1000_pl1000","Graph Draw Options",200,10,600,400);
  gr_120_LLP_mh1000_pl1000_1GeV->SetLineColor(4); // blue 
  gr_120_LLP_mh1000_pl1000_1GeV->GetHistogram()->SetMinimum(-5.);
  gr_120_LLP_mh1000_pl1000_1GeV->Draw("AL*");
  gr_120_LLP_mh1000_pl1000_2GeV->SetLineColor(3);
  gr_120_LLP_mh1000_pl1000_3GeV->SetLineColor(2);
  gr_120_LLP_mh1000_pl1000_4GeV->SetLineColor(1);
  gr_120_LLP_mh1000_pl1000_2GeV->Draw("L*");
  gr_120_LLP_mh1000_pl1000_3GeV->Draw("L*");
  gr_120_LLP_mh1000_pl1000_4GeV->Draw("L*");
  gr_120_LLP_mh1000_pl1000_1GeV->SetTitle("htSum Rate vs. Signal Efficiency for One Jet Delayed Hits and Hit Fraction;LLP, mh=1TeV, c#scale[1.2]{#tau}=1m Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  m_LLPMh1000pl1000_htSum430->SetMarkerStyle(21);
  m_LLPMh1000pl1000_htSum430->SetMarkerColor(3);
  m_LLPMh1000pl1000_htSum430->Draw();
  auto legend_htSum_mh1000_pl1000 = new TLegend(0.45,0.15,0.9,0.35);
  legend_htSum_mh1000_pl1000->AddEntry(gr_120_LLP_mh1000_pl1000_1GeV,"H_{T}>120 GeV, with timing cuts, 1GeV");
  legend_htSum_mh1000_pl1000->AddEntry(gr_120_LLP_mh1000_pl1000_2GeV,"H_{T}>120 GeV, with timing cuts, 2GeV");
  legend_htSum_mh1000_pl1000->AddEntry(gr_120_LLP_mh1000_pl1000_3GeV,"H_{T}>120 GeV, with timing cuts, 3GeV");
  legend_htSum_mh1000_pl1000->AddEntry(gr_120_LLP_mh1000_pl1000_4GeV,"H_{T}>120 GeV, with timing cuts, 4GeV");
  legend_htSum_mh1000_pl1000->AddEntry(m_LLPMh1000pl1000_htSum430,"H_{T}>430 GeV, no timing cuts");
  legend_htSum_mh1000_pl1000->Draw();
  gr_120_LLP_mh1000_pl1000_1GeV->GetXaxis()->SetLimits(0.,1.);
  gr_120_LLP_mh1000_pl1000_1GeV->GetHistogram()->SetMinimum(1.);
  gr_120_LLP_mh1000_pl1000_1GeV->GetHistogram()->SetMaximum(100000.);
  c1_htSum_mh1000_pl1000->SetLogy();
  c1_htSum_mh1000_pl1000->SetGrid();
  c1_htSum_mh1000_pl1000->SaveAs("plots/NuGun_singleJetTrigger_htSumRates_vs_SignalEff_Pl1000Mh1000.pdf");
  c1_htSum_mh1000_pl1000->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/NuGun_singleJetTrigger_htSumRates_vs_SignalEff_Pl1000Mh1000.pdf");
  // ct = 10m
  TCanvas *c1_htSum_mh1000_pl10000 = new TCanvas("c1_htSum_mh1000_pl10000","Graph Draw Options",200,10,600,400);
  gr_120_LLP_mh1000_pl10000_1GeV->SetLineColor(4); // blue
  gr_120_LLP_mh1000_pl10000_1GeV->GetHistogram()->SetMinimum(-5.);
  gr_120_LLP_mh1000_pl10000_1GeV->Draw("AL*");
  gr_120_LLP_mh1000_pl10000_2GeV->SetLineColor(3); 
  gr_120_LLP_mh1000_pl10000_3GeV->SetLineColor(2);
  gr_120_LLP_mh1000_pl10000_4GeV->SetLineColor(1);
  gr_120_LLP_mh1000_pl10000_2GeV->Draw("L*");
  gr_120_LLP_mh1000_pl10000_3GeV->Draw("L*");
  gr_120_LLP_mh1000_pl10000_4GeV->Draw("L*");
  gr_120_LLP_mh1000_pl10000_1GeV->SetTitle("htSum Rate vs. Signal Efficiency for One Jet Delayed Hits and Hit Fraction;LLP, mh=1TeV, c#scale[1.2]{#tau}=10m Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  m_LLPMh1000pl10000_htSum430->SetMarkerStyle(21);
  m_LLPMh1000pl10000_htSum430->SetMarkerColor(3);
  m_LLPMh1000pl10000_htSum430->Draw();
  auto legend_htSum_mh1000_pl10000 = new TLegend(0.45,0.15,0.9,0.35);
  legend_htSum_mh1000_pl10000->AddEntry(gr_120_LLP_mh1000_pl10000_1GeV,"H_{T}>120 GeV, with timing cuts, 1GeV");
  legend_htSum_mh1000_pl10000->AddEntry(gr_120_LLP_mh1000_pl10000_2GeV,"H_{T}>120 GeV, with timing cuts, 2GeV");
  legend_htSum_mh1000_pl10000->AddEntry(gr_120_LLP_mh1000_pl10000_3GeV,"H_{T}>120 GeV, with timing cuts, 3GeV");
  legend_htSum_mh1000_pl10000->AddEntry(gr_120_LLP_mh1000_pl10000_4GeV,"H_{T}>120 GeV, with timing cuts, 4GeV");
  legend_htSum_mh1000_pl10000->AddEntry(m_LLPMh1000pl10000_htSum430,"H_{T}>430 GeV, no timing cuts");
  legend_htSum_mh1000_pl10000->Draw();
  gr_120_LLP_mh1000_pl10000_1GeV->GetXaxis()->SetLimits(0.,1.);
  gr_120_LLP_mh1000_pl10000_1GeV->GetHistogram()->SetMinimum(1.);
  gr_120_LLP_mh1000_pl10000_1GeV->GetHistogram()->SetMaximum(100000.);
  c1_htSum_mh1000_pl10000->SetLogy();
  c1_htSum_mh1000_pl10000->SetGrid();
  c1_htSum_mh1000_pl10000->SaveAs("plots/NuGun_singleJetTrigger_htSumRates_vs_SignalEff_Pl10000Mh1000.pdf");
  c1_htSum_mh1000_pl10000->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/NuGun_singleJetTrigger_htSumRates_vs_SignalEff_Pl10000Mh1000.pdf");

  // htSum rate, LLP mh=350GeV, ct=0.5m efficiency, regional
  TCanvas *c1_htSum_mh350_pl500 = new TCanvas("c1_htSum_mh350_pl500","Graph Draw Options",200,10,600,400);
  gr_120_LLP_mh350_pl500_1GeV->SetLineColor(4); // blue
  gr_120_LLP_mh350_pl500_1GeV->GetHistogram()->SetMinimum(-5.);
  gr_120_LLP_mh350_pl500_1GeV->Draw("AL*");
  gr_120_LLP_mh350_pl500_2GeV->SetLineColor(3);
  gr_120_LLP_mh350_pl500_3GeV->SetLineColor(2);
  gr_120_LLP_mh350_pl500_4GeV->SetLineColor(1);
  gr_120_LLP_mh350_pl500_2GeV->Draw("L*");
  gr_120_LLP_mh350_pl500_3GeV->Draw("L*");
  gr_120_LLP_mh350_pl500_4GeV->Draw("L*");
  gr_120_LLP_mh350_pl500_1GeV->SetTitle("htSum Rate vs. Signal Efficiency for One Jet Delayed Hits and Hit Fraction;LLP, mh=350GeV, c#scale[1.2]{#tau}=0.5m Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  m_LLPMh350pl500_htSum430->SetMarkerStyle(21);
  m_LLPMh350pl500_htSum430->SetMarkerColor(3);
  m_LLPMh350pl500_htSum430->Draw();
  auto legend_htSum_mh350_pl500 = new TLegend(0.45,0.15,0.9,0.35);
  legend_htSum_mh350_pl500->AddEntry(gr_120_LLP_mh350_pl500_1GeV,"H_{T}>120 GeV, with timing cuts, 1GeV");
  legend_htSum_mh350_pl500->AddEntry(gr_120_LLP_mh350_pl500_2GeV,"H_{T}>120 GeV, with timing cuts, 2GeV");
  legend_htSum_mh350_pl500->AddEntry(gr_120_LLP_mh350_pl500_3GeV,"H_{T}>120 GeV, with timing cuts, 3GeV");
  legend_htSum_mh350_pl500->AddEntry(gr_120_LLP_mh350_pl500_4GeV,"H_{T}>120 GeV, with timing cuts, 4GeV");
  legend_htSum_mh350_pl500->AddEntry(m_LLPMh350pl500_htSum430,"H_{T}>430 GeV, no timing cuts");
  legend_htSum_mh350_pl500->Draw();
  gr_120_LLP_mh350_pl500_1GeV->GetXaxis()->SetLimits(0.,1.);
  gr_120_LLP_mh350_pl500_1GeV->GetHistogram()->SetMinimum(1.);
  gr_120_LLP_mh350_pl500_1GeV->GetHistogram()->SetMaximum(100000.);
  c1_htSum_mh350_pl500->SetLogy();
  c1_htSum_mh350_pl500->SetGrid();
  c1_htSum_mh350_pl500->SaveAs("plots/NuGun_singleJetTrigger_htSumRates_vs_SignalEff_Pl500Mh350.pdf");
  c1_htSum_mh350_pl500->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/NuGun_singleJetTrigger_htSumRates_vs_SignalEff_Pl500Mh350.pdf");
  // ct=1m
  TCanvas *c1_htSum_mh350_pl1000 = new TCanvas("c1_htSum_mh350_pl1000","Graph Draw Options",200,10,600,400);
  gr_120_LLP_mh350_pl1000_1GeV->SetLineColor(4); // blue 
  gr_120_LLP_mh350_pl1000_1GeV->GetHistogram()->SetMinimum(-5.);
  gr_120_LLP_mh350_pl1000_1GeV->Draw("AL*");
  gr_120_LLP_mh350_pl1000_2GeV->SetLineColor(3);
  gr_120_LLP_mh350_pl1000_3GeV->SetLineColor(2);
  gr_120_LLP_mh350_pl1000_4GeV->SetLineColor(1);
  gr_120_LLP_mh350_pl1000_2GeV->Draw("L*");
  gr_120_LLP_mh350_pl1000_3GeV->Draw("L*");
  gr_120_LLP_mh350_pl1000_4GeV->Draw("L*");
  gr_120_LLP_mh350_pl1000_1GeV->SetTitle("htSum Rate vs. Signal Efficiency for One Jet Delayed Hits and Hit Fraction;LLP, mh=350GeV, c#scale[1.2]{#tau}=1m Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  m_LLPMh350pl1000_htSum430->SetMarkerStyle(21);
  m_LLPMh350pl1000_htSum430->SetMarkerColor(3);
  m_LLPMh350pl1000_htSum430->Draw();
  auto legend_htSum_mh350_pl1000 = new TLegend(0.45,0.15,0.9,0.35);
  legend_htSum_mh350_pl1000->AddEntry(gr_120_LLP_mh350_pl1000_1GeV,"H_{T}>120 GeV, with timing cuts, 1GeV");
  legend_htSum_mh350_pl1000->AddEntry(gr_120_LLP_mh350_pl1000_2GeV,"H_{T}>120 GeV, with timing cuts, 2GeV");
  legend_htSum_mh350_pl1000->AddEntry(gr_120_LLP_mh350_pl1000_3GeV,"H_{T}>120 GeV, with timing cuts, 3GeV");
  legend_htSum_mh350_pl1000->AddEntry(gr_120_LLP_mh350_pl1000_4GeV,"H_{T}>120 GeV, with timing cuts, 4GeV");
  legend_htSum_mh350_pl1000->AddEntry(m_LLPMh350pl1000_htSum430,"H_{T}>430 GeV, no timing cuts");
  legend_htSum_mh350_pl1000->Draw();
  gr_120_LLP_mh350_pl1000_1GeV->GetXaxis()->SetLimits(0.,1.);
  gr_120_LLP_mh350_pl1000_1GeV->GetHistogram()->SetMinimum(1.);
  gr_120_LLP_mh350_pl1000_1GeV->GetHistogram()->SetMaximum(100000.);
  c1_htSum_mh350_pl1000->SetLogy();
  c1_htSum_mh350_pl1000->SetGrid();
  c1_htSum_mh350_pl1000->SaveAs("plots/NuGun_singleJetTrigger_htSumRates_vs_SignalEff_Pl1000Mh350.pdf");
  c1_htSum_mh350_pl1000->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/NuGun_singleJetTrigger_htSumRates_vs_SignalEff_Pl1000Mh350.pdf");
  // ct=10m
  TCanvas *c1_htSum_mh350_pl10000 = new TCanvas("c1_htSum_mh350_pl10000","Graph Draw Options",200,10,600,400);
  gr_120_LLP_mh350_pl10000_1GeV->SetLineColor(4); // blue                                                                                       
  gr_120_LLP_mh350_pl10000_1GeV->GetHistogram()->SetMinimum(-5.);
  gr_120_LLP_mh350_pl10000_1GeV->Draw("AL*");
  gr_120_LLP_mh350_pl10000_2GeV->SetLineColor(3);
  gr_120_LLP_mh350_pl10000_3GeV->SetLineColor(2);
  gr_120_LLP_mh350_pl10000_4GeV->SetLineColor(1);
  gr_120_LLP_mh350_pl10000_2GeV->Draw("L*");
  gr_120_LLP_mh350_pl10000_3GeV->Draw("L*");
  gr_120_LLP_mh350_pl10000_4GeV->Draw("L*");
  gr_120_LLP_mh350_pl10000_1GeV->SetTitle("htSum Rate vs. Signal Efficiency for One Jet Delayed Hits and Hit Fraction;LLP, mh=350GeV, c#scale[1.2]{#tau}=10m Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  m_LLPMh350pl10000_htSum430->SetMarkerStyle(21);
  m_LLPMh350pl10000_htSum430->SetMarkerColor(3);
  m_LLPMh350pl10000_htSum430->Draw();
  auto legend_htSum_mh350_pl10000 = new TLegend(0.45,0.15,0.9,0.35);
  legend_htSum_mh350_pl10000->AddEntry(gr_120_LLP_mh350_pl10000_1GeV,"H_{T}>120 GeV, with timing cuts, 1GeV");
  legend_htSum_mh350_pl10000->AddEntry(gr_120_LLP_mh350_pl10000_2GeV,"H_{T}>120 GeV, with timing cuts, 2GeV");
  legend_htSum_mh350_pl10000->AddEntry(gr_120_LLP_mh350_pl10000_3GeV,"H_{T}>120 GeV, with timing cuts, 3GeV");
  legend_htSum_mh350_pl10000->AddEntry(gr_120_LLP_mh350_pl10000_4GeV,"H_{T}>120 GeV, with timing cuts, 4GeV");
  legend_htSum_mh350_pl10000->AddEntry(m_LLPMh350pl10000_htSum430,"H_{T}>430 GeV, no timing cuts");
  legend_htSum_mh350_pl10000->Draw();
  gr_120_LLP_mh350_pl10000_1GeV->GetXaxis()->SetLimits(0.,1.);
  gr_120_LLP_mh350_pl10000_1GeV->GetHistogram()->SetMinimum(1.);
  gr_120_LLP_mh350_pl10000_1GeV->GetHistogram()->SetMaximum(100000.);
  c1_htSum_mh350_pl10000->SetLogy();
  c1_htSum_mh350_pl10000->SetGrid();
  c1_htSum_mh350_pl10000->SaveAs("plots/NuGun_singleJetTrigger_htSumRates_vs_SignalEff_Pl10000Mh350.pdf");
  c1_htSum_mh350_pl10000->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/NuGun_singleJetTrigger_htSumRates_vs_SignalEff_Pl10000Mh350.pdf");
  // mh = 250 GeV
  // htSum rate, LLP mh=250GeV, ct=0.5m efficiency, regional                      
  TCanvas *c1_htSum_mh250_pl500 = new TCanvas("c1_htSum_mh250_pl500","Graph Draw Options",200,10,600,400);
  gr_120_LLP_mh250_pl500_1GeV->SetLineColor(4); // blue                                
  gr_120_LLP_mh250_pl500_1GeV->GetHistogram()->SetMinimum(-5.);
  gr_120_LLP_mh250_pl500_1GeV->Draw("AL*");
  gr_120_LLP_mh250_pl500_2GeV->SetLineColor(3);
  gr_120_LLP_mh250_pl500_3GeV->SetLineColor(2);
  gr_120_LLP_mh250_pl500_4GeV->SetLineColor(1);
  gr_120_LLP_mh250_pl500_2GeV->Draw("L*");
  gr_120_LLP_mh250_pl500_3GeV->Draw("L*");
  gr_120_LLP_mh250_pl500_4GeV->Draw("L*");
  gr_120_LLP_mh250_pl500_1GeV->SetTitle("htSum Rate vs. Signal Efficiency for One Jet Delayed Hits and Hit Fraction;LLP, mh=250GeV, c#scale[1.2]{#tau}=0.5m Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  m_LLPMh250pl500_htSum430->SetMarkerStyle(21);
  m_LLPMh250pl500_htSum430->SetMarkerColor(3);
  m_LLPMh250pl500_htSum430->Draw();
  auto legend_htSum_mh250_pl500 = new TLegend(0.45,0.15,0.9,0.35);
  legend_htSum_mh250_pl500->AddEntry(gr_120_LLP_mh250_pl500_1GeV,"H_{T}>120 GeV, with timing cuts, 1GeV");
  legend_htSum_mh250_pl500->AddEntry(gr_120_LLP_mh250_pl500_2GeV,"H_{T}>120 GeV, with timing cuts, 2GeV");
  legend_htSum_mh250_pl500->AddEntry(gr_120_LLP_mh250_pl500_3GeV,"H_{T}>120 GeV, with timing cuts, 3GeV");
  legend_htSum_mh250_pl500->AddEntry(gr_120_LLP_mh250_pl500_4GeV,"H_{T}>120 GeV, with timing cuts, 4GeV");
  legend_htSum_mh250_pl500->AddEntry(m_LLPMh250pl500_htSum430,"H_{T}>430 GeV, no timing cuts");
  legend_htSum_mh250_pl500->Draw();
  gr_120_LLP_mh250_pl500_1GeV->GetXaxis()->SetLimits(0.,1.);
  gr_120_LLP_mh250_pl500_1GeV->GetHistogram()->SetMinimum(1.);
  gr_120_LLP_mh250_pl500_1GeV->GetHistogram()->SetMaximum(100000.);
  c1_htSum_mh250_pl500->SetLogy();
  c1_htSum_mh250_pl500->SetGrid();
  c1_htSum_mh250_pl500->SaveAs("plots/NuGun_singleJetTrigger_htSumRates_vs_SignalEff_Pl500Mh250.pdf");
  c1_htSum_mh250_pl500->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/NuGun_singleJetTrigger_htSumRates_vs_SignalEff_Pl500Mh250.pdf");
  // ct=1m
  TCanvas *c1_htSum_mh250_pl1000 = new TCanvas("c1_htSum_mh250_pl1000","Graph Draw Options",200,10,600,400);
  gr_120_LLP_mh250_pl1000_1GeV->SetLineColor(4); // blue  
  gr_120_LLP_mh250_pl1000_1GeV->GetHistogram()->SetMinimum(-5.);
  gr_120_LLP_mh250_pl1000_1GeV->Draw("AL*");
  gr_120_LLP_mh250_pl1000_2GeV->SetLineColor(3);
  gr_120_LLP_mh250_pl1000_3GeV->SetLineColor(2);
  gr_120_LLP_mh250_pl1000_4GeV->SetLineColor(1);
  gr_120_LLP_mh250_pl1000_2GeV->Draw("L*");
  gr_120_LLP_mh250_pl1000_3GeV->Draw("L*");
  gr_120_LLP_mh250_pl1000_4GeV->Draw("L*");
  gr_120_LLP_mh250_pl1000_1GeV->SetTitle("htSum Rate vs. Signal Efficiency for One Jet Delayed Hits and Hit Fraction;LLP, mh=250GeV, c#scale[1.2]{#tau}=1m Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  m_LLPMh250pl1000_htSum430->SetMarkerStyle(21);
  m_LLPMh250pl1000_htSum430->SetMarkerColor(3);
  m_LLPMh250pl1000_htSum430->Draw();
  auto legend_htSum_mh250_pl1000 = new TLegend(0.45,0.15,0.9,0.35);
  legend_htSum_mh250_pl1000->AddEntry(gr_120_LLP_mh250_pl1000_1GeV,"H_{T}>120 GeV, with timing cuts, 1GeV");
  legend_htSum_mh250_pl1000->AddEntry(gr_120_LLP_mh250_pl1000_2GeV,"H_{T}>120 GeV, with timing cuts, 2GeV");
  legend_htSum_mh250_pl1000->AddEntry(gr_120_LLP_mh250_pl1000_3GeV,"H_{T}>120 GeV, with timing cuts, 3GeV");
  legend_htSum_mh250_pl1000->AddEntry(gr_120_LLP_mh250_pl1000_4GeV,"H_{T}>120 GeV, with timing cuts, 4GeV");
  legend_htSum_mh250_pl1000->AddEntry(m_LLPMh250pl1000_htSum430,"H_{T}>430 GeV, no timing cuts");
  legend_htSum_mh250_pl1000->Draw();
  gr_120_LLP_mh250_pl1000_1GeV->GetXaxis()->SetLimits(0.,1.);
  gr_120_LLP_mh250_pl1000_1GeV->GetHistogram()->SetMinimum(1.);
  gr_120_LLP_mh250_pl1000_1GeV->GetHistogram()->SetMaximum(100000.);
  c1_htSum_mh250_pl1000->SetLogy();
  c1_htSum_mh250_pl1000->SetGrid();
  c1_htSum_mh250_pl1000->SaveAs("plots/NuGun_singleJetTrigger_htSumRates_vs_SignalEff_Pl1000Mh250.pdf");
  c1_htSum_mh250_pl1000->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/NuGun_singleJetTrigger_htSumRates_vs_SignalEff_Pl1000Mh250.pdf");
  // ct=10m      
  TCanvas *c1_htSum_mh250_pl10000 = new TCanvas("c1_htSum_mh250_pl10000","Graph Draw Options",200,10,600,400);
  gr_120_LLP_mh250_pl10000_1GeV->SetLineColor(4); // blue
  gr_120_LLP_mh250_pl10000_1GeV->GetHistogram()->SetMinimum(-5.);
  gr_120_LLP_mh250_pl10000_1GeV->Draw("AL*");
  gr_120_LLP_mh250_pl10000_2GeV->SetLineColor(3);
  gr_120_LLP_mh250_pl10000_3GeV->SetLineColor(2);
  gr_120_LLP_mh250_pl10000_4GeV->SetLineColor(1);
  gr_120_LLP_mh250_pl10000_2GeV->Draw("L*");
  gr_120_LLP_mh250_pl10000_3GeV->Draw("L*");
  gr_120_LLP_mh250_pl10000_4GeV->Draw("L*");
  gr_120_LLP_mh250_pl10000_1GeV->SetTitle("htSum Rate vs. Signal Efficiency for One Jet Delayed Hits and Hit Fraction;LLP, mh=250GeV, c#scale[1.2]{#tau}=10m Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  m_LLPMh250pl10000_htSum430->SetMarkerStyle(21);
  m_LLPMh250pl10000_htSum430->SetMarkerColor(3);
  m_LLPMh250pl10000_htSum430->Draw();
  auto legend_htSum_mh250_pl10000 = new TLegend(0.45,0.15,0.9,0.35);
  legend_htSum_mh250_pl10000->AddEntry(gr_120_LLP_mh250_pl10000_1GeV,"H_{T}>120 GeV, with timing cuts, 1GeV");
  legend_htSum_mh250_pl10000->AddEntry(gr_120_LLP_mh250_pl10000_2GeV,"H_{T}>120 GeV, with timing cuts, 2GeV");
  legend_htSum_mh250_pl10000->AddEntry(gr_120_LLP_mh250_pl10000_3GeV,"H_{T}>120 GeV, with timing cuts, 3GeV");
  legend_htSum_mh250_pl10000->AddEntry(gr_120_LLP_mh250_pl10000_4GeV,"H_{T}>120 GeV, with timing cuts, 4GeV");
  legend_htSum_mh250_pl10000->AddEntry(m_LLPMh250pl10000_htSum430,"H_{T}>430 GeV, no timing cuts");
  legend_htSum_mh250_pl10000->Draw();
  gr_120_LLP_mh250_pl10000_1GeV->GetXaxis()->SetLimits(0.,1.);
  gr_120_LLP_mh250_pl10000_1GeV->GetHistogram()->SetMinimum(1.);
  gr_120_LLP_mh250_pl10000_1GeV->GetHistogram()->SetMaximum(100000.);
  c1_htSum_mh250_pl10000->SetLogy();
  c1_htSum_mh250_pl10000->SetGrid();
  c1_htSum_mh250_pl10000->SaveAs("plots/NuGun_singleJetTrigger_htSumRates_vs_SignalEff_Pl10000Mh250.pdf");
  c1_htSum_mh250_pl10000->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/NuGun_singleJetTrigger_htSumRates_vs_SignalEff_Pl10000Mh250.pdf");
  // mh = 125 GeV
  // htSum rate, LLP mh=125GeV, ct=0.5m efficiency, regional                         
  TCanvas *c1_htSum_mh125_pl500 = new TCanvas("c1_htSum_mh125_pl500","Graph Draw Options",200,10,600,400);
  gr_120_LLP_mh125_pl500_1GeV->SetLineColor(4); // blue                                   
  gr_120_LLP_mh125_pl500_1GeV->GetHistogram()->SetMinimum(-5.);
  gr_120_LLP_mh125_pl500_1GeV->Draw("AL*");
  gr_120_LLP_mh125_pl500_2GeV->SetLineColor(3);
  gr_120_LLP_mh125_pl500_3GeV->SetLineColor(2);
  gr_120_LLP_mh125_pl500_4GeV->SetLineColor(1);
  gr_120_LLP_mh125_pl500_2GeV->Draw("L*");
  gr_120_LLP_mh125_pl500_3GeV->Draw("L*");
  gr_120_LLP_mh125_pl500_4GeV->Draw("L*");
  gr_120_LLP_mh125_pl500_1GeV->SetTitle("htSum Rate vs. Signal Efficiency for One Jet Delayed Hits and Hit Fraction;LLP, mh=125GeV, c#scale[1.2]{#tau}=0.5m Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  m_LLPMh125pl500_htSum430->SetMarkerStyle(21);
  m_LLPMh125pl500_htSum430->SetMarkerColor(3);
  m_LLPMh125pl500_htSum430->Draw();
  auto legend_htSum_mh125_pl500 = new TLegend(0.45,0.15,0.9,0.35);
  legend_htSum_mh125_pl500->AddEntry(gr_120_LLP_mh125_pl500_1GeV,"H_{T}>120 GeV, with timing cuts, 1GeV");
  legend_htSum_mh125_pl500->AddEntry(gr_120_LLP_mh125_pl500_2GeV,"H_{T}>120 GeV, with timing cuts, 2GeV");
  legend_htSum_mh125_pl500->AddEntry(gr_120_LLP_mh125_pl500_3GeV,"H_{T}>120 GeV, with timing cuts, 3GeV");
  legend_htSum_mh125_pl500->AddEntry(gr_120_LLP_mh125_pl500_4GeV,"H_{T}>120 GeV, with timing cuts, 4GeV");
  legend_htSum_mh125_pl500->AddEntry(m_LLPMh125pl500_htSum430,"H_{T}>430 GeV, no timing cuts");
  legend_htSum_mh125_pl500->Draw();
  gr_120_LLP_mh125_pl500_1GeV->GetXaxis()->SetLimits(0.,1.);
  gr_120_LLP_mh125_pl500_1GeV->GetHistogram()->SetMinimum(1.);
  gr_120_LLP_mh125_pl500_1GeV->GetHistogram()->SetMaximum(100000.);
  c1_htSum_mh125_pl500->SetLogy();
  c1_htSum_mh125_pl500->SetGrid();
  c1_htSum_mh125_pl500->SaveAs("plots/NuGun_singleJetTrigger_htSumRates_vs_SignalEff_Pl500Mh125.pdf");
  c1_htSum_mh125_pl500->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/NuGun_singleJetTrigger_htSumRates_vs_SignalEff_Pl500Mh125.pdf");
  // ct=1m                                           
  TCanvas *c1_htSum_mh125_pl1000 = new TCanvas("c1_htSum_mh125_pl1000","Graph Draw Options",200,10,600,400);
  gr_120_LLP_mh125_pl1000_1GeV->SetLineColor(4); // blue  
  gr_120_LLP_mh125_pl1000_1GeV->GetHistogram()->SetMinimum(-5.);
  gr_120_LLP_mh125_pl1000_1GeV->Draw("AL*");
  gr_120_LLP_mh125_pl1000_2GeV->SetLineColor(3);
  gr_120_LLP_mh125_pl1000_3GeV->SetLineColor(2);
  gr_120_LLP_mh125_pl1000_4GeV->SetLineColor(1);
  gr_120_LLP_mh125_pl1000_2GeV->Draw("L*");
  gr_120_LLP_mh125_pl1000_3GeV->Draw("L*");
  gr_120_LLP_mh125_pl1000_4GeV->Draw("L*");
  gr_120_LLP_mh125_pl1000_1GeV->SetTitle("htSum Rate vs. Signal Efficiency for One Jet Delayed Hits and Hit Fraction;LLP, mh=125GeV, c#scale[1.2]{#tau}=1m Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  m_LLPMh125pl1000_htSum430->SetMarkerStyle(21);
  m_LLPMh125pl1000_htSum430->SetMarkerColor(3);
  m_LLPMh125pl1000_htSum430->Draw();
  auto legend_htSum_mh125_pl1000 = new TLegend(0.45,0.15,0.9,0.35);
  legend_htSum_mh125_pl1000->AddEntry(gr_120_LLP_mh125_pl1000_1GeV,"H_{T}>120 GeV, with timing cuts, 1GeV");
  legend_htSum_mh125_pl1000->AddEntry(gr_120_LLP_mh125_pl1000_2GeV,"H_{T}>120 GeV, with timing cuts, 2GeV");
  legend_htSum_mh125_pl1000->AddEntry(gr_120_LLP_mh125_pl1000_3GeV,"H_{T}>120 GeV, with timing cuts, 3GeV");
  legend_htSum_mh125_pl1000->AddEntry(gr_120_LLP_mh125_pl1000_4GeV,"H_{T}>120 GeV, with timing cuts, 4GeV");
  legend_htSum_mh125_pl1000->AddEntry(m_LLPMh125pl1000_htSum430,"H_{T}>430 GeV, no timing cuts");
  legend_htSum_mh125_pl1000->Draw();
  gr_120_LLP_mh125_pl1000_1GeV->GetXaxis()->SetLimits(0.,1.);
  gr_120_LLP_mh125_pl1000_1GeV->GetHistogram()->SetMinimum(1.);
  gr_120_LLP_mh125_pl1000_1GeV->GetHistogram()->SetMaximum(100000.);
  c1_htSum_mh125_pl1000->SetLogy();
  c1_htSum_mh125_pl1000->SetGrid();
  c1_htSum_mh125_pl1000->SaveAs("plots/NuGun_singleJetTrigger_htSumRates_vs_SignalEff_Pl1000Mh125.pdf");
  c1_htSum_mh125_pl1000->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/NuGun_singleJetTrigger_htSumRates_vs_SignalEff_Pl1000Mh125.pdf");
  // ct=10m 
  TCanvas *c1_htSum_mh125_pl10000 = new TCanvas("c1_htSum_mh125_pl10000","Graph Draw Options",200,10,600,400);
  gr_120_LLP_mh125_pl10000_1GeV->SetLineColor(4); // blue 
  gr_120_LLP_mh125_pl10000_1GeV->GetHistogram()->SetMinimum(-5.);
  gr_120_LLP_mh125_pl10000_1GeV->Draw("AL*");
  gr_120_LLP_mh125_pl10000_2GeV->SetLineColor(3);
  gr_120_LLP_mh125_pl10000_3GeV->SetLineColor(2);
  gr_120_LLP_mh125_pl10000_4GeV->SetLineColor(1);
  gr_120_LLP_mh125_pl10000_2GeV->Draw("L*");
  gr_120_LLP_mh125_pl10000_3GeV->Draw("L*");
  gr_120_LLP_mh125_pl10000_4GeV->Draw("L*");
  gr_120_LLP_mh125_pl10000_1GeV->SetTitle("htSum Rate vs. Signal Efficiency for One Jet Delayed Hits and Hit Fraction;LLP, mh=125GeV, c#scale[1.2]{#tau}=10m Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  m_LLPMh125pl10000_htSum430->SetMarkerStyle(21);
  m_LLPMh125pl10000_htSum430->SetMarkerColor(3);
  m_LLPMh125pl10000_htSum430->Draw();
  auto legend_htSum_mh125_pl10000 = new TLegend(0.45,0.15,0.9,0.35);
  legend_htSum_mh125_pl10000->AddEntry(gr_120_LLP_mh125_pl10000_1GeV,"H_{T}>120 GeV, with timing cuts, 1GeV");
  legend_htSum_mh125_pl10000->AddEntry(gr_120_LLP_mh125_pl10000_2GeV,"H_{T}>120 GeV, with timing cuts, 2GeV");
  legend_htSum_mh125_pl10000->AddEntry(gr_120_LLP_mh125_pl10000_3GeV,"H_{T}>120 GeV, with timing cuts, 3GeV");
  legend_htSum_mh125_pl10000->AddEntry(gr_120_LLP_mh125_pl10000_4GeV,"H_{T}>120 GeV, with timing cuts, 4GeV");
  legend_htSum_mh125_pl10000->AddEntry(m_LLPMh125pl10000_htSum430,"H_{T}>430 GeV, no timing cuts");
  legend_htSum_mh125_pl10000->Draw();
  gr_120_LLP_mh125_pl10000_1GeV->GetXaxis()->SetLimits(0.,1.);
  gr_120_LLP_mh125_pl10000_1GeV->GetHistogram()->SetMinimum(1.);
  gr_120_LLP_mh125_pl10000_1GeV->GetHistogram()->SetMaximum(100000.);
  c1_htSum_mh125_pl10000->SetLogy();
  c1_htSum_mh125_pl10000->SetGrid();
  c1_htSum_mh125_pl10000->SaveAs("plots/NuGun_singleJetTrigger_htSumRates_vs_SignalEff_Pl10000Mh125.pdf");
  c1_htSum_mh125_pl10000->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/NuGun_singleJetTrigger_htSumRates_vs_SignalEff_Pl10000Mh125.pdf");

  // htSum rate, QCD, regional
  TCanvas *c2_htSum = new TCanvas("c2_htSum","Graph Draw Options",200,10,600,400);
  gr_120_QCD_1GeV->SetLineColor(4); // blue  
  gr_120_QCD_1GeV->GetHistogram()->SetMinimum(-5.);
  gr_120_QCD_1GeV->Draw("AL*");
  gr_120_QCD_2GeV->SetLineColor(3);
  gr_120_QCD_3GeV->SetLineColor(2);
  gr_120_QCD_4GeV->SetLineColor(1);
  gr_120_QCD_2GeV->Draw("L*");
  gr_120_QCD_3GeV->Draw("L*");
  gr_120_QCD_4GeV->Draw("L*");
  gr_120_QCD_1GeV->SetTitle("htSum Rate vs. Background Efficiency for One Jet Delayed Hits and Hit Fraction;QCD Efficiency;Neutrino Gun Rate (kHz, unnorm)   ");
  m_QCD_htSum430->SetMarkerStyle(21);
  m_QCD_htSum430->SetMarkerColor(3);
  m_QCD_htSum430->Draw();
  gr_120_QCD_1GeV->GetXaxis()->SetLimits(0.,1.);
  c2_htSum->SetGrid();
  auto legend2_htSum = new TLegend(0.45,0.15,0.9,0.35);
  legend2_htSum->AddEntry(gr_120_QCD_1GeV,"H_{T}>120 GeV, with timing cuts, 1GeV");
  legend2_htSum->AddEntry(gr_120_QCD_2GeV,"H_{T}>120 GeV, with timing cuts, 2GeV");
  legend2_htSum->AddEntry(gr_120_QCD_3GeV,"H_{T}>120 GeV, with timing cuts, 3GeV");
  legend2_htSum->AddEntry(gr_120_QCD_4GeV,"H_{T}>120 GeV, with timing cuts, 4GeV");
  legend2_htSum->AddEntry(m_QCD_htSum430,"H_{T}>430 GeV, no timing cuts");
  legend2_htSum->Draw();
  c2_htSum->SetLogy();
  gr_120_QCD_1GeV->GetHistogram()->SetMinimum(1.);
  gr_120_QCD_1GeV->GetHistogram()->SetMaximum(100000.);
  c2_htSum->SaveAs("plots/NuGun_singleJetTrigger_htSumRates_vs_BackgroundEff.pdf");
  c2_htSum->SaveAs("/eos/user/g/gkopp/www/HCAL_LLP/NuGun_singleJetTrigger_htSumRates_vs_BackgroundEff.pdf");


  return 0;
}
