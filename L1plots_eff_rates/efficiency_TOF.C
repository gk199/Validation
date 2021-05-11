#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <new>
#include<climits>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TEfficiency.h"
#include "TF1.h"
#include "TH2.h"
#include "TProfile.h"

#include "TStyle.h"
#include "TAttFill.h"
#include "TPaveStats.h"
#include "TMinuit.h"
#include "TPostScript.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TRandom.h"
#include "TPaletteAxis.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TCanvas.h"

using namespace std;

Double_t EMG(Double_t *x, Double_t *par){
  double x_prime = par[1]*(x[0] - par[0]);
  double term_1 = 0.5 * (1 + TMath::Erf(x_prime));
  double exp_modif = TMath::Exp(-par[2]*(x_prime - par[2]*0.5) );
  double term_2 = 0.5 * (1 + TMath::Erf(x_prime - par[2]));
  double func = term_1 - exp_modif*term_2;
  return func;
}

void efficiency_TOF(){
  TCanvas *c1 =new TCanvas("c1", " ", 0, 0,1000,1000);
  c1->Range(0,0,1,1);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetFrameBorderMode(0);
  c1->Draw();
  c1->SetGrid();
  gStyle->SetOptStat(0);

  string file_type = "LLP_mh1000_pl10000";

  TFile *g1 =TFile::Open("../rates_new_cond_LLP_mh1000_pl10000.root");
  TH1F *h1 = (TH1F*)g1->Get("TOFdelay");
  TH1F *h2 = (TH1F*)g1->Get("TOFdelay_trigger");
  TH1F *h3 = (TH1F*)g1->Get("TOFdelay_120trigger");

  h1->SetLineColorAlpha(kWhite, 1.);
  h1->SetTitle("Delayed Jet TOF Delay Efficiency for LLP_mh1000_pl10000");
  h1->GetXaxis()->SetTitle("LLP TOF + b TOF - TOF expected (ns)");
  h1->GetXaxis()->SetTitleSize(0.045);
  h1->GetXaxis()->SetTitleOffset(1.1);
  h1->GetYaxis()->SetTitle("Delayed Jet Efficiency");
  h1->GetYaxis()->SetTitleSize(0.045);
  h1->GetYaxis()->SetTitleOffset(1.1);
  h1->GetYaxis()->SetRangeUser(0.,1.5);
  if (file_type == "QCD") h1->GetYaxis()->SetRangeUser(0.,0.2);
  h1->GetXaxis()->SetRangeUser(0.,200);
  h1->Draw();

  c1->Update();
  c1->RedrawAxis();

  TEfficiency* pEff1 = 0;
  TEfficiency* pEff2 = 0;

  if(TEfficiency::CheckConsistency(*h2,*h1) && TEfficiency::CheckConsistency(*h3,*h1))
    {
      pEff1 = new TEfficiency(*h2,*h1);
      pEff1->SetLineWidth(3.);
      pEff1->SetLineColor(kGray);
      pEff1->Draw("same");
      pEff2 = new TEfficiency(*h3,*h1);
      pEff2->SetLineWidth(3.);
      pEff2->SetLineColor(kBlack);
      pEff2->Draw("same");
    }

  TLegend *legend1 = new TLegend(0.15, 0.7, 0.6, 0.8);
  legend1->SetTextFont(42);
  legend1->SetLineColor(0);
  legend1->SetTextSize(0.04);
  legend1->SetFillColor(0);
  legend1->AddEntry(pEff1, "Delayed L1 jet, p_{T}>40 GeV", "l");
  legend1->AddEntry(pEff2, "Delayed L1 jet, p_{T}>40 GeV, event HT>120", "l");

  legend1->Draw("same");

  char saveFile[100];
  sprintf(saveFile,"/eos/user/g/gkopp/www/HCAL_LLP/TimingBit/112X_TDCsim_DelayedJet/TOFefficiency_LLP_mh1000_pl10000.pdf");
  c1->SaveAs(saveFile);

}

