// Match HCAL TP with L1em. Check if fg[0] is set..

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TChain.h"
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include "boost/program_options.hpp"

// #include "L1Trigger/L1TNtuples/interface/L1AnalysisEventDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoVertexDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisCaloTPDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisGeneratorDataFormat.h"

void fgBitAnalysis(const std::string &inputFileDirectory, double hcal_l1_dR);
// void energyRatioAnalysis(const std::string &inputFileDirectory, int energy, double genThresh, double shortThresh, double longThresh, double offset, double slope, double coneSize);
void energyRatioAnalysis(const std::string &inputFileDirectory, double genThresh, double shortThresh, double longThresh, double offset, double slope, double coneSize);

namespace po = boost::program_options;

int main(int argc, char *argv[])
{
  // int energy;
  double genThresh, bitShortThresh, bitLongthresh, bitOffset, bitSlope, coneSize;
  std::string ntuplePath;

  po::options_description desc("Allowed Program Options");
  desc.add_options()
                  ("help", "produce help messages")
                  ("input", po::value<std::string>(&ntuplePath)->default_value("l1ntuples/"), "Path with input files")
                  ("genThresh", po::value<double>(&genThresh)->default_value(3.0), "Generator Particle Threshold")
                  ("bitShortThresh", po::value<double>(&bitShortThresh)->default_value(10.0), "Bit Short Threshold")
                  ("bitLongThresh", po::value<double>(&bitLongthresh)->default_value(10.0), "Bit Long Threshold")
                  ("bitOffset", po::value<double>(&bitOffset)->default_value(10.1), "Bit Offset")
                  ("bitSlope", po::value<double>(&bitSlope)->default_value(100.2), "Bit Slope")
                  ("coneSize", po::value<double>(&coneSize)->default_value(0.20), "matching coneSize");

                  // ("energy", po::value<int>(&energy)->default_value(30), "Energy to use")
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help"))
  {
    std::cout << desc << std::endl;
    return 1;
  }

  energyRatioAnalysis(ntuplePath, genThresh, bitShortThresh, bitLongthresh, bitOffset, bitSlope, coneSize);
  // energyRatioAnalysis(ntuplePath, energy, genThresh, bitShortThresh, bitLongthresh, bitOffset, bitSlope, coneSize);

  return 0;
}

double deltaPhi(double phi1, double phi2)
{
  double result = phi1 - phi2;
  if (fabs(result) > 9999)
    return result;
  while (result > TMath::Pi())
    result -= 2 * TMath::Pi();
  while (result <= -TMath::Pi())
    result += 2 * TMath::Pi();
  return result;
}

double deltaR(double eta1, double phi1, double eta2, double phi2)
{
  double deta = eta1 - eta2;
  double dphi = deltaPhi(phi1, phi2);
  return sqrt(deta * deta + dphi * dphi);
}

// void energyRatioAnalysis(const std::string &inputFileDirectory, int energy, double genThresh, double shortThresh, double longThresh, double offset, double slope, double coneSize)
void energyRatioAnalysis(const std::string &inputFileDirectory, double genThresh, double shortThresh, double longThresh, double offset, double slope, double coneSize)
{
  std::cout << "generator threshold set to: " << genThresh << std::endl;
  std::cout << "short threshold set to: " << shortThresh << std::endl;
  std::cout << "long threshold set to: " << longThresh << std::endl;
  std::cout << "offset set to: " << offset << std::endl;
  std::cout << "slope set to: " << slope << std::endl;

  std::string inputFileZ(inputFileDirectory);
  // inputFileZ += "/L1Ntuple_zee.root";
  inputFileZ += "/singleElectronE";
  // inputFileZ += std::to_string(energy);
  inputFileZ += ".root";

  // std::string inputFileQ(inputFileDirectory);
  // inputFileQ += "/L1Ntuple_QCDPU.root";

  std::string inputFilePi(inputFileDirectory);
  // inputFilePi += "/L1Ntuple_singlepi.root";
  inputFilePi += "/singlePionE";
  // inputFilePi += std::to_string(energy);
  inputFilePi += ".root";

  std::string outputFilename = "energyRatioPlots";
  // outputFilename += std::to_string(energy);
  outputFilename += ".root";

  TFile *outfile = TFile::Open(outputFilename.c_str(), "recreate");

  std::cout << "loading up the TChain" << std::endl;

  // Z->ee Trees
  TChain *treeL1emuZ = new TChain("l1UpgradeEmuTree/L1UpgradeTree");
  treeL1emuZ->Add(inputFileZ.c_str());

  TChain *treeL1TPemuZ = new TChain("l1CaloTowerEmuTree/L1CaloTowerTree");
  treeL1TPemuZ->Add(inputFileZ.c_str());

  TChain *treeL1GenZ = new TChain("l1GeneratorTree/L1GenTree");
  treeL1GenZ->Add(inputFileZ.c_str());

  // Single Pi Trees
  TChain *treeL1emuPi = new TChain("l1UpgradeEmuTree/L1UpgradeTree");
  treeL1emuPi->Add(inputFilePi.c_str());

  TChain *treeL1TPemuPi = new TChain("l1CaloTowerEmuTree/L1CaloTowerTree");
  treeL1TPemuPi->Add(inputFilePi.c_str());

  TChain *treeL1GenPi = new TChain("l1GeneratorTree/L1GenTree");
  treeL1GenPi->Add(inputFilePi.c_str());

  // Z->ee Branches
  // L1Analysis::L1AnalysisL1UpgradeDataFormat *l1emuZ_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  // treeL1emuZ->SetBranchAddress("L1Upgrade", &l1emuZ_);

  L1Analysis::L1AnalysisCaloTPDataFormat *l1TPemuZ_ = new L1Analysis::L1AnalysisCaloTPDataFormat();
  treeL1TPemuZ->SetBranchAddress("CaloTP", &l1TPemuZ_);

  L1Analysis::L1AnalysisGeneratorDataFormat *l1GenZ_ = new L1Analysis::L1AnalysisGeneratorDataFormat();
  treeL1GenZ->SetBranchAddress("Generator", &l1GenZ_);

  // Single Pi Branches
  // L1Analysis::L1AnalysisL1UpgradeDataFormat *l1emuPi_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  // treeL1emuPi->SetBranchAddress("L1Upgrade", &l1emuPi_);

  L1Analysis::L1AnalysisCaloTPDataFormat *l1TPemuPi_ = new L1Analysis::L1AnalysisCaloTPDataFormat();
  treeL1TPemuPi->SetBranchAddress("CaloTP", &l1TPemuPi_);

  L1Analysis::L1AnalysisGeneratorDataFormat *l1GenPi_ = new L1Analysis::L1AnalysisGeneratorDataFormat();
  treeL1GenPi->SetBranchAddress("Generator", &l1GenPi_);

  TH1F *minGenETPdR = new TH1F("minTPGenE", "", 60, 0, 6);
  TH1F *minGenPiTPdR = new TH1F("minTPGenPi", "", 60, 0, 6);

  TH1F *shortEnergyE = new TH1F("shortEnE", "", 125, -10, 240);
  TH1F *longEnergyE = new TH1F("longEnE", "", 125, -10, 240);
  TH1F *energyRatioE = new TH1F("nrgRatioE", "", 40, 0, 20);

  TH1F *singleTPshortEnergyE = new TH1F("sTPshortE","",100,0,100);
  TH1F *singleTPlongEnergyE = new TH1F("sTPlongE","",100,0,100);

  TH1F *shortEnergyPi = new TH1F("shortEnPi", "", 125, -10, 240);
  TH1F *longEnergyPi = new TH1F("longEnPi", "", 125, -10, 240);
  TH1F *energyRatioPi = new TH1F("nrgRatioPi", "", 40, 0, 20);

  TH1F *singleTPshortEnergyPi = new TH1F("sTPshortPi","",100,0,100);
  TH1F *singleTPlongEnergyPi = new TH1F("sTPlongPi","",100,0,100);

  TH1F *hcalTPMultipilictyE = new TH1F("hcalTPMulE", "", 20, 0, 20);
  TH1F *hcalTPMultipilictyPi = new TH1F("hcalTPMulPi", "", 20, 0, 20);

  TH1F *singleExprE = new TH1F("singleTPexprE", "", 250, 0, 1000);
  singleExprE->SetCanExtend(TH1::kAllAxes);
  TH1F *singleExprPi = new TH1F("singleTPexprPi", "", 250, 0, 1000);
  singleExprPi->SetCanExtend(TH1::kAllAxes);

  TH1F *eff_num = new TH1F("numEff", "", 100, 0, 200);
  TH1F *eff_den = new TH1F("denEff", "", 100, 0, 200);

  TH1F *mistag_num = new TH1F("numMiss","",50,0,50);
  TH1F *mistag_den = new TH1F("denMiss","",50,0,50);

  Long64_t nentries;
  nentries = treeL1GenZ->GetEntries();

  int ntupleFineGrain = 0;
  int newFineGrain = 0;
  int nParticles = 0;
  int nMatchedWithBit = 0;

  for (Long64_t jentry = 0; jentry < nentries; jentry++)
  {
    // if ((jentry % 1000) == 0)
    // std::cout << "Done " << jentry << " events of Zee " << nentries << std::endl;

    treeL1TPemuZ->GetEntry(jentry);
    // treeL1emuZ->GetEntry(jentry);
    treeL1GenZ->GetEntry(jentry);

    for (int g = 0; g < l1GenZ_->nPart; g++)
    {
      if (fabs(l1GenZ_->partEta[g]) <= 2.8)
        continue;
      if (fabs(l1GenZ_->partEta[g]) >= 5.3)
        continue;
      if (abs(l1GenZ_->partId[g]) != 11) // Electron
        continue;
      // if (abs(l1GenZ_->partParent[g]) != 23) // Parent is Z
      //   continue;
      // if (!l1GenZ_->partFromHard[g]) // From Hard Process
      //   continue;
      if (l1GenZ_->partE[g] < genThresh)
        continue;

      eff_den->Fill(l1GenZ_->partPt[g]);
      nParticles++; // At this point we are sure this is an electron..
      int matchedTPs = 0;
      float shortSum = 0.0;
      float longSum = 0.0;
      double minGenTPdR = 999.;
      bool emfbit = false;

      for (int t = 0; t < l1TPemuZ_->nHCALTP; t++)
      {
        if (l1TPemuZ_->hcalTPet[t] <= 0.5 or abs(l1TPemuZ_->hcalTPieta[t]) <= 28)
          continue;

        double dR = deltaR(l1GenZ_->partEta[g], l1GenZ_->partPhi[g], l1TPemuZ_->hcalTPeta[t], l1TPemuZ_->hcalTPphiFunc[t]);

        if (dR < coneSize)
        {
          matchedTPs++;
          shortSum += l1TPemuZ_->hcalTPshortFiberE[t];
          longSum += l1TPemuZ_->hcalTPlongFiberE[t];
          
          singleTPshortEnergyE->Fill(l1TPemuZ_->hcalTPshortFiberE[t]);
          singleTPlongEnergyE->Fill(l1TPemuZ_->hcalTPlongFiberE[t]);
          if (l1TPemuZ_->hcalTPfineGrain[t])
            ntupleFineGrain++;

          if (l1TPemuZ_->hcalTPshortFiberE[t] > shortThresh and l1TPemuZ_->hcalTPlongFiberE[t] > longThresh)
          {
            singleExprE->Fill((l1TPemuZ_->hcalTPlongFiberE[t] - offset) * slope - l1TPemuZ_->hcalTPshortFiberE[t]);

            if (l1TPemuZ_->hcalTPshortFiberE[t] < (l1TPemuZ_->hcalTPlongFiberE[t] - offset) * slope) // fineGrainBit is set in this case!!
            {
              emfbit = true; // once this is true it should remain true..
              newFineGrain++;
            }
          }
        }
        if (dR < minGenTPdR)
        {
          minGenTPdR = dR;
        }

      } // End HCAL TP loop

      if (emfbit) // was the bit set by any TP?
      {
        nMatchedWithBit++; // increase the number of particles matched with a bit !!!
        eff_num->Fill(l1GenZ_->partPt[g]);
      }

      minGenETPdR->Fill(minGenTPdR);
      if (matchedTPs > 0)
      {
        energyRatioE->Fill(longSum / shortSum);
        shortEnergyE->Fill(shortSum);
        longEnergyE->Fill(longSum);
      }
      hcalTPMultipilictyE->Fill(matchedTPs);

    } // End Gen Particle Loop
  }   // End Z -> ee Loop

  std::cout << "ntuple FineGrains : " << ntupleFineGrain << std::endl;
  std::cout << "electron calculated finegrains: " << newFineGrain << std::endl;
  std::cout << "number of generator electrons : " << nParticles << std::endl;
  std::cout << "number of electrons matched with a set bit : " << nMatchedWithBit << std::endl;

  // reset counters for pion sample..
  nentries = treeL1GenPi->GetEntries();
  ntupleFineGrain = 0;
  newFineGrain = 0;
  nParticles = 0;
  nMatchedWithBit = 0;

  // run over Pion Sample
  for (Long64_t jentry = 0; jentry < nentries; jentry++)
  {
    // if ((jentry % 1000) == 0)
    // std::cout << "Done " << jentry << " events of Single Pi " << nentries << std::endl;

    treeL1TPemuPi->GetEntry(jentry);
    // treeL1emuZ->GetEntry(jentry);
    treeL1GenPi->GetEntry(jentry);

    for (int g = 0; g < l1GenPi_->nPart; g++)
    {
      if (fabs(l1GenPi_->partEta[g]) <= 2.8)
        continue;
      if (fabs(l1GenPi_->partEta[g]) >= 5.3)
        continue;
      if (abs(l1GenPi_->partId[g]) != 211) // Pion
        continue;
      // if (!l1GenPi_->partFromHard[g]) // From Hard Process
      // continue;
      if (l1GenPi_->partE[g] < genThresh)
        continue;

      mistag_den->Fill(l1GenPi_->partPt[g]);
      nParticles++; // generator pion "Accepted"

      int matchedTPs = 0;
      float shortSum = 0.0;
      float longSum = 0.0;
      double minGenTPdR = 999.;
      bool emfbit = false;

      for (int t = 0; t < l1TPemuPi_->nHCALTP; t++)
      {
        if (l1TPemuPi_->hcalTPet[t] <= 0.5 or abs(l1TPemuPi_->hcalTPieta[t]) <= 28)
          continue;

        double dR = deltaR(l1GenPi_->partEta[g], l1GenPi_->partPhi[g], l1TPemuPi_->hcalTPeta[t], l1TPemuPi_->hcalTPphiFunc[t]);

        if (dR < coneSize)
        {
          matchedTPs++;
          shortSum += l1TPemuPi_->hcalTPshortFiberE[t];
          longSum += l1TPemuPi_->hcalTPlongFiberE[t];

          singleTPshortEnergyPi->Fill(l1TPemuPi_->hcalTPshortFiberE[t]);
          singleTPlongEnergyPi->Fill(l1TPemuPi_->hcalTPlongFiberE[t]);

          if (l1TPemuPi_->hcalTPfineGrain[t])
            ntupleFineGrain++;

          if (l1TPemuPi_->hcalTPshortFiberE[t] > shortThresh and l1TPemuPi_->hcalTPlongFiberE[t] > longThresh)
          {
            singleExprPi->Fill((l1TPemuPi_->hcalTPlongFiberE[t] - offset) * slope - l1TPemuPi_->hcalTPshortFiberE[t]);

            if (l1TPemuPi_->hcalTPshortFiberE[t] < (l1TPemuPi_->hcalTPlongFiberE[t] - offset) * slope)
            {
              emfbit = true; // once this is true it should remain true..
              newFineGrain++;
            }
          }
        }
        if (dR < minGenTPdR)
        {
          minGenTPdR = dR;
        }

      } // End Loop over HCAL TPs

      if (emfbit)
      {
        nMatchedWithBit++;
        mistag_num->Fill(l1GenPi_->partPt[g]);
      }

      minGenPiTPdR->Fill(minGenTPdR);

      if (matchedTPs > 0)
      {
        energyRatioPi->Fill(longSum / shortSum);
        shortEnergyPi->Fill(shortSum);
        longEnergyPi->Fill(longSum);
      }

      hcalTPMultipilictyPi->Fill(matchedTPs);

    } // End Loop over Gen particles
  }   // End Single Pi Loop

  std::cout << "ntuple FineGrains : " << ntupleFineGrain << std::endl;
  std::cout << "pion calculated finegrains: " << newFineGrain << std::endl;
  std::cout << "number of generator pions : " << nParticles << std::endl;
  std::cout << "number of pions matched with a set bit : " << nMatchedWithBit << std::endl;

  outfile->cd();

  minGenETPdR->Write();
  minGenPiTPdR->Write();

  shortEnergyE->Write();
  longEnergyE->Write();
  energyRatioE->Write();

  shortEnergyPi->Write();
  longEnergyPi->Write();
  energyRatioPi->Write();

  hcalTPMultipilictyE->Write();
  hcalTPMultipilictyPi->Write();

  singleExprE->Write();
  singleExprPi->Write();

  eff_num->Write();
  eff_den->Write();

  mistag_num->Write();
  mistag_den->Write();

  singleTPshortEnergyE->Write();
  singleTPlongEnergyE->Write();
  singleTPshortEnergyPi->Write();
  singleTPlongEnergyPi->Write();
}

void fgBitAnalysis(const std::string &inputFileDirectory, double hcal_l1_dR)
{
  std::string inputFileZ(inputFileDirectory);
  inputFileZ += "/L1Ntuple_zee.root";

  std::string inputFileQ(inputFileDirectory);
  inputFileQ += "/L1Ntuple_qcd.root";

  std::string outputFilename = "l1egfinegrain_";
  outputFilename += std::to_string(int(hcal_l1_dR));
  outputFilename += ".root";

  TFile *outfile = TFile::Open(outputFilename.c_str(), "recreate");

  std::cout << "loading up the TChain" << std::endl;

  TChain *treeL1emuZ = new TChain("l1UpgradeEmuTree/L1UpgradeTree");
  treeL1emuZ->Add(inputFileZ.c_str());

  TChain *treeL1TPemuZ = new TChain("l1CaloTowerEmuTree/L1CaloTowerTree");
  treeL1TPemuZ->Add(inputFileZ.c_str());

  TChain *treeL1GenZ = new TChain("l1GeneratorTree/L1GenTree");
  treeL1GenZ->Add(inputFileZ.c_str());

  TChain *treeL1emuQ = new TChain("l1UpgradeEmuTree/L1UpgradeTree");
  treeL1emuQ->Add(inputFileQ.c_str());

  TChain *treeL1TPemuQ = new TChain("l1CaloTowerEmuTree/L1CaloTowerTree");
  treeL1TPemuQ->Add(inputFileQ.c_str());

  TChain *treeL1GenQ = new TChain("l1GeneratorTree/L1GenTree");
  treeL1GenQ->Add(inputFileQ.c_str());

  L1Analysis::L1AnalysisL1UpgradeDataFormat *l1emuZ_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  treeL1emuZ->SetBranchAddress("L1Upgrade", &l1emuZ_);

  L1Analysis::L1AnalysisCaloTPDataFormat *l1TPemuZ_ = new L1Analysis::L1AnalysisCaloTPDataFormat();
  treeL1TPemuZ->SetBranchAddress("CaloTP", &l1TPemuZ_);

  // add generator level tree to check mc truth
  L1Analysis::L1AnalysisGeneratorDataFormat *l1GenZ_ = new L1Analysis::L1AnalysisGeneratorDataFormat();
  treeL1GenZ->SetBranchAddress("Generator", &l1GenZ_);

  L1Analysis::L1AnalysisL1UpgradeDataFormat *l1emuQ_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  treeL1emuQ->SetBranchAddress("L1Upgrade", &l1emuQ_);

  L1Analysis::L1AnalysisCaloTPDataFormat *l1TPemuQ_ = new L1Analysis::L1AnalysisCaloTPDataFormat();
  treeL1TPemuQ->SetBranchAddress("CaloTP", &l1TPemuQ_);

  // add generator level tree to check mc truth
  // L1Analysis::L1AnalysisGeneratorDataFormat *l1GenQ_ = new L1Analysis::L1AnalysisGeneratorDataFormat();
  // treeL1GenQ->SetBranchAddress("Generator", &l1GenQ_);

  Long64_t nentries;
  nentries = treeL1emuZ->GetEntries();

  TH1F *minDrL1eHcalTP = new TH1F("minddrL1HcalTP", ";min dR(l1eg,hcalTP);", 100, 0, 10);
  TH1F *minDrL1eGene = new TH1F("mindrl1eggen", ";;min dR(l1eg,genEle)", 60, 0, 6);

  TH1F *eff_num = new TH1F("numEff", "", 100, 0, 200);
  TH1F *eff_den = new TH1F("denEff", "", 100, 0, 200);

  TH1F *purity_num = new TH1F("numPur", "", 100, 0, 200);
  TH1F *purity_den = new TH1F("denPur", "", 100, 0, 200);

  // zee sample loop
  for (Long64_t jentry = 0; jentry < nentries; jentry++)
  {
    if ((jentry % 1000) == 0)
      std::cout << "Done " << jentry << " events of Zee " << nentries << std::endl;

    treeL1TPemuZ->GetEntry(jentry);
    treeL1emuZ->GetEntry(jentry);
    treeL1GenZ->GetEntry(jentry);

    for (int i = 0; i < l1emuZ_->nEGs; i++)
    {
      // only look at eg with abs(ieta) > 28 and et > 3 GeV
      if (abs(l1emuZ_->egIEta[i]) <= 28 or l1emuZ_->egEt[i] <= 3.0)
        continue;
      double minL1GenDR = 999;
      int genIdx = -1;
      // Loop over gen particles to find a gen level electron
      // with a Z parent
      for (int g = 0; g < l1GenZ_->nPart; g++)
      {
        if (abs(l1GenZ_->partId[g]) != 11) // make sure it's an electron..
          continue;
        if (abs(l1GenZ_->partParent[g]) != 23) // make sure parent is Z
          continue;
        if (!l1GenZ_->partFromHard[g]) // from hard scatter process (only for efficiencies?)
          continue;

        double dR = deltaR(l1emuZ_->egEta[i], l1emuZ_->egPhi[i], l1GenZ_->partEta[g], l1GenZ_->partPhi[g]);
        if (dR < minL1GenDR)
        {
          minL1GenDR = dR;
          if (dR <= 0.2)
            genIdx = g; // id of the generator particle that is matched to the L1eg
        }
      }

      if (genIdx >= 0) // matched with a gen E from Zee..
      {
        eff_den->Fill(l1emuZ_->egEt[i]);
      }

      minDrL1eGene->Fill(minL1GenDR);

      double minL1hcalDR = 999;
      bool emfbit = false;
      for (int j = 0; j < l1TPemuZ_->nHCALTP; j++)
      {
        if (l1TPemuZ_->hcalTPet[j] <= 0.5 or abs(l1TPemuZ_->hcalTPieta[j]) <= 28)
          continue;

        double dR = deltaR(l1emuZ_->egEta[i], l1emuZ_->egPhi[i], l1TPemuZ_->hcalTPeta[j], l1TPemuZ_->hcalTPphiFunc[j]);
        // double dR = deltaR(l1emuZ_->egEta[i], l1emuZ_->egPhi[i], hcalEta, hcalPhi);
        if (dR <= hcal_l1_dR)
        {
          emfbit = emfbit or l1TPemuZ_->hcalTPfineGrain[j];
        }
        if (dR < minL1hcalDR)
        {
          minL1hcalDR = dR;
        }
      }

      if (genIdx >= 0) // relax to see what is going on?
      {
        minDrL1eHcalTP->Fill(minL1hcalDR);
      }
      if (emfbit)
      {
        purity_den->Fill(l1emuZ_->egEt[i]);
      }
      if (genIdx >= 0 and emfbit) // matched to hcalTP & emfbit is "on"
      {
        eff_num->Fill(l1emuZ_->egEt[i]);
        purity_num->Fill(l1emuZ_->egEt[i]);
      }
    }
  } // finish zee sample loop

  TH1F *fake_num = new TH1F("num_fake", "", 100, 0, 200);
  TH1F *fake_den = new TH1F("den_fake", "", 100, 0, 200);

  // qcd sample loop
  nentries = treeL1emuQ->GetEntries();
  for (Long64_t jentry = 0; jentry < nentries; jentry++)
  {
    if ((jentry % 1000) == 0)
      std::cout << "Done " << jentry << " events of qcd " << nentries << std::endl;

    treeL1TPemuQ->GetEntry(jentry);
    treeL1emuQ->GetEntry(jentry);
    // treeL1GenQ->GetEntry(jentry);

    for (int i = 0; i < l1emuQ_->nEGs; i++)
    {
      if (abs(l1emuQ_->egIEta[i]) <= 28 or l1emuQ_->egEt[i] <= 3.0)
        continue;

      fake_den->Fill(l1emuQ_->egEt[i]);

      bool emfbitQ = false;

      for (int j = 0; j < l1TPemuQ_->nHCALTP; j++)
      {
        if (l1TPemuQ_->hcalTPet[j] <= 0.5 or abs(l1TPemuQ_->hcalTPieta[j]) <= 28)
          continue;

        double dR = deltaR(l1emuQ_->egEta[i], l1emuQ_->egPhi[i], l1TPemuQ_->hcalTPeta[j], l1TPemuQ_->hcalTPphiFunc[j]);
        if (dR < hcal_l1_dR)
        {
          emfbitQ = emfbitQ or l1TPemuQ_->hcalTPfineGrain[j]; // if any are true, emfbit stays true
        }
      }
      if (emfbitQ)
      {
        fake_num->Fill(l1emuQ_->egEt[i]);
      }
    }
  }

  //
  outfile->cd();
  minDrL1eHcalTP->Write();
  minDrL1eGene->Write();
  eff_num->Write();
  eff_den->Write();
  purity_num->Write();
  purity_den->Write();
  fake_num->Write();
  fake_den->Write();
}
