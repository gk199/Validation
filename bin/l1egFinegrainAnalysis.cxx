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
#include <iomanip> // std::setprecision
#include <fstream>
#include <string>
#include <stdlib.h>
#include "boost/program_options.hpp"
#include <TNtuple.h>
#include <TLorentzVector.h>

// #include "L1Trigger/L1TNtuples/interface/L1AnalysisEventDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoVertexDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisCaloTPDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisGeneratorDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/HcalTPCoordMapUtils.h"

void fgBitAnalysis(const std::string &inputFileDirectory);
void fgBitGenAnalysis(const std::string &inputFileDirectory);
// void energyRatioAnalysis(const std::string &inputFileDirectory, int energy, double genThresh, double energyFracCut, double shortThresh, double longThresh, double offset, double slope, double coneSize, bool pileup);
// void energyRatioAnalysis(const std::string &inputFileDirectory, double genThresh, double shortThresh, double longThresh, double offset, double slope, double coneSize);

namespace po = boost::program_options;

int main(int argc, char *argv[])
{
  std::vector<int> energies;
  double genThresh, eFC, bitShortThresh, bitLongthresh, bitOffset, bitSlope, coneSize;
  std::string ntuplePath;
  int mode;
  std::vector<int> defaultEnergies{30, 50, 70, 100, 150, 300};
  bool runZee, pileup;

  po::options_description desc("Allowed Program Options");
  desc.add_options()("help", "produce help messages")
 ("input", po::value<std::string>(&ntuplePath)->default_value(""), "Path with input files")
    ("energy", po::value<std::vector<int>>(&energies)->multitoken()->default_value(defaultEnergies, "30 50 70 100 150 300"), "vector of sample energies")
    ("genThresh", po::value<double>(&genThresh)->default_value(3.0), "Generator Particle Threshold")
    ("energyFracCut", po::value<double>(&eFC)->default_value(99.0), "maximum dE/E")
    ("bitShortThresh", po::value<double>(&bitShortThresh)->default_value(10.0), "Bit Short Threshold")
    ("bitLongThresh", po::value<double>(&bitLongthresh)->default_value(10.0), "Bit Long Threshold")
    ("bitOffset", po::value<double>(&bitOffset)->default_value(10.1), "Bit Offset")
    ("bitSlope", po::value<double>(&bitSlope)->default_value(100.2), "Bit Slope")
    ("coneSize", po::value<double>(&coneSize)->default_value(0.20), "matching coneSize")
    ("mode", po::value<int>(&mode)->default_value(0), "output optimization mode flag")
    ("runZee", po::value<bool>(&runZee)->default_value(false), "run on Zee sample")
    ("pileup", po::value<bool>(&pileup)->default_value(false), "run on pileup samples");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help"))
  {
    std::cout << desc << std::endl;
    return 1;
  }

  // std::cout << ntuplePath << std::endl;
  std::map<int, std::string> modeFlags = {{0, ""}, {1, "short,"}, {2, "long,"}, {3, "offset,"}, {4, "slope,"}};

  if (runZee)
    {
      fgBitGenAnalysis(ntuplePath);
      fgBitAnalysis(ntuplePath);
    }
  // else
  // {
  //   for (int energy : energies)
  //   {
  //     std::cout << modeFlags[mode]
  //               << energy << ","
  //               << genThresh << ","
  //               << eFC << ","
  //               << bitShortThresh << ","
  //               << bitLongthresh << ","
  //               << bitOffset << ","
  //               << bitSlope << ","
  //               << coneSize << ",";
  //     energyRatioAnalysis(ntuplePath, energy, genThresh, eFC, bitShortThresh, bitLongthresh, bitOffset, bitSlope, coneSize, pileup);
  //   }
  // }
  // energyRatioAnalysis(ntuplePath, genThresh, bitShortThresh, bitLongthresh, bitOffset, bitSlope, coneSize);

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

// void energyRatioAnalysis(const std::string &inputFileDirectory, int energy, double genThresh, double energyFracCut, double shortThresh, double longThresh, double offset, double slope, double coneSize, bool pileup)
// // void energyRatioAnalysis(const std::string &inputFileDirectory, double genThresh, double shortThresh, double longThresh, double offset, double slope, double coneSize)
// {

//   std::string inputFileZ(inputFileDirectory);
//   std::string inputFilePi(inputFileDirectory);
//   std::string outputFilename("");

//   inputFileZ += "/singleElectronE";
//   inputFileZ += std::to_string(energy);
//   if (pileup)
//     inputFileZ += "PU";
//   inputFileZ += "Run3.root";

//   // std::string inputFileQ(inputFileDirectory);
//   // inputFileQ += "/L1Ntuple_QCDPU.root";

//   // inputFilePi += "/L1Ntuple_singlepi.root";
//   inputFilePi += "/singlePionE";
//   inputFilePi += std::to_string(energy);
//   if (pileup)
//     inputFilePi += "PU";
//   inputFilePi += "Run3.root";

//   outputFilename += "energyPlotsRun3/energyRatioPlots";
//   outputFilename += std::to_string(energy);
//   if (pileup)
//     outputFilename += "PU";
//   outputFilename += ".root";

//   // std::cout << inputFileZ << std::endl;
//   // std::cout << inputFilePi << std::endl;
//   // std::cout << outputFilename << std::endl;

//   TFile *outfile = TFile::Open(outputFilename.c_str(), "recreate");

//   // std::cout << "loading up the TChain" << std::endl;

//   // Z->ee Trees
//   // TChain *treeL1emuZ = new TChain("l1UpgradeEmuTree/L1UpgradeTree");
//   // treeL1emuZ->Add(inputFileZ.c_str());

//   TChain *treeL1TPemuZ = new TChain("l1CaloTowerEmuTree/L1CaloTowerTree");
//   treeL1TPemuZ->Add(inputFileZ.c_str());

//   TChain *treeL1GenZ = new TChain("l1GeneratorTree/L1GenTree");
//   treeL1GenZ->Add(inputFileZ.c_str());

//   // Single Pi Trees
//   // TChain *treeL1emuPi = new TChain("l1UpgradeEmuTree/L1UpgradeTree");
//   // treeL1emuPi->Add(inputFilePi.c_str());

//   TChain *treeL1TPemuPi = new TChain("l1CaloTowerEmuTree/L1CaloTowerTree");
//   treeL1TPemuPi->Add(inputFilePi.c_str());

//   TChain *treeL1GenPi = new TChain("l1GeneratorTree/L1GenTree");
//   treeL1GenPi->Add(inputFilePi.c_str());

//   // Z->ee Branches
//   // L1Analysis::L1AnalysisL1UpgradeDataFormat *l1emuZ_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
//   // treeL1emuZ->SetBranchAddress("L1Upgrade", &l1emuZ_);

//   L1Analysis::L1AnalysisCaloTPDataFormat *l1TPemuZ_ = new L1Analysis::L1AnalysisCaloTPDataFormat();
//   treeL1TPemuZ->SetBranchAddress("CaloTP", &l1TPemuZ_);

//   L1Analysis::L1AnalysisGeneratorDataFormat *l1GenZ_ = new L1Analysis::L1AnalysisGeneratorDataFormat();
//   treeL1GenZ->SetBranchAddress("Generator", &l1GenZ_);

//   // Single Pi Branches
//   // L1Analysis::L1AnalysisL1UpgradeDataFormat *l1emuPi_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
//   // treeL1emuPi->SetBranchAddress("L1Upgrade", &l1emuPi_);

//   L1Analysis::L1AnalysisCaloTPDataFormat *l1TPemuPi_ = new L1Analysis::L1AnalysisCaloTPDataFormat();
//   treeL1TPemuPi->SetBranchAddress("CaloTP", &l1TPemuPi_);

//   L1Analysis::L1AnalysisGeneratorDataFormat *l1GenPi_ = new L1Analysis::L1AnalysisGeneratorDataFormat();
//   treeL1GenPi->SetBranchAddress("Generator", &l1GenPi_);

//   TH1F *minGenETPdR = new TH1F("minTPGenE", "", 50, 0, 5);
//   TH1F *minGenPiTPdR = new TH1F("minTPGenPi", "", 50, 0, 5);

//   TH1F *minGenETPdEt = new TH1F("deltaEtE", "", 70, -2, 5);
//   TH1F *minGenPiTPdEt = new TH1F("deltaEtPi", "", 70, -2, 5);

//   TH1F *genEleEt = new TH1F("genEleEt", "", 150, 0, 150);
//   TH1F *genPionEt = new TH1F("genPionEt", "", 150, 0 , 150);

//   TH1F *shortEnergyE = new TH1F("shortEnE", "", 125, -10, 240);
//   TH1F *longEnergyE = new TH1F("longEnE", "", 125, -10, 240);
//   TH1F *energyRatioE = new TH1F("nrgRatioE", "", 40, 0, 20);

//   TH1F *singleTPshortEnergyE = new TH1F("sTPshortE", "", 140, 0, 140);
//   TH1F *singleTPlongEnergyE = new TH1F("sTPlongE", "", 140, 0, 140);

//   TH1F *shortEnergyPi = new TH1F("shortEnPi", "", 125, -10, 240);
//   TH1F *longEnergyPi = new TH1F("longEnPi", "", 125, -10, 240);
//   TH1F *energyRatioPi = new TH1F("nrgRatioPi", "", 40, 0, 20);

//   TH1F *singleTPshortEnergyPi = new TH1F("sTPshortPi", "", 140, 0, 140);
//   TH1F *singleTPlongEnergyPi = new TH1F("sTPlongPi", "", 140, 0, 140);

//   TH1F *hcalTPMultipilictyE = new TH1F("hcalTPMulE", "", 20, 0, 20);
//   TH1F *hcalTPMultipilictyPi = new TH1F("hcalTPMulPi", "", 20, 0, 20);

//   TH1F *singleExprE = new TH1F("singleTPexprE", "", 250, 0, 1000);
//   singleExprE->SetCanExtend(TH1::kAllAxes);
//   TH1F *singleExprPi = new TH1F("singleTPexprPi", "", 250, 0, 1000);
//   singleExprPi->SetCanExtend(TH1::kAllAxes);

//   TH1F *singleTPRatioE = new TH1F("singleTPRatioE", "", 40, 0, 20);
//   TH1F *singleTPRatioPi = new TH1F("singleTPRatioPi", "", 40, 0, 20);

//   TH1F *eff_num = new TH1F("numEff", "", 50, 0, 50);
//   TH1F *eff_den = new TH1F("denEff", "", 50, 0, 50);

//   TH1F *mistag_num = new TH1F("numMiss", "", 50, 0, 50);
//   TH1F *mistag_den = new TH1F("denMiss", "", 50, 0, 50);

//   TH1F *matched_numE = new TH1F("matchNumE", "", 50, 0, 50);
//   TH1F *matched_denE = new TH1F("matchDenE", "", 50, 0, 50);

//   TH1F *matched_numPi = new TH1F("matchNumPi", "", 50, 0, 50);
//   TH1F *matched_denPi = new TH1F("matchDenPi", "", 50, 0, 50);

//   TNtuple *ntuple = new TNtuple("ntupledData", "Energy Data Ntuple", "label:short:long");
//   TNtuple *ntupleSum = new TNtuple("ntupledSum", "Energy Sum Data Ntuple", "label:shortSum:longSum");

//   Long64_t nentries;
//   nentries = treeL1GenZ->GetEntries();

//   int ntupleFineGrain = 0;
//   int newFineGrain = 0;
//   int nParticles = 0;
//   int nMatchedWithTP = 0;
//   int nMatchedWithBit = 0;

//   for (Long64_t jentry = 0; jentry < nentries; jentry++)
//   {
//     // if ((jentry % 1000) == 0)
//     // std::cout << "Done " << jentry << " events of Zee " << nentries << std::endl;

//     treeL1TPemuZ->GetEntry(jentry);
//     // treeL1emuZ->GetEntry(jentry);
//     treeL1GenZ->GetEntry(jentry);

//     for (int g = 0; g < l1GenZ_->nPart; g++)
//     {
//       if (fabs(l1GenZ_->partEta[g]) < 2.9)
//         continue;
//       if (fabs(l1GenZ_->partEta[g]) > 5.1)
//         continue;
//       if (abs(l1GenZ_->partId[g]) != 11) // Electron
//         continue;
//       // if (abs(l1GenZ_->partParent[g]) != 23) // Parent is Z
//       //   continue;
//       // if (!l1GenZ_->partFromHard[g]) // From Hard Process
//       //   continue;
//       if (l1GenZ_->partE[g] < genThresh)
//         continue;

//       TLorentzVector genV;
//       genV.SetPtEtaPhiE(l1GenZ_->partPt[g], l1GenZ_->partEta[g], l1GenZ_->partPhi[g], l1GenZ_->partE[g]);

//       genEleEt->Fill(genV.Et());
//       matched_denE->Fill(l1GenZ_->partPt[g]);
//       nParticles++; // At this point we are sure this is an electron..
//       int matchedTPs = 0;
//       float shortSum = 0.0;
//       float longSum = 0.0;
//       double minGenTPdR = 999.;
//       bool emfbit = false;

//       //      double minGenTPdEt = 999.;

//       for (int t = 0; t < l1TPemuZ_->nHCALTP; t++)
//       {
//         if (l1TPemuZ_->hcalTPet[t] <= 0.5 or abs(l1TPemuZ_->hcalTPieta[t]) <= 28)
//           continue;

//         double dR = deltaR(l1GenZ_->partEta[g], l1GenZ_->partPhi[g], l1TPemuZ_->hcalTPeta[t], phiVal(l1TPemuZ_->hcalTPCaliphi[t]));

//         if (dR <= coneSize)
//         {
//           minGenETPdEt->Fill((l1TPemuZ_->hcalTPet[t] - genV.Et()) / genV.Et());
//         }

//         if (fabs((l1TPemuZ_->hcalTPet[t] - genV.Et()) / genV.Et()) > energyFracCut)
//           continue;
        
//         if (dR <= coneSize)
//         {
//           matchedTPs++;
//           shortSum += l1TPemuZ_->hcalTPshortFiberE[t];
//           longSum += l1TPemuZ_->hcalTPlongFiberE[t];


//           singleTPshortEnergyE->Fill(l1TPemuZ_->hcalTPshortFiberE[t]);
//           singleTPlongEnergyE->Fill(l1TPemuZ_->hcalTPlongFiberE[t]);
//           ntuple->Fill(1.00, l1TPemuZ_->hcalTPshortFiberE[t], l1TPemuZ_->hcalTPlongFiberE[t]);
//           if (l1TPemuZ_->hcalTPfineGrain[t])
//             ntupleFineGrain++;

//           if (l1TPemuZ_->hcalTPshortFiberE[t] > shortThresh and l1TPemuZ_->hcalTPlongFiberE[t] > longThresh)
//           {
//             singleExprE->Fill((l1TPemuZ_->hcalTPlongFiberE[t] - offset) * slope - l1TPemuZ_->hcalTPshortFiberE[t]);
//             singleTPRatioE->Fill((l1TPemuZ_->hcalTPlongFiberE[t] / l1TPemuZ_->hcalTPshortFiberE[t]));
//             if (l1TPemuZ_->hcalTPshortFiberE[t] < (l1TPemuZ_->hcalTPlongFiberE[t] - offset) * slope) // fineGrainBit is set in this case!!
//             {
//               emfbit = true; // once this is true it should remain true..
//               newFineGrain++;
//             }
//           }
//           // if (fabs((l1TPemuZ_->hcalTPet[t] - genV.Et()) / genV.Et()) < fabs(minGenTPdEt))
//           // minGenTPdEt = (l1TPemuZ_->hcalTPet[t] - genV.Et()) / genV.Et();
//         }
//         if (dR < minGenTPdR)
//         {
//           minGenTPdR = dR;
//           //          minGenTPdEt = (l1TPemuZ_->hcalTPet[t] - genV.Et()) / genV.Et();
//         }

//       } // End HCAL TP loop

//       if (emfbit) // was the bit set by any TP?
//       {
//         nMatchedWithBit++; // increase the number of particles matched with a bit !!!
//         eff_num->Fill(l1GenZ_->partPt[g]);
//       }

//       //      minGenETPdEt->Fill(minGenTPdEt);
//       minGenETPdR->Fill(minGenTPdR);
//       if (matchedTPs > 0)
//       {
//         eff_den->Fill(l1GenZ_->partPt[g]);
//         matched_numE->Fill(l1GenZ_->partPt[g]);
//         nMatchedWithTP++;
//         energyRatioE->Fill(longSum / shortSum);
//         shortEnergyE->Fill(shortSum);
//         longEnergyE->Fill(longSum);
//         ntupleSum->Fill(1.00, shortSum, longSum);
//       }
//       hcalTPMultipilictyE->Fill(matchedTPs);

//     } // End Gen Particle Loop
//   }   // End single Electron Loop

//   // std::cout << "ntuple FineGrains : " << ntupleFineGrain << std::endl;
//   // std::cout << "electron calculated finegrains: " << newFineGrain << std::endl;
//   // std::cout << "number of generator electrons : " << nParticles << std::endl;
//   // std::cout << "number of electrons matched with a set bit : " << nMatchedWithBit << std::endl;
//   // std::cout << " Electrons Tagged : " <<  std::fixed << std::setprecision(2) << (double(nMatchedWithBit)/nMatchedWithTP) << std::endl;
//   std::cout << std::fixed << std::setprecision(2) << (double(nMatchedWithBit) / nMatchedWithTP) << ",";
//   std::cout << std::fixed << std::setprecision(2) << (double(nMatchedWithBit) / nParticles) << ",";

//   // reset counters for pion sample..
//   nentries = treeL1GenPi->GetEntries();
//   ntupleFineGrain = 0;
//   newFineGrain = 0;
//   nParticles = 0;
//   nMatchedWithTP = 0;
//   nMatchedWithBit = 0;

//   // run over Pion Sample
//   for (Long64_t jentry = 0; jentry < nentries; jentry++)
//   {
//     // if ((jentry % 1000) == 0)
//     // std::cout << "Done " << jentry << " events of Single Pi " << nentries << std::endl;

//     treeL1TPemuPi->GetEntry(jentry);
//     // treeL1emuZ->GetEntry(jentry);
//     treeL1GenPi->GetEntry(jentry);

//     for (int g = 0; g < l1GenPi_->nPart; g++)
//     {
//       if (fabs(l1GenPi_->partEta[g]) < 2.9)
//         continue;
//       if (fabs(l1GenPi_->partEta[g]) > 5.1)
//         continue;
//       if (abs(l1GenPi_->partId[g]) != 211) // Pion
//         continue;
//       // if (!l1GenPi_->partFromHard[g]) // From Hard Process
//       // continue;
//       if (l1GenPi_->partE[g] < genThresh)
//         continue;

//       TLorentzVector genV;
//       genV.SetPtEtaPhiE(l1GenPi_->partPt[g], l1GenPi_->partEta[g], l1GenPi_->partPhi[g], l1GenPi_->partE[g]);
      
//       genPionEt->Fill(genV.Et());
//       matched_denPi->Fill(l1GenPi_->partPt[g]);
//       nParticles++; // generator pion "Accepted"

//       int matchedTPs = 0;
//       float shortSum = 0.0;
//       float longSum = 0.0;
//       double minGenTPdR = 999.;
//       bool emfbit = false;

//       //      double minGenTPdEt = 999.;

//       for (int t = 0; t < l1TPemuPi_->nHCALTP; t++)
//       {
//         if (l1TPemuPi_->hcalTPet[t] <= 0.5 or abs(l1TPemuPi_->hcalTPieta[t]) <= 28)
//           continue;

//         double dR = deltaR(l1GenPi_->partEta[g], l1GenPi_->partPhi[g], l1TPemuPi_->hcalTPeta[t], phiVal(l1TPemuPi_->hcalTPCaliphi[t]));

//         if (dR <= coneSize)
//         {
//           minGenPiTPdEt->Fill((l1TPemuPi_->hcalTPet[t] - genV.Et()) / genV.Et());          
//         }

//         if (fabs((l1TPemuPi_->hcalTPet[t] - genV.Et()) / genV.Et()) > energyFracCut)
//           continue;

//         if (dR <= coneSize)
//         {
//           matchedTPs++;
//           shortSum += l1TPemuPi_->hcalTPshortFiberE[t];
//           longSum += l1TPemuPi_->hcalTPlongFiberE[t];
//           singleTPshortEnergyPi->Fill(l1TPemuPi_->hcalTPshortFiberE[t]);
//           singleTPlongEnergyPi->Fill(l1TPemuPi_->hcalTPlongFiberE[t]);
//           ntuple->Fill(0.00, l1TPemuPi_->hcalTPshortFiberE[t], l1TPemuPi_->hcalTPlongFiberE[t]);
//           if (l1TPemuPi_->hcalTPfineGrain[t])
//             ntupleFineGrain++;

//           if (l1TPemuPi_->hcalTPshortFiberE[t] > shortThresh and l1TPemuPi_->hcalTPlongFiberE[t] > longThresh)
//           {
//             singleExprPi->Fill((l1TPemuPi_->hcalTPlongFiberE[t] - offset) * slope - l1TPemuPi_->hcalTPshortFiberE[t]);
//             singleTPRatioPi->Fill((l1TPemuPi_->hcalTPlongFiberE[t] / l1TPemuPi_->hcalTPshortFiberE[t]));

//             if (l1TPemuPi_->hcalTPshortFiberE[t] < (l1TPemuPi_->hcalTPlongFiberE[t] - offset) * slope)
//             {
//               emfbit = true; // once this is true it should remain true..
//               newFineGrain++;
//             }
//           }
//           // if (fabs((l1TPemuPi_->hcalTPet[t] - genV.Et()) / genV.Et()) < fabs(minGenTPdEt))
//           // minGenTPdEt = (l1TPemuPi_->hcalTPet[t] - genV.Et()) / genV.Et();
//           // minGenPiTPdEt->Fill(((l1TPemuPi_->hcalTPet[t] - genV.Et()) / genV.Et()));
//         }
//         if (dR < minGenTPdR)
//         {
//           minGenTPdR = dR;
//           //          minGenTPdEt = (l1TPemuPi_->hcalTPet[t] - genV.Et()) / genV.Et();
//         }

//       } // End Loop over HCAL TPs

//       if (emfbit)
//       {
//         nMatchedWithBit++;
//         mistag_num->Fill(l1GenPi_->partPt[g]);
//       }

//       //      minGenPiTPdEt->Fill(minGenTPdEt);

//       minGenPiTPdR->Fill(minGenTPdR);

//       if (matchedTPs > 0)
//       {
//         mistag_den->Fill(l1GenPi_->partPt[g]);
//         matched_numPi->Fill(l1GenPi_->partPt[g]);
//         nMatchedWithTP++;
//         energyRatioPi->Fill(longSum / shortSum);
//         shortEnergyPi->Fill(shortSum);
//         longEnergyPi->Fill(longSum);
//         ntupleSum->Fill(0.00, shortSum, longSum);
//       }

//       hcalTPMultipilictyPi->Fill(matchedTPs);

//     } // End Loop over Gen particles
//   }   // End Single Pi Loop

//   // std::cout << "ntuple FineGrains : " << ntupleFineGrain << std::endl;
//   // std::cout << "pion calculated finegrains: " << newFineGrain << std::endl;
//   // std::cout << "number of generator pions : " << nParticles << std::endl;
//   // std::cout << "number of pions matched with a set bit : " << nMatchedWithBit << std::endl;
//   // std::cout << " Pions Mistagged: " <<  std::fixed << std::setprecision(2) << (double(nMatchedWithBit)/nMatchedWithTP) << std::endl;
//   std::cout << std::fixed << std::setprecision(2) << (double(nMatchedWithBit) / nMatchedWithTP) << ",";
//   std::cout << std::fixed << std::setprecision(2) << (double(nMatchedWithBit) / nParticles) << std::endl;

//   outfile->cd();

//   minGenETPdR->Write();
//   minGenPiTPdR->Write();

//   shortEnergyE->Write();
//   longEnergyE->Write();
//   energyRatioE->Write();

//   shortEnergyPi->Write();
//   longEnergyPi->Write();
//   energyRatioPi->Write();

//   hcalTPMultipilictyE->Write();
//   hcalTPMultipilictyPi->Write();

//   singleExprE->Write();
//   singleExprPi->Write();

//   eff_num->Write();
//   eff_den->Write();

//   mistag_num->Write();
//   mistag_den->Write();

//   singleTPshortEnergyE->Write();
//   singleTPlongEnergyE->Write();
//   singleTPshortEnergyPi->Write();
//   singleTPlongEnergyPi->Write();
//   singleTPRatioE->Write();
//   singleTPRatioPi->Write();
//   matched_numE->Write();
//   matched_denE->Write();
//   matched_numPi->Write();
//   matched_denPi->Write();
//   minGenETPdEt->Write();
//   minGenPiTPdEt->Write();
//   ntuple->Write();
//   ntupleSum->Write();

//   genEleEt->Write();
//   genPionEt->Write();

//   outfile->Close();
// }

void fgBitGenAnalysis(const std::string &inputFileDirectory)
{
  std::string inputFileZ(inputFileDirectory);
  std::string outputFilename("l1egfinegrainGen.root");

  inputFileZ += "L1Ntuple_Zee.root";
  TFile *outfile = TFile::Open(outputFilename.c_str(), "recreate");

  TChain *treeL1TPemuZ = new TChain("l1CaloTowerEmuTree/L1CaloTowerTree");
  treeL1TPemuZ->Add(inputFileZ.c_str());

  TChain *treeL1GenZ = new TChain("l1GeneratorTree/L1GenTree");
  treeL1GenZ->Add(inputFileZ.c_str());

  L1Analysis::L1AnalysisCaloTPDataFormat *l1TPemuZ_ = new L1Analysis::L1AnalysisCaloTPDataFormat();
  treeL1TPemuZ->SetBranchAddress("CaloTP", &l1TPemuZ_);

  L1Analysis::L1AnalysisGeneratorDataFormat *l1GenZ_ = new L1Analysis::L1AnalysisGeneratorDataFormat();
  treeL1GenZ->SetBranchAddress("Generator", &l1GenZ_);

  TH1F *minGenETPdR = new TH1F("minTPGenZee", "", 60, 0, 6);
  // TH1F *shortEnergyE = new TH1F("shortEnZee", "", 125, -10, 240); 
  // TH1F *longEnergyE = new TH1F("longEnZee", "", 125, -10, 240);
  // TH1F *energyRatioE = new TH1F("nrgRatioZee", "", 40, 0, 20);
  TH1F *hcalTPMultipilictyE = new TH1F("hcalTPMulZee", "", 20, 0, 20);
  TH1F *eff_num = new TH1F("numEffZee", "", 150, 0, 150);
  TH1F *eff_den = new TH1F("denEffZee", "", 150, 0, 150);

  // TH1F *eff_numE = new TH1F("numEffZeeE", "", 100, 0, 1000);
  // TH1F *eff_denE = new TH1F("denEffZeeE", "", 100, 0, 1000);

  Long64_t nentries;
  nentries = treeL1GenZ->GetEntries();

  int ntupleFineGrain = 0;
  // int newFineGrain = 0;	    
  // int nParticles = 0;
  int nMatchedWithBit = 0;
  int nMatchedWithTP = 0;

  for (Long64_t jentry = 0; jentry < nentries; jentry++)
  {
    treeL1TPemuZ->GetEntry(jentry);
    treeL1GenZ->GetEntry(jentry);

    for (int g = 0; g < l1GenZ_->nPart; g++)
    {
      if (fabs(l1GenZ_->partEta[g]) < 2.9)
        continue;
      if (fabs(l1GenZ_->partEta[g]) > 5.1)
        continue;
      if (abs(l1GenZ_->partId[g]) != 11) // Electron
        continue;
      if (abs(l1GenZ_->partParent[g]) != 23) // Parent is Z
        continue;
      if (!l1GenZ_->partIsHard[g]) // From Hard Process
        continue;
      if (l1GenZ_->partE[g] < 3.0)
        continue;

      eff_den->Fill(l1GenZ_->partPt[g]);
      int matchedTPs = 0;
      double minGenTPdR = 999.;
      bool emfbit = false;

      for (int t = 0; t < l1TPemuZ_->nHCALTP; t++)
      {
        if (l1TPemuZ_->hcalTPet[t] <= 0.5 or abs(l1TPemuZ_->hcalTPieta[t]) <= 28)
          continue;

        double dR = deltaR(l1GenZ_->partEta[g], l1GenZ_->partPhi[g], etaVal(l1TPemuZ_->hcalTPieta[t]), phiVal(l1TPemuZ_->hcalTPCaliphi[t]));

        if (dR <= 0.2)
        {
          matchedTPs++;
          // shortSum += l1TPemuZ_->hcalTPshortFiberE[t];
          // longSum += l1TPemuZ_->hcalTPlongFiberE[t];

          if (l1TPemuZ_->hcalTPfineGrain[t])
            {
	      ntupleFineGrain++;
	      emfbit = true;
	    }

          // if (l1TPemuZ_->hcalTPshortFiberE[t] > shortThresh and l1TPemuZ_->hcalTPlongFiberE[t] > longThresh)
          // {

          //   if (l1TPemuZ_->hcalTPshortFiberE[t] < (l1TPemuZ_->hcalTPlongFiberE[t] - offset) * slope) // fineGrainBit is set in this case!!
          //   {
          //     emfbit = true; // once this is true it should remain true..
          //     newFineGrain++;
          //   }
          // }
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
        // eff_numE->Fill(l1GenZ_->partE[g]);
      }

      minGenETPdR->Fill(minGenTPdR);
      if (matchedTPs > 0)
      {
        nMatchedWithTP++;
        // energyRatioE->Fill(longSum / shortSum);
        // shortEnergyE->Fill(shortSum);
        // longEnergyE->Fill(longSum);
      }
      hcalTPMultipilictyE->Fill(matchedTPs);

    } // End Gen Particle Loop
  }   // End z->ee loop

  // std::cout << "Z->ee sample results:" << std::endl;
  // std::cout << "ntuple FineGrains : " << ntupleFineGrain << std::endl;
  // std::cout << "electron calculated finegrains: " << newFineGrain << std::endl;
  // std::cout << "number of generator electrons : " << nParticles << std::endl;
  // std::cout << "number of electrons matched with a set bit : " << nMatchedWithBit << std::endl;
  // std::cout << "Zee Electrons Tagged : " <<  std::fixed << std::setprecision(2) << (double(nMatchedWithBit)/nMatchedWithTP) << std::endl;

  outfile->cd();
  minGenETPdR->Write();
  // shortEnergyE->Write();
  // longEnergyE->Write();
  // energyRatioE->Write();
  hcalTPMultipilictyE->Write();
  eff_num->Write();
  // eff_numE->Write();
  eff_den->Write();
  // eff_denE->Write();
  outfile->Close();
}
 
void fgBitAnalysis(const std::string &inputFileDirectory)
{
  std::string inputFileZ(inputFileDirectory);
  std::string inputFileQ(inputFileDirectory);
  inputFileZ += "L1Ntuple_Zee.root";
  inputFileQ += "L1Ntuple_Qcd.root";

  std::string outputFilename = "l1egfinegrain.root";

  TFile *outfile = TFile::Open(outputFilename.c_str(), "recreate");
  
  std::cout << "loading up the TChain" << std::endl;

  TChain *treeL1emuZ = new TChain("l1UpgradeEmuTree/L1UpgradeTree");
  treeL1emuZ->Add(inputFileZ.c_str());

  TChain *treeL1TPemuZ = new TChain("l1CaloTowerEmuTree/L1CaloTowerTree");
  treeL1TPemuZ->Add(inputFileZ.c_str());

  TChain *treeL1GenZ = new TChain("l1GeneratorTree/L1GenTree");
  treeL1GenZ->Add(inputFileZ.c_str());

  TChain *treeL1emuQ = new TChain("l1UpgradeEmuTree/L1UpgradeTree");
  treeL1emuZ->Add(inputFileQ.c_str());

  TChain *treeL1TPemuQ = new TChain("l1CaloTowerEmuTree/L1CaloTowerTree");
  treeL1TPemuZ->Add(inputFileQ.c_str());

  // TChain *treeL1GenZ = new TChain("l1GeneratorTree/L1GenTree");
  // treeL1GenZ->Add(inputFileZ.c_str());

  L1Analysis::L1AnalysisL1UpgradeDataFormat *l1emuZ_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  treeL1emuZ->SetBranchAddress("L1Upgrade", &l1emuZ_);

  L1Analysis::L1AnalysisCaloTPDataFormat *l1TPemuZ_ = new L1Analysis::L1AnalysisCaloTPDataFormat();
  treeL1TPemuZ->SetBranchAddress("CaloTP", &l1TPemuZ_);

  // add generator level tree to check mc truth
  L1Analysis::L1AnalysisGeneratorDataFormat *l1GenZ_ = new L1Analysis::L1AnalysisGeneratorDataFormat();
  treeL1GenZ->SetBranchAddress("Generator", &l1GenZ_);

  // QCD Branches
  L1Analysis::L1AnalysisL1UpgradeDataFormat *l1emuQ_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  treeL1emuQ->SetBranchAddress("L1Upgrade", &l1emuQ_);

  L1Analysis::L1AnalysisCaloTPDataFormat *l1TPemuQ_ = new L1Analysis::L1AnalysisCaloTPDataFormat();
  treeL1TPemuQ->SetBranchAddress("CaloTP", &l1TPemuQ_);

  Long64_t nentries;
  nentries = treeL1emuZ->GetEntries();

  TH1F *minDrL1eHcalTP = new TH1F("minddrL1HcalTP", ";min dR(l1eg,hcalTP);", 100, 0, 10);
  TH1F *minDrL1eGene = new TH1F("mindrl1eggen", ";;min dR(l1eg,genEle)", 60, 0, 6);

  TH1F *eff_num = new TH1F("numEff", "", 100, 0, 200);
  TH1F *eff_den = new TH1F("denEff", "", 100, 0, 200);

  TH1F *purity_num = new TH1F("numPur", "", 100, 0, 200);
  TH1F *purity_den = new TH1F("denPur", "", 100, 0, 200);

  // TH1F *mistag_num = new TH1F("numMiss", "", 100, 0 ,200);
  // TH1F *mistag_den = new TH1F("denMiss", "", 100, 0 , 200);
  
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
      if (abs(l1emuZ_->egTowerIEta[i]) <= 28 or l1emuZ_->egEt[i] <= 3.0)
        continue;

      double minL1GenDR = 999;
      int genIdx = -1;
      // Loop over gen particles to find a gen level electron
      // with a Z parent matched to the L1 particle
      for (int g = 0; g < l1GenZ_->nPart; g++)
      {
	if (fabs(l1GenZ_->partEta[g]) < 2.9)
	  continue;
	if (fabs(l1GenZ_->partEta[g]) > 5.1)
	  continue;
        if (abs(l1GenZ_->partId[g]) != 11) // make sure it's an electron..
          continue;
        if (abs(l1GenZ_->partParent[g]) != 23) // make sure parent is Z
          continue;
        if (!l1GenZ_->partIsHard[g]) // from hard scatter process (only for efficiencies?)
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
        if (l1TPemuZ_->hcalTPet[j] <= 5.0 or abs(l1TPemuZ_->hcalTPieta[j]) <= 28)
          continue;

        double dR = deltaR(l1emuZ_->egEta[i], l1emuZ_->egPhi[i], etaVal(l1TPemuZ_->hcalTPieta[j] ), phiVal(l1TPemuZ_->hcalTPCaliphi[j] ) );
        // double dR = deltaR(l1emuZ_->egEta[i], l1emuZ_->egPhi[i], hcalEta, hcalPhi);
        if (dR <= 0.2)
        {
          emfbit = emfbit or l1TPemuZ_->hcalTPfineGrain[j];
        }
        if (dR < minL1hcalDR)
        {
          minL1hcalDR = dR;
        }
      }

      if (genIdx >= 0) // if the particle is matched to a gen particle??
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
        if (l1TPemuQ_->hcalTPet[j] <= 5.0 or abs(l1TPemuQ_->hcalTPieta[j]) <= 28)
          continue;

        double dR = deltaR(l1emuQ_->egEta[i], l1emuQ_->egPhi[i], l1TPemuQ_->hcalTPieta[j], l1TPemuQ_->hcalTPCaliphi[j]);
        if (dR < hcal_l1_dR)
        {
          emfbitQ = emfbitQ or l1TPemuQ_->hcalTPfineGrain[j]; // if any are true, emfbit stays true
        }
      }
      if (emfbitQ)
      {
        fake_num->Fill(l1emuQ_->egEt[i]);
      }
    } // finish qcd l1egamma loop
  } // finish qcd sample loop

  outfile->cd();
  minDrL1eHcalTP->Write();
  minDrL1eGene->Write();
  eff_num->Write();
  eff_den->Write();
  purity_num->Write();
  purity_den->Write();
  // fake_num->Write();
  // fake_den->Write();
}
