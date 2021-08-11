
#include "TChain.h"
#include "TString.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include <iostream>
#include <set>
#include <array>
#include <algorithm>
#include <cmath>

/* map centrality to Ncoll for MC event weight */
float findNcoll(int);

void mc(TString fname = "root/HiForestAOD_ZS_8-2.root",
        TString outName = "mcmatch.root") {
  TChain *tch = new TChain("hiEvtAnalyzer/HiTree");
  tch->Add(fname);
  TChain *tdata = new TChain("ggHiNtuplizer/EventTree");
  tdata->Add(fname);

  TTreeReader reader(tch);
  TTreeReaderValue<Int_t> centrality(reader, "hiBin");
  TTreeReader dreader(tdata);
  TTreeReaderValue<Int_t> nPho(dreader, "nPho");
  TTreeReaderArray<Int_t> mcPID(dreader, "mcPID");
  TTreeReaderArray<Int_t> genId(dreader, "pho_genMatchedIndex");
  // uncorrected supercluster energy
  TTreeReaderArray<Float_t> recoE(dreader, "phoSCRawE");
  // gen energy
  TTreeReaderArray<Float_t> mcE(dreader, "mcE");
  // reading pT into a vector since we need its length
  TTreeReaderValue<std::vector<Float_t> > mcPt(dreader, "mcPt");
  // For spike rejection
  TTreeReaderArray<Float_t> phoSigmaIEtaIEta_2012(dreader, "phoSigmaIEtaIEta_2012");
  TTreeReaderArray<Float_t> pho_swissCrx(dreader, "pho_swissCrx");
  TTreeReaderArray<Float_t> pho_seedTime(dreader, "pho_seedTime");
  auto passedPhoSpikeRejection = [&] (int i) {
    return ( phoSigmaIEtaIEta_2012[i] > 0.002 &&
             pho_swissCrx[i] < 0.9 && std::fabs(pho_seedTime[i]) < 3 );
  };
  // For HI18HEM failure
  TTreeReaderArray<Float_t> phoSCEta(dreader, "phoSCEta");
  TTreeReaderArray<Float_t> phoSCPhi(dreader, "phoSCPhi");
  auto passedHI18HEMfailurePho = [&] (int i) {
    return !(phoSCEta[i] < -1.39 && phoSCPhi[i] < -0.9 && phoSCPhi[i] > -1.6);
  };

  // For photon ID selections
  TTreeReaderArray<Float_t> phoHoverE(dreader, "phoHoverE");
  TTreeReaderArray<Float_t> mcCalIsoDR04(dreader, "mcCalIsoDR04");

  TFile fout(outName, "recreate");
  // histograms for efficiency study
  std::array<float, 6> ptbins = {15, 30, 40, 50, 80, 120};
  // total number of gen photon
  TH1F total("tot", "total;pT/GeV", ptbins.size() - 1, ptbins.data());
  // gen photon that's matched to a RECO photon
  TH1F matched("mat", "matched;pT/GeV", ptbins.size() - 1, ptbins.data());
  // histograms for energy scale study
  std::array<TH1F *, ptbins.size()> eScale;
  for (int i = 0; i < ptbins.size() - 1; ++i) {
    eScale[i] =
        new TH1F(TString::Format("esc%.0f_%.0f", ptbins[i], ptbins[i + 1]),
                 "Energy scale", 30, 0.4, 1.6);
  }
  while (reader.Next()) {
    double centWeight = findNcoll(*centrality);
    dreader.Next();
    // gen matching
    std::set<int> matchedGamma;
    for (int iPho = 0; iPho < *nPho; ++iPho) {
      int genMatchedIndex = genId[iPho];
      bool isMatched = (genMatchedIndex >= 0);
      bool isMatched2GenPhoton = (isMatched && mcPID[genMatchedIndex] == 22);
      // store matched mc photon
      if (isMatched2GenPhoton) {
        matchedGamma.insert(genMatchedIndex);
        float genPt = (*mcPt)[genMatchedIndex];
        // skip photons with out-of-range pT
        if (genPt < ptbins.front() || genPt > ptbins.back()) {
          continue;
          // reject spike photon
        } else if (!passedPhoSpikeRejection(iPho)) {
          continue;
          // reject photons that failed HEM modules
        } else if (!passedHI18HEMfailurePho(genMatchedIndex)) {
          continue;
          // photon selection
        } else if (phoHoverE[iPho] > 0.1) {
          continue;
          // gen photon selection
        } else if (mcCalIsoDR04[genMatchedIndex] > 5) {
          continue;
        }
        // Fill  energy scale / resolution
        int ptbin = std::lower_bound(ptbins.cbegin(), ptbins.cend(), genPt) -
                    ptbins.cbegin() - 1;
        eScale[ptbin]->Fill(recoE[iPho]/mcE[genMatchedIndex], centWeight);
        // store matched gen photon index for later use
        if (isMatched2GenPhoton) {
          matchedGamma.insert(genMatchedIndex);
        }
      }
      // for every gen photon, fill hists based on its matching status
      for (int iGen = 0; iGen < (*mcPt).size(); ++iGen) {
        float genPt = (*mcPt)[iGen];
        if (genPt < ptbins.front() || genPt > ptbins.back()) {
          continue;
          // reject photons that failed HEM modules
        } else if (!passedHI18HEMfailurePho(iGen)) {
          continue;
          // gen photon selection
        } else if (mcCalIsoDR04[iGen] > 5) {
          continue;
        }
        if (mcPID[iGen] != 22)
          continue;
        total.Fill(genPt, centWeight);
        if (matchedGamma.count(iGen)) {
          matched.Fill(genPt, centWeight);
        }
      }
    }
  }
  fout.Write();
  fout.Close();
}

float findNcoll(int hiBin) {
  const int nbins = 200;
  const float Ncoll[nbins] = {
      1976.95, 1944.02, 1927.29, 1891.9,  1845.3,  1807.2,  1760.45, 1729.18,
      1674.8,  1630.3,  1590.52, 1561.72, 1516.1,  1486.5,  1444.68, 1410.88,
      1376.4,  1347.32, 1309.71, 1279.98, 1255.31, 1219.89, 1195.13, 1165.96,
      1138.92, 1113.37, 1082.26, 1062.42, 1030.6,  1009.96, 980.229, 955.443,
      936.501, 915.97,  892.063, 871.289, 847.364, 825.127, 806.584, 789.163,
      765.42,  751.187, 733.001, 708.31,  690.972, 677.711, 660.682, 640.431,
      623.839, 607.456, 593.307, 576.364, 560.967, 548.909, 530.475, 519.575,
      505.105, 490.027, 478.133, 462.372, 451.115, 442.642, 425.76,  416.364,
      405.154, 392.688, 380.565, 371.167, 360.28,  348.239, 340.587, 328.746,
      320.268, 311.752, 300.742, 292.172, 281.361, 274.249, 267.025, 258.625,
      249.931, 240.497, 235.423, 228.63,  219.854, 214.004, 205.425, 199.114,
      193.618, 185.644, 180.923, 174.289, 169.641, 161.016, 157.398, 152.151,
      147.425, 140.933, 135.924, 132.365, 127.017, 122.127, 117.817, 113.076,
      109.055, 105.16,  101.323, 98.098,  95.0548, 90.729,  87.6495, 84.0899,
      80.2237, 77.2201, 74.8848, 71.3554, 68.7745, 65.9911, 63.4136, 61.3859,
      58.1903, 56.4155, 53.8486, 52.0196, 49.2921, 47.0735, 45.4345, 43.8434,
      41.7181, 39.8988, 38.2262, 36.4435, 34.8984, 33.4664, 31.8056, 30.351,
      29.2074, 27.6924, 26.7754, 25.4965, 24.2802, 22.9651, 22.0059, 21.0915,
      19.9129, 19.1041, 18.1487, 17.3218, 16.5957, 15.5323, 14.8035, 14.2514,
      13.3782, 12.8667, 12.2891, 11.61,   11.0026, 10.3747, 9.90294, 9.42648,
      8.85324, 8.50121, 7.89834, 7.65197, 7.22768, 6.7755,  6.34855, 5.98336,
      5.76555, 5.38056, 5.11024, 4.7748,  4.59117, 4.23247, 4.00814, 3.79607,
      3.68702, 3.3767,  3.16309, 2.98282, 2.8095,  2.65875, 2.50561, 2.32516,
      2.16357, 2.03235, 1.84061, 1.72628, 1.62305, 1.48916, 1.38784, 1.28366,
      1.24693, 1.18552, 1.16085, 1.12596, 1.09298, 1.07402, 1.06105, 1.02954};
  return Ncoll[hiBin];
}

