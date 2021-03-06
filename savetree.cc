
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
  TFile* f = TFile::Open(fname);
  TTree *evt;
  TTree *ntuple;
  f->GetObject("hiEvtAnalyzer/HiTree", evt);
  f->GetObject("ggHiNtuplizer/EventTree", ntuple);
  evt->SetBranchStatus("*", 0);
  ntuple->SetBranchStatus("*", 0);

  for (auto activeBranchName : {"nPho", "mcPID", "pho_genMatchedIndex",
                                "phoSCRawE", "mcE", "mcPt", "phoSigmaIEtaIEta_2012",
                                "pho_swissCrx", "pho_seedTime",
                                "phoHoverE", "mcCalIsoDR04",
                                "phoSCEta", "phoSCPhi", "mcEta", "mcPhi"}) {
    ntuple->SetBranchStatus(activeBranchName, 1);
  }
  evt->SetBranchStatus("hiBin", 1);

  TFile fout(outName, "recreate");
  auto newevt = evt->CloneTree();
  auto newntu = ntuple->CloneTree();

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

void savetree() {
TString zs4 = "root://xrootd.cmsaf.mit.edu//store/user/pchou/run3/ECAL/QCDPhoton_pThat-15_TuneCP5_HydjetDrumMB_5p02TeV_Pythia8/HINPbPbAutumn18GSHIMix-103X_upgrade2018_realistic_HI_v11-v1_FOREST_ECAL_ZS_EB_EE_4_SR_HI_4_MI_2/210629_221847/0000/HiForestAOD_%d.root";
TString zs8 = "root://xrootd.cmsaf.mit.edu//store/user/pchou/run3/ECAL/QCDPhoton_pThat-15_TuneCP5_HydjetDrumMB_5p02TeV_Pythia8/HINPbPbAutumn18GSHIMix-103X_upgrade2018_realistic_HI_v11-v1_FOREST_ECAL_ZS_EB_EE_8_SR_HI_4_MI_2/210629_221354/0000/HiForestAOD_%d.root";
TString zs10 = "root://xrootd.cmsaf.mit.edu//store/user/pchou/run3/ECAL/QCDPhoton_pThat-15_TuneCP5_HydjetDrumMB_5p02TeV_Pythia8/HINPbPbAutumn18GSHIMix-103X_upgrade2018_realistic_HI_v11-v1_FOREST_ECAL_ZS_EB_EE_10_SR_HI_4_MI_2/210629_221929/0000/HiForestAOD_%d.root";
TString official = "root://xrootd.cmsaf.mit.edu//store/user/katatar/official/HIRun2018PbPb/QCDPhoton_pThat-15_TuneCP5_HydjetDrumMB_5p02TeV_Pythia8/HINPbPbAutumn18DR-mva98_103X_upgrade2018_realistic_HI_v11-v1-FOREST/190918_070650/0000/HiForestAOD_%d.root";
for (auto i = 1; i <= 110; ++i) {
	cout << "saving " << i << "\n";
	mc(TString::Format(zs4, i), TString::Format("pruned/zs04_%d.root", i));
	mc(TString::Format(zs8, i), TString::Format("pruned/zs08_%d.root", i));
	mc(TString::Format(zs10, i), TString::Format("pruned/zs10_%d.root", i));
}
 for (auto i = 1; i <= 328; ++i) {
   mc(TString::Format(official, i), TString::Format("pruned/official_%d.root", i));
 }
}
