#include "HtoAA/Utilities/interface/Config.h"
#include "HtoAA/Utilities/src/Config.cc"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#include "CondFormats/BTauObjects/interface/BTagEntry.h"

using namespace std;

int main(int argc, char * argv[]) {

  if (argc!=2) {
    std::cout << "Specify btag file" << std::endl;
    exit(-1);
  }

  string bTagAlgorithm("bTagAlgorithm");
  string fileName = string(argv[1]);
  string BtagSfFile = string("/nfs/dust/cms/user/rasp/BTagReshaping/")+fileName;
  BTagCalibration calib = BTagCalibration(bTagAlgorithm, BtagSfFile);
  BTagCalibrationReader reader_B = BTagCalibrationReader(BTagEntry::OP_TIGHT, "central",{"up","down"});
  BTagCalibrationReader reader_C = BTagCalibrationReader(BTagEntry::OP_TIGHT, "central",{"up","down"});
  BTagCalibrationReader reader_Light = BTagCalibrationReader(BTagEntry::OP_TIGHT, "central",{"up","down"});
  reader_B.load(calib, BTagEntry::FLAV_B, "comb");
  reader_C.load(calib, BTagEntry::FLAV_C, "comb");
  reader_Light.load(calib, BTagEntry::FLAV_UDSG, "comb");
  
  std::vector<double> ptVec = {25,50,100};
  std::vector<double> etaVec = {0.4,2.0};

  // UDSG-Jets
  std::cout << std::endl;
  std::cout << "udsg jets..." << std::endl; 
  for (auto eta : etaVec) {
    for (auto pt : ptVec) {
      double sf = reader_Light.eval_auto_bounds("central",BTagEntry::FLAV_UDSG,eta,pt);
      double sf_up = reader_Light.eval_auto_bounds("up",BTagEntry::FLAV_UDSG,eta,pt);
      double sf_down = reader_Light.eval_auto_bounds("down",BTagEntry::FLAV_UDSG,eta,pt);
      printf(" pt = %5.1f   eta = %4.1f   SF = %4.1f  + %4.1f - %4.1f\n",pt,eta,sf,sf_up/sf,sf_down/sf);
    }
  }
  // C-Jets
  std::cout << std::endl;
  std::cout << "c jets..." << std::endl;
  for (auto eta : etaVec) {
    for (auto pt : ptVec) {
      double sf = reader_Light.eval_auto_bounds("central",BTagEntry::FLAV_C,eta,pt);
      double sf_up = reader_Light.eval_auto_bounds("up",BTagEntry::FLAV_C,eta,pt);
      double sf_down = reader_Light.eval_auto_bounds("down",BTagEntry::FLAV_C,eta,pt);
      printf(" pt = %5.1f   eta = %4.1f   SF = %4.1f  + %4.1f - %4.1f\n",pt,eta,sf,sf_up/sf,sf_down/sf);
    }
  }
  // B-Jets
  std::cout << std::endl;
  std::cout << "b jets..." << std::endl;
  for (auto eta : etaVec) {
    for (auto pt : ptVec) {
      double sf = reader_Light.eval_auto_bounds("central",BTagEntry::FLAV_B,eta,pt);
      double sf_up = reader_Light.eval_auto_bounds("up",BTagEntry::FLAV_B,eta,pt);
      double sf_down = reader_Light.eval_auto_bounds("down",BTagEntry::FLAV_B,eta,pt);
      printf(" pt = %5.1f   eta = %4.1f   SF = %4.1f  + %4.1f - %4.1f\n",pt,eta,sf,sf_up/sf,sf_down/sf);
    }
  }

}
