#ifndef SSL_SRP_H
#define SSL_SRP_H

#include "kaldi-common.h"
#include "spicax-eigen.h"
#include "array-geometry.h"
#include "beamform-das.h"

namespace spicax {
enum WeightType {
  SRPPHAT = 1,
};

struct SrpOptions {
  int fs;
  int fftlen;
  int angle_range;
  int resolution;
  int lowerfreq;
  int upperfreq;
  int block_frames;
  float mini_height;
  int mini_peak_dist;
  int Npeaks;
  WeightType weight_type;
  SrpOptions() {
    fs = 16000;
    fftlen = 512;
    angle_range = 360;
    resolution = 5;
    lowerfreq = 500;
    upperfreq = 4000;
    block_frames = 100;
    weight_type = WeightType::SRPPHAT;
    mini_height = 0.9;
    mini_peak_dist = 30;
    Npeaks = 2;
  };
};

class SrpComputer {
 public:
  typedef SrpOptions Options;
  SrpComputer() = default;
  SrpComputer(const ArrayGeometry & array, const SrpOptions &opts) {Initialize(array, opts);};
  SrpComputer(const SrpComputer & other) = delete;
  SrpComputer&operator=(const SrpComputer & other) = delete;

  void Initialize(const ArrayGeometry & array, const SrpOptions &opts);
  int Accumlate(const Eigen::MatrixXcf &array_spec);
  int Localize()const {return pks_loc_(0);};
  const Eigen::VectorXi& GetMultiDoa() {return pks_loc_;};
  void Reset() {
    as_.resize(0);
    as_matrix_.resize(0, 0);
    frame_index_ = 0;
  };

 private:
  void FindPeaks(const Eigen::VectorXf& v);
  SrpOptions opts_;
  int frame_index_;
  int num_angles_;
  int num_bins_;
  int num_chan_;
  int block_frames_plus1_;
  Eigen::VectorXi pks_loc_;
  Eigen::VectorXf as_;// angular spectrum
  Eigen::VectorXf as_normalized_;// angular spectrum
  Eigen::MatrixXf as_matrix_;//angular spectrum matrix, store about 1 seconds.
  Eigen::MatrixXcf array_weighted_;
  DasComputer das_beam_;
};
}
#endif
