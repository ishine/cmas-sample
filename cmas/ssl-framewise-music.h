#ifndef SSL_FRAMEWISE_MUSIC_H
#define SSL_FRAMEWISE_MUSIC_H

#include "spicax-eigen.h"
#include "array-geometry.h"
#include <vector>

namespace spicax {
struct  FramewiseMusicOptions {
  int lowerfreq;
  int upperfreq;
  int fs;
  int fftlen;
  int num_angles;
  int max_angle;
  float mini_height;
  int mini_peak_dist;
  int Npeaks;
  FramewiseMusicOptions(): lowerfreq(100), upperfreq(7900), fs(16000),
    fftlen(1024), num_angles(180), max_angle(180), mini_height(0.9), mini_peak_dist(30), Npeaks(2) {};
};

class FramewiseMusicComputer {
 public:
  typedef FramewiseMusicOptions Options;
  FramewiseMusicComputer(): ssl_angle_(-1) {};
  FramewiseMusicComputer(const ArrayGeometry & array, const FramewiseMusicOptions &opts) {Init(array, opts);};

  FramewiseMusicComputer(const FramewiseMusicComputer & other) = delete;
  FramewiseMusicComputer& operator=(const FramewiseMusicComputer &other) = delete;

  void Init(const ArrayGeometry & array, const FramewiseMusicOptions &opts);

  // int Accumlate(const Eigen::MatrixXcf &array_spec);

  // int Localize()const {return ssl_angle_;};

  // Eigen::VectorXi& LocalizeMultiSources() {return pks_loc_;};

  // Simplest method for Framewise Cgmm Beamform
  Eigen::VectorXi &LocalizeMultiSources(const std::vector<Eigen::MatrixXcf> &noisy_cov);

 private:
  void FindPeaks();
  ArrayGeometry array_;
  FramewiseMusicOptions opts_;
  int ssl_angle_;
  int lowerbin_;
  int upperbin_;
  int num_bins_;
  int num_chan_;
  int num_angles_;
  int resolution_;
  // int block_index_;

  Eigen::VectorXi pks_loc_;
  std::vector<Eigen::MatrixXcf> sv_matrix_;//num_bins * (num_angles * num_chan)
  Eigen::MatrixXf music_spec_;//num_angles * num_bins
  Eigen::VectorXf music_;//num_angles*1
};
}


#endif