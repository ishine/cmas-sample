#ifndef BEMFORM_DAS_H
#define BEMFORM_DAS_H

#include "array-geometry.h"
#include "spicax-eigen.h"
#include <vector>
#include <atomic>

namespace spicax {
struct DasOptions {
  int fs;
  int fftlen;
  float lowerfreq;
  float upperfreq;
  Eigen::VectorXi expected_angles;
  DasOptions() {
    fs = 16000;
    fftlen = 512;
    lowerfreq = 100;
    upperfreq = 7950;
    expected_angles.setZero(1);
    expected_angles(0) = 90;
  };
};

class DasComputer {
 public:
  typedef DasOptions Options;
  DasComputer() = default;
  DasComputer(const ArrayGeometry &array, const DasOptions &opts) {Initialize(array, opts);};
  // Initialize Delay-and-sum Beamformer
  bool Initialize(const ArrayGeometry &array, const DasOptions &opts);
  // Set expected angles
  void SetExpectedAngles(const Eigen::VectorXi &expected_angles) {
    if (expected_angles.size() == num_output_) {
      if (opts_.expected_angles.cwiseEqual(expected_angles).count() == num_output_) {
        update_weight_ = false;
        return;
      }
    }
    opts_.expected_angles = expected_angles;
    update_weight_ = true;
  };
  // Beamform
  void Beamform(const Eigen::MatrixXcf &array_spec, Eigen::MatrixXcf &out_spec, const Eigen::VectorXi &est_doa = Eigen::VectorXi::Zero(0));
  // Beamform the spectra and get the steered response power
  void Beamform(const Eigen::MatrixXcf &array_spec, Eigen::VectorXf &out_amp2);
  // Reset delay-and-sum Beamformer
  bool Reset() {
    weight_.resize(0);
    update_weight_ = false;
    num_output_ = 0;
    return true;
  };
 private:
  void InitializeWeight(const Eigen::VectorXi & weight_angles);
  DasOptions opts_;
  ArrayGeometry array_;
  int num_bins_;
  int num_chan_;
  std::atomic_int num_output_;
  int lowerbin_;
  int upperbin_;
  std::atomic_bool update_weight_;
  std::vector<Eigen::MatrixXcf> weight_;       // num_output * (num_bins * num_chan)
};
}
#endif