#include "beamform-das.h"
#include "kaldi-common.h"

namespace spicax {
// Initialize the Delay-and-Sum Beamformer
bool DasComputer::Initialize(ArrayGeometry &array, const DasOptions &opts) {
  array_ = array;
  opts_ = opts;
  num_bins_ = opts_.fftlen / 2 + 1;
  num_chan_ = array_.GetNumMic();
  if (num_chan_ < 1) {return false;}
  lowerbin_ = opts_.lowerfreq * opts_.fftlen / opts_.fs;
  upperbin_ = opts_.upperfreq * opts_.fftlen / opts_.fs;
  InitializeWeight();
  update_weight_ = false;
  return true;
}
// Initialize weight based on Far Field propagation model
void DasComputer::InitializeWeight() {
  num_output_ = opts_.expected_angles.size();
  weight_.resize(num_output_);
  KALDI_LOG << "opts_.expected_angles= " << opts_.expected_angles.transpose() << "\n";
  for (int iout = 0; iout < num_output_; iout++) {
    array_.ComputeSteeringVector(opts_.expected_angles(iout), opts_.fs, opts_.fftlen, lowerbin_, upperbin_, weight_[iout]);
    weight_[iout] = weight_[iout].conjugate() / num_chan_;
  }
}
// Beamform
void DasComputer::Beamform(const Eigen::MatrixXcf &array_spec, Eigen::MatrixXcf &out_spec, const Eigen::VectorXi &est_doa) {
  if (update_weight_) {
    InitializeWeight();
    update_weight_ = false;
  }
  // Check if it is initialized
  if (weight_.size() == 0) {KALDI_LOG << "Beamformer is not initialized\n"; return ;}

  int num_frames = array_spec.cols();
  out_spec.resize(num_bins_ * num_output_, num_frames);
  for (int iframe = 0; iframe < num_frames; iframe++) {
    Eigen::Map<const Eigen::MatrixXcf> array_frame(&(array_spec(0, iframe)), num_bins_, num_chan_);
    for (int iout = 0; iout < num_output_; iout++) {
      out_spec.col(iframe).middleRows(iout * num_bins_, num_bins_) = array_frame.cwiseProduct(weight_[iout]).rowwise().sum();
    }
  }
}
}