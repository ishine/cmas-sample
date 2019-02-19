#include "beamform-das.h"
#include "kaldi-common.h"

namespace spicax {
// Initialize the Delay-and-Sum Beamformer
bool DasComputer::Initialize(const ArrayGeometry &array, const DasOptions &opts) {
  array_ = array;
  opts_ = opts;
  num_bins_ = opts_.fftlen / 2 + 1;
  num_chan_ = array_.GetNumMic();
  if (num_chan_ < 1) {return false;}
  lowerbin_ = opts_.lowerfreq * opts_.fftlen / opts_.fs;
  upperbin_ = opts_.upperfreq * opts_.fftlen / opts_.fs;
  InitializeWeight(opts_.expected_angles);
  update_weight_ = false;
  return true;
}
// Initialize weight based on Far Field propagation model
void DasComputer::InitializeWeight(const Eigen::VectorXi & weight_angles) {
  num_output_ = weight_angles.size();//opts_.expected_angles.size();
  weight_.resize(num_output_);
  // KALDI_LOG << "opts_.expected_angles= " << opts_.expected_angles.transpose() << "\n";
  for (int iout = 0; iout < num_output_; iout++) {
    array_.ComputeSteeringVector(weight_angles(iout), opts_.fs, opts_.fftlen, lowerbin_, upperbin_, weight_[iout]);
    weight_[iout] = weight_[iout].conjugate() / num_chan_;
    // KALDI_LOG << iout << ",weight_(100)=" << weight_[iout].row(100) << "\n";
  }
}
// Beamform
void DasComputer::Beamform(const Eigen::MatrixXcf &array_spec, Eigen::MatrixXcf &out_spec, const Eigen::VectorXi &est_doa) {
  if (update_weight_) {
    InitializeWeight(opts_.expected_angles);
    update_weight_ = false;
  }
  // Check if it is initialized
  if (weight_.size() == 0) {KALDI_LOG << "Beamformer is not initialized\n"; return ;}

  if (est_doa.size() == 0) {
    int num_frames = array_spec.cols();
    out_spec.resize(num_bins_ * num_output_, num_frames);
    for (int iframe = 0; iframe < num_frames; iframe++) {
      Eigen::Map<const Eigen::MatrixXcf> array_frame(&(array_spec(0, iframe)), num_bins_, num_chan_);
      for (int iout = 0; iout < num_output_; iout++) {
        out_spec.col(iframe).middleRows(iout * num_bins_, num_bins_) = array_frame.cwiseProduct(weight_[iout]).rowwise().sum();
      }
    }
  } else {
    int num_frames = array_spec.cols();
    if ((est_doa.size() == 1) && (est_doa(0) < 0)) {
      out_spec.setZero(num_bins_ * num_output_, num_frames);
      return ;
    }
    out_spec.resize(num_bins_ * num_output_, num_frames);
    Eigen::VectorXi delta;
    for (int iframe = 0; iframe < num_frames; iframe++) {
      Eigen::Map<const Eigen::MatrixXcf> array_frame(&(array_spec(0, iframe)), num_bins_, num_chan_);
      for (int iout = 0; iout < num_output_; iout++) {
        delta = est_doa.array() - opts_.expected_angles(iout);
        delta = delta.cwiseAbs();
        float delta_min = delta.minCoeff();
        if (delta_min < 60) {
          out_spec.col(iframe).middleRows(iout * num_bins_, num_bins_) = array_frame.cwiseProduct(weight_[iout]).rowwise().sum();
        } else {
          out_spec.col(iframe).middleRows(iout * num_bins_, num_bins_).setZero();
        }
      }
    }
  }
}
// Compute the steered response power for Sound source localization
void DasComputer::Beamform(const Eigen::MatrixXcf &array_spec, Eigen::VectorXf &out_amp2) {
  int num_frames = array_spec.cols();
  out_amp2.setZero(num_output_);
  for (int iframe = 0; iframe < num_frames; iframe++) {
    Eigen::Map<const Eigen::MatrixXcf> array_frame(&(array_spec(0, iframe)), num_bins_, num_chan_);
    for (int iout = 0; iout < num_output_; iout++) {
      float value = array_frame.cwiseProduct(weight_[iout]).rowwise().sum().cwiseAbs2().sum();
      out_amp2(iout) += value;
    }
  }
}
}