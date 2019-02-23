#ifndef BEAMFORM_ONLINE_CGMM_H
#define BEAMFORM_ONLINE_CGMM_H

#include "kaldi-common.h"
#include "spicax-eigen.h"
#include "array-geometry.h"
#include "ssl-framewise-music.h"
#include <atomic>

namespace spicax {
struct OnlineCgmmMvdrOptions {
  int fs;           // sample rate
  int fftlen;       // fft samples
  int num_iter;     // when set to 1, probably online process
  int num_mix;      // num of mix components
  int lowerfreq;    // lower frequency in Hz
  int upperfreq;    // upper frequency in Hz
  int block_frames;  // number of frames used to compute noise covariance matrix
  int delay_frames;  // algorithm process delay in frames
  Eigen::VectorXi preset_angles;// preset angles
  Eigen::VectorXi expected_angles;
  OnlineCgmmMvdrOptions() {
    fs = 16000;
    fftlen = 512;
    num_iter = 1;
    num_mix = 3;
    lowerfreq = 100;
    upperfreq = 7800;
    block_frames = 31;
    delay_frames = block_frames / 2; // selected index of the output
    // linear microphone array design
    preset_angles.setZero(num_mix);
    preset_angles(0) = 90;
    preset_angles(1) = 30;
    preset_angles(2) = 150;
    expected_angles.setZero(1);
    expected_angles(0) = 90;
  }
};

class OnlineCgmmMvdrComputer {
 public:
  typedef OnlineCgmmMvdrOptions Options;
  OnlineCgmmMvdrComputer() = default;
  OnlineCgmmMvdrComputer(ArrayGeometry & array, const OnlineCgmmMvdrOptions & opts)
    : array_(array), opts_(opts) {
    Initialize(array, opts);
  };
  // Intialize CGMM model based on the preset angles
  void Initialize(ArrayGeometry& array, const OnlineCgmmMvdrOptions& opts);
  // Set the output expected angles
  void SetExpectedAngles(Eigen::VectorXi &expected_angles) {
    if (expected_angles.size() == num_output_) {
      if (opts_.expected_angles.cwiseEqual(expected_angles).count() == num_output_) {
        update_weight_ = false;
        return;
      }
    }
    expected_angles_ = expected_angles;
    update_weight_ = true;
  };
  // All in one interface
  void Beamform(const Eigen::MatrixXcf &array_spec, Eigen::MatrixXcf &out_spec, const Eigen::VectorXi &est_doa = Eigen::VectorXi::Zero(0));
  // Reset all the intermediate parameters
  bool Reset() {
    if (orginal_spec_block_.size() > 0) {orginal_spec_block_.resize(0);}
    if (reshuff_spec_block_.size() > 0) {reshuff_spec_block_.resize(0);}
    if (inv_sigma_.size() > 0) {inv_sigma_.resize(0);}
    if (steer_vector_.size() > 0) {steer_vector_.resize(0);}
    if (cov_accumulator_.size() > 0) {cov_accumulator_.resize(0);}
    if (lambda_block_.size() > 0) {lambda_block_.resize(0);}
    if (phi_block_.size() > 0) {phi_block_.resize(0);}
    if (weight_.size() > 0) {weight_.resize(0);}
    if (noisy_cov_.size() > 0) {noisy_cov_.resize(0);}
    absolute_index_ = 0;
    selected_mix_ = 0;
    Q_ = 0.0f;
    return true;
  };
 private:
  // Initialize the CGMM model based on the preset angles,
  // Note only one parameter: covariance matrix \Simga
  void InitializeCgmmModel();
  // Prepare
  void PrepareSpectraBlock(const Eigen::MatrixXcf &array_spec);
  // Compute lambda(posteriori probablity) of each mixture for eacho t-f unit
  void ComputeLambda();
  // Compute phi
  void ComputePhi();
  // Accumulate covariance matrix
  void AccumulateSigma();
  // Compute MVDR weights
  void ComputeMvdrWeight();

  // Train the CGMM model
  void Train(const Eigen::MatrixXcf &array_spec);
  // Expectation Step: compute the posteriori probablity of each Component
  void Expectation();
  // Maximization Step: update the CGMM model
  void Maximization();

  // Compute loglikes
  void ComputeLoglikes(Eigen::VectorXf &loglikes) {
    for (int imix = 0; imix < num_mix_; imix++) {
      loglikes(imix) = -num_chan_ * std::log(loglikes(imix));
    }
  }
  // Apply SoftMax for loglike
  void ApplySoftMax(Eigen::VectorXf &v) {
    float max = v.maxCoeff();
    v = (v.array() - max).exp();
    float sum = v.sum();
    v /= sum;
  }

  ArrayGeometry array_;
  OnlineCgmmMvdrOptions opts_;

  FramewiseMusicOptions frame_music_opts_;
  FramewiseMusicComputer frame_music_;

  int num_bins_;
  int lowerbin_;
  int upperbin_;
  int num_mix_;
  int num_chan_;
  int num_output_;

  int block_frames_plus1_; // equals opts_.block_frames + 1, for iterative computation
  int frame_index_; // indicate where is the newest frame
  int absolute_index_;
  int selected_mix_;
  std::atomic_bool update_weight_;

  Eigen::VectorXi expected_angles_;
  // Loss value
  float Q_;

  // Original array spectra block
  std::vector<Eigen::MatrixXcf> orginal_spec_block_;       // block_frames_plus1_ * (num_bins * num_chan)
  // Reshuffled input array spectra
  std::vector<Eigen::MatrixXcf> reshuff_spec_block_;       // num_bins * (num_chan * block_frames_plus1_)
  // CGMM model paramters for each mixture and bin i.e. (\Sigma)^(-1)
  std::vector<Eigen::MatrixXcf> inv_sigma_;            // num_bins * (num_chan * num_chan * num_mix)
  // Steering vector for each preset angles
  std::vector<Eigen::MatrixXcf> steer_vector_;        // num_mix * (num_bins * num_chan)

  // Accumulated covariance matrix of mixture
  std::vector<Eigen::MatrixXcf> cov_accumulator_;       // num_bins * (num_chan * num_chan * num_mix)
  // Posterioris probablity of mixtures for each time-frequency unit
  std::vector<Eigen::MatrixXf> lambda_block_;           // num_bins * (num_mix * block_frames_plus1)
  // Phi of mixtures for each time-frequency unit
  std::vector<Eigen::MatrixXf> phi_block_;           // num_bins * (num_mix * block_frames_plus1)

  // Noisy covariance matrix
  std::vector<Eigen::MatrixXcf> noisy_cov_;
  // Noise covariance matrix
  std::vector<Eigen::MatrixXcf> noise_cov_;
  // Localization results
  Eigen::VectorXi est_doa_;
  // Beamform weight matrix, will be computed for every input frame
  std::vector<Eigen::MatrixXcf> weight_;       // num_mix * (num_bins * num_chan)
};
}
#endif