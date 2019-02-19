#ifndef ONLINE_WEIGHTED_OVERLAP_ADD_H
#define ONLINE_WEIGHTED_OVERLAP_ADD_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>
#include "spicax-eigen.h"
#include "kaldi-common.h"
#include "emphasis.h"
#include "r2fft.h"
#include <vector>

namespace spicax {
// all the options in decomposition and reconstruction;
struct OnlineWolaOptions {
  int num_chan_;
  int out_chan_;
  int block_;
  int fs_;
  int frame_shift_;
  int frame_len_;
  int fft_len_;
  int num_bins_;
  int overlap_len_;
  int padding_zero_len_;
  float wola_cons_;
  std::string window_type;

  OnlineWolaOptions(int num_chan = 1, int out_chan = 1, int block = 10, int fs = 16000, float frame_shift_ms = 0.016,
                    float frame_len_ms = 0.032, int fft_len = 512)
    : num_chan_(num_chan),
      out_chan_(out_chan),
      block_(block),
      fs_(fs),
      frame_shift_(fs * frame_shift_ms),
      frame_len_(fs * frame_len_ms),
      fft_len_(fft_len),
      num_bins_(fft_len_ / 2 + 1),
      overlap_len_(frame_len_ - frame_shift_),
      padding_zero_len_(fft_len_ - frame_len_),
      wola_cons_(1.0f),
      window_type("sqrt_hanning") {
    // 75% overlap rate
    if (frame_len_ == 4 * frame_shift_) {wola_cons_ = 1.0f / 2.38f;}
  };
};

class OnlineWola {
 public:
  // Default constructor
  OnlineWola() = default;
  explicit OnlineWola(const OnlineWolaOptions &opts) {Init(opts);};

  bool Init(const OnlineWolaOptions &opts) {
    opts_ = opts;
    // initilize window
    if (win_.size() != opts_.frame_len_) {
      win_.setZero(opts_.frame_len_);
      if (opts_.window_type == "sqrt_hanning") {
        for (int i = 0; i < opts_.frame_len_; i++) {win_[i] = std::sin(M_PI * i / opts_.frame_len_);}
      }
    }
    in_saved_frame_ = 0;
    in_saved_samples_ = opts_.frame_len_ - opts_.frame_shift_;

    in_time_.setZero(opts_.num_chan_ * opts_.frame_len_, opts_.block_);

    pcm_in_float_.setZero(opts_.block_ * opts_.frame_shift_ * opts_.num_chan_ );
    pcm_out_float_.setZero(opts_.block_ * opts_.frame_shift_ * opts_.out_chan_);

    // output overlap buffer
    out_overlap_.setZero(opts_.out_chan_, opts_.frame_len_);
    pre_emph_.resize(opts_.num_chan_);
    return true;
  };
  // Reset
  bool Reset() {
    in_saved_frame_ = 0;
    in_saved_samples_ = opts_.frame_len_ - opts_.frame_shift_;
    if (in_time_.size() > 0)in_time_.setZero();
    if (pcm_in_float_.size() > 0)pcm_in_float_.setZero();
    if (pcm_out_float_.size() > 0)pcm_out_float_.setZero();
    if (out_overlap_.size() > 0)out_overlap_.setZero();
    return true;
  };

  // 1. slice pcm_in into overlap parts based on the frame_shift and frame_len
  // 2. transform all the frames into frequency domain
  // 3. pcm_in is interleaved<currently>, i.e. chanels are interleaved sample by sample
  // 4. return null matrix if data less than the block
  // 5. make sure samples_per_chan is less than the block_size
  const Eigen::MatrixXcf& Decompose(const short* pcm_in, int samples_per_chan, bool interleaved = true);
  const Eigen::MatrixXcf& Decompose(const float* pcm_in, int samples_per_chan, bool interleaved = true);
  // pcm_out is allocated outside, returns how many samples are reconstructed
  int Reconstruct(short* pcm_out);
  int Reconstruct(float* pcm_out);

  Eigen::MatrixXcf& GetInSpec() { return in_spec_; }
  Eigen::MatrixXcf& GetOutSpec() { return out_spec_; }

 private:
  void Analyze();
  OnlineWolaOptions opts_;

  Eigen::VectorXf pcm_in_float_;
  Eigen::VectorXf pcm_out_float_;

  int in_saved_samples_;
  int in_saved_frame_;

  Eigen::MatrixXf in_time_;       // (frame_len * num_chan_) * num_frames
  Eigen::MatrixXf out_overlap_;   // out_chan_ * frame_len
  Eigen::MatrixXcf in_spec_;      // (num_bins_ * num_chan_ ) * num_frames
  Eigen::MatrixXcf out_spec_;     // (num_bins * out_chan_) * num_frames, out_chan_ = 1

  Eigen::VectorXf win_;
  Radix2RealFft<float> normal_rfft_;

  std::vector<PreEmphasis> pre_emph_;
  DeEmphasis de_emph_;
};
}  // namespace spicax
#endif
