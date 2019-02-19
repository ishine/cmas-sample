#ifndef WEIGHTED_OVERLAP_ADD_H
#define WEIGHTED_OVERLAP_ADD_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include "spicax-eigen.h"
#include "emphasis.h"
#include "r2fft.h"

namespace spicax {
// all the options in decomposition and reconstruction;
struct WolaOptions {
  int num_chan_;
  int out_chan_;
  int fs_;
  int frame_shift_;
  int frame_len_;
  int fft_len_;
  int num_bins_;
  int overlap_len_;
  int padding_zero_len_;
  float wola_cons_;
  std::string win_type_;

  WolaOptions(int num_chan = 1, int out_chan = 1, int fs = 16000, float frame_shift_ms = 0.016,
              float frame_len_ms = 0.032, int fft_len = 512)
    : num_chan_(num_chan),
      out_chan_(out_chan),
      fs_(fs),
      frame_shift_(fs * frame_shift_ms),
      frame_len_(fs * frame_len_ms),
      fft_len_(fft_len),
      num_bins_(fft_len_ / 2 + 1),
      overlap_len_(frame_len_ - frame_shift_),
      padding_zero_len_(fft_len_ - frame_len_),
      wola_cons_(1.0f),
      win_type_("sqrt_hanning") {
    // 75% overlap rate
    if (frame_len_ == 4 * frame_shift_) {wola_cons_ = 1.0f / 2.38f;}
    if (frame_len_ == 3 * frame_shift_) {wola_cons_ = 1.0f / 1.26f;}
  };
  void SetParameters(int num_chan = 6, int fs = 16000, float frame_shift_ms = 0.016,
                     float frame_len_ms = 0.032, int fft_len = 512) {
    num_chan_ = num_chan;
    fs_ = fs;
    frame_shift_ = fs * frame_shift_ms;
    frame_len_ = fs * frame_len_ms;
    fft_len_ = fft_len;
    num_bins_ = fft_len_ / 2 + 1;
    overlap_len_ = frame_len_ - frame_shift_;
    padding_zero_len_ = fft_len_ - frame_len_;
    wola_cons_ = 1.0f;
    if (frame_len_ == 4 * frame_shift_) {wola_cons_ = 1.0f / 2.38f;}
    if (frame_len_ == 3 * frame_shift_) {wola_cons_ = 1.0f / 1.26f;}
  };
};

// WOLA(weighted overlap-add) is a method that could decompose and reconstruct speech stream perfectly.
class Wola {
 public:
  Wola() = default;
  explicit Wola(WolaOptions opts) {Init(opts);};
  bool Init(WolaOptions opts) {
    opts_ = opts;
    win_.setZero(opts_.frame_len_);
    if (opts_.win_type_ == "sqrt_hanning") {
      for (int i = 0; i < opts_.frame_len_; i++) {
        win_(i) = std::sin(M_PI * i / opts_.frame_len_);
      }
    }
    out_time_.setZero(opts_.out_chan_, opts_.frame_len_);
    return true;
  };
  // decompose a multichannel block, pcm_in is interleaved. samples is the length of each channel
  const Eigen::MatrixXcf& Decompose(short* pcm_in, int samples_each_chan, bool interleaved = true);
  // return how many samples are reconstructed
  int Reconstruct(short* pcm_out);

  // Eigen::MatrixXcf& GetInSpec() { return in_spec_; };
  Eigen::MatrixXcf& GetOutSpec() { return out_spec_; };

 private:
  WolaOptions opts_;
  Eigen::MatrixXcf in_spec_;   // (num_chan_ * num_bins_) * num_frames
  // std::vector<Eigen::MatrixXcf> array_spec_;
  Eigen::MatrixXcf out_spec_;  // (num_bins ) * num_frames
  Eigen::MatrixXf out_time_;

  Eigen::VectorXf win_;
  Radix2RealFft<float> normal_rfft_;

  PreEmphasis pre_emph_;
  DeEmphasis de_emph_;
};
}  // namespace spicax
#endif
