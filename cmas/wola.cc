#include "wola.h"
#include <algorithm>
#include <complex>
#include <cstddef>
#include <string>
#include "kaldi-common.h"

namespace spicax {
// Note pcm_in is interleaved
const Eigen::MatrixXcf& Wola::Decompose(short* pcm_in, int samples_each_chan, bool interleaved) {
  if (pcm_in == nullptr) {
    in_spec_.resize(0, 0);
    return in_spec_;
  }
  // slice frames according to the overlap_rate
  int num_frames_available = (samples_each_chan - opts_.overlap_len_) / opts_.frame_shift_;
  Eigen::MatrixXf wave_frames;
  wave_frames.setZero(opts_.num_chan_ * opts_.frame_len_, num_frames_available);
  for (int i = 0; i < num_frames_available; i++) {
    if (interleaved) {
      int offset = i * opts_.frame_shift_ * opts_.num_chan_;
      // channels  sample by sample
      for (int j = 0; j < opts_.num_chan_; j++) {
        for (int k = 0; k < opts_.frame_len_; k++) {
          wave_frames(j * opts_.frame_len_ + k, i) = (float)pcm_in[offset + k * opts_.num_chan_ + j];
        }
      }
    } else {
      // channels block by block
      for (int j = 0; j < opts_.num_chan_; j++) {
        for (int k = 0; k < opts_.frame_len_; k++) {
          wave_frames(j * opts_.frame_len_ + k, i) =
            (float)pcm_in[j * samples_each_chan + i * opts_.frame_shift_ + k];
        }
      }
    }
  }

  // transform to frequency domain
  in_spec_.setZero(opts_.num_bins_ * opts_.num_chan_, num_frames_available);
  float* pf = reinterpret_cast<float*>(in_spec_.data());
  for (int i = 0; i < num_frames_available; i++) {
    for (int j = 0; j < opts_.num_chan_; j++) {
      for (int k = 0; k < opts_.frame_len_; k++) {
        pf[k] = wave_frames(j * opts_.frame_len_ + k, i) * win_(k);
      }
      normal_rfft_.Compute(pf, opts_.fft_len_, true);
      pf += (opts_.num_bins_ * 2);
    }
  }
  return in_spec_;
}

int Wola::Reconstruct(short* pcm_out) {
  if (pcm_out == nullptr) {
    KALDI_LOG << "pcm_out is nullptr\n";
    return 0;
  }
  int num_frames_available = out_spec_.cols();
  if (num_frames_available == 0) {
    KALDI_LOG << "no available frames\n";
    return 0;
  }

  opts_.out_chan_ = out_spec_.rows() / opts_.num_bins_;
  if (out_spec_.rows() % opts_.num_bins_ != 0) {
    KALDI_LOG << "mismatch: out_spec_.rows()=" << out_spec_.rows() << " is not an intergral multiple of " << opts_.num_bins_ << "\n";
    return 0;
  }
  out_time_.setZero(opts_.out_chan_, opts_.frame_len_);

  // transform each column of out_spec back to time domain
  float* pf = reinterpret_cast<float*>(out_spec_.data());
  for (int i = 0; i < num_frames_available; i++) {
    for (int ichan = 0; ichan < opts_.out_chan_; ichan++) {
      normal_rfft_.Compute(pf, opts_.fft_len_, false);
      for (int j = 0; j < opts_.frame_len_; j++) {
        out_time_(ichan, j) += pf[j] * win_[j];
      }
      pf += (opts_.num_bins_ * 2);
    }

    // the head part of out_time_ is the output
    short* p_out = pcm_out + i * opts_.out_chan_ * opts_.frame_shift_;
    float* p_overlap = out_time_.data();
    for (int j = 0; j < opts_.frame_shift_ * opts_.out_chan_; j++) {
      p_out[j] = (short)(std::max(std::min((p_overlap[j] * opts_.wola_cons_), 32767.0f), -32767.0f));
    }

    // move ahead
    out_time_.topLeftCorner(opts_.out_chan_, opts_.overlap_len_) =
      out_time_.bottomRightCorner(opts_.out_chan_, opts_.overlap_len_).eval();
    out_time_.bottomRightCorner(opts_.out_chan_, opts_.frame_shift_).setZero();
  }
  return (num_frames_available * opts_.frame_shift_);
}

}  // namespace spicax