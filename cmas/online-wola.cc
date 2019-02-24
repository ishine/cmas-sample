#include "online-wola.h"
#include <algorithm>
#include <cstddef>
#include <string>
#include "kaldi-common.h"

namespace spicax {

const Eigen::MatrixXcf& OnlineWola::Decompose(const float* pcm_in, int samples_per_chan, bool interleaved) {
  if ((pcm_in == nullptr) || (samples_per_chan == 0)) {in_spec_.resize(0, 0); return in_spec_;}
  //
  Eigen::Map<const Eigen::MatrixXf> data_in(pcm_in, opts_.num_chan_, samples_per_chan);
  //calculate how many frames could be gotten
  int samples_available = samples_per_chan + in_saved_samples_;
  int num_frames_in = (samples_available - opts_.overlap_len_) / opts_.frame_shift_;
  num_frames_in = std::max(num_frames_in, 0);
  int num_blocks = (num_frames_in + in_saved_frame_) / opts_.block_;
  if (num_blocks == 0) {in_spec_.resize(0, 0);}

  int part1_samples = std::min(opts_.frame_len_ - in_saved_samples_, samples_per_chan);
  //array_frame: frame_len_ * num_chan_
  Eigen::Map<Eigen::MatrixXf> array_frame(&in_time_(0, in_saved_frame_), opts_.frame_len_, opts_.num_chan_);
  //new frame
  array_frame.middleRows(in_saved_samples_, part1_samples) = data_in.leftCols(part1_samples).transpose();

  if (num_frames_in == 0) {
    in_saved_samples_ = (in_saved_samples_ + part1_samples) % opts_.frame_len_;
    return in_spec_;
  }

  int last_frame = in_saved_frame_;
  in_saved_frame_++;

  if (in_saved_frame_ == opts_.block_) {
    Analyze();
    in_saved_frame_ = 0;
  } else {
    Eigen::Map<Eigen::MatrixXf> array_frame_last(&in_time_(0, last_frame), opts_.frame_len_, opts_.num_chan_);
    Eigen::Map<Eigen::MatrixXf> array_frame(&in_time_(0, in_saved_frame_), opts_.frame_len_, opts_.num_chan_);
    array_frame.topRows(opts_.overlap_len_) = array_frame_last.bottomRows(opts_.overlap_len_);
  }
  int pos = part1_samples;

  for (int iframe = 1; iframe < num_frames_in; iframe++) {
    Eigen::Map<Eigen::MatrixXf> array_frame(&in_time_(0, in_saved_frame_), opts_.frame_len_, opts_.num_chan_);
    //new frame
    if (pos < samples_per_chan)
      array_frame.bottomRows(opts_.frame_shift_) = data_in.middleCols(pos, opts_.frame_shift_).transpose();
    pos += opts_.frame_shift_;

    last_frame = (last_frame + 1) % opts_.block_;
    in_saved_frame_++;
    if (in_saved_frame_ == opts_.block_) {
      Analyze();
      in_saved_frame_ = 0;
    } else {
      Eigen::Map<Eigen::MatrixXf> array_frame_new(&in_time_(0, in_saved_frame_), opts_.frame_len_, opts_.num_chan_);
      array_frame_new.topRows(opts_.overlap_len_) = array_frame.bottomRows(opts_.overlap_len_);
    }
  }

  pos = part1_samples + (num_frames_in - 1) * opts_.frame_shift_;
  int samples_leftover = samples_per_chan - pos;
  if (samples_leftover > 0) {
    Eigen::Map<Eigen::MatrixXf> array_frame(&in_time_(0, in_saved_frame_), opts_.frame_len_, opts_.num_chan_);
    array_frame.middleRows(opts_.overlap_len_, samples_leftover) = data_in.middleCols(pos, samples_leftover).transpose();
  }
  in_saved_samples_ = opts_.overlap_len_ + samples_leftover;
  return in_spec_;
}

const Eigen::MatrixXcf& OnlineWola::Decompose(const short* pcm_in, int samples_per_chan, bool interleaved) {
  int num_samples = samples_per_chan * opts_.num_chan_;
  if(num_samples > pcm_in_float_.size()){
    pcm_in_float_.setZero(num_samples);
    pcm_out_float_.setZero(num_samples);
  }
  for (int isample = 0; isample < num_samples; isample++ ) {
    pcm_in_float_(isample) = (float)pcm_in[isample];
  }
  return Decompose(pcm_in_float_.data(), samples_per_chan, interleaved);
}

void OnlineWola::Analyze() {
  in_spec_.setZero(opts_.num_bins_ * opts_.num_chan_, opts_.block_);
  Eigen::Map<Eigen::MatrixXf> array_frame(&(in_time_(0, 0)), opts_.frame_len_, opts_.num_chan_);
  Eigen::Map<Eigen::MatrixXf> array_frame_last(&(in_time_(0, opts_.block_ - 1)), opts_.frame_len_, opts_.num_chan_);

  float *pf = reinterpret_cast<float*>(in_spec_.data());
  Eigen::Map<Eigen::MatrixXf> array_spec(pf, opts_.num_bins_ * 2, opts_.num_chan_);

  for (int ichan = 0; ichan < opts_.num_chan_; ichan++) {
    // pre_emph_[ichan].Compute(&array_frame(0, ichan), opts_.frame_len_);
    array_spec.col(ichan).head(opts_.frame_len_) = array_frame.col(ichan).cwiseProduct(win_);
    normal_rfft_.Compute(pf + ichan * opts_.num_bins_ * 2, opts_.fft_len_, true);
  }

  array_frame.topRows(opts_.overlap_len_) = array_frame_last.bottomRows(opts_.overlap_len_);
  array_frame.bottomRows(opts_.frame_shift_).setZero();

  for (int iframe = 1; iframe < opts_.block_; iframe++) {
    Eigen::Map<Eigen::MatrixXf> array_frame(&(in_time_(0, iframe)), opts_.frame_len_, opts_.num_chan_);
    pf += opts_.num_bins_ * 2 * opts_.num_chan_;
    Eigen::Map<Eigen::MatrixXf> array_spec(pf, opts_.num_bins_ * 2, opts_.num_chan_);
    for (int ichan = 0; ichan < opts_.num_chan_; ichan++) {
      // pre_emph_[ichan].Compute(&array_frame(0, ichan), opts_.frame_len_);
      array_spec.col(ichan).head(opts_.frame_len_) = array_frame.col(ichan).cwiseProduct(win_);
      normal_rfft_.Compute(pf + ichan * opts_.num_bins_ * 2, opts_.fft_len_, true);
    }
  }
  in_time_.rightCols(opts_.block_ - 1).setZero();
}

// pcm_out is allocated outside
int OnlineWola::Reconstruct(short* pcm_out) {
  int ret = Reconstruct(pcm_out_float_.data());
  for (int isample = 0 ; isample < ret * opts_.out_chan_; isample++) {
    pcm_out[isample] = pcm_out_float_(isample);
  }
  return ret;
}

int OnlineWola::Reconstruct(float* pcm_out) {
  int num_frames_available = out_spec_.cols();
  if (num_frames_available == 0) {return 0;}

  // opts_.out_chan_ = out_spec_.rows() / opts_.num_bins_;
  if (opts_.num_bins_ * opts_.out_chan_ != out_spec_.rows()) {return 0;}

  int samples_per_chan = num_frames_available * opts_.frame_shift_;
  // transform each col of out_spec back to time domain
  float* pf = reinterpret_cast<float*>(out_spec_.data());
  float* p_out = pcm_out;

  for (int i = 0; i < num_frames_available; i++) {
    for (int ichan = 0; ichan < opts_.out_chan_; ichan++) {
      normal_rfft_.Compute(pf, opts_.fft_len_, false);
      // de_emph_.Compute(pf, opts_.frame_len_);
      for (int j = 0; j < opts_.frame_len_; j++) {
        out_overlap_.row(ichan)(j) += pf[j] * win_[j];
      }
      pf += (opts_.num_bins_ * 2);
    }

    float* poverlap = out_overlap_.data();
    // the head block[num_chan * frame_shift] of out_overlap_ is the output
    for (int j = 0; j < opts_.frame_shift_ * opts_.out_chan_; j++) {
      p_out[j] = std::max(std::min((poverlap[j] * opts_.wola_cons_), 32767.0f), -32767.0f);
    }
    // de_emph_.Compute(p_out, opts_.frame_len_);
    p_out += opts_.frame_shift_ * opts_.out_chan_;

    // move ahead
    out_overlap_.leftCols(opts_.overlap_len_) = out_overlap_.rightCols(opts_.overlap_len_).eval();
    out_overlap_.rightCols(opts_.frame_shift_).setZero();
  }
  return samples_per_chan;
}
}  // namespace spicax