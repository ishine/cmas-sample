// feat/resample.cc

// Copyright    2013  Pegah Ghahremani
//              2014  IMSL, PKU-HKUST (author: Wei Shi)
//              2014  Yanqing Sun, Junjie Wang
//              2014  Johns Hopkins University (author: Daniel Povey)

// See ../../COPYING for clarification regarding multiple authors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//  http://www.apache.org/licenses/LICENSE-2.0
//
// THIS CODE IS PROVIDED *AS IS* BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, EITHER EXPRESS OR IMPLIED, INCLUDING WITHOUT LIMITATION ANY IMPLIED
// WARRANTIES OR CONDITIONS OF TITLE, FITNESS FOR A PARTICULAR PURPOSE,
// MERCHANTABLITY OR NON-INFRINGEMENT.
// See the Apache 2 License for the specific language governing permissions and
// limitations under the License.


#include <algorithm>
#include <limits>
// #include "feat/feature-functions.h"
// #include "matrix/matrix-functions.h"
#include "resample.h"

// namespace kaldi {
LinearResample::LinearResample(int samp_rate_in_hz,
                               int samp_rate_out_hz,
                               float filter_cutoff_hz,
                               int num_zeros):
  samp_rate_in_(samp_rate_in_hz),
  samp_rate_out_(samp_rate_out_hz),
  filter_cutoff_(filter_cutoff_hz),
  num_zeros_(num_zeros) {
  KALDI_ASSERT(samp_rate_in_hz > 0.0 &&
               samp_rate_out_hz > 0.0 &&
               filter_cutoff_hz > 0.0 &&
               filter_cutoff_hz * 2 <= samp_rate_in_hz &&
               filter_cutoff_hz * 2 <= samp_rate_out_hz &&
               num_zeros > 0);

  // base_freq is the frequency of the repeating unit, which is the gcd
  // of the input frequencies.
  int base_freq = kaldi::Gcd(samp_rate_in_, samp_rate_out_);
  input_samples_in_unit_ = samp_rate_in_ / base_freq;
  output_samples_in_unit_ = samp_rate_out_ / base_freq;

  SetIndexesAndWeights();
  Reset();
}

void LinearResample::Init(int samp_rate_in_hz,
                          int samp_rate_out_hz,
                          float filter_cutoff_hz,
                          int num_zeros) {
  samp_rate_in_ = (samp_rate_in_hz);
  samp_rate_out_ = (samp_rate_out_hz);
  filter_cutoff_ = (filter_cutoff_hz);
  num_zeros_ = (num_zeros);
  KALDI_ASSERT(samp_rate_in_hz > 0.0 &&
               samp_rate_out_hz > 0.0 &&
               filter_cutoff_hz > 0.0 &&
               filter_cutoff_hz * 2 <= samp_rate_in_hz &&
               filter_cutoff_hz * 2 <= samp_rate_out_hz &&
               num_zeros > 0);

  // base_freq is the frequency of the repeating unit, which is the gcd
  // of the input frequencies.
  int base_freq = kaldi::Gcd(samp_rate_in_, samp_rate_out_);
  input_samples_in_unit_ = samp_rate_in_ / base_freq;
  output_samples_in_unit_ = samp_rate_out_ / base_freq;
  SetIndexesAndWeights();
  Reset();
}

int LinearResample::GetNumOutputSamples(int input_num_samp,
                                        bool flush) const {
  // For exact computation, we measure time in "ticks" of 1.0 / tick_freq,
  // where tick_freq is the least common multiple of samp_rate_in_ and
  // samp_rate_out_.
  int tick_freq = kaldi::Lcm(samp_rate_in_, samp_rate_out_);
  int ticks_per_input_period = tick_freq / samp_rate_in_;

  // work out the number of ticks in the time interval
  // [ 0, input_num_samp/samp_rate_in_ ).
  int interval_length_in_ticks = input_num_samp * ticks_per_input_period;
  if (!flush) {
    double window_width = num_zeros_ / (2.0 * filter_cutoff_);
    // To count the window-width in ticks we take the floor.  This
    // is because since we're looking for the largest integer num-out-samp
    // that fits in the interval, which is open on the right, a reduction
    // in interval length of less than a tick will never make a difference.
    // For example, the largest integer in the interval [ 0, 2 ) and the
    // largest integer in the interval [ 0, 2 - 0.9 ) are the same (both one).
    // So when we're subtracting the window-width we can ignore the fractional
    // part.
    int window_width_ticks = floor(window_width * tick_freq);
    // The time-period of the output that we can sample gets reduced
    // by the window-width (which is actually the distance from the
    // center to the edge of the windowing function) if we're not
    // "flushing the output".
    // interval_length_in_ticks -= window_width_ticks;
  }
  if (interval_length_in_ticks <= 0)
    return 0;
  int ticks_per_output_period = tick_freq / samp_rate_out_;
  // Get the last output-sample in the closed interval, i.e. replacing [ ) with
  // [ ].  Note: integer division rounds down.  See
  // http://en.wikipedia.org/wiki/Interval_(mathematics) for an explanation of
  // the notation.
  int last_output_samp = interval_length_in_ticks / ticks_per_output_period;
  // We need the last output-sample in the open interval, so if it takes us to
  // the end of the interval exactly, subtract one.
  if (last_output_samp * ticks_per_output_period == interval_length_in_ticks)
    last_output_samp--;
  // First output-sample index is zero, so the number of output samples
  // is the last output-sample plus one.
  int num_output_samp = last_output_samp + 1;
  return num_output_samp;
}

void LinearResample::SetIndexesAndWeights() {
  first_index_.resize(output_samples_in_unit_);
  weights_.resize(output_samples_in_unit_);

  double window_width = num_zeros_ / (2.0 * filter_cutoff_);

  int max_len = 0, min_len = 256;
  for (int i = 0; i < output_samples_in_unit_; i++) {
    double output_t = i / static_cast<double>(samp_rate_out_);
    double min_t = output_t - window_width, max_t = output_t + window_width;
    // we do ceil on the min and floor on the max, because if we did it
    // the other way around we would unnecessarily include indexes just
    // outside the window, with zero coefficients.  It's possible
    // if the arguments to the ceil and floor expressions are integers
    // (e.g. if filter_cutoff_ has an exact ratio with the sample rates),
    // that we unnecessarily include something with a zero coefficient,
    // but this is only a slight efficiency issue.
    int min_input_index = ceil(min_t * samp_rate_in_),
        max_input_index = floor(max_t * samp_rate_in_),
        num_indices = max_input_index - min_input_index + 1;
    first_index_[i] = min_input_index;
    weights_[i].resize(num_indices);
    for (int j = 0; j < num_indices; j++) {
      int input_index = min_input_index + j;
      double input_t = input_index / static_cast<double>(samp_rate_in_),
             delta_t = input_t - output_t;
      // sign of delta_t doesn't matter.
      weights_[i](j) = FilterFunc(delta_t) / samp_rate_in_;
    }
    if (max_len < num_indices) {max_len = num_indices;}
    if (min_len > num_indices) {min_len = num_indices;}
  }
  // KALDI_LOG << "max_len=" << max_len << ",min_len=" << min_len << "\n";
  input_remainder_.setZero(max_len);
  causal_offset_ = min_len / 2;
}


// inline
void LinearResample::GetIndexes(int samp_out, int *first_samp_in, int *samp_out_wrapped) const {
  // A unit is the smallest nonzero amount of time that is an exact
  // multiple of the input and output sample periods.  The unit index
  // is the answer to "which numbered unit we are in".
  int unit_index = samp_out / output_samples_in_unit_;
  // samp_out_wrapped is equal to samp_out % output_samples_in_unit_
  *samp_out_wrapped = static_cast<int>(samp_out - unit_index * output_samples_in_unit_);
  *first_samp_in = first_index_[*samp_out_wrapped] + unit_index * input_samples_in_unit_;
}


void LinearResample::Resample(const Eigen::VectorXf &input, bool flush, Eigen::VectorXf &output) {
  int input_dim = input.size();
  int tot_input_samp = input_sample_offset_ + input_dim,
      tot_output_samp = GetNumOutputSamples(tot_input_samp, flush);

  KALDI_ASSERT(tot_output_samp >= output_sample_offset_);

  output.setZero(tot_output_samp - output_sample_offset_);
  // samp_out is the index into the total output signal, not just the part
  // of it we are producing here.
  for (int samp_out = output_sample_offset_; samp_out < tot_output_samp; samp_out++) {
    int first_samp_in;
    int samp_out_wrapped;
    GetIndexes(samp_out, &first_samp_in, &samp_out_wrapped);
    const Eigen::VectorXf &weights = weights_[samp_out_wrapped];
    // first_input_index is the first index into "input" that we have a weight
    // for.
    int first_input_index = static_cast<int>(first_samp_in - input_sample_offset_);
    first_input_index -= causal_offset_;
    float this_output;
    if (first_input_index >= 0 && first_input_index + weights.size() <= input_dim) {
      Eigen::VectorXf input_part = input.segment(first_input_index, weights.size());
      this_output = input_part.dot(weights);
    } else {  // Handle edge cases.
      this_output = 0.0;
      for (int i = 0; i < weights.size(); i++) {
        float weight = weights(i);
        int input_index = first_input_index + i;
        if (input_index < 0 && input_remainder_.size() + input_index >= 0) {
          this_output += weight * input_remainder_(input_remainder_.size() + input_index);
        } else if (input_index >= 0 && input_index < input_dim) {
          this_output += weight * input(input_index);
        } else if (input_index >= input_dim) {
          // We're past the end of the input and are adding zero; should only
          // happen if the user specified flush == true, or else we would not
          // be trying to output this sample.
          KALDI_ASSERT(flush);
        }
      }
    }
    int output_index = static_cast<int>(samp_out - output_sample_offset_);
    output(output_index) = this_output;
  }

  if (flush) {
    Reset();  // Reset the internal state.
  } else {
    SetRemainder(input);
    if (tot_input_samp > samp_rate_in_) tot_input_samp -= samp_rate_in_;
    if (tot_output_samp > samp_rate_out_) tot_output_samp -= samp_rate_out_;
    input_sample_offset_ = tot_input_samp;
    output_sample_offset_ = tot_output_samp;
  }
}

void LinearResample::SetRemainder(const Eigen::VectorXf &input) {
  Eigen::VectorXf old_remainder(input_remainder_);
  // max_remainder_needed is the width of the filter from side to side,
  // measured in input samples.  you might think it should be half that,
  // but you have to consider that you might be wanting to output samples
  // that are "in the past" relative to the beginning of the latest
  // input... anyway, storing more remainder than needed is not harmful.
  int max_remainder_needed = ceil(samp_rate_in_ * num_zeros_ / filter_cutoff_);
  input_remainder_.resize(max_remainder_needed);
  for (int index = - input_remainder_.size(); index < 0; index++) {
    // we interpret "index" as an offset from the end of "input" and
    // from the end of input_remainder_.
    int input_index = index + input.size();
    if (input_index >= 0)
      input_remainder_(index + input_remainder_.size()) = input(input_index);
    else if (input_index + old_remainder.size() >= 0)
      input_remainder_(index + input_remainder_.size()) =
        old_remainder(input_index + old_remainder.size());
    // else leave it at zero.
  }
}

void LinearResample::Reset() {
  input_sample_offset_ = 0;
  output_sample_offset_ = 0;
  input_remainder_.resize(0);
}

/** Here, t is a time in seconds representing an offset from
    the center of the windowed filter function, and FilterFunction(t)
    returns the windowed filter function, described
    in the header as h(t) = f(t)g(t), evaluated at t.
*/
float LinearResample::FilterFunc(float t) const {
  float window,  // raised-cosine (Hanning) window of width num_zeros_/2*filter_cutoff_
        filter;  // sinc filter function
  if (fabs(t) < num_zeros_ / (2.0 * filter_cutoff_))
    window = 0.5 * (1 + cos(M_2PI * filter_cutoff_ / num_zeros_ * t));
  else
    window = 0.0;  // outside support of window function
  if (t != 0)
    filter = sin(M_2PI * filter_cutoff_ * t) / (M_PI * t);
  else
    filter = 2 * filter_cutoff_;  // limit of the function at t = 0
  return filter * window;
}

// void DownsampleWaveForm(float orig_freq, const Eigen::VectorXf &wave,
//                         float new_freq, Eigen::VectorXf &new_wave) {
//   KALDI_ASSERT(new_freq < orig_freq);
//   float lowpass_cutoff = 0.99 * 0.5 * new_freq;
//   int lowpass_filter_width = 6;
//   LinearResample signal_downsampler(orig_freq, new_freq,
//                                     lowpass_cutoff, lowpass_filter_width);
//   signal_downsampler.Resample(wave, true, new_wave);
// }
// }  // namespace kaldi
