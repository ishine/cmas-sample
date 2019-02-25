#include "online-agc.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>

namespace spicax {
bool OnlineAgc::Compute(short *xin, short *xout, int N) {
  if (N > preset_len_) {
    if (xin_float_ != nullptr) {delete [] xin_float_; xin_float_ = nullptr;}
    if (xout_float_ != nullptr) {delete[]xout_float_; xout_float_ = nullptr;}
    xin_float_ = new float [N];
    xout_float_ = new float[N];
  }

  if (N <= 0 || xin == nullptr || xout == nullptr) {return false;}

  if (max_limit_ < 1.0f) {max_limit_ *= 32767;}
  for (int isample = 0; isample < N ; isample++) {
    xin_float_[isample] = xin[isample];
  }
  Compute(xin_float_, xout_float_, N);
  for (int isample = 0 ; isample < N; isample++) {
    xout[isample] = (short)floor(xout_float_[isample]);
  }
  return true;
}

bool OnlineAgc::Compute(float *xin, float *xout, int N) {
  if (N == 0 || xin == nullptr || xout == nullptr) {return false;}

  float max_val = 0.0;

  //find the maximum absolute value in this block
  for (int i = 0; i < N; i++) {
    float absvalue = (std::abs(xin[i]));
    if (max_val < absvalue) max_val = absvalue;
  }

  //find the maxium value in the preset search window
  if (frame_index_ == win_len_) {
    Smax_ = std::max(max_val, Stmp_);
    Stmp_ = max_val;
    frame_index_ = 1;
  } else {
    Stmp_ = std::max(max_val, Stmp_);
    Smax_ = std::max(max_val, Smax_);
    frame_index_++;
  }

  float max_output = (float)(Smax_ * gain_pre_);

  //avoid clipping
  float gain_target;
  if (max_output >= max_limit_) {
    gain_target = (float)max_limit_ / Smax_;
  } else {
    gain_target = gain_pre_;
  }

  //noise pieces
  if (Smax_ < gate_) {
    gain_target = gain_pre_ * Smax_ / gate_;
  }

  float alpha, alpha1;
  alpha = 0.95f;
  alpha1 = 1 - alpha;
  float gain = gain_real_;

  for (int i = 0; i < N; i++) {
    gain = alpha * gain + alpha1 * gain_target;
    xout[i] = std::max(std::min((xin[i] * gain), max_limit_), -max_limit_);
  }
  gain_real_ = gain;
  return true;
}

}