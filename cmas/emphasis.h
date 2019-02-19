#ifndef EMPHASIS_H
#define EMPHASIS_H
#include <algorithm>
#include <cmath>
#include <iostream>
namespace spicax {
class PreEmphasis {
 public:
  PreEmphasis(float alpha = 0.9) : alpha_(alpha), last_(0), lastf_(0.0f) {}
  void Compute(short *in, int len) {
    short tmp = in[len - 1];
    for (int i = len - 1; i > 0; i--) {
      in[i] -= (short)(std::floor((in[i - 1]) * alpha_));
    }
    in[0] -= (short)(std::floor((last_ * alpha_)));
    last_ = tmp;
  }
  void Compute(float *in, int len) {
    float tmp = in[len - 1];
    for (int i = len - 1; i > 0; i--) {
      in[i] -= in[i - 1] * alpha_;
    }
    in[0] -= (std::floor((lastf_ * alpha_)));
    lastf_ = tmp;
  }

 private:
  float alpha_;
  short last_;
  float lastf_;
};

class DeEmphasis {
 public:
  DeEmphasis(float alpha = 0.9) : alpha_(alpha), last_(0), lastf_(0.0f) {}
  void Compute(short *in, int len) {
    int tmp = (int)in[0] + (std::floor((last_ * alpha_)));
    // in[0] += (short)(std::floor((last_ * alpha_)));
    in[0] = (short)(std::max(std::min(tmp, 32767), -32768));
    for (int i = 1; i < len; i++) {
      tmp = (int)in[i] + (std::floor((tmp * alpha_)));
      in[i] = (short)(std::max(std::min(tmp, 32767), -32768));
      // in[i] += (short)(std::floor((in[i - 1] * alpha_)));
    }
    last_ = in[len - 1];
  }

  void Compute(float *in, int len) {
    float tmp = in[0] + last_ * alpha_;
    in[0] = (std::max(std::min(tmp, 1.0f), -1.0f));
    for (int i = 1; i < len; i++) {
      tmp = in[i] + tmp * alpha_;
      in[i] = std::max(std::min(tmp, 1.0f), -1.0f);
      // in[i] += (short)(std::floor((in[i - 1] * alpha_)));
    }
    lastf_ = in[len - 1];
  }

 private:
  float alpha_;
  int last_;
  float lastf_;
};
}  // namespace spicax
#endif