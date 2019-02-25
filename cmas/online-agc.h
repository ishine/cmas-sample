#ifndef ONLINE_AGC_H
#define ONLINE_AGC_H

namespace spicax {
class OnlineAgc {
 public:
  // max_limit: Sarturation Level
  // noise_gate: not amplify when signal is lower than noise_gate
  // gain_pre: default amplification factor
  // win_len : searching the maximum amplitude during last win_len frames
  OnlineAgc(float max_limit = 30000.0f, float noise_gate = 0.0f, float gain_pre = 1.0f, int win_len = 1) {
    Stmp_ = 0.0f;
    Smax_ = 0.0f;
    Init(max_limit, noise_gate, gain_pre, win_len);
    frame_index_ = 0;
    preset_len_ = 16000;
    xin_float_ = new float [preset_len_];
    xout_float_ = new float [preset_len_];
  };

  ~OnlineAgc() {
    if (xin_float_ != nullptr) {delete[]xin_float_; xin_float_ = nullptr;}
    if (xout_float_ != nullptr) {delete[]xout_float_; xout_float_ = nullptr;}
  }

  void Init(float max_limit, float noise_gate, float gain_pre, int win_len) {
    gate_ = noise_gate;
    gain_pre_ = gain_pre;
    max_limit_ = max_limit;
    win_len_ = win_len;
  };
  //for short data
  bool Compute(short *xin, short *xout, int N);

  //for float data, the main body, if max_limit is a large int(e.g. 32767) mostly due to not nomarlized to 1.0
  bool Compute(float *xin, float *xout, int N);

 private:
  float Stmp_;    //
  float Smax_;    //
  float max_limit_;  // MaxValue limit
  float gate_;
  float gain_pre_;   //
  float gain_real_;  //
  int win_len_;  // Window Size
  int frame_index_;   // frame index in Window Size
  int preset_len_;
  float *xin_float_;
  float *xout_float_;
};
}
#endif
