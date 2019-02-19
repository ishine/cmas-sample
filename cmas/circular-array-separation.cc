#include "cmas-api.h"
#include "beamform-das.h"
#include "online-wola.h"
#include "string-split.h"
#include "spicax-eigen.h"
#include "kaldi-common.h"
#include <sstream>

class CmasSeparator {
 public:
  CmasSeparator() = default;
  CmasSeparator(const char * expected_beam_index) {Initialize(expected_beam_index);};
  // Initialize array geometry, wola and online cgmm.
  bool Initialize(const char * expected_beam_index) {
    // get the expected angle
    std::string expect_angle_string(expected_beam_index);
    Eigen::VectorXi expected_angles;
    std::vector<std::string> expect_vect = spicax::split(expect_angle_string, "|");
    int out_chan = expect_vect.size();
    expected_angles.setZero(out_chan);
    std::stringstream ss;
    for (int imix = 0; imix < out_chan; imix++) {
      ss.clear();
      ss.str("");
      ss << expect_vect[imix];
      ss >> expected_angles(imix);
    }
    expected_angles *= 60;
    //
    int num_chan = 4;
    int array_type = 2;
    Eigen::MatrixXf mic_pos_cart;
    //
    if (array_type == 1) {
      expected_angles.array() += 30;
      // linear microphone array
      mic_pos_cart.setZero(3, num_chan);
      mic_pos_cart.row(0) = Eigen::ArrayXf::LinSpaced(4, -0.0495, 0.0495);
    }
    if (array_type == 2) {
      // circular microphone array
      num_chan = 6;
      mic_pos_cart.setZero(3, num_chan);
      float radius = 0.043f;
      mic_pos_cart(0, 0) = radius;
      mic_pos_cart(1, 0) = 0.0f;
      mic_pos_cart(0, 1) = radius * 0.5f;
      mic_pos_cart(1, 1) = -radius * std::sqrt(3) * 0.5f;
      mic_pos_cart(0, 2) = -radius * 0.5f;
      mic_pos_cart(1, 2) = -radius * std::sqrt(3) * 0.5f;
      mic_pos_cart(0, 3) = -radius;
      mic_pos_cart(1, 3) = 0.0f;
      mic_pos_cart(0, 4) = -radius * 0.5f;
      mic_pos_cart(1, 4) = radius * std::sqrt(3) * 0.5f;
      mic_pos_cart(0, 5) = radius * 0.5f;
      mic_pos_cart(1, 5) = radius * std::sqrt(3) * 0.5f;
    }
    bool array_ret = array_.Init(mic_pos_cart);
    if (!array_ret) {return false;}

    // Weighted Overlap Add
    int fs = 16000;
    float frame_shift_ms = 0.016f;
    // int frame_shift_samples = frame_shift_ms * fs;
    int overlap_rate = 2;
    float frame_len_ms = frame_shift_ms * overlap_rate;
    int fft_len = 1;
    while (fft_len < (frame_len_ms * fs)) {fft_len <<= 1;}
    int frames_per_block = 1;
    spicax::OnlineWolaOptions wola_opts(num_chan, out_chan, frames_per_block, fs, frame_shift_ms, frame_len_ms, fft_len);
    bool wola_ret = online_wola_.Init(wola_opts);
    if (!wola_ret) {return false;}

    spicax::DasOptions online_das_opts;

    online_das_opts.expected_angles = expected_angles;
    bool ret2 = online_das_.Initialize(array_, online_das_opts);
    if (!ret2) {return false;}
    return true;
  };

  int Separate(const char * pcm_in, int samples_per_chan, char *pcm_out) {
    const short *pcm_in_short = (const short *)pcm_in;
    short *pcm_out_short = (short*)pcm_out;
    const Eigen::MatrixXcf& in_spec = online_wola_.Decompose(pcm_in_short, samples_per_chan);
    if (in_spec.cols() == 0 ) {return 0;}
    Eigen::MatrixXcf& out_spec = online_wola_.GetOutSpec();
    online_das_.Beamform(in_spec, out_spec);
    int num_samples = online_wola_.Reconstruct(pcm_out_short);
    return num_samples;
  };

  bool Release() {
    online_wola_.Reset();
    online_das_.Reset();
    return true;
  };
 private:
  ArrayGeometry array_;
  spicax::OnlineWola online_wola_;
  spicax::DasComputer online_das_;
};

CmasSeparator separator;

bool CmasInitialize(const char * expected_beam_index) {
  return separator.Initialize(expected_beam_index);
}

int CmasSeparate(const char *pcm_in, int samples_per_chan, char *pcm_out) {
  return separator.Separate(pcm_in, samples_per_chan, pcm_out);
}

bool CmasRelease() {
  return separator.Release();
}
