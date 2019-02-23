#include "cmas-api.h"
#include "beamform-das.h"
#include "ssl-srp.h"
#include "online-wola.h"
#include "string-split.h"
#include "spicax-eigen.h"
#include "kaldi-common.h"
#include "beamform-online-cgmm.h"
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
    out_chan_ = expect_vect.size();
    if (out_chan_ < 1)return false;
    expected_angles.setZero(out_chan_);
    std::stringstream ss;
    for (int imix = 0; imix < out_chan_; imix++) {
      ss.clear();
      ss.str("");
      ss << expect_vect[imix];
      ss >> expected_angles(imix);
    }
    expected_angles *= 60;
    //
    num_chan_ = 6;
    int array_type = 2;
    Eigen::MatrixXf mic_pos_cart;
    // //
    // if (array_type == 1) {
    //   expected_angles.array() += 30;
    //   // linear microphone array
    //   mic_pos_cart.setZero(3, num_chan_);
    //   mic_pos_cart.row(0) = Eigen::ArrayXf::LinSpaced(4, -0.0495, 0.0495);
    // }
    if (array_type == 2) {
      // circular microphone array
      num_chan_ = 6;
      mic_pos_cart.setZero(3, num_chan_);
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
    int fftlen = 1;
    while (fftlen < (frame_len_ms * fs)) {fftlen <<= 1;}
    int frames_per_block = 1;
    spicax::OnlineWolaOptions wola_opts(num_chan_, out_chan_, frames_per_block, fs, frame_shift_ms, frame_len_ms, fftlen);
    bool wola_ret = online_wola_.Init(wola_opts);
    if (!wola_ret) {return false;}

    spicax::DasOptions online_das_opts;
    online_das_opts.expected_angles = expected_angles;
    bool ret2 = online_das_.Initialize(array_, online_das_opts);
    if (!ret2) {return false;}

    spicax::SrpOptions srp_opts;
    srp_opts.block_frames = 31;
    online_srp_.Initialize(array_, srp_opts);

    frame_index_ = 0;
    delay_frames_ = 10;
    block_frames_ = 50;
    delay_index_ = block_frames_ - delay_frames_;
    num_bins_ = fftlen / 2 + 1;
    block_spec_.setZero((num_chan_ * num_bins_), block_frames_);

    spicax::OnlineCgmmMvdrOptions online_cgmm_opts;
    if (array_type == 2) {
      int num_mix = 6;
      Eigen::VectorXi preset_angles;
      preset_angles.setZero(num_mix);
      preset_angles = Eigen::ArrayXi::LinSpaced(num_mix, 0, 300);
      online_cgmm_opts.num_mix = num_mix;
      online_cgmm_opts.preset_angles = preset_angles;
    }
    online_cgmm_opts.expected_angles = expected_angles;
    online_cgmm_opts.block_frames = 101;
    online_cgmm_.Initialize(array_, online_cgmm_opts);

    return true;
  };

  int Separate(const char * pcm_in, int samples_per_chan, char *pcm_out) {
    const short *pcm_in_short = (const short *)pcm_in;
    short *pcm_out_short = (short*)pcm_out;
    const Eigen::MatrixXcf& in_spec = online_wola_.Decompose(pcm_in_short, samples_per_chan);
    if (in_spec.cols() == 0 ) {return 0;}
    // in_spec.rows() != (num_chan_ * num_bins_);
    int num_frames = in_spec.cols();
    Eigen::MatrixXcf& out_spec = online_wola_.GetOutSpec();
    if (0) {
      out_spec.setZero(num_bins_ * out_chan_, num_frames);
      Eigen::MatrixXcf in_spec_frame, out_spec_frame;
      for (int iframe = 0; iframe < num_frames; iframe++) {
        block_spec_.col(frame_index_) = in_spec.col(iframe);
        in_spec_frame = in_spec.col(iframe);
        online_srp_.Accumlate(in_spec_frame);
        const Eigen::VectorXi & est_doa = online_srp_.GetMultiDoa();
        in_spec_frame = block_spec_.col(delay_index_);
        online_das_.Beamform(in_spec_frame, out_spec_frame, est_doa);
        out_spec.col(iframe) = out_spec_frame;
        frame_index_ = (frame_index_ + 1) % block_frames_;
        delay_index_ = (delay_index_ + 1) % block_frames_;
      }
    } else {
      Eigen::MatrixXcf& out_spec = online_wola_.GetOutSpec();
      online_cgmm_.Beamform(in_spec, out_spec);
    }

    // online_das_.Beamform(in_spec, out_spec);
    int num_samples = online_wola_.Reconstruct(pcm_out_short);
    return num_samples;
  };

  bool Release() {
    online_wola_.Reset();
    online_das_.Reset();
    online_srp_.Reset();
    frame_index_ = 0;
    delay_frames_ = 10;
    block_frames_ = 50;
    delay_index_ = block_frames_ - delay_frames_;
    block_spec_.resize(0, 0);
    return true;
  };
 private:
  int delay_frames_;
  int block_frames_;
  int frame_index_;
  int delay_index_;
  int num_chan_;
  int out_chan_;
  int num_bins_;
  Eigen::MatrixXcf block_spec_;
  ArrayGeometry array_;
  spicax::OnlineWola online_wola_;
  spicax::DasComputer online_das_;
  spicax::SrpComputer online_srp_;
  spicax::OnlineCgmmMvdrComputer online_cgmm_;
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
