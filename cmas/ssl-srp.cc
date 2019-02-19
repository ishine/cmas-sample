#include "ssl-srp.h"
#include <algorithm>
#include <vector>
#include <functional>
#include <numeric>

namespace spicax {
template<typename T>
int sort_index(const std::vector<T> &v, std::vector<int> &idx) {
  idx.resize(v.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&v](int i1, int i2) { return v[i1] > v[i2]; });
  return 1;
}

void SrpComputer::Initialize(const ArrayGeometry & array, const SrpOptions &opts) {
  opts_ = opts;
  num_bins_ = opts_.fftlen / 2 + 1;
  num_chan_ = array.GetNumMic();
  frame_index_ = 0;
  pks_loc_.setZero(1);
  num_angles_ = opts_.angle_range / opts_.resolution;

  as_.setZero(num_angles_);
  as_normalized_.setZero(num_angles_);

  block_frames_plus1_ = opts_.block_frames + 1;
  as_matrix_.setZero(num_angles_, block_frames_plus1_);

  DasOptions das_opts;
  das_opts.fs = opts_.fs;
  das_opts.fftlen = opts_.fftlen;
  das_opts.lowerfreq = opts_.lowerfreq;
  das_opts.upperfreq = opts_.upperfreq;
  das_opts.expected_angles = Eigen::VectorXi::LinSpaced(num_angles_, 0, opts_.angle_range);
  das_beam_.Initialize(array, das_opts);
}

int SrpComputer::Accumlate(const Eigen::MatrixXcf &array_spec) {
  int num_rows = array_spec.rows();
  int num_cols = array_spec.cols();
  if ((num_cols < 1) || (num_rows < 1)) {return 0;}
  if (((num_rows != num_bins_) && (num_cols != num_chan_)) && (num_rows != (num_bins_ * num_chan_))) {
    KALDI_LOG << "Expected " << num_bins_ << " or " << (num_bins_ * num_chan_)
              << " rows. but get" << array_spec.rows() << " rows.\n";
    return 0;
  }
  // Single frame in matrix form
  if (num_rows == num_bins_) {
    Eigen::MatrixXcf array_tmp;
    array_tmp = array_spec;
    Eigen::Map<Eigen::MatrixXcf> array_tmp_map(array_tmp.data(), num_chan_ * num_bins_, 1);
    if (opts_.weight_type == WeightType::SRPPHAT) {
      array_weighted_ = array_tmp_map.cwiseAbs().cwiseMax(1e-3f);
      array_weighted_ = array_tmp_map.cwiseQuotient(array_weighted_);
    } else {
      array_weighted_ = array_tmp_map;
    }
  } else {
    if (opts_.weight_type == WeightType::SRPPHAT) {
      array_weighted_ = array_spec.cwiseAbs().cwiseMax(1e-3f);
      array_weighted_ = array_spec.cwiseQuotient(array_weighted_);
    } else {
      array_weighted_ = array_spec;
    }
  }
  int num_frames = array_weighted_.cols();
  Eigen::VectorXf das_amp2, delta;
  Eigen::MatrixXcf array_frame;

  for (int iframe = 0; iframe < num_frames; iframe++) {
    array_frame = array_weighted_.col(iframe);
    das_beam_.Beamform(array_frame, das_amp2);
    // KALDI_LOG << "array_frame(0)=" << array_frame(100) << ",das_amp2=" << das_amp2.transpose() << "\n";
    as_matrix_.col(frame_index_) = das_amp2;
    int oldest_index = (frame_index_ + 1) % block_frames_plus1_;
    delta = (das_amp2 - as_matrix_.col(oldest_index));
    // KALDI_LOG << "delta=" << delta.transpose() << "\n";
    as_ += (das_amp2 - as_matrix_.col(oldest_index));
    frame_index_ = (frame_index_ + 1) % block_frames_plus1_;
  }

  float max_value = as_.maxCoeff();
  float min_value = as_.minCoeff();
  float diff = std::abs(min_value - max_value) / max_value;
  if (diff < 0.1) {
    pks_loc_.resize(1);
    pks_loc_(0) = -1000;
  } else {
    as_normalized_ = as_ / max_value;
    FindPeaks(as_normalized_);
  }
  return 0;
}

void SrpComputer::FindPeaks(const Eigen::VectorXf& v) {
  std::vector<int> locs;
  std::vector<float> peaks;
  // float mini_height = 0.7;
  // int mini_peak_dist = 30;
  float mini_height = opts_.mini_height;
  int mini_peak_dist = opts_.mini_peak_dist;
  int N_Peaks = opts_.Npeaks;
  // if (((v(0) > v(1)) && (v(0) > v(num_angles_ - 1))) && (v(0) >= mini_height)) {
  //   locs.push_back(0);
  //   peaks.push_back(v(0));
  // }
  // if ((v(num_angles_ - 1) > v(0)) && (v(num_angles_ - 1) > v(num_angles_ - 2)) && (v(num_angles_ - 1) >= mini_height)) {
  //   locs.push_back(0);
  //   peaks.push_back(v(0));
  // }
  // KALDI_LOG << "v=" << v.transpose() << "\n";
  for (int iangle = 1; iangle < num_angles_ - 1; iangle++) {
    if ((v(iangle) > v(iangle - 1)) && (v(iangle) > v(iangle + 1)) && (v(iangle) >= mini_height)) {
      locs.push_back(iangle);
      peaks.push_back(v(iangle));
    }
  }

  int num_pks = locs.size();

  if (num_pks == 0 ) {
    pks_loc_.resize(1);
    pks_loc_(0) = 0;
    return ;
  }
  if (num_pks == 1) {
    pks_loc_.resize(1);
    pks_loc_(0) = locs[0] * opts_.resolution;
    return;
  }

  std::vector<int> index;
  sort_index(peaks, index);
  std::vector<int> locs_final;
  locs_final.push_back(locs[index[0]]);
  for (int ipks = 1; ipks < num_pks; ipks++) {
    bool dist_test = true;
    int num_locs = locs_final.size();
    for (int iloc = 0; iloc < num_locs; iloc++) {
      if (std::abs(locs[index[ipks]] - locs_final[iloc]) <= mini_peak_dist) dist_test = false;
    }
    if (dist_test) {locs_final.push_back(locs[index[ipks]]);}
    if (locs_final.size() == N_Peaks)break;
  }

  int num_peaks_final = locs_final.size();
  pks_loc_.resize(num_peaks_final);
  for (int ipks = 0; ipks < num_peaks_final; ipks++) {
    pks_loc_(ipks) = locs_final[ipks] * opts_.resolution;
  }
}
}
