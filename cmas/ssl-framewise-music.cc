#include "ssl-framewise-music.h"
#include "kaldi-common.h"
#include <algorithm>
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

void FramewiseMusicComputer::Init(const ArrayGeometry & array, const FramewiseMusicOptions &opts) {
  array_ = array;
  opts_ = opts;
  lowerbin_ = opts_.lowerfreq * opts_.fftlen / opts_.fs;
  upperbin_ = opts_.upperfreq * opts_.fftlen / opts_.fs;
  num_bins_ = opts_.fftlen / 2 + 1;
  num_chan_ = array.GetNumMic();
  num_angles_ = opts_.num_angles;

  resolution_  = opts_.max_angle / num_angles_;
  sv_matrix_.resize(num_bins_);

  for (int ibin = lowerbin_; ibin < upperbin_; ibin++) {
    sv_matrix_[ibin].setZero(num_angles_, num_chan_);
  }
  Eigen::MatrixXcf sv;//num_bins * num_chan
  for (int iangle = 0; iangle < num_angles_; iangle++) {
    int angle = resolution_ * iangle;
    array_.ComputeSteeringVector(angle, opts_.fs, opts_.fftlen, lowerbin_, upperbin_, sv);
    for (int ibin = lowerbin_; ibin < upperbin_; ibin++) {
      sv_matrix_[ibin].row(iangle) = sv.row(ibin).conjugate();
    }
  }
  music_spec_.setZero(num_angles_, num_bins_);
}

// int FramewiseMusicComputer::Accumlate(const Eigen::MatrixXcf &array_spec) {
//   Eigen::MatrixXcf noisy_cov;
//   Eigen::VectorXcf obs;
//   int num_frames = array_spec.cols();

//   Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcf> es(num_chan_);
//   for (int ibin = lowerbin_; ibin < upperbin_; ibin++) {
//     noisy_cov.setZero(num_chan_, num_chan_);
//     for (int iframe = 0; iframe < num_frames; iframe++) {
//       Eigen::Map<const Eigen::MatrixXcf> array_frame(&(array_spec(0, iframe)), num_bins_, num_chan_);
//       obs = array_frame.row(ibin) * 0.001f;
//       noisy_cov += (obs * obs.adjoint());
//     }
//     es.compute(noisy_cov);
//     Eigen::MatrixXcf speech_ev = es.eigenvectors();
//     Eigen::MatrixXcf P_spec = sv_matrix_[ibin] * speech_ev.leftCols(num_chan_ - 1);
//     music_spec_.col(ibin) = P_spec.cwiseAbs().rowwise().sum();
//   }
//   music_ = music_spec_.rowwise().sum();
//   music_ = music_.cwiseInverse();

//   float max_value = music_.maxCoeff(&ssl_angle_);
//   music_ = music_ / max_value;
//   FindPeaks();

//   ssl_angle_ *= resolution_;

//   block_index_++;
//   return 0;
// }


void FramewiseMusicComputer::FindPeaks() {
  std::vector<int> locs;
  std::vector<float> peaks;
  // float mini_height = 0.7;
  // int mini_peak_dist = 30;
  float mini_height = opts_.mini_height;
  int mini_peak_dist = opts_.mini_peak_dist;
  int N_Peaks = opts_.Npeaks;

  for (int iangle = 1; iangle < num_angles_ - 1; iangle++) {
    if ((music_(iangle) > music_(iangle - 1)) && (music_(iangle) > music_(iangle + 1)) && (music_(iangle) >= mini_height)) {
      locs.push_back(iangle);
      peaks.push_back(music_(iangle));
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
    pks_loc_(0) = locs[0] * resolution_;
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
      if (std::abs(locs[index[ipks]] - locs_final[0]) <= mini_peak_dist) dist_test = false;
    }
    if (dist_test) {locs_final.push_back(locs[index[ipks]]);}
    if (locs_final.size() == N_Peaks)break;
  }

  int num_peaks_final = locs_final.size();
  pks_loc_.resize(num_peaks_final);
  for (int ipks = 0; ipks < num_peaks_final; ipks++) {
    pks_loc_(ipks) = locs_final[ipks] * resolution_;
  }
}

Eigen::VectorXi &FramewiseMusicComputer::LocalizeMultiSources(const std::vector<Eigen::MatrixXcf> &noisy_cov) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcf> es(num_chan_);
  for (int ibin = lowerbin_; ibin < upperbin_; ibin++) {
    es.compute(noisy_cov[ibin]);
    Eigen::MatrixXcf speech_ev = es.eigenvectors();
    Eigen::MatrixXcf P_spec = sv_matrix_[ibin] * speech_ev.leftCols(num_chan_ - 1);
    music_spec_.col(ibin) = P_spec.cwiseAbs().rowwise().sum();
  }
  music_ = music_spec_.rowwise().sum();
  music_ = music_.cwiseInverse();

  float max_value = music_.maxCoeff(&ssl_angle_);
  music_ = music_ / max_value;
  FindPeaks();
  return pks_loc_;
}
}
