#include "beamform-online-cgmm.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <cstring>
#include <cstddef>
#include <complex>

#ifndef M_LOG_2PI
#define M_LOG_2PI 1.8378770664093454835606594728112
#endif

#ifndef M_LOG_PI
#define M_LOG_PI 1.144729885849400
#endif

namespace spicax {
// Initialize
void OnlineCgmmMvdrComputer::Initialize(ArrayGeometry& array, const OnlineCgmmMvdrOptions& opts) {
  array_ = array;
  opts_ = opts;
  num_chan_ = array.GetNumMic();
  num_bins_ = opts_.fftlen / 2 + 1;
  lowerbin_ = opts_.lowerfreq * opts_.fftlen / opts_.fs;
  upperbin_ = opts_.upperfreq * opts_.fftlen / opts_.fs;
  num_mix_  = opts_.num_mix;

  block_frames_plus1_ = opts_.block_frames + 1;
  frame_index_ = 0;
  absolute_index_ = 0;

  // original array spectra block
  orginal_spec_block_.resize(block_frames_plus1_);
  for (int iframe = 0 ; iframe < block_frames_plus1_; iframe++) {
    orginal_spec_block_[iframe].setZero(num_bins_, num_chan_);
  }
  // reshuffle the array spectra according to the frequency bin
  reshuff_spec_block_.resize(num_bins_);
  noisy_cov_.resize(num_bins_);
  noise_cov_.resize(num_bins_);
  for (int ibin = lowerbin_; ibin < upperbin_; ibin++) {
    // Each column represents a multi-channel signal vector
    reshuff_spec_block_[ibin].setZero(num_chan_, block_frames_plus1_);
    noisy_cov_[ibin].setZero(num_chan_, num_chan_);
    noise_cov_[ibin].setZero(num_chan_, num_chan_);
  }
  if (num_chan_ == 6) {
    std::ifstream is;
    is.open("noise-6-256-512-16000.model", std::ios::in | std::ios::binary);
    if (is.is_open()) {
      for (int ibin = 0; ibin < num_bins_; ibin++) {
        noise_cov_[ibin].resize(num_chan_, num_chan_);
        is.read((char*)(noise_cov_[ibin].data()), num_chan_ * num_chan_ * sizeof(float) * 2);
      }
      is.close();
    } else {
      KALDI_LOG << "Cannot find noise-6-256-512-16000.model\n";
    }
  }

  selected_mix_ = 0;
  InitializeCgmmModel();

  frame_music_opts_.lowerfreq = std::max(300, opts_.lowerfreq);
  frame_music_opts_.upperfreq = std::min(4000, opts_.upperfreq);
  frame_music_opts_.fs = opts_.fs;
  frame_music_opts_.fftlen = opts_.fftlen;
  if (num_chan_ == 4) {frame_music_opts_.max_angle = 180;}
  if (num_chan_ == 6 || num_chan_ == 13) {frame_music_opts_.max_angle = 360;}
  frame_music_opts_.num_angles = frame_music_opts_.max_angle / 2;
  frame_music_opts_.mini_height = 0.9;
  frame_music_.Init(array_, frame_music_opts_);
}
// Initialize the CGMM model
void OnlineCgmmMvdrComputer::InitializeCgmmModel() {
  num_output_ = opts_.expected_angles.size();
  // (\Sigma)^(-1) and cov_accumulator
  inv_sigma_.resize(num_bins_);
  cov_accumulator_.resize(num_bins_);
  phi_block_.resize(num_bins_);
  lambda_block_.resize(num_bins_);
  for (int ibin = lowerbin_; ibin < upperbin_; ibin++) {
    inv_sigma_[ibin].setZero(num_chan_, num_chan_ * num_mix_);
    cov_accumulator_[ibin].setZero(num_chan_, num_chan_ * num_mix_);
    phi_block_[ibin].setZero(num_mix_, block_frames_plus1_);
    lambda_block_[ibin].setZero(num_mix_, block_frames_plus1_);
  }

  // Beamform weights matrix
  weight_.resize(num_mix_);
  steer_vector_.resize(num_mix_);
  for (int imix = 0; imix < num_mix_; imix++) {
    weight_[imix].setZero(num_bins_, num_chan_);
    steer_vector_[imix].setZero(num_bins_, num_chan_);
  }

  Eigen::MatrixXcf sv, mix_cov;
  Eigen::VectorXcf sv_bin;
  for (int imix = 0 ; imix < num_mix_; imix++) {
    array_.ComputeSteeringVector(opts_.preset_angles(imix), opts_.fs, opts_.fftlen, lowerbin_, upperbin_, steer_vector_[imix]);
    weight_[imix] = steer_vector_[imix].conjugate() / num_chan_;
    for (int ibin = lowerbin_; ibin < upperbin_; ibin++) {
      sv_bin = steer_vector_[imix].row(ibin);
      mix_cov = sv_bin * sv_bin.adjoint();
      float trace = mix_cov.trace().real() * 1e-4f;
      for (int ichan = 0 ; ichan < num_chan_; ichan++) { mix_cov(ichan, ichan) += trace;}
      inv_sigma_[ibin].middleCols(imix * num_chan_, num_chan_) = mix_cov.inverse();
    }
  }
}
// // Estimate the beamform weights note target_angle is preset thus the target_angle will be ignored here
// void OnlineCgmmMvdrComputer::Estimate(const Eigen::MatrixXcf &array_spec, int target_angle) {
//   // int num_frames = array_spec.cols();
//   Eigen::MatrixXcf array_frame;
//   PrepareSpectraBlock(array_spec);
//   ComputePhi();
//   ComputeLambda();
//   AccumulateSigma();
//   if (absolute_index_ >= opts_.delay_frames) {
//     est_doa_ = frame_music_.LocalizeMultiSources(noisy_cov_);
//   }
//   ComputeMvdrWeight();
// }
// // Beamform the input
// void OnlineCgmmMvdrComputer::Beamform(const Eigen::MatrixXcf &array_spec, Eigen::MatrixXcf &out_spec) {
//   out_spec.resize(num_bins_ * num_output_, 1);
//   int selected_frame;
//   if (frame_index_ >= opts_.delay_frames)
//     selected_frame = frame_index_ - opts_.delay_frames;
//   else
//     selected_frame = frame_index_ + block_frames_plus1_ - opts_.delay_frames;

//   if (absolute_index_ < opts_.delay_frames) {
//     selected_frame = frame_index_;
//     absolute_index_ ++;
//   }
//   int selected_mix;
//   Eigen::VectorXi delta;
//   for (int iout = 0; iout < num_output_; iout++) {
//     delta = opts_.preset_angles.array() - opts_.expected_angles(iout);
//     delta.cwiseAbs().minCoeff(&selected_mix);
//     out_spec.middleRows(num_bins_ * iout, num_bins_) = orginal_spec_block_[selected_frame].cwiseProduct(weight_[selected_mix]).rowwise().sum();
//   }
//   frame_index_ = (frame_index_ + 1) % block_frames_plus1_;
// }
//
void OnlineCgmmMvdrComputer::Beamform(const Eigen::MatrixXcf &array_spec, Eigen::MatrixXcf &out_spec, const Eigen::VectorXi &est_doa) {
  if (update_weight_) {
    opts_.expected_angles = expected_angles_;
    InitializeCgmmModel();
    update_weight_ = false;
  }
  if (weight_.size() == 0) {KALDI_LOG << "Beamformer is not initialized\n"; return ;}
  int num_frames = array_spec.cols();
  out_spec.resize(num_bins_ * num_output_, num_frames);
  Eigen::MatrixXcf array_frame;
  for (int iframe = 0; iframe < num_frames; iframe++) {
    array_frame = array_spec.col(iframe);
    PrepareSpectraBlock(array_frame);
    ComputePhi();
    ComputeLambda();
    AccumulateSigma();
    if (absolute_index_ >= opts_.delay_frames) {
      est_doa_ = frame_music_.LocalizeMultiSources(noisy_cov_);
    }
    ComputeMvdrWeight();

    int selected_frame;
    if (frame_index_ >= opts_.delay_frames)
      selected_frame = frame_index_ - opts_.delay_frames;
    else
      selected_frame = frame_index_ + block_frames_plus1_ - opts_.delay_frames;

    if (absolute_index_ < opts_.delay_frames) {
      selected_frame = frame_index_;
      absolute_index_ ++;
    }
    int selected_mix;
    Eigen::VectorXi delta;
    for (int iout = 0; iout < num_output_; iout++) {
      delta = opts_.preset_angles.array() - opts_.expected_angles(iout);
      delta.cwiseAbs().minCoeff(&selected_mix);
      out_spec.block(num_bins_ * iout, iframe, num_bins_, 1) =
        orginal_spec_block_[selected_frame].cwiseProduct(weight_[selected_mix]).rowwise().sum();
    }
    frame_index_ = (frame_index_ + 1) % block_frames_plus1_;
  }
}

// Prepare the block matrix
void OnlineCgmmMvdrComputer::PrepareSpectraBlock(const Eigen::MatrixXcf &array_spec) {
  // Store the orignal spectrum
  Eigen::Map<const Eigen::MatrixXcf> array_frame(&(array_spec(0, 0)), num_bins_, num_chan_);
  orginal_spec_block_[frame_index_] = array_frame;
  // Reshuffle the spectrum
  for (int ibin = lowerbin_; ibin < upperbin_; ibin++) {
    reshuff_spec_block_[ibin].col(frame_index_) = array_frame.row(ibin) / 32767;
  }
}
// Compute phi
void OnlineCgmmMvdrComputer::ComputePhi() {
  Eigen::VectorXcf obs;
  Eigen::MatrixXcf inv_sigma;
  for (int ibin = lowerbin_; ibin < upperbin_; ibin++) {
    obs = reshuff_spec_block_[ibin].col(frame_index_);
    for (int imix = 0 ; imix < num_mix_; imix++) {
      inv_sigma = inv_sigma_[ibin].middleCols(imix * num_chan_, num_chan_);
      phi_block_[ibin](imix, frame_index_) = (obs.adjoint() * inv_sigma * obs).sum().real() / num_chan_;
    }
  }
}
// Compute lambda
void OnlineCgmmMvdrComputer::ComputeLambda() {
  Eigen::VectorXf phi;
  for (int ibin = lowerbin_; ibin < upperbin_; ibin++) {
    phi = phi_block_[ibin].col(frame_index_);
    // phi = phi.cwiseMax(1e-12f);
    ComputeLoglikes(phi);
    ApplySoftMax(phi);
    lambda_block_[ibin].col(frame_index_) = phi;
  }
}
// Accumulate all the convariance matrix of each mixture for every t-f unit
void OnlineCgmmMvdrComputer::AccumulateSigma() {
  Eigen::MatrixXcf delta;
  Eigen::VectorXcf obs_newest, obs_lastest;
  Eigen::MatrixXcf cov_newest, cov_lastest;
  for (int ibin = lowerbin_; ibin < upperbin_; ibin++) {
    obs_newest = reshuff_spec_block_[ibin].col(frame_index_);
    cov_newest = obs_newest * obs_newest.adjoint();
    int lastest_index = (frame_index_ + 1) % block_frames_plus1_;
    obs_lastest = reshuff_spec_block_[ibin].col(lastest_index);
    cov_lastest = obs_lastest * obs_lastest.adjoint();

    noisy_cov_[ibin] += (cov_newest - cov_lastest);

    for (int imix = 0; imix < num_mix_; imix++) {
      delta = cov_newest * lambda_block_[ibin](imix, frame_index_)
              - cov_lastest * lambda_block_[ibin](imix, lastest_index);
      cov_accumulator_[ibin].middleCols(imix * num_chan_, num_chan_) += delta;
    }
  }
}
// Compute MVDR beamforming weights
void OnlineCgmmMvdrComputer::ComputeMvdrWeight() {
  float cons = std::sqrt((float)(num_chan_));
  Eigen::VectorXcf steer_vector;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcf> es(num_chan_);
  Eigen::MatrixXcf speech_cov, noise_cov, noisy_cov;
  float expect_norm = 20 * num_chan_;
  if (est_doa_.size() == 1) {
    // Only one sources, special process
    int nearest_mix;
    Eigen::VectorXi delta = opts_.preset_angles.array() - est_doa_(0);
    // KALDI_LOG
    delta.cwiseAbs().minCoeff(&nearest_mix);
    for (int ibin = lowerbin_; ibin < upperbin_; ibin++) {
      noisy_cov = noisy_cov_[ibin];
      for (int imix = 0; imix < num_mix_; imix++) {
        if (imix == nearest_mix) {
          // speech_cov = noisy_cov;
          speech_cov = cov_accumulator_[ibin].middleCols(imix * num_chan_, num_chan_);
          if (num_chan_ == 6) {
            noise_cov = noise_cov_[ibin];
            es.compute(speech_cov);
            Eigen::MatrixXcf speech_ev = es.eigenvectors();
            steer_vector = speech_ev.col(num_chan_ - 1) * cons;
            Eigen::VectorXcf num = noise_cov.lu().solve(steer_vector);
            float deno = (steer_vector.adjoint() * num).sum().real();
            weight_[imix].row(ibin) = num.conjugate() / deno;
          } else {
            Eigen::MatrixXcf speech_ev = es.eigenvectors();
            steer_vector = speech_ev.col(num_chan_ - 1) * cons;
            weight_[imix].row(ibin) = steer_vector.adjoint() / num_chan_;
          }
        } else {
          noise_cov = noisy_cov;
          steer_vector = steer_vector_[imix].row(ibin);
          // float trace = noise_cov.trace().real() * 1e-4f;
          // for (int ichan = 0; ichan < num_chan_; ichan++) {noise_cov(ichan, ichan) += trace;}
          Eigen::VectorXcf num = noise_cov.lu().solve(steer_vector);
          float deno = (steer_vector.adjoint() * num).sum().real();
          weight_[imix].row(ibin) = num.conjugate() / deno;
        }
        Eigen::VectorXcf weight_bin = weight_[imix].row(ibin);
        float p1norm = weight_bin.cwiseAbs().sum();
        if ( p1norm > expect_norm) {
          weight_[imix].row(ibin) *= expect_norm;
          weight_[imix].row(ibin) /= p1norm;
        }
      }
    }
    return ;
  }

  for (int ibin = lowerbin_; ibin < upperbin_; ibin++) {
    noisy_cov = noisy_cov_[ibin];
    for (int imix = 0; imix < num_mix_; imix++) {
      speech_cov = cov_accumulator_[ibin].middleCols(imix * num_chan_, num_chan_);
      noise_cov = noisy_cov - speech_cov + noise_cov_[ibin] * opts_.block_frames * 2;
      es.compute(speech_cov);
      Eigen::MatrixXcf speech_ev = es.eigenvectors();
      steer_vector = speech_ev.col(num_chan_ - 1) * cons;
      // float trace = noise_cov.trace().real() * 1e-4f;
      // for (int ichan = 0; ichan < num_chan_; ichan++) {noise_cov(ichan, ichan) += trace;}
      Eigen::VectorXcf num = noise_cov.lu().solve(steer_vector);
      float deno = (steer_vector.adjoint() * num).sum().real();
      weight_[imix].row(ibin) = num.conjugate() / deno;
      Eigen::VectorXcf weight_bin = weight_[imix].row(ibin);
      float p1norm = weight_bin.cwiseAbs().sum();

      if ( p1norm > expect_norm) {
        weight_[imix].row(ibin) *= expect_norm;
        weight_[imix].row(ibin) /= p1norm;
      }
    }
  }
}
// void OnlineCgmmMvdrComputer::Train(const Eigen::MatrixXcf &array_spec) {
//   for (int iter = 0; iter < opts_.num_iter_; iter++) {
//     Q_ = 0.0f;
//     Expectation(array_spec);
//     Maximization();
//     int num_frames = array_spec.cols();
//     Q_ /= (num_frames * (upperbin_ - lowerbin_));
//     KALDI_LOG << "iter " << iter << ":Q_=" << Q_ << "\n";
//   }
// }

// void OnlineCgmmMvdrComputer::Expectation(const Eigen::MatrixXcf &array_spec) {
//   Eigen::VectorXcf obs;
//   Eigen::VectorXf posterior;
//   Eigen::MatrixXcf obs_cor;
//   int num_frames = array_spec.cols();
//   for (int ibin = lowerbin_; ibin < upperbin_; ibin++) {
//     cov_accumulator_[ibin].setZero();
//   }
//   tfmask_.setZero();

//   for (int iframe = 0; iframe < num_frames; iframe++) {
//     Eigen::Map<const Eigen::MatrixXcf> array_frame(&(array_spec(0, iframe)), num_bins_, num_chan_);
//     for (int ibin = lowerbin_; ibin < upperbin_; ibin++) {
//       obs = array_frame.row(ibin);
//       ComponentPosteriors(obs, ibin, posterior);
//       speech_mask_(ibin, iframe) = posterior(selected_mix_);
//       obs_cor = obs * obs.adjoint();
//       for (int imix = 0; imix < num_mix_; imix++) {
//         cov_accumulator_[ibin].middleCols(imix * num_chan_, num_chan_) += posterior(imix) * obs_cor;
//       }
//     }
//   }
// }

// void OnlineCgmmMvdrComputer::Maximization() {
//   Eigen::MatrixXcf mix_cov;
//   for (int ibin = lowerbin_; ibin < upperbin_; ibin++) {
//     for (int imix = 0; imix < num_mix_; imix++) {
//       mix_cov = cov_accumulator_[ibin].middleCols(imix * num_chan_, num_chan_) / tfmask_(imix, ibin);
//       float trace_mix_cov = mix_cov.trace().real();
//       for (int ichan = 0; ichan < num_chan_; ichan++) {mix_cov(ichan, ichan) += trace_mix_cov * 1e-4f;}
//       cov_[ibin].middleCols(imix * num_chan_, num_chan_) = mix_cov;
//       gconsts_(imix, ibin) = -std::log(mix_cov.determinant().real());
//       inv_cov_[ibin].middleCols(imix * num_chan_, num_chan_) = mix_cov.inverse();
//     }
//   }
// }

// void OnlineCgmmMvdrComputer::AccumComponentCov(const Eigen::MatrixXcf &array_spec) {
//   Eigen::VectorXcf obs;
//   Eigen::VectorXf posterior;
//   Eigen::MatrixXcf obs_cor;
//   int num_frames = array_spec.cols();
//   for (int iter = 0 ; iter < opts_.num_iter_; iter++) {
//     // KALDI_LOG << "iter=" << iter << "-1\n";
//     for (int ibin = lowerbin_; ibin < upperbin_; ibin++) {cov_accumulator_[ibin].setZero();}
//     tfmask_.setZero();
//     // KALDI_LOG << "iter=" << iter << "-2\n";
//     //compute cgmm posteriori probablities
//     for (int iframe = 0; iframe < num_frames; iframe++) {
//       Eigen::Map<const Eigen::MatrixXcf> array_frame(&(array_spec(0, iframe)), num_bins_, num_chan_);
//       // KALDI_LOG << "iter=" << iter << "-3\n";
//       for (int ibin = lowerbin_; ibin < upperbin_; ibin++) {
//         obs = array_frame.row(ibin) * 1e-3f;
//         ComponentPosteriors(obs, ibin, posterior);
//         // KALDI_LOG << "iter=" << iter << "-4\n";
//         obs_cor = obs * obs.adjoint();
//         for (int imix = 0; imix < num_mix_; imix++) {
//           cov_accumulator_[ibin].middleCols(imix * num_chan_, num_chan_).noalias() += posterior(imix) * obs_cor;
//         }
//       }
//     }
//   }
// }


// void OnlineCgmmMvdrComputer::ComputeWeightAllNoise() {
//   Eigen::VectorXcf steer_vector;
//   Eigen::MatrixXcf noise_cov;
//   for (int ibin = lowerbin_; ibin < upperbin_; ibin++) {
//     steer_vector = steer_vector_.row(ibin);
//     noise_cov = noisy_cov_[ibin];
//     Eigen::VectorXcf num = noise_cov.lu().solve(steer_vector);
//     float deno = (steer_vector.adjoint() * num).sum().real();
//     weight_.row(ibin) = num.conjugate() / deno;
//   }
// }

// void OnlineCgmmMvdrComputer::ComputeWeightMaxEigen() {
//   Eigen::VectorXcf steer_vector;
//   float cons = std::sqrt((float)(num_chan_));
//   Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcf> es(num_chan_);
//   Eigen::MatrixXcf speech_cov;
//   for (int ibin = lowerbin_; ibin < upperbin_; ibin++) {
//     speech_cov = noisy_cov_[ibin];
//     es.compute(speech_cov);
//     Eigen::MatrixXcf speech_ev = es.eigenvectors();
//     steer_vector = speech_ev.col(num_chan_ - 1) * cons;
//     weight_.row(ibin) = steer_vector.adjoint() / num_chan_;
//   }
// }

// void OnlineCgmmMvdrComputer::EstimateNoisyCov(const Eigen::MatrixXcf &array_spec) {
//   Eigen::VectorXcf obs;
//   Eigen::VectorXf posterior;
//   Eigen::MatrixXcf obs_cor;
//   // float alpha = 0.5f;
//   // float _alpha = 1 - alpha;
//   int num_frames = array_spec.cols();
//   for (int ibin = lowerbin_; ibin < upperbin_; ibin++) {
//     noisy_cov_[ibin].setZero(num_chan_, num_chan_);
//   }
//   for (int iframe = 0; iframe < num_frames; iframe++) {
//     Eigen::Map<const Eigen::MatrixXcf> array_frame(&(array_spec(0, iframe)), num_bins_, num_chan_);
//     for (int ibin = lowerbin_; ibin < upperbin_; ibin++) {
//       obs = array_frame.row(ibin);
//       obs_cor = obs * obs.adjoint();
//       noisy_cov_[ibin] += obs_cor;
//     }
//   }
//   for (int ibin = lowerbin_; ibin < upperbin_; ibin++) {
//     noisy_cov_[ibin] /= num_frames;
//   }
// }


// //calculate mvdr weights
// void OnlineCgmmMvdrComputer::ComputeMvdrWeight() {
//   weight_.setZero(num_bins_, num_chan_);
//   float cons = std::sqrt((float)(num_chan_));
//   Eigen::VectorXcf steer_vector;
//   Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcf> es(num_chan_);
//   Eigen::MatrixXcf speech_cov, noise_cov;
//   for (int ibin = lowerbin_; ibin < upperbin_; ibin++) {
//     speech_cov = cov_accumulator_[ibin].middleCols(selected_mix_ * num_chan_, num_chan_);
//     noise_cov.setZero(num_chan_, num_chan_);
//     for (int imix = 0; imix < selected_mix_; imix++) {
//       noise_cov += cov_accumulator_[ibin].middleCols(imix * num_chan_, num_chan_);
//     }
//     for (int imix = selected_mix_ + 1; imix < num_mix_; imix++) {
//       noise_cov += cov_accumulator_[ibin].middleCols(imix * num_chan_, num_chan_);
//     }

//     es.compute(speech_cov);
//     Eigen::MatrixXcf speech_ev = es.eigenvectors();
//     steer_vector = speech_ev.col(num_chan_ - 1) * cons;

//     float trace = noise_cov.trace().real() * 1e-4f;
//     for (int ichan = 0; ichan < num_chan_; ichan++) {noise_cov(ichan, ichan) += trace;}
//     Eigen::VectorXcf num = noise_cov.lu().solve(steer_vector);
//     float deno = (steer_vector.adjoint() * num).sum().real();
//     weight_.row(ibin) = num.conjugate() / deno;
//   }
// }
// //apply softmax
// void OnlineCgmmMvdrComputer::ApplySoftMax(Eigen::VectorXf &v) const {
//   // exp(y)/sum(exp(y))
//   float max = v.maxCoeff();
//   v = (v.array() - max).exp();
//   float sum = v.sum();
//   v /= sum;
// }
// //calculate log likelihood
// void OnlineCgmmMvdrComputer::LogLikelihoods(const Eigen::VectorXcf &obs, int ibin, Eigen::VectorXf &loglikes) {
//   loglikes = gconsts_.col(ibin);
//   Eigen::MatrixXcf inv_cov;
// #if 1
//   ComputePhi(obs, ibin);
//   for (int imix = 0; imix < num_mix_; imix++) {
//     // loglike = gconst - num_chan * log(phi)
//     loglikes(imix) -= num_chan_ * std::log(phi_(imix, ibin));
//   }
// #else
//   for (int imix = 0; imix < num_mix_; imix++) {
//     // loglike = gconst - y^H * inv(R) * y
//     inv_cov = inv_cov_[ibin].middleCols(imix * num_chan_, num_chan_);
//     loglikes(imix) -= (obs.adjoint() * inv_cov * obs).sum().real();
//   }
// #endif
// }
// //expectation step, compute the posterioris of each component
// void OnlineCgmmMvdrComputer::ComponentPosteriors(const Eigen::VectorXcf &obs, int ibin, Eigen::VectorXf &posterior) {
//   // KALDI_LOG << "ComponentPosteriors-1\n";
//   // KALDI_LOG << "obs.rows()=" << obs.rows() << ",ibin=" << ibin << ",num_mix=" << num_mix_ << "\n";
//   LogLikelihoods(obs, ibin, posterior);
//   // KALDI_LOG << "ComponentPosteriors-2\n";
//   Eigen::VectorXf loglike = posterior;
//   // KALDI_LOG << "posterior.rows()=" << posterior.rows() << "\n";
//   ApplySoftMax(posterior);
//   // KALDI_LOG << "ComponentPosteriors-3\n";
//   Q_ += loglike.cwiseProduct(posterior).sum();
//   tfmask_.col(ibin) += posterior;
// }
// //compute the energy of each components
// void OnlineCgmmMvdrComputer::ComputePhi(const Eigen::VectorXcf &obs, int ibin) {
//   Eigen::MatrixXcf inv_cov;
//   for (int imix = 0; imix < num_mix_; imix++) {
//     inv_cov = inv_cov_[ibin].middleCols(imix * num_chan_, num_chan_);
//     phi_(imix, ibin) = (obs.adjoint() * inv_cov * obs).sum().real() / num_chan_;
//     // phi_(imix, ibin) = std::max(phi_(imix, ibin), 1e-12f);
//   }
// }
}