#ifndef ARRAY_GEOMETRY_H
#define ARRAY_GEOMETRY_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include "kaldi-common.h"
#include "spicax-eigen.h"

//array geometry class
class ArrayGeometry {
 public:
  // explicit default constructor
  ArrayGeometry() : num_chan_(0) {};
  // unary constructor
  explicit ArrayGeometry(const Eigen::MatrixXf &mic_pos_cart) {Init(mic_pos_cart);};
  // directly config array geometry, mic_pos_cart: 3(x,y,z) * num_mic
  bool Init(const Eigen::MatrixXf &mic_pos_cart) {
    if ((mic_pos_cart.rows() != 3) && (mic_pos_cart.cols() != 3)) {
      KALDI_LOG << "Expected 3-D cartesian coordinates(x, y, z) for all the microphones\n";
      return false;
    }
    if (mic_pos_cart.rows() == 3) {
      mic_pos_cart_ = mic_pos_cart;
      num_chan_ = mic_pos_cart_.cols();
    }
    if (mic_pos_cart.cols() == 3) {
      mic_pos_cart_ = mic_pos_cart.transpose();
      num_chan_ = mic_pos_cart_.rows();
    }
    return true;
  }
  // config array geometry from a file, deprecated method, too many bugs may appear.
  bool Init(const char* array_config_file) {
    try {
      std::ifstream is;
      is.open(array_config_file, std::ios::in);
      if (!is.is_open()) {
        KALDI_LOG << "can not find file: " << array_config_file << "." << std::endl;
        return false;
      }
      std::string line;
      // parse all the configs, set the number of microphones first
      std::getline(is, line);
      int pos = line.find("num_mic=");
      std::string value = line.substr(pos + 8);
      num_chan_ = std::stoi(value);
      mic_pos_cart_.setZero(3, num_chan_);
      // the second line
      std::getline(is, line);
      pos = line.find("radius=");
      if (pos != std::string::npos) {
        // polar form: radius + angle, probably circular microphone array
        value = line.substr(pos + 7);
        float radius = std::stof(value);  // radius
        Eigen::VectorXf mic_pos_polar;
        mic_pos_polar.setZero(num_chan_);
        for (int i = 0; i < num_chan_; i++) {
          std::string expect_token = "mic" + std::to_string(i) + "=";
          std::getline(is, line);
          pos = line.find(expect_token);
          if (pos == std::string::npos) {
            KALDI_LOG << "Error: can not find " << expect_token << ".\n";
            is.close();
            return false;
          }
          value = line.substr(pos + expect_token.size());
          mic_pos_polar(i) = std::stof(value);
          double rad = mic_pos_polar(i) * M_PI / 180;
          mic_pos_cart_(0, i) = radius * std::cos(rad);
          mic_pos_cart_(1, i) = radius * std::sin(rad);
        }
      } else {
        pos = line.find("mic0=(");
        if (pos != std::string::npos) {
          // cartesian form
          for (int i = 0; i < num_chan_; i++) {
            if (!ParseLine(line, i)) {
              is.close();
              return false;
            }
            std::getline(is, line);
          }
        } else {
          KALDI_LOG << "Error, config file must either be in polar form or in cartesian form.\n";
          is.close();
          return false;
        }
      }
      is.close();
      return true;
    } catch (const std::exception& e) {
      std::cerr << e.what() << std::endl;
      return false;
    }
  }
  // compute steeering vector based on far-field plane wave propagation
  void ComputeSteeringVector(int target_angle, int fs, int fftlen, int lowerbin, int upperbin, Eigen::MatrixXcf &sv) const {
    if (mic_pos_cart_.cols() == 0) {
      KALDI_LOG << "Array geometry is not initialized\n";
      return ;
    }
    int num_bins = fftlen / 2 + 1;
    sv.setZero(num_bins, num_chan_);
    float rad_s = target_angle * M_PI / 180;
    Eigen::RowVector3f direction;
    direction(0) = std::cos(rad_s);
    direction(1) = std::sin(rad_s);
    direction(2) = 0;
    Eigen::RowVectorXf dist = direction * mic_pos_cart_;
    float vs = 340.0f;
    float cons = 2 * M_PI * fs / (vs * fftlen);
    for (int ibin = lowerbin; ibin < upperbin; ibin++) {
      for (int ichan = 0; ichan < num_chan_; ichan++) {
        sv(ibin, ichan) = std::polar(1.0f, (cons * dist(ichan) * ibin));
      }
    }
  };
  //compute steeering vector based on far-field plane wave propagation
  void ComputeSteeringVector(int target_angle, int fs, int fftlen, float lowerfreq, float upperfreq, Eigen::MatrixXcf &sv) const {
    if (mic_pos_cart_.cols() == 0) {
      KALDI_LOG << "Array geometry is not initialized\n";
      return ;
    }
    int num_bins = fftlen / 2 + 1;
    int lowerbin = (int)lowerfreq * fftlen / fs;
    int upperbin = (int)upperfreq * fftlen / fs;
    sv.setZero(num_bins, num_chan_);
    float rad_s = target_angle * M_PI / 180;
    Eigen::RowVector3f direction;
    direction(0) = std::cos(rad_s);
    direction(1) = std::sin(rad_s);
    direction(2) = 0;
    Eigen::RowVectorXf dist = direction * mic_pos_cart_;
    float vs = 340.0f;
    float cons = 2 * M_PI * fs / (vs * fftlen);
    for (int ibin = lowerbin; ibin < upperbin; ibin++) {
      for (int ichan = 0; ichan < num_chan_; ichan++) {
        sv(ibin, ichan) = std::polar(1.0f, (cons * dist(ichan) * ibin));
      }
    }
  };
  //
  void ComputeDiffuseMatrix(int fs, int fftlen, float lowerfreq, float upperfreq, std::vector<Eigen::MatrixXf> &diffuse_cov)const {
    if (mic_pos_cart_.cols() == 0) {
      KALDI_LOG << "Array geometry is not initialized\n";
      return ;
    }
    Eigen::MatrixXf dist_mic;
    dist_mic.setZero(num_chan_, num_chan_);

    Eigen::MatrixXf dist_mic_x;
    for (int ichan = 0; ichan < num_chan_; ichan++) {
      dist_mic_x = mic_pos_cart_;
      dist_mic_x.colwise() -= mic_pos_cart_.col(ichan);
      dist_mic.row(ichan) = dist_mic_x.colwise().norm();
    }
    float freq;
    float freq_cons = (float)fs / fftlen * 2 * M_PI / 340.0f;
    int num_bins = fftlen / 2 + 1;
    if (diffuse_cov.size() == 0) {diffuse_cov.resize(num_bins);}
    int lowerbin = (int)lowerfreq * fftlen / fs;
    int upperbin = (int)upperfreq * fftlen / fs;
    for (int ibin = lowerbin; ibin < upperbin; ibin++) {
      freq = (float)ibin * freq_cons;
      dist_mic_x = dist_mic * freq;
      diffuse_cov[ibin] = dist_mic_x.array().sin();
      diffuse_cov[ibin] = diffuse_cov[ibin].array() / (dist_mic_x.array() + 1e-12f);
      for (int ichan = 0; ichan < num_chan_; ichan++) {diffuse_cov[ibin](ichan, ichan) = 1.001f;}
    }
  }
  //get num_chan
  int GetNumMic() const { return num_chan_; };

 private:
  bool ParseLine(std::string& line, int i) {
    std::string expect_token = "mic" + std::to_string(i) + "=(";
    std::size_t pos = line.find(expect_token);
    std::size_t pos1 = line.find_first_of(')');
    if (pos == std::string::npos || pos1 == std::string::npos) {
      KALDI_LOG << "Error: can not find the coordinates of mic " << i;
      return false;
    }
    pos += expect_token.size();
    std::string pointxyz = line.substr(pos, pos1 - pos);
    // need to check if pointxyz is legal

    // parse point in (x, y, z) form, separator is ','
    pos = pointxyz.find_first_of(',');
    std::string pointx = pointxyz.substr(0, pos);
    pos1 = pointxyz.find_last_of(',');
    pos += 1;
    std::string pointy = pointxyz.substr(pos, pos1 - pos);
    std::string pointz = pointxyz.substr(pos1 + 1);
    mic_pos_cart_(0, i) = std::stof(pointx);
    mic_pos_cart_(1, i) = std::stof(pointy);
    mic_pos_cart_(2, i) = std::stof(pointz);
    return true;
  };

  int num_chan_;                    // numbers of microphones
  Eigen::MatrixXf mic_pos_cart_;    // in cartesian coordinate, [x, y, z] 3 * num_chan_
};

#endif
