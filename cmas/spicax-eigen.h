#ifndef SPICAX_EIGEN_H
#define SPICAX_EIGEN_H
#include "Eigen/Dense"
#include "Eigen/SVD"
#include "Eigen/LU"

// #define CIRCULAR_ARRAY 2

#define LINEAR_ARRAY 1
// #define CLOUND_MINDS_PEPPER 3

namespace spicax {

typedef Eigen::Matrix<short, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;
typedef Eigen::Matrix<short, Eigen::Dynamic, 1> VectorXs;

enum BeamformType {
  BLOCK_MAX_EIGEN = 1,
  BLOCK_CGMM = 2
};

enum SteerVectorType {
  FarField = 1,
  Noisy = 10,
  Speech = 20
};

enum LocalizationType {
  BLOCK_SRP = 1,
  BLOCK_MUSIC = 2
};
}

#endif