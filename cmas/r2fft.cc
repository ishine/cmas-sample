#include "r2fft.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <cstring>
#include <iostream>

#ifndef M_2PI
#define M_2PI 6.283185307179586476925286766559005
#endif

template<typename Real>
void Radix2ComplexFft<Real>::Init(int N) {
  Release();
  N_ = N;
  ComputeTwiddle();
}

template<typename Real>
void Radix2ComplexFft<Real>::Release() {
  N_ = 0;
  if (NULL != butter_) {
    delete[] butter_;
    butter_ = NULL;
  }
  if (NULL != twiddle_) {
    delete[] twiddle_;
    twiddle_ = NULL;
  }
  if (NULL != bitswap_) {
    delete[] bitswap_;
    bitswap_ = NULL;
  }
}

template<typename Real>
void Radix2ComplexFft<Real>::ComputeTwiddle() {
  int C, L, K;
  Real theta;
  int N_2 = N_ >> 1;
  butter_ = new int[N_2];
  bitswap_ = new int[N_];
  twiddle_ = new Real[N_];

  for (C = 0; C < N_2; C++) {
    theta = (M_2PI * C) / N_;
    twiddle_[2 * C] = (Real)cos(theta);
    twiddle_[2 * C + 1] = (Real)sin(theta);
  }

  butter_[0] = 0;
  L = 1;
  K = N_ >> 2;
  while (K >= 1) {
    for (C = 0; C < L; C++) butter_[C + L] = butter_[C] + K;
    L <<= 1;
    K >>= 1;
  }
}
template<typename Real>
Radix2ComplexFft<Real>::~Radix2ComplexFft() {
  Release();
}
template<typename Real>
Radix2ComplexFft<Real>::Radix2ComplexFft(const Radix2ComplexFft<Real> &other) {
  if (NULL != butter_) delete[] butter_;
  butter_ = NULL;
  if (NULL != twiddle_) delete[] twiddle_;
  twiddle_ = NULL;
  if (NULL != bitswap_) delete[] bitswap_;
  bitswap_ = NULL;
  N_ = other.N_;

  ComputeTwiddle();
}

template<typename Real>
void Radix2ComplexFft<Real>::Compute(Real *x, int N, bool forward) {
  if (N_ != N) Init(N);

  int Cycle, C, S, NC;
  int Step = N_ >> 1;
  int K1, K2;
  Real R1, I1, R2, I2;
  Real ReFFTPhi, ImFFTPhi;

  for (Cycle = 1; Cycle < N_; Cycle <<= 1, Step >>= 1) {
    K1 = 0;
    K2 = Step << 1;

    for (C = 0; C < Cycle; C++) {
      NC = butter_[C] << 1;
      ReFFTPhi = twiddle_[NC];
      ImFFTPhi = twiddle_[NC + 1];
      if (!forward) ImFFTPhi = -ImFFTPhi;
      for (S = 0; S < Step; S++) {
        R1 = x[K1];
        I1 = x[K1 + 1];
        R2 = x[K2];
        I2 = x[K2 + 1];

        x[K1++] = R1 + ReFFTPhi * R2 + ImFFTPhi * I2;
        x[K1++] = I1 - ImFFTPhi * R2 + ReFFTPhi * I2;
        x[K2++] = R1 - ReFFTPhi * R2 - ImFFTPhi * I2;
        x[K2++] = I1 + ImFFTPhi * R2 - ReFFTPhi * I2;
      }
      K1 = K2;
      K2 = K1 + (Step << 1);
    }
  }

  NC = N_ >> 1;
  for (C = 0; C < NC; C++) {
    bitswap_[C] = butter_[C] << 1;
    bitswap_[C + NC] = 1 + bitswap_[C];
  }
  for (C = 0; C < N_; C++) {
    if ((S = bitswap_[C]) != C) {
      bitswap_[S] = S;
      K1 = C << 1;
      K2 = S << 1;
      R1 = x[K1];
      x[K1++] = x[K2];
      x[K2++] = R1;
      R1 = x[K1];
      x[K1] = x[K2];
      x[K2] = R1;
    }
  }
  if (!forward) {
    NC = N_ << 1;
    for (C = 0; C < NC;) x[C++] /= N_;
  }
}

template<typename Real>
void Radix2RealFft<Real>::Init(int N) {
  if (N_ == N) return;
  if ((N & (N - 1)) != 0)
    return ;
  Release();

  N_ = N;
  r2cfft_ = new Radix2ComplexFft<Real>;
  Wnk_ = new Real[N_];
  y_ = new Real[N_];

  Real *p;
  p = Wnk_;
  Real theta;
  for (int i = 0; i < N_ / 2; i++) {
    theta = (M_2PI * i) / N_;
    *p++ = cos(theta);
    *p++ = sin(theta);
  }
}

template<typename Real>
void Radix2RealFft<Real>::Release() {
  N_ = 0;
  if (Wnk_ != NULL) delete[] Wnk_;
  Wnk_ = NULL;
  if (y_ != NULL) delete[] y_;
  y_ = NULL;
  if (r2cfft_ != NULL) delete r2cfft_;
  r2cfft_ = NULL;
}

template<typename Real>
Radix2RealFft<Real>::~Radix2RealFft() {
  Release();
}

template<typename Real>
Radix2RealFft<Real>::Radix2RealFft(const Radix2RealFft &other) : N_(other.N_) {
  Init(N_);
}

template<typename Real>
void Radix2RealFft<Real>::Compute(Real *x, int N, bool forward) {
  if (N_ != N) Init(N);

  Real forward_backward;
  if (forward) {
    forward_backward = 1;
    r2cfft_->Compute(x, N / 2, true);
  } else {
    forward_backward = -1;
  }

  Real R1, I1, R2, I2, Re1, Im1, Re2, Im2;

  Real Wre, Wim;
  for (int i = 1; i < N_ / 2; i++) {
    R1 = x[2 * i];
    I1 = x[2 * i + 1];

    I2 = x[N_ - 2 * i + 1];
    R2 = x[N_ - 2 * i];

    Re1 = (R1 + R2);
    Im1 = (I1 - I2);

    Re2 = (I1 + I2);
    Im2 = (R2 - R1);

    Wre = Wnk_[2 * i] * forward_backward;
    Wim = Wnk_[2 * i + 1];

    y_[2 * i] = Re1 + Re2 * Wre + Im2 * Wim;
    y_[2 * i] *= 0.5;
    y_[2 * i + 1] = Im1 - Re2 * Wim + Im2 * Wre;
    y_[2 * i + 1] *= 0.5;
  }
  if (forward) {
    y_[0] = x[0] + x[1];
    y_[1] = x[0] - x[1];
  } else {
    y_[0] = (x[0] + x[1]) / 2;
    y_[1] = (x[0] - x[1]) / 2;
  }

  std::memcpy(x, y_, N_ * sizeof(Real));

  if (!forward) r2cfft_->Compute(x, N / 2, false);
}

template class Radix2ComplexFft<float>;
template class Radix2ComplexFft<double>;
template class Radix2RealFft<float>;
template class Radix2RealFft<double>;
