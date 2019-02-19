#ifndef RADIX2FFT_H
#define RADIX2FFT_H
#include <complex>

using std::complex;
template<typename Real>
class Radix2ComplexFft {
 public:
  // N is the number of complex points (must be a power of two, or this
  // will crash).  Note that the constructor does some work so it's best to
  // initialize the object once and do the computation many times.
  Radix2ComplexFft() : N_(0), twiddle_(NULL), butter_(NULL), bitswap_(NULL) {}

  // destructor
  ~Radix2ComplexFft();
  // Copy constructor
  Radix2ComplexFft(const Radix2ComplexFft &other);

  // This version of Compute takes a single array of size N*2,
  // containing [ r0 im0 r1 im1 ... ].  Otherwise its behavior is  the
  // same as the version above.
  void Compute(Real *x, int N, bool forward);

 private:
  void ComputeTwiddle();
  void Init(int N);
  void Release();

  int N_;          // fft len
  Real *twiddle_;      // twiddler factor : same as equation expression Wnk
  int *butter_;   // butterfly order
  int *bitswap_;  // for final permutation

  // Disallow assignment.
  Radix2ComplexFft &operator=(const Radix2ComplexFft &other);
};

template<typename Real>
class Radix2RealFft {
 public:
  // default constructor
  Radix2RealFft() : N_(0), Wnk_(NULL), y_(NULL), r2cfft_(NULL) {}

  // desctructor
  ~Radix2RealFft();
  // Copy constructor
  Radix2RealFft(const Radix2RealFft &other);

  /// If forward == true, this function transforms from a sequence of N float
  /// points to its complex fourier
  /// transform; otherwise it goes in the reverse direction.  If you call it
  /// in the forward and then reverse direction and multiply by 1.0/N, you
  /// will get back the original data.
  /// The interpretation of the complex-FFT data is as follows: the array
  /// is a sequence of complex numbers C_n of length N/2 with (float, im)
  /// format,
  /// i.e. [real0, real_{N/2}, real1, im1, real2, im2, real3, im3, ...].
  void Compute(Real *x, int N, bool forward);

  /// This is as the other Compute() function, but it is a const version that
  /// uses a user-supplied buffer.
  // void Compute(float *x, bool forward, std::vector<float> *temp_buffer)
  // const;

 private:
  // Disallow assignment.
  Radix2RealFft &operator=(const Radix2RealFft &other);
  // initialize Wnk_ and r2cfft_
  void Init(int N);
  // Release all the allocated buffer
  void Release();
  int N_;
  Real *Wnk_;
  Real *y_;  // a temperal buffer to avoid endless new-delete
  Radix2ComplexFft<Real> *r2cfft_;
};

#endif