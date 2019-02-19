/*
// THE WAVE FORMAT IS SPECIFIED IN:
// https:// ccrma.stanford.edu/courses/422/projects/WaveFormat/
//
//
//
//  RIFF
//  |
//  WAVE
//  |    \    \   \
//  fmt_ data ... data
//
//
//  Riff is a general container, which usually contains one WAVE chunk
//  each WAVE chunk has header sub-chunk 'fmt_'
//  and one or more data sub-chunks 'data'
//
//  [Note from Dan: to say that the wave format was ever "specified" anywhere is
//   not quite right.  The guy who invented the wave format attempted to create
//   a formal specification but it did not completely make sense.  And there
//   doesn't seem to be a consensus on what makes a valid wave file,
//   particularly where the accuracy of header information is concerned.]
*/

#ifndef WAVE_HEADER_H_
#define WAVE_HEADER_H_

#include <cstring>
#include <iostream>

/// For historical reasons, we scale waveforms to the range
/// (2^15-1)*[-1, 1], not the usual default DSP range [-1, 1].
// const float kWaveSampleMax = 32768.0;

/// This class reads and hold wave file header information.
class WaveHead {
 public:
  WaveHead(): samp_freq_(0), samp_count_(0), num_channels_(0), header_size_(0),
    bits_per_sample_(0), reverse_bytes_(0) {}

  /// Is stream size unknown? Duration and SampleCount not valid if true.
  // bool IsStreamed() const { return samp_count_ < 0; }

  /// Sample frequency, Hz.
  float SampFreq() const { return samp_freq_; }

  /// Number of samples in stream. Invalid if IsStreamed() is true.
  unsigned int SampleCount() const { return samp_count_; }

  /// Approximate duration, seconds. Invalid if IsStreamed() is true.
  float Duration() const { return samp_count_ / samp_freq_; }

  /// Number of channels, 1 to 16.
  int NumChannels() const { return num_channels_; }

  /// Bytes per sample.
  size_t BlockAlign() const { return 2 * num_channels_; }

  /// Wave data bytes. Invalid if IsStreamed() is true.
  size_t DataBytes() const { return samp_count_ * BlockAlign(); }

  //return header bytes
  int HeaderSize() const {return header_size_;}

  int BitsPerSample()const {return bits_per_sample_;}

  /// Is data file byte order different from machine byte order?
  bool ReverseBytes() const { return reverse_bytes_; }

  /// 'is' should be opened in binary mode. Read() will throw on error.
  /// On success 'is' will be positioned at the beginning of wave data.
  void Read(std::istream &is);
  void Read(char * buffer);
  void Write(std::ostream &os, unsigned int num_samp, int sample_rate = 16000,
             int num_chan = 1, int bytes_per_samp = 2) const;
  void Write(char * buffer, unsigned int num_samp,
             int samplerate = 16000, int num_chan = 1, int bytes_per_samp = 2) const;
//   void Write(std::ostream &os, unsigned int num_samp, int sample_rate);

 private:
  float samp_freq_;
  unsigned int samp_count_;  // 0 if empty, -1 if undefined length.
  // unsigned int datalen_; // in samples
  short num_channels_;
  int header_size_;
  int bits_per_sample_;
  bool reverse_bytes_;  // File endianness differs from host.
};

#endif
