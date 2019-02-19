#ifndef AUDIO_DEVICE_H
#define AUDIO_DEVICE_H

#include "portaudio.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <chrono>
#include <thread>
#include <atomic>
#include "wave-header.h"
#include "spicax-eigen.h"
#include "kaldi-common.h"
#include "circular-buffer.h"
#include "resample.h"

class AudioDevicePa {
 public:
  // Default constructor
  AudioDevicePa():
    record_stream_(nullptr), default_samplerate_(16000),
    frames_per_buffer_(default_samplerate_ / 100), target_samplerate_(16000),
    target_channels_(1), num_record_device_(0), save_file_(true) {};
  // Set requierd sample rate and channels
  void SetRequiredRateChannels(int required_sample_rate, int required_channels) {
    target_samplerate_ = required_sample_rate;
    target_channels_ = required_channels;
  };
  // Destructor, handle all the errors
  ~AudioDevicePa() {StopRecording();};
  // Scan to get information of available recording devices
  const char* ScanAvailableDevice();
  // Startrecording, downsampling if necessary
  bool StartRecording(std::string device_key_name, int required_sample_rate, int required_channels);
  // stop recording
  void StopRecording();
  // pause recording
  void PauseRecording() {}
  // resume recording
  void ResumeRecording() {}
  // read pcm data
  int ReadPcmData(short *buffer, int required_samples) {
    int num_samples = circular_buffer_.PopWithTimeOut(buffer, required_samples * target_channels_, std::chrono::seconds(2));
    return num_samples;
  }
  // actual callback function
  int PushQueue(const void *input, int frames);

 private:
  // Copy contructor
  AudioDevicePa(const AudioDevicePa & other) = delete;
  // Operator=
  AudioDevicePa &operator=(const AudioDevicePa & other) = delete;
  // Prepare to record
  bool PrepareStream(std::string device_key_name);

  PaStreamParameters record_paras_;
  PaStream *record_stream_;

  int default_samplerate_;
  int frames_per_buffer_;
  int target_samplerate_;
  int target_channels_;

  std::atomic_bool record_flag_;

  spicax::MatrixXs orig_data_;
  Eigen::MatrixXf input_;
  Eigen::MatrixXf output_;
  spicax::MatrixXs output_short_;

  // information of all the record devices found
  std::string record_device_info_;
  // number of input devices found
  int num_record_device_;

  bool save_file_;
  std::ofstream os_;

  spicax::CircularBuffer circular_buffer_;
  DownSampler downsampler_;
};

#endif