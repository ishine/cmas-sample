#include "audio-device-pa.h"
// A global callback function due to PortAudio
int RecordCallback(const void *input, void *output, unsigned long frames,
                   const PaStreamCallbackTimeInfo *timeInfo,
                   PaStreamCallbackFlags statusFlags, void *userData) {
  AudioDevicePa *audio_device = reinterpret_cast<AudioDevicePa *>(userData);
  return audio_device->PushQueue(input, frames);
}
// Scan all the available device
const char* AudioDevicePa::ScanAvailableDevice() {
  // Case 1: scan device during recording, return the pre-saved information
  if (record_stream_ != nullptr) {
    if (Pa_IsStreamActive(record_stream_)) {return record_device_info_.c_str();}
    Pa_Terminate();
  }
  // Make sure only one Port Audio
  PaError paerr = Pa_Initialize();
  if (paerr != paNoError) {
    KALDI_LOG << "Can not find ALSA\n";
    return record_device_info_.c_str();
  }
  int num_devices = Pa_GetDeviceCount();
  num_record_device_ = 0;
  for (int i = 0; i < num_devices; i++) {
    const PaDeviceInfo * device_info = Pa_GetDeviceInfo(i);
    if (device_info->maxInputChannels > 0) {num_record_device_++;}
    std::cout << "[" << i << "]: name = " << device_info->name
              << ", hostApi=" << Pa_GetHostApiInfo(device_info->hostApi)->name
              << ", maxInputChannels = " << device_info->maxInputChannels
              << ", maxOutputChannels=" << device_info->maxOutputChannels << ".\n";
  }
  //Close PortAudio
  std::stringstream ss;
  // separate by '|'
  ss << num_record_device_; ss << "|";

  for (int i = 0; i < num_devices; i++) {
    const PaDeviceInfo * device_info = Pa_GetDeviceInfo(i);
    if (device_info->maxInputChannels > 0) {
      ss << device_info->name; ss << "|";
      ss << device_info->maxInputChannels; ss << "|";
      ss << (int)device_info->defaultSampleRate; ss << "|";
    }
  }
  getline(ss, record_device_info_);
  Pa_Terminate();
  return record_device_info_.c_str();
}
// Start recording
bool AudioDevicePa::StartRecording(std::string device_key_name, int required_sample_rate, int required_channels) {
  // Case 1: record_stream_ is not null
  if ((record_stream_ != nullptr) ) {
    // Case 1-1: stream is recording, do nothing
    if (Pa_IsStreamActive(record_stream_)) {
      KALDI_LOG << "Stream is recording\n";
      return true;
    }
    // Case 1-2: some errors Occurs, reset first
    StopRecording();
  }
  // Case 2: list all the record device if necessary
  if (num_record_device_ == 0) {ScanAvailableDevice();}
  // Set required sample rate and recording channels
  SetRequiredRateChannels(required_sample_rate, required_channels);
  Pa_Terminate();
  // Prepare recording stream
  PaError paerr = Pa_Initialize();
  if (paerr != paNoError) {KALDI_LOG << Pa_GetErrorText(paerr) << "\n";}
  if (!PrepareStream(device_key_name)) {return false;}
  // Open a stream for recording
  paerr = Pa_OpenStream(&record_stream_, &record_paras_, NULL, (double)default_samplerate_, frames_per_buffer_, paClipOff, RecordCallback, this);
  if (paerr != paNoError) {
    KALDI_LOG << "Error occurs:" << Pa_GetErrorText(paerr) << "\n";
    return false;
  }
  // Start the recording stream
  paerr = Pa_StartStream(record_stream_);
  if (paerr != paNoError) {
    KALDI_LOG << "Error occurs:" << Pa_GetErrorText(paerr) << "\n";
    return false;
  }
  record_flag_ = true;
  KALDI_LOG << "Start recording ...\n";
  return true;
};
// Stop recording
void AudioDevicePa::StopRecording() {
  // Case 1: stream is already stopped
  if (record_stream_ == nullptr) {return ;}
  // Case 2: stream is recording, let's stop it
  if (Pa_IsStreamActive(record_stream_) == 1) {Pa_StopStream(record_stream_);}
  Pa_CloseStream(record_stream_);
  Pa_Terminate();
  downsampler_.Reset();
  circular_buffer_.Resize(0);
  orig_data_.resize(0, 0);
  output_short_.resize(0, 0);
  input_.resize(0, 0);
  output_.resize(0, 0);
  num_record_device_ = 0;
  record_flag_ = false;
  record_stream_ = nullptr;
  if (os_.is_open()) {os_.close();}
  KALDI_LOG << "Stop recording\n";
};
// Main Callback function
int AudioDevicePa::PushQueue(const void *input, int frames) {
  std::memcpy(orig_data_.data(), input, frames_per_buffer_ * record_paras_.channelCount * sizeof(short));
  if (os_.is_open()) {os_.write((const char*)input, sizeof(short) * frames * record_paras_.channelCount);}
  if (default_samplerate_ == target_samplerate_) {
    if (target_channels_ != record_paras_.channelCount) {
      output_short_ = orig_data_.topRows(target_channels_);
    }
    circular_buffer_.PushWithOverflow(output_short_.data(), frames_per_buffer_ * target_channels_);
  } else {
    //Allow Downsampling
    for (int i = 0 ; i < frames_per_buffer_; i++) {
      for (int j = 0; j < target_channels_; j++) {input_(j, i) = (float)orig_data_(j, i);}
    }
    downsampler_.Resample(input_, output_);
    int num_samples = output_.cols();
    for (int i = 0; i < num_samples; i++) {
      for (int j = 0; j < target_channels_; j++) {output_short_(j, i) = (short)output_(j, i);}
    }
    circular_buffer_.PushWithOverflow(output_short_.data(), num_samples * target_channels_);
  }
  return paContinue;
};
// Prepare and initialize all the variable to record audio stream
bool AudioDevicePa::PrepareStream(std::string device_key_name) {
  // Let's
  int num_devices = Pa_GetDeviceCount();
  if (device_key_name.empty()) {
    record_paras_.device = Pa_GetDefaultInputDevice();
  } else {
    record_paras_.device = -1;
    for (int i = 0; i < num_devices; i++) {
      std::string device_name(Pa_GetDeviceInfo(i)->name);
      if (device_name.find(device_key_name) != std::string::npos) {
        record_paras_.device = i;
        break;
      }
    }
  }
  if (record_paras_.device == -1) {return false;}
  record_paras_.channelCount = Pa_GetDeviceInfo(record_paras_.device)->maxInputChannels;
  if (record_paras_.channelCount == 0) {
    KALDI_LOG << "Error, maxInputChannels of " << device_key_name << "=" << Pa_GetDeviceInfo(record_paras_.device)->maxInputChannels << " \n";
    return false;
  }
  record_paras_.sampleFormat = paInt16;
  record_paras_.suggestedLatency = Pa_GetDeviceInfo(record_paras_.device)->defaultLowInputLatency;
  record_paras_.hostApiSpecificStreamInfo = NULL;

  default_samplerate_ = Pa_GetDeviceInfo(record_paras_.device)->defaultSampleRate;
  frames_per_buffer_ = default_samplerate_ / 100;
  if (target_channels_ > record_paras_.channelCount) {target_channels_ = record_paras_.channelCount;}
  downsampler_.Init(default_samplerate_, target_samplerate_, target_channels_);
  circular_buffer_.Resize(target_samplerate_ * target_channels_ * 20);
  orig_data_.setZero(record_paras_.channelCount, frames_per_buffer_);
  input_.setZero(target_channels_, frames_per_buffer_);
  output_.setZero(target_channels_, frames_per_buffer_);
  int num_samples = frames_per_buffer_ * target_samplerate_ / default_samplerate_;
  output_short_.setZero(target_channels_, num_samples);
  if (save_file_) {os_.open("org.pcm", std::ios::out | std::ios::binary);}
  return true;
}
