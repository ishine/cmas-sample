#include "audio-device-pa.h"
#include <iostream>
#include <fstream>
#include <string>
#include "kaldi-common.h"
#include "parse-options.h"

int main(int argc, char **argv) {

  const char *usage =
    "test-audio-devie-pa \n"
    "Usage: test-audio-devie-pa [options] <wav-file> \n"
    " e.g.: test-audio-devie-pa saved-file-name.wav\n";
  kaldi::ParseOptions po(usage);

  int fs = 16000;
  int channels = 1;
  std::string key_name = "";
  po.Register("expected-ch", &channels, "expected channels");
  po.Register("expected-fs", &fs, "expected sample rate");
  po.Register("expected-device", &key_name, "some key words in the name of expected device");
  po.Read(argc, argv);
  if (po.NumArgs() != 1) {
    po.PrintUsage();
    return -1;
  }

  AudioDevicePa device;
  const char * info = device.ScanAvailableDevice();
  KALDI_LOG << "device_information:\n" << info << "\n";

  device.StartRecording(key_name, fs, channels);

  WaveHead wave_head;
  std::ofstream os;
  std::string out_name = po.GetArg(1);
  os.open(out_name, std::ios::out | std::ios::binary);
  wave_head.Write(os, 3 * fs , fs, channels);
  short * buffer = new short[fs * channels];
  std::memset(buffer, 0x0, sizeof(short) * fs * channels);
  int all_samples = 0;

  int required_samples = fs;
  KALDI_LOG << "start to recording 3 seconds, target_samplerate=" << fs << ",channels=" << channels << "\n";
  for (int i = 0; i < 3; i++) {
    int num_samples = device.ReadPcmData(buffer, required_samples);
    KALDI_LOG << "read " << num_samples << " samples\n";
    KALDI_LOG << "sleep for 1 second to simulate some other operations\n";
    std::this_thread::sleep_for(std::chrono::seconds(1));
    all_samples += num_samples;
    os.write((const char*)buffer, num_samples * channels * sizeof(short));
  }
  os.seekp(0, std::ios::beg);
  wave_head.Write(os, all_samples, fs, channels);
  os.close();
  if (buffer != nullptr) {delete [] buffer; buffer = nullptr;}
  device.StopRecording();
  return 0;
}