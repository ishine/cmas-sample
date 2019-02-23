#include "kaldi-common.h"
#include "spicax-eigen.h"
#include "wave-header.h"
#include "cmas-api.h"
#include "timer.h"
#include "parse-options.h"
#include "string-split.h"

using namespace std;
using namespace kaldi;

int main(int argc, char** argv) {
  const char *usage =
    "test-cmas-api. framewise complex Gaussian mixture model MVDR beamformer, a data-driven beamformer.\n"
    " We first compute the probablity of each mixture, and then estimate the speech and noise covariance matrix.\n"
    " Finally, we compute the beamformer weight based on the MVDR/GEV beamformer. Mostly MVDR beamformer is better.\n"
    "Usage: test-cmas-api [options] <wav-in-file> <wav-out-file>\n"
    " e.g.: test-cmas-api in.wav out.wav\n"
    " e.g.: test-cmas-api in.wav\n";
  kaldi::ParseOptions po(usage);

  std::string expect_beam_string = "1|2";
  po.Register("expect-beams", &expect_beam_string, "specify the expected angles, separated by '|', "
              "this match how many channels will be output");
  po.Read(argc, argv);

  if ((po.NumArgs() != 2) && (po.NumArgs() != 1)) {
    po.PrintUsage();
    exit(1);
  }

  std::string input_filename = po.GetArg(1);
  std::ifstream is;
  is.open(input_filename, std::ios::in | std::ios::binary);
  if (!is.is_open()) {
    KALDI_LOG << "Cannot open file " << input_filename << "\n";
    is.close();
    return -1;
  }

  // get the actual bytes of input wave file
  is.seekg(0, is.end);
  size_t in_wave_bytes = is.tellg();
  is.seekg(0, is.beg);

  // check the header of the input wave file
  WaveHead header;
  header.Read(is);
  int num_chan = 6;
  if (header.NumChannels() != num_chan) {
    KALDI_LOG << "Expected " << num_chan << " channels, but actual " << header.NumChannels()
              << " channels for wave file " << input_filename << "\n";
    is.close();
    return -1;
  }
  if ((int)header.SampFreq() != 16000) {
    KALDI_LOG << "Expected sample rate is " << 16000 << " but get " << (int)header.SampFreq() << "\n";
    is.close();
    return -1;
  }

  // compute the samples per channel
  int samples_all = (in_wave_bytes - header.HeaderSize()) / (sizeof(short) * num_chan);

  int frame_shift_samples = 256;
  int fs = 16000;
  // int num_frames = samples_all / frame_shift_samples;
  // build input and output buffer
  const int samples_per_chan = frame_shift_samples;
  int num_blocks = samples_all / samples_per_chan;
  spicax::MatrixXs pcm_in, pcm_out;
  pcm_in.setZero(num_chan, samples_per_chan);
  std::vector<std::string> expec_angles = spicax::split(expect_beam_string, "|");
  int out_chan = expec_angles.size();
  pcm_out.setZero(out_chan, samples_per_chan);

  std::string out_name;
  if (po.NumArgs() == 1) {
    std::string subsuffix = "cmas";
    out_name = input_filename.substr(0, input_filename.size() - 4) + "-separated-" + subsuffix + ".wav";
  }
  if (po.NumArgs() == 2)out_name = po.GetArg(2);
  std::ofstream os;
  os.open(out_name, ios::out | ios::binary);
  header.Write(os, samples_all, fs, out_chan);

  kaldi::Timer ATimer;
  ATimer.Reset();
  int out_samples = 0;
  // KALDI_LOG << "CmasInitialize-1\n";
  CmasInitialize(expect_beam_string.data());
  // KALDI_LOG << "CmasInitialize-2\n";
  // the main for loop function
  for (int i = 0; i < num_blocks; i++) {
    is.read((char*)pcm_in.data(), samples_per_chan * num_chan * sizeof(short));
    int num_samples = CmasSeparate((const char*)pcm_in.data(), samples_per_chan, (char*)pcm_out.data());
    // KALDI_LOG << "num_samples=" << num_samples << "\n";
    os.write(reinterpret_cast<const char*>(pcm_out.data()), sizeof(short) * num_samples * out_chan);
    out_samples += num_samples;
  }
  // KALDI_LOG<<"1\n";
  CmasRelease();
  {
    // zero padding for a better visualization
    // out_samples /= out_chan;
    int zero_padding_len = (samples_all - out_samples);
    // KALDI_LOG << "input_samples=" << samples_all << ", output_samples=" << out_samples
    //           << ",padding_samples=" << zero_padding_len << "\n";
    spicax::VectorXs zero_padding;
    zero_padding.setZero(zero_padding_len * out_chan);
    os.write((const char*)zero_padding.data(), sizeof(short) * zero_padding_len * out_chan);
  }
  float time_elapsed = ATimer.Elapsed();
  KALDI_LOG << "RTF:" <<  time_elapsed / header.Duration() << "\n";
  os.close();
  is.close();
  return 0;
}