#include "kaldi-common.h"
#include "spicax-eigen.h"
#include "parse-options.h"
#include "wave-header.h"
#include "wola.h"

int main(int argc, char **argv) {
  const char *usage =
    "generate-stationary-noise-matrix \n"
    "Usage: generate-stationary-noise-matrix [options] <wav-in-file> <noise-model>\n"
    " e.g.: generate-stationary-noise-matrix in.wav noise.model\n"
    " e.g.: generate-stationary-noise-matrix in.wav\n";
  kaldi::ParseOptions po(usage);

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
  int num_chan = header.NumChannels();
  if ((int)header.SampFreq() != 16000) {
    KALDI_LOG << "Expected sample rate is " << 16000 << " but get " << (int)header.SampFreq() << "\n";
    is.close();
    return -1;
  }

  // compute the samples per channel
  int samples_all = (in_wave_bytes - header.HeaderSize()) / (sizeof(short) * num_chan);
  int fs = 16000;
  float frame_shift_ms = 0.016f;
  int frame_shift_samples = frame_shift_ms * fs;
  int overlap_rate = 2;
  float frame_len_ms = frame_shift_ms * overlap_rate;
  int fft_len = 1;
  while (fft_len < (frame_len_ms * fs)) {fft_len <<= 1;}

  int out_chan = 1;
  spicax::WolaOptions wola_opts(num_chan, out_chan, fs, frame_shift_ms, frame_len_ms, fft_len);
  spicax::Wola wola;
  wola.Init(wola_opts);

  spicax::MatrixXs pcm_in;
  pcm_in.setZero(num_chan, samples_all);
  is.read((char*)pcm_in.data(), samples_all * num_chan * sizeof(short));
  const Eigen::MatrixXcf& in_spec = wola.Decompose(pcm_in.data(), samples_all);
  is.close();

  std::string out_name;
  if (po.NumArgs() == 1) {
    std::string subsuffix = std::to_string(num_chan) + "-" + std::to_string(frame_shift_samples) + "-" +
                            std::to_string(fft_len) + "-" + std::to_string(fs);
    out_name = "noise-" + subsuffix + ".model";
  }
  if (po.NumArgs() == 2)out_name = po.GetArg(2);

  std::ofstream os;
  os.open(out_name, std::ios::out | std::ios::binary);

  std::vector<Eigen::MatrixXcf> noise_cov;
  int num_bins = fft_len / 2 + 1;
  noise_cov.resize(num_bins);
  for (int ibin = 0; ibin < num_bins; ibin++) {
    noise_cov[ibin].setZero(num_chan, num_chan);
  }
  int num_frames = in_spec.cols();
  KALDI_LOG << "num_frames=" << num_frames << "\n";
  Eigen::VectorXcf obs;
  for (int iframe = 0; iframe < num_frames; iframe++) {
    Eigen::Map<const Eigen::MatrixXcf> array_frame(&(in_spec(0, iframe)), num_bins, num_chan);
    for (int ibin = 0; ibin < num_bins; ibin++) {
      obs = array_frame.row(ibin) / 32767.0f;
      noise_cov[ibin] += obs * obs.adjoint();
    }
  }
  for (int ibin = 0; ibin < num_bins; ibin++) {
    noise_cov[ibin] /= num_frames;
    os.write((const char*)(noise_cov[ibin].data()), num_chan * num_chan * sizeof(float) * 2);
  }
  os.close();
  return -1;
}