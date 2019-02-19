#include "wave-header.h"
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <limits>
#include <sstream>
#include <cstring>
#include <string>
#include <cstddef>
#include <vector>

// A utility class for reading wave header.
struct WaveHeaderReadGofer {
  std::istream &is;
  bool swap;
  char tag[5];

  WaveHeaderReadGofer(std::istream &is) : is(is), swap(false) { memset(tag, '\0', sizeof tag); }

  void Expect4ByteTag(const char *expected) {
    is.read(tag, 4);
    if (is.fail()) std::cout << "WaveData: expected " << expected << ", failed to read anything";
    if (strcmp(tag, expected)) std::cout << "WaveData: expected " << expected << ", got " << tag;
  }

  void Read4ByteTag() {
    is.read(tag, 4);
    if (is.fail()) std::cout << "WaveData: expected 4-byte chunk-name, got read error";
  }

  unsigned int ReadUint32() {
    union {
      char result[4];
      unsigned int ans;
    } u;
    is.read(u.result, 4);
    // std::cout << "(int)u.ans =" << u.ans << std::endl;
    if (is.fail()) std::cout << "WaveData: unexpected end of file or read error";
    return u.ans;
  }

  unsigned short ReadUint16() {
    union {
      char result[2];
      unsigned short ans;
    } u;
    is.read(u.result, 2);

    if (is.fail()) std::cout << "WaveData: unexpected end of file or read error";
    return u.ans;
  }
};

static void WriteUint32(std::ostream &os, unsigned int i) {
  union {
    char buf[4];
    unsigned int i;
  } u;
  u.i = i;
  os.write(u.buf, 4);
  if (os.fail()) std::cout << "WaveData: error writing to stream.";
}

static void WriteUint16(std::ostream &os, unsigned short i) {
  union {
    char buf[2];
    unsigned short i;
  } u;
  u.i = i;

  os.write(u.buf, 2);
  if (os.fail()) std::cout << "WaveData: error writing to stream.";
}

void WaveHead::Read(std::istream &is) {
  WaveHeaderReadGofer reader(is);
  reader.Read4ByteTag();
  if (strcmp(reader.tag, "RIFF") == 0)
    reverse_bytes_ = false;
  else if (strcmp(reader.tag, "RIFX") == 0)
    reverse_bytes_ = true;
  else
    std::cout << "WaveData: expected RIFF or RIFX, got " << reader.tag;

  reader.swap = reverse_bytes_;

  unsigned int riff_chunk_size = reader.ReadUint32();
  reader.Expect4ByteTag("WAVE");

  // std::cout << "1\n";
  unsigned int riff_chunk_read = 0;
  riff_chunk_read += 4;  // WAVE included in riff_chunk_size.

  reader.Expect4ByteTag("fmt ");
  unsigned int subchunk1_size = reader.ReadUint32();
  unsigned short audio_format = reader.ReadUint16();
  num_channels_ = reader.ReadUint16();
  unsigned int sample_rate = reader.ReadUint32();
  unsigned int byte_rate = reader.ReadUint32();
  unsigned int block_align = reader.ReadUint16();
  unsigned int bits_per_sample = reader.ReadUint16();
  samp_freq_ = static_cast<float>(sample_rate);

  // std::cout << "2\n";

  unsigned int fmt_chunk_read = 16;
  if (audio_format == 1) {
    if (subchunk1_size < 16) {
      std::cout << "WaveData: expect PCM format data to have fmt chunk "
                << "of at least size 16.";
    }
  } else if (audio_format == 0xFFFE) {  // WAVE_FORMAT_EXTENSIBLE
    unsigned short extra_size = reader.ReadUint16();
    if (subchunk1_size < 40 || extra_size < 22) {
      std::cout << "WaveData: malformed WAVE_FORMAT_EXTENSIBLE format data.";
    }
    reader.ReadUint16();  // Unused for PCM.
    reader.ReadUint32();  // Channel map: we do not care.
    unsigned int guid1 = reader.ReadUint32(), guid2 = reader.ReadUint32(), guid3 = reader.ReadUint32(),
                 guid4 = reader.ReadUint32();
    fmt_chunk_read = 40;

    // Support only KSDATAFORMAT_SUBTYPE_PCM for now. Interesting formats:
    // ("00000001-0000-0010-8000-00aa00389b71", KSDATAFORMAT_SUBTYPE_PCM)
    // ("00000003-0000-0010-8000-00aa00389b71", KSDATAFORMAT_SUBTYPE_IEEE_FLOAT)
    // ("00000006-0000-0010-8000-00aa00389b71", KSDATAFORMAT_SUBTYPE_ALAW)
    // ("00000007-0000-0010-8000-00aa00389b71", KSDATAFORMAT_SUBTYPE_MULAW)
    if (guid1 != 0x00000001 || guid2 != 0x00100000 || guid3 != 0xAA000080 || guid4 != 0x719B3800) {
      std::cout << "WaveData: unsupported WAVE_FORMAT_EXTENSIBLE format.";
    }
  } else {
    std::cout << "WaveData: can read only PCM data, format id in file is: " << audio_format;
  }
  // std::cout << "3\n";
  for (unsigned int i = fmt_chunk_read; i < subchunk1_size; ++i) is.get();  // use up extra data.

  if (num_channels_ == 0) std::cout << "WaveData: no channels present";
  if (bits_per_sample != 16) std::cout << "WaveData: unsupported bits_per_sample = " << bits_per_sample;
  bits_per_sample_ = bits_per_sample;
  if (byte_rate != sample_rate * bits_per_sample / 8 * num_channels_)
    std::cout << "Unexpected byte rate " << byte_rate << " vs. " << sample_rate << " * " << (bits_per_sample / 8)
              << " * " << num_channels_;
  if (block_align != num_channels_ * bits_per_sample / 8)
    std::cout << "Unexpected block_align: " << block_align << " vs. " << num_channels_ << " * "
              << (bits_per_sample / 8);// << std::endl;

  riff_chunk_read += 8 + subchunk1_size;
  // size of what we just read, 4 bytes for "fmt " + 4
  // for subchunk1_size + subchunk1_size itself.

  // We support an optional "fact" chunk (which is useless but which
  // we encountered), and then a single "data" chunk.

  reader.Read4ByteTag();
  riff_chunk_read += 4;

  // Skip any subchunks between "fmt" and "data".  Usually there will
  // be a single "fact" subchunk, but on Windows there can also be a
  // "list" subchunk.
  while (strcmp(reader.tag, "data") != 0) {
    // We will just ignore the data in these chunks.
    unsigned int chunk_sz = reader.ReadUint32();
    if (chunk_sz != 4 && strcmp(reader.tag, "fact") == 0) std::cout << "Expected fact chunk to be 4 bytes long.";
    for (unsigned int i = 0; i < chunk_sz; i++) is.get();
    riff_chunk_read += 4 + chunk_sz;  // for chunk_sz (4) + chunk contents (chunk-sz)

    // Now read the next chunk name.
    reader.Read4ByteTag();
    riff_chunk_read += 4;
  }
  // std::cout << "header size = " << riff_chunk_read - 4 << "\n";
  // header_size_ = riff_chunk_read + 4;
  if (strcmp(reader.tag, "data")) std::cout << "WaveData: expected data chunk, got instead " << reader.tag;

  unsigned int data_chunk_size = reader.ReadUint32();
  // std::cout << "data_chunk_size=" << data_chunk_size << std::endl;
  // datalen_ = data_chunk_size / (bits_per_sample/8);
  riff_chunk_read += 4;
  header_size_ = riff_chunk_read + 8;
  // Figure out if the file is going to be read to the end. Values as
  // observed in the wild:
  bool is_stream_mode = riff_chunk_size == 0 || riff_chunk_size == 0xFFFFFFFF || data_chunk_size == 0 ||
                        data_chunk_size == 0xFFFFFFFF || data_chunk_size == 0x7FFFF000;  // This value is used by SoX.

  if (is_stream_mode)
    std::cout << "Read in RIFF chunk size: " << riff_chunk_size << ", data chunk size: " << data_chunk_size
              << ". Assume 'stream mode' (reading data to EOF).";

  if (!is_stream_mode && std::abs(static_cast<long long>(riff_chunk_read) + static_cast<long long>(data_chunk_size) -
                                  static_cast<long long>(riff_chunk_size)) > 1) {
    // We allow the size to be off by one without warning, because there is a
    // weirdness in the format of RIFF files that means that the input may
    // sometimes be padded with 1 unused byte to make the total size even.
    std::cout << "Expected " << riff_chunk_size << " bytes in RIFF chunk, but "
              << "after first data block there will be " << riff_chunk_read << " + " << data_chunk_size << " bytes "
              << "(we do not support reading multiple data chunks)." ;//<< std::endl;
    std::cout << std::endl;
  }

  if (is_stream_mode)
    samp_count_ = -1;
  else {
    // samp_count_ = data_chunk_size / block_align;
    samp_count_ = data_chunk_size / (num_channels_ * bits_per_sample / 8);
  }
}

void WaveHead::Read(char *buffer) {
  if (buffer == nullptr) {
    std::cout << "buffer is nullptr\n";
    return ;
  }
  //RIFF chunk
  //"RIFF" = 0x52494646, big end == 0x46464952 little end
  int riff_chunk_read = 0;
  unsigned int tag = *((int *)(buffer));
  if (tag != 0x46464952) {
    std::cout << "Expected RIFF tag.\n";
    return ;
  }
  riff_chunk_read += 4;
  unsigned int riff_chunk_size = *((int *)(buffer + riff_chunk_read));
  riff_chunk_read += 4;


  //"WAVE" 0x57415645 = 0x45564157
  tag = *((int *)(buffer + riff_chunk_read));
  if (tag != 0x45564157) {
    std::cout << "Expected WAVE tag\n";
    return ;
  }
  riff_chunk_read += 4;

  //fmt chunk
  //"fmt " = 0x666d7420 = 0x20746d66
  tag = *((int *)(buffer + riff_chunk_read));
  if (tag != 0x20746d66) {
    std::cout << "Expected fmt tag\n";
    return ;
  }
  riff_chunk_read += 4;

  int fmt_chunk_size = *((int *)(buffer + riff_chunk_read));
  // std::cout << "fmt_chunk_size=" << fmt_chunk_size << "\n";
  if ( fmt_chunk_size != 16 && fmt_chunk_size != 18 && fmt_chunk_size != 40) {
    std::cout << "fmt_chunk_size=" << fmt_chunk_size << "\n";
    return ;
  }
  riff_chunk_read += 4;
  int fmt_chunk_read = riff_chunk_read;

  short fmt_type = *((short *)(buffer + riff_chunk_read));
  // std::cout << "fmt_type=" << fmt_type << "\n";
  if (fmt_type != 1) {
    std::cout << "not pcm format, currently we support pcm format only\n";
    return ;
  }
  riff_chunk_read += 2;

  num_channels_ = *((short *)(buffer + riff_chunk_read));
  riff_chunk_read += 2;

  int sample_rate = *((int *)(buffer + riff_chunk_read));
  samp_freq_ = (float)sample_rate;
  riff_chunk_read += 4;

  int byte_rate = *((int *)(buffer + riff_chunk_read));
  riff_chunk_read += 4;

  short block_bytes = *((short *)(buffer + riff_chunk_read));
  if (byte_rate != block_bytes * sample_rate) {
    std::cout << "byte_rate != block_bytes * sample_rate \n";
    return ;
  }
  riff_chunk_read += 2;

  short bits_per_sample = *((short *)(buffer + riff_chunk_read));
  std::cout << "bits_per_sample=" << bits_per_sample << "\n";
  riff_chunk_read += 2;

  riff_chunk_read = fmt_chunk_read + fmt_chunk_size;
  std::cout << "riff_chunk_read=" << riff_chunk_read << "\n";

  //skip all chunk not data
  //"data" = 0x64617461 = 0x61746164
  int chunk_name = *((int *)(buffer + riff_chunk_read));
  while (chunk_name != 0x61746164) {
    std::cout << "chunk_name is not data, skip this chunk\n";
    riff_chunk_read += 4;
    int other_size = *((int *)(buffer + riff_chunk_read));
    riff_chunk_read += other_size;
    chunk_name = *((int *)(buffer + riff_chunk_read));
  }

  std::cout << "riff_chunk_read=" << riff_chunk_read << "\n";
  riff_chunk_read += 4;
  size_t data_chunk_size = *((int *)(buffer + riff_chunk_read));
  samp_count_ = data_chunk_size / (num_channels_ * bits_per_sample / 8);
  if ((data_chunk_size + riff_chunk_read - 4) != riff_chunk_size) {
    std::cout << "Expected " << riff_chunk_size << " bytes in RIFF chunk, but after first data block there will be "
              << riff_chunk_read << " + " << data_chunk_size << " bytes\n";
  }
}

void WaveHead::Write(std::ostream &os, unsigned int num_samp, int sample_rate, int num_chan, int bytes_per_samp) const {
  if (num_samp == 0) std::cout << "Error: attempting to write empty WAVE file";

  os << "RIFF";
  // int num_chan = 1, bytes_per_samp = 2;
  unsigned int subchunk2size = num_samp * bytes_per_samp * num_chan;
  unsigned int chunk_size = 36 + subchunk2size;
  WriteUint32(os, chunk_size);
  os << "WAVE";
  os << "fmt ";
  WriteUint32(os, 16);
  WriteUint16(os, 1);
  WriteUint16(os, num_chan);
  WriteUint32(os, sample_rate);
  WriteUint32(os, sample_rate * num_chan * bytes_per_samp);
  WriteUint16(os, num_chan * bytes_per_samp);
  WriteUint16(os, 8 * bytes_per_samp);
  os << "data";
  WriteUint32(os, subchunk2size);
}

void WaveHead::Write(char * buffer, unsigned int num_samp, int sample_rate_, int num_chan, int bytes_per_samp)const {
  // RIFF/WAVE header
  buffer[0] = 'R';  buffer[1] = 'I'; buffer[2] = 'F'; buffer[3] = 'F';
  unsigned int total_len = num_samp * num_chan + 36;
  buffer[4] = (total_len & 0xff); buffer[5] = ((total_len >> 8) & 0xff);
  buffer[6] = ((total_len >> 16) & 0xff); buffer[7] = ((total_len >> 24) & 0xff);
  buffer[8] = 'W'; buffer[9] = 'A'; buffer[10] = 'V'; buffer[11] = 'E';
  // 'fmt ' chunk
  buffer[12] = 'f'; buffer[13] = 'm'; buffer[14] = 't'; buffer[15] = ' ';
  // 4 bytes: size of 'fmt ' chunk
  buffer[16] = 16; buffer[17] = 0; buffer[18] = 0; buffer[19] = 0;
  // format = 1
  buffer[20] = 1; buffer[21] = 0;
  buffer[22] = num_chan; buffer[23] = 0;
  buffer[24] = (sample_rate_ & 0xff); buffer[25] = ((sample_rate_ >> 8) & 0xff);
  buffer[26] = ((sample_rate_ >> 16) & 0xff); buffer[27] = ((sample_rate_ >> 24) & 0xff);
  int byteRate = num_chan * sample_rate_ * sizeof(short);
  buffer[28] = (byteRate & 0xff); buffer[29] = ((byteRate >> 8) & 0xff);
  buffer[30] = ((byteRate >> 16) & 0xff); buffer[31] = ((byteRate >> 24) & 0xff);
  // block align
  buffer[32] = (num_chan * sizeof(short)); buffer[33] = 0;
  // bits per sample
  buffer[34] = 16; buffer[35] = 0;
  buffer[36] = 'd'; buffer[37] = 'a'; buffer[38] = 't'; buffer[39] = 'a';
  buffer[40] = (num_samp & 0xff); buffer[41] = ((num_samp >> 8) & 0xff);
  buffer[42] = ((num_samp >> 16) & 0xff); buffer[43] = ((num_samp >> 24) & 0xff);
}

