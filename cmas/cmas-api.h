// Six Circular Microphone Array, microphone are arranged in a counterclockwise way
//            mic2--mic1
//          /           \
//      mic3             mic0
//          \           /
//            mic4--mic5
//        Angular positions    beam index
// mic0 :        0 degree        Beam 0
// mic1 :       60 degree        Beam 1
// mic2 :      120 degree        Beam 2
// mic3 :      180 degree        Beam 3
// mic4 :      240 degree        Beam 4
// mic5 :      300 degree        Beam 5

// Initialize circular separator using the required ouput beam index, separated by '|'
// e.g. CsInitialize("0|4") means CsSeparate will get Beam0 and Beam4 at the same time
bool CmasInitialize(const char * expected_beam_index);

// // Initialize circular separator using the required ouput beam index, separated by '|'
// // e.g. CmasSetExpectedBeam("0|4") means CsSeparate will get Beam0 and Beam4 at the same time
// bool CmasSetExpectedBeam(const char * expected_beam_index);

// Do the separation, output the expected beams
// pcm_in, unprocessed audio stream, samples are interleaved,
//         i.e. sample sequence: mic0,mic1,mic2,mic3,mic4,mic5,mic0,mic1,mic2,mic3,mic4,mic5,mic0,...
//                              |<----------(t-1)------------>|<-----------(t)------------->|
// samples_per_chan, the input length of each microphone, required at least a frame-shift.
// pcm_out, interleaved, the output buffer, allocated outside.
// return, how many samples are output, note we process frame-by-frame,
// thus this value must be multiple times of frame-shift(0.016s, 256 samples for 16 kHz)
int CmasSeparate(const char *pcm_in, int samples_per_chan, char *pcm_out);

// Release circular separator
bool CmasRelease();