project(spicax-cmas)
include_directories(../)
aux_source_directory(. base-src)
set(base-src audio-device-pa.cc beamform-das.cc circular-array-separation.cc
io-funcs.cc kaldi-error.cc kaldi-math.cc kaldi-utils.cc online-wola.cc
parse-options.cc r2fft.cc resample.cc srfft.cc string-split.cc
text-utils.cc wave-header.cc wola.cc)
add_library(cmas-api ${base-src})

find_library(corefound CoreFoundation)
find_library(iokit IOKit)
find_library(coreaudio CoreAudio)
find_library(audiotoolbox AudioToolbox)
find_library(coreservices CoreServices)
find_library(audiounit AudioUnit)
find_library(carbon Carbon)

link_directories("/Users/danhui/project/07-voiceAI/circular-microphone-array-separator/cmas/")
target_link_libraries(cmas-api libportaudio.a ${iokit} ${corefound} ${coreaudio} ${audiotoolbox} ${coreservices} ${carbon} ${audiounit})

add_executable(test-cmas-api test-cmas-api.cc)
target_link_libraries(test-cmas-api cmas-api)

add_executable(test-audio-device-pa test-audio-device-pa.cc)
target_link_libraries(test-audio-device-pa cmas-api)