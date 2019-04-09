cmake_minimum_required (VERSION 2.8)
project(spicax-cmas)

set(CMAKE_F "-Wall -std=c++11 -fPIC -O3")
set(CMAKE_CXX_FLAGS ${CMAKE_F})

aux_source_directory(${CMAKE_CURRENT_LIST_DIR}/cmas spciax-cmas-src)
include_directories(${CMAKE_CURRENT_LIST_DIR}/ ${CMAKE_CURRENT_LIST_DIR}/cmas /Users/danhui/project/03-opensource/eigen-git-mirror/)

add_library(cmas-api ${spciax-cmas-src})

find_library(corefound CoreFoundation)
find_library(iokit IOKit)
find_library(coreaudio CoreAudio)
find_library(audiotoolbox AudioToolbox)
find_library(coreservices CoreServices)
find_library(audiounit AudioUnit)
find_library(carbon Carbon)

link_directories(${CMAKE_CURRENT_LIST_DIR}/cmas)
target_link_libraries(cmas-api libportaudio.a ${iokit} ${corefound} ${coreaudio} ${audiotoolbox} ${coreservices} ${carbon} ${audiounit})

add_executable(test-cmas-api test-cmas-api.cc)
target_link_libraries(test-cmas-api cmas-api)

add_executable(test-audio-device-pa test-audio-device-pa.cc)
target_link_libraries(test-audio-device-pa cmas-api)

add_executable(generate-stationary-noise-matrix generate-stationary-noise-matrix.cc)
target_link_libraries(generate-stationary-noise-matrix cmas-api)