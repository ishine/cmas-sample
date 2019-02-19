// base/kaldi-common.h

// Copyright 2009-2011 Microsoft Corporation

// See ../../COPYING for clarification regarding multiple authors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//  http://www.apache.org/licenses/LICENSE-2.0
//
// THIS CODE IS PROVIDED *AS IS* BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, EITHER EXPRESS OR IMPLIED, INCLUDING WITHOUT LIMITATION ANY IMPLIED
// WARRANTIES OR CONDITIONS OF TITLE, FITNESS FOR A PARTICULAR PURPOSE,
// MERCHANTABLITY OR NON-INFRINGEMENT.
// See the Apache 2 License for the specific language governing permissions and
// limitations under the License.

#ifndef KALDI_BASE_KALDI_COMMON_H_
#define KALDI_BASE_KALDI_COMMON_H_ 1

#include <cstddef>
#include <cstdlib>
#include <cstring>  // C string stuff like strcpy
#include <string>
#include <sstream>
#include <stdexcept>
#include <cassert>
#include <vector>
#include <iostream>
#include <fstream>

#ifdef __ANDROID__
#include <android/log.h>

#define LOGANDROID "spicax" // ������Զ����LOG�ı�ʶ
#define LOGD(...)  __android_log_print(ANDROID_LOG_DEBUG,LOGANDROID,__VA_ARGS__) // ����LOGD����
#define LOGI(...)  __android_log_print(ANDROID_LOG_INFO,LOGANDROID,__VA_ARGS__) // ����LOGI����
#define LOGW(...)  __android_log_print(ANDROID_LOG_WARN,LOGANDROID,__VA_ARGS__) // ����LOGW����
#define LOGE(...)  __android_log_print(ANDROID_LOG_ERROR,LOGANDROID,__VA_ARGS__) // ����LOGE����
#define LOGF(...)  __android_log_print(ANDROID_LOG_FATAL,LOGANDROID,__VA_ARGS__) // ����LOGF����

#endif

#include "kaldi-utils.h"
#include "kaldi-error.h"
#include "kaldi-types.h"
#include "io-funcs.h"
#include "kaldi-math.h"

#endif  // KALDI_BASE_KALDI_COMMON_H_
