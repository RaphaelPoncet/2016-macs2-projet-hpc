// Copyright 2016 Raphael Poncet.

// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at

//      http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef COMMON_HPP
#define COMMON_HPP

#include <vector>

#include "plog/Log.h"
#include "plog/Appenders/ColorConsoleAppender.h"
#include "plog/Appenders/ConsoleAppender.h"
#include "plog/Appenders/RollingFileAppender.h"

template <typename T>
T* VectorRawData(std::vector<T>& vec) { 
 
 return (vec.size() == 0 ? NULL : &vec.at(0));

}

template <typename T>
const T* ConstVectorRawData(const std::vector<T>& vec) { 
 
 return (vec.size() == 0 ? NULL : &vec.at(0));

}

const int MAX_NB_CHARS = 8;

typedef enum _DataType {NONE, FLOAT32, FLOAT64, NB_DATATYPES} DataType;

// In bytes.
static const int DATATYPE_SIZE[NB_DATATYPES] = {0, 4, 8};

// VTK names.
static const char DATATYPE_NAME[NB_DATATYPES][MAX_NB_CHARS] = {"None", "Float32", "Float64"};

#define UNUSED(expr) do { (void)(expr); } while (0)

#endif // COMMON_HPP
