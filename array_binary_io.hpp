// Copyright 2016 Raphael Poncet.

// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at

//      http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language goveroning permissions and
// limitations under the License.

#ifndef ARRAY_BINARY_IO_HPP
#define ARRAY_BINARY_IO_HPP

#include <iosfwd>
#include <string>

#include "common.hpp"

typedef enum _GriddedDataFormat {GRIDDED_ASCII, GRIDDED_BINARY} GriddedDataFormat;

void ReadBinaryVariableHeader(std::istream* is_ptr, 
                              std::string* name_ptr, 
                              DataType* data_type_ptr,
                              int* nb_components_ptr, 
                              int* n_fast_ptr,
                              int* n_medium_ptr,
                              int* n_slow_ptr,
                              RealT* x_fast_min_ptr,
                              RealT* x_fast_max_ptr,
                              RealT* x_medium_min_ptr,
                              RealT* x_medium_max_ptr,
                              RealT* x_slow_min_ptr,
                              RealT* x_slow_max_ptr);


void ReadBinaryVariable(std::istream& input_stream, 
                        int n_fast, int n_fast_padding, int n_medium, 
                        int n_slow, int nb_components, RealT* data);

void WriteBinaryVariableHeader(const std::string& variable_name, 
                               DataType data_type, 
                               GriddedDataFormat gridded_format,
                               int nb_components, 
                               int n_fast, 
                               int n_medium, 
                               int n_slow, 
                               RealT x_fast_min, 
                               RealT x_fast_max,
                               RealT x_medium_min, 
                               RealT x_medium_max,
                               RealT x_slow_min, 
                               RealT x_slow_max,
                               std::ostream* os_ptr);

void WriteBinaryVariableSliceSlowDimension(DataType data_type,
                                           GriddedDataFormat gridded_format,
                                           int n_fast, int n_fast_padding, 
                                           int n_medium, int islow,
                                           const RealT* data, std::ostream* os_ptr);

void WriteBinaryVariable(DataType data_type,
                         GriddedDataFormat gridded_format,
                         int nb_components,
                         int n_fast, int n_fast_padding, 
                         int n_medium, int n_slow,
                         const RealT* data, std::ostream* os_ptr);

#endif // ARRAY_BINARY_IO_HPP
