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

#ifndef GRID_OUTPUT_HPP
#define GRID_OUTPUT_HPP

#include <iosfwd>
#include <string>

#include "multidimensional_storage.hpp"

class RectilinearGrid3D;

void OutputGridAndData(const std::string& format,
                       const std::string& ascii_or_binary,
                       const RectilinearGrid3D& grid, 
                       const MultiDimensionalStorage4D& data,
                       std::ofstream* os_ptr);

#endif // GRID_OUTPUT_HPP
