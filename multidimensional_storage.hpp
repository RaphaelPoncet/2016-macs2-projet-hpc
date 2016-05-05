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

#ifndef MULTIDIMENSIONAL_STORAGE_HPP
#define MULTIDIMENSIONAL_STORAGE_HPP

// Intended to hold a family of variables defined on a 2D cartesian
// grid (aka an image). Index for family corresponds to the slowest
// dimension. Fastest dimension can be padded.
class MultiDimensionalStorage3D {
public:
  MultiDimensionalStorage3D(int n1, int n2, int n3);
  MultiDimensionalStorage3D(int n1, int n2, int n3, int padding);
  void Allocate();
  void DeAllocate();
  void Validate();
  RealT* m_data;
private:
  int m_n_fast;
  int m_n2;
  int m_n_slow;
  int m_padding;

};

#endif // MULTIDIMENSIONAL_STORAGE_HPP
