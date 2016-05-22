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

// Intended to hold a family of variables defined on a 3D rectilinear
// grid (aka a volume). Index for family corresponds to the slowest
// dimension. Fastest dimension can be padded.
class MultiDimensionalStorage4D {
public:
  MultiDimensionalStorage4D();
  MultiDimensionalStorage4D(int n1, int n2, int n3, int n4);
  MultiDimensionalStorage4D(int n1, int n2, int n3, int n4, int padding);
  int n_fast() { return m_n_fast; }
  int n_fast() const { return m_n_fast; }
  int n2() { return m_n2; }
  int n2() const { return m_n2; }
  int n3() { return m_n3; }
  int n3() const { return m_n3; }
  int n_slow() { return m_n_slow; }
  int n_slow() const { return m_n_slow; }
  int n_fast_padding() { return m_n_fast_padding; }
  int n_fast_padding() const { return m_n_fast_padding; }
  void Allocate();
  void DeAllocate();
  void Validate();
  void Validate() const;
  RealT* RawDataSlowDimension(int i);
  const RealT* RawDataSlowDimension(int i) const;
private:
  int m_n_fast;
  int m_n2;
  int m_n3;
  int m_n_slow;
  int m_n_fast_padding;
  RealT* m_data;
};


#endif // MULTIDIMENSIONAL_STORAGE_HPP
