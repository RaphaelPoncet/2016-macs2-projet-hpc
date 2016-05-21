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

#ifndef RECTILINEAR_GRID_HPP
#define RECTILINEAR_GRID_HPP

#include <fstream>
#include <vector>

// Rectilinear 3D grid (also handles 2D as a particular case).
class RectilinearGrid3D {
public:
  RectilinearGrid3D(RealT x_fast_min, RealT x_fast_max, int n_fast,
                    RealT x_medium_min, RealT x_medium_max, int n_medium,
                    RealT x_slow_min, RealT x_slow_max, int n_slow);
  int n_fast() const { return m_n_fast; }
  int n_medium() const { return m_n_medium; }
  int n_slow() const { return m_n_slow; }
  int x_fast_min() const {return m_x_fast_min; }
  int x_fast_max() const {return m_x_fast_max; }
  int x_medium_min() const {return m_x_medium_min; }
  int x_medium_max() const {return m_x_medium_max; }
  int x_slow_min() const {return m_x_slow_min; }
  int x_slow_max() const {return m_x_slow_max; }
  void WriteHeaderVTKXml(std::ofstream* os_ptr) const;
  void WriteHeaderVTKXml(std::ofstream* os_ptr);
  void WriteFooterVTKXml(std::ofstream* os_ptr) const;
  void WriteFooterVTKXml(std::ofstream* os_ptr);
  void WriteVTKXmlAscii(std::ofstream* os_ptr) const;
  void WriteVTKXmlAscii(std::ofstream* os_ptr);
private:
  int m_n_fast;
  int m_n_medium;
  int m_n_slow;
  RealT m_x_fast_min;
  RealT m_x_fast_max;
  RealT m_x_medium_min;
  RealT m_x_medium_max;
  RealT m_x_slow_min;
  RealT m_x_slow_max;
  std::vector<RealT> m_fast_coordinates;
  std::vector<RealT> m_medium_coordinates;
  std::vector<RealT> m_slow_coordinates;
};

#endif // RECTILINEAR_GRID_HPP
