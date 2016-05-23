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

#ifndef WAVE_SOLVER_OPTIONS_HPP
#define WAVE_SOLVER_OPTIONS_HPP

#include <iosfwd>
#include <vector>

class RectilinearGrid3D;

class WaveSolverOptions {
public:
  WaveSolverOptions();
  WaveSolverOptions(int sponge_width, RealT sponge_base_value, RealT sponge_slope,
                    const RectilinearGrid3D& grid);
  void InitSponge(const RectilinearGrid3D& grid);
  void SetSponge();
  void Set(int sponge_width, RealT sponge_base_value, RealT sponge_slope,
           const RectilinearGrid3D& grid);
  void DumpAllSponges(std::ostream* os_fast_ptr,
                      std::ostream* os_medium_ptr,
                      std::ostream* os_slow_ptr);
  std::vector<RealT>& sponge_slow() {return m_sponge_slow; }
  const std::vector<RealT>& sponge_slow() const {return m_sponge_slow; }
  std::vector<RealT>& sponge_medium() {return m_sponge_medium; }
  const std::vector<RealT>& sponge_medium() const {return m_sponge_medium; }
  std::vector<RealT>& sponge_fast() {return m_sponge_fast; }
  const std::vector<RealT>& sponge_fast() const {return m_sponge_fast; }
private:
  int m_sponge_width;
  RealT m_sponge_base_value;
  RealT m_sponge_slope;
  std::vector<RealT> m_sponge_fast;
  std::vector<RealT> m_sponge_medium;
  std::vector<RealT> m_sponge_slow;
};

#endif // WAVE_SOLVER_OPTIONS_HPP
