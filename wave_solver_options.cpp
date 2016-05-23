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

#include "wave_solver_options.hpp"

#include <cassert>
#include <fstream>
#include <vector>

#include "common.hpp"
#include "rectilinear_grid.hpp"

static void DumpArrayAscii(int n, RealT* sponge_array, std::ostream* os_ptr) {

  assert(0 < n);
  assert(sponge_array != NULL);

  if (!*os_ptr) {

    LOG_ERROR << "Invalid output stream";
    std::abort();

  }

  for (int i = 0; i < n; ++i)
    *os_ptr << sponge_array[i] << "\n";

}

static void InitSpongeArray(int sponge_width, RealT base_value, RealT slope, int n, RealT* out) {

  assert(0 <= sponge_width);
  assert(0 < n);
  assert(out != NULL);

  // We implement the method of Cerjan et al. (1985): we multiply the
  // wave amplitudes by a smooth decaying factor. Choice of parameters
  // is empirical.

  for (int i = 0; i < n; ++i)
    out[i] = 1.0;

  if (sponge_width < n - sponge_width) {

    for (int i = 0; i < sponge_width; ++i) {
     
      const RealT x = (RealT)((sponge_width - i) * 20 / sponge_width);
      out[i] = base_value + (1.0 - base_value) * expf(- slope * x * x);

    }

    for (int i = n - sponge_width; i < n; ++i) {
     
      const RealT x = (RealT)((n - sponge_width - i) * 20 / sponge_width);
      out[i] = base_value + (1.0 - base_value) * expf(- slope * x * x);

    }

  } else {

    // The array is too short to support a sponge. Do not attenuate
    // waves in this case.
    LOG_INFO << "Will not initalise sponge layer: "
             << "sponge_width ("
             << sponge_width << ") larger than half length (" << n << ")";

  }
}

WaveSolverOptions::WaveSolverOptions():
  m_sponge_width(0), m_sponge_base_value(1.0), m_sponge_slope(0.0) {};

WaveSolverOptions::WaveSolverOptions(int sponge_width, RealT sponge_base_value, 
                                     RealT sponge_slope, const RectilinearGrid3D& grid):

  m_sponge_width(sponge_width), m_sponge_base_value(sponge_base_value), 
  m_sponge_slope(sponge_slope) {

  UNUSED(grid);

};

void WaveSolverOptions::Set(int sponge_width, RealT sponge_base_value, 
                            RealT sponge_slope, const RectilinearGrid3D& grid) {

  m_sponge_width = sponge_width;
  m_sponge_base_value = sponge_base_value;
  m_sponge_slope = sponge_slope;

  UNUSED(grid);

}


void WaveSolverOptions::InitSponge(const RectilinearGrid3D& grid) {

  const int n_fast = grid.n_fast();
  const int n_medium = grid.n_medium();
  const int n_slow = grid.n_slow();
  
  m_sponge_fast.assign(n_fast, 0.0);
  m_sponge_medium.assign(n_medium, 0.0);
  m_sponge_slow.assign(n_slow, 0.0);
  
}

void WaveSolverOptions::SetSponge() {

  InitSpongeArray(m_sponge_width, m_sponge_base_value, m_sponge_slope, 
                  m_sponge_fast.size(), VectorRawData<RealT>(m_sponge_fast));

  InitSpongeArray(m_sponge_width, m_sponge_base_value, m_sponge_slope, 
                  m_sponge_medium.size(), VectorRawData<RealT>(m_sponge_medium));

  InitSpongeArray(m_sponge_width, m_sponge_base_value, m_sponge_slope, 
                  m_sponge_slow.size(), VectorRawData<RealT>(m_sponge_slow));

}

void WaveSolverOptions::DumpAllSponges(std::ostream* os_fast_ptr,
                                       std::ostream* os_medium_ptr,
                                       std::ostream* os_slow_ptr) {

  DumpArrayAscii(m_sponge_fast.size(), VectorRawData<RealT>(m_sponge_fast), os_fast_ptr);
  DumpArrayAscii(m_sponge_medium.size(), VectorRawData<RealT>(m_sponge_medium), os_medium_ptr);
  DumpArrayAscii(m_sponge_slow.size(), VectorRawData<RealT>(m_sponge_slow), os_slow_ptr);

}
