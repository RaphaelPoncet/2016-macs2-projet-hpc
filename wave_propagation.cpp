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

#include <cmath>
#include <fstream>

#include "common.hpp"
#include "wave_propagation.hpp"


// Basic absorbing boundary conditions.
void ApplySpongeLayer_0(int n_fast, int n_fast_padding, 
                        int n_medium, int n_slow,
                        const RealT* sponge_fast,
                        const RealT* sponge_medium,
                        const RealT* sponge_slow,
                        RealT* data) {

  const int n_fast_pad = n_fast + n_fast_padding;

  // Sponge layer.
  for (int islow = 0; islow < n_slow; ++islow) {
    for (int imedium = 0; imedium < n_medium; ++imedium) {
      for (int ifast = 0; ifast < n_fast; ++ifast) {
        
        const size_t index = 
          n_medium * n_fast_pad * islow + n_fast_pad * imedium + ifast;
        
        const RealT scaling = sponge_slow[islow] * sponge_medium[imedium] * sponge_fast[ifast];

        data[index] *= scaling;

      }
    }
  }
}

void ComputeLaplacian_0(int n_fast, int n_fast_padding, 
                        int n_medium, int n_slow,
                        int stencil_radius,
                        RealT dx, RealT dy, RealT dz, 
                        const RealT* p, RealT* laplace_p) {
  
  UNUSED(dy);
  assert(stencil_radius == 1);
  
  const int n_fast_pad = n_fast + n_fast_padding;

  const bool two_dimensional_mode = (n_medium == 1);

  if (two_dimensional_mode) {

    for (int islow = 0; islow < n_slow; ++islow) {
      for (int ifast = 0; ifast < n_fast; ++ifast) {
          
        const size_t index = n_fast_pad * islow + ifast;

        const RealT pressure_po = (ifast < n_fast - 1 ? p[index + 1] : p[index - (n_fast - 1)]); 
        const RealT pressure_mo = (ifast >= 1 ? p[index - 1] : p[index + (n_fast - 1)]); 

        const RealT pressure_op = (islow < n_slow - 1 ? p[index + n_fast_pad] : p[index - (n_slow - 1) * n_fast_pad]); 
        const RealT pressure_om = (islow >= 1 ? p[index - n_fast_pad] : p[index + (n_slow - 1) * n_fast_pad]); 

        laplace_p[index] = (pressure_po - 2.0 * p[index] + pressure_mo) / (dx * dx);
        laplace_p[index] += (pressure_op - 2.0 * p[index] + pressure_om) / (dz * dz);

      }
    }
    
  } else {

    LOG_ERROR << "3D wave propagation not implemented";
    std::abort();

  }

}

void AdvanceWavePressure_0(int n_fast, int n_fast_padding, 
                           int n_medium, int n_slow,
                           int stencil_radius,
                           RealT dt, 
                           const RealT* velocity,
                           const RealT* laplace_p, 
                           const RealT* p0, 
                           RealT* p1) {
  
  assert(stencil_radius == 1);

  const int n_fast_pad = n_fast + n_fast_padding;

  const bool two_dimensional_mode = (n_medium == 1);

  if (two_dimensional_mode) {

    for (int islow = 0; islow < n_slow; ++islow) {
      for (int ifast = 0; ifast < n_fast; ++ifast) {

        const size_t index = n_fast_pad * islow + ifast;
        const RealT s = velocity[index] * velocity[index] * dt * dt;
                
        p1[index] = 2.0 * p0[index] - p1[index] + s * laplace_p[index];

      }
    }
    
  } else {

    LOG_ERROR << "3D wave propagation not implemented";
    std::abort();

  }
}
