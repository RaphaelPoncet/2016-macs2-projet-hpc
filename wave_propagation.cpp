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

    const size_t n_fast_medium = n_fast_pad * n_medium;

    // for (int islow = 1; islow < n_slow - 1; ++islow) {
    //   for (int imedium = 1; imedium < n_medium - 1; ++imedium) {
    //     for (int ifast = 1; ifast < n_fast - 1; ++ifast) {
          
    //       const size_t index = n_fast_medium * islow + n_fast_pad * imedium + ifast;
          
    //       const RealT pressure_poo = p[index + 1];
    //       const RealT pressure_moo = p[index - 1];

    //       const RealT pressure_opo = p[index + n_fast_pad];
    //       const RealT pressure_omo = p[index - n_fast_pad];

    //       const RealT pressure_oop = p[index + n_fast_medium];
    //       const RealT pressure_oom = p[index - n_fast_medium];

    //       UNUSED(pressure_poo);
    //       UNUSED(pressure_moo);
    //       UNUSED(pressure_opo);
    //       UNUSED(pressure_omo);
    //       UNUSED(pressure_oop);
    //       UNUSED(pressure_oom);
    //       laplace_p[index] = (pressure_poo - 2.0 * p[index] + pressure_moo) / (dx * dx);
    //       laplace_p[index] += (pressure_opo - 2.0 * p[index] + pressure_omo) / (dz * dz);
    //       laplace_p[index] += (pressure_oop - 2.0 * p[index] + pressure_oom) / (dy * dy);
    //       // laplace_p[index] = 0.0;
    //     }
    //   }
    // }

    for (int islow = 0; islow < n_slow; ++islow) {
      for (int imedium = 0; imedium < n_medium; ++imedium) {
        for (int ifast = 0; ifast < n_fast; ++ifast) {
          
          const size_t index = n_medium * n_fast_pad * islow + n_fast_pad * imedium + ifast;

          const RealT pressure_poo = (ifast < n_fast - 1 ? p[index + 1] : p[index - (n_fast - 1)]); 
          const RealT pressure_moo = (ifast >= 1 ? p[index - 1] : p[index + (n_fast - 1)]); 

          const RealT pressure_opo = (islow < n_slow - 1 ? p[index + n_fast_pad] : p[index - (n_slow - 1) * n_fast_pad]); 
          const RealT pressure_omo = (islow >= 1 ? p[index - n_fast_pad] : p[index + (n_slow - 1) * n_fast_pad]); 

          const RealT pressure_oop = (islow < n_slow - 1 ? p[index + n_fast_medium] : p[index - (n_slow - 1) * n_fast_medium]); 
          const RealT pressure_oom = (islow >= 1 ? p[index - n_fast_medium] : p[index + (n_slow - 1) * n_fast_medium]); 

          laplace_p[index] = (pressure_poo - 2.0 * p[index] + pressure_moo) / (dx * dx);
          laplace_p[index] += (pressure_opo - 2.0 * p[index] + pressure_omo) / (dz * dz);
          laplace_p[index] += (pressure_oop - 2.0 * p[index] + pressure_oom) / (dy * dy);

        }
      }
    }

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

    for (int islow = 0; islow < n_slow; ++islow) {
      for (int imedium = 0; imedium < n_medium; ++imedium) {
        for (int ifast = 0; ifast < n_fast; ++ifast) {
          
          const size_t index = n_medium * n_fast_pad * islow + n_fast_pad * imedium + ifast;
          const RealT s = velocity[index] * velocity[index] * dt * dt;
          
          p1[index] = 2.0 * p0[index] - p1[index] + s * laplace_p[index];

        }
      }
    }

  }
}
