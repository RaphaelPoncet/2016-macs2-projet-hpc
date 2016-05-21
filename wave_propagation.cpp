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

void InitSpongeArray(int sponge_width, int n, RealT* out) {

  assert(0 < sponge_width);
  assert(0 < n);
  assert(out != NULL);

  // We implement the method of Cerjan et al. (1985): we multiply the
  // wave amplitudes by a smooth factor.

  for (int i = 0; i < n; ++i)
    out[i] = 1.0;

  if (sponge_width < n - sponge_width) {

    for (int i = 0; i < sponge_width; ++i) {
     
      const RealT x = (RealT)((sponge_width - i) * 20 / sponge_width);
      out[i] = expf(- 0.015 * x * x);

    }

    for (int i = n - sponge_width; i < n; ++i) {
     
      const RealT x = (RealT)((n - sponge_width - i) * 20 / sponge_width);
      out[i] = expf(- 0.015 * x * x);

    }

  } else {

    // The array is too short to support a sponge. Do not attenuate
    // waves in this case.
    LOG_INFO << "Will not initalise sponge layer: "
             << "sponge_width ("
             << sponge_width << ") larger than half length (" << n << ")";

  }
}

void DumpSpongeArray(int n, RealT* sponge_array, std::ostream* os_ptr) {

  DumpArrayAscii(n, sponge_array, os_ptr);

}

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
        // const RealT scaling = 1.0;

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
  
  const int n_fast_pad = n_fast + n_fast_padding;

  const bool two_dimensional_mode = (n_medium == 1);

  if (two_dimensional_mode) {

    for (int islow = 0; islow < n_slow; ++islow) {
      if ((islow >= stencil_radius) && (islow < n_slow - stencil_radius)) {
        for (int ifast = 0; ifast < n_fast; ++ifast) {
          if ((ifast >= stencil_radius) && (ifast < n_fast - stencil_radius)) {

            const size_t index = n_fast_pad * islow + ifast;

            laplace_p[index] = (p[index + 1] - 2.0 * p[index] + p[index - 1]) / (dx * dx);
            laplace_p[index] += (p[index + n_fast_pad] - 2.0 * p[index] + p[index - n_fast_pad]) / (dz * dz);
              
          }
        }
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
  
  const int n_fast_pad = n_fast + n_fast_padding;

  const bool two_dimensional_mode = (n_medium == 1);

  if (two_dimensional_mode) {

    for (int islow = 0; islow < n_slow; ++islow) {
      if ((islow >= stencil_radius) && (islow < n_slow - stencil_radius)) {
        for (int ifast = 0; ifast < n_fast; ++ifast) {
          if ((ifast >= stencil_radius) && (ifast < n_fast - stencil_radius)) {

                const size_t index = n_fast_pad * islow + ifast;
                
                const RealT s = velocity[index] * velocity[index] * dt * dt;
                
                p1[index] = 2.0 * p0[index] - p1[index] + s * laplace_p[index];;
              
          }
        }
      }
    }

  } else {

    LOG_ERROR << "3D wave propagation not implemented";
    std::abort();

  }

}
