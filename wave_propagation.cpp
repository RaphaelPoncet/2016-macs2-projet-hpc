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
