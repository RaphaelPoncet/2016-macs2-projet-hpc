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

#include "special_functions.hpp"

#include <cmath>

// The Ricker function is the second derivative of a Gaussian.
void Ricker(RealT amplitude, RealT frequency, RealT x0, int n, const RealT* x, RealT* out) {

  const RealT c = M_PI * M_PI * frequency * frequency;

  for (int i = 0; i < n; ++i)
    out[i] = amplitude * (1.0 - 2.0 * c * (x[i] - x0) * (x[i] - x0)) * exp(- c * (x[i] - x0) * (x[i] - x0));
  
}
