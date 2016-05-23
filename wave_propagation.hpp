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

#ifndef WAVE_PROPAGATION_HPP
#define WAVE_PROPAGATION_HPP

#include <iosfwd>

#include "common.hpp"

// Baseline non optimized kernel.
void ApplySpongeLayer_0(int n_fast, int n_fast_padding, 
                        int n_medium, int n_slow,
                        const RealT* sponge_fast,
                        const RealT* sponge_medium,
                        const RealT* sponge_slow,
                        RealT* data);

// Baseline non optimized kernel.
void ComputeLaplacian_0(int n_fast, int n_fast_padding, 
                        int n_medium, int n_slow,
                        int stencil_radius,
                        RealT dx, RealT dy, RealT dz, 
                        const RealT* p, RealT* laplace_p);

// Baseline non optimized kernel.
void AdvanceWavePressure_0(int n_fast, int n_fast_padding, 
                           int n_medium, int n_slow,
                           int stencil_radius,
                           RealT dt, 
                           const RealT* velocity,
                           const RealT* laplace_p, 
                           const RealT* p0, 
                           RealT* p1);

#endif // WAVE_PROPAGATION_HPP
