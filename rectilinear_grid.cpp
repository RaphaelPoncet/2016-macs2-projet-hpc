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

#include "common.hpp"
#include "rectilinear_grid.hpp"

static void ValidateInput(int n1, int n2, int n3,
			  RealT x1_min, RealT x1_max,
			  RealT x2_min, RealT x2_max,
			  RealT x3_min, RealT x3_max) { 

  assert(0 < n1);
  assert(0 < n2);
  assert(0 < n3);

  assert(x1_min < x1_max);
  assert(x2_min < x2_max);
  assert(x3_min < x3_max);

}

RectilinearGrid3D::RectilinearGrid3D(RealT x_fast_min, RealT x_fast_max, int n_fast,
				     RealT x_medium_min, RealT x_medium_max, int n_medium,
				     RealT x_slow_min, RealT x_slow_max, int n_slow):
  m_n_fast(n_fast), m_n_medium(n_medium), m_n_slow(n_slow), 
  m_x_fast_min(x_fast_min), m_x_fast_max(x_fast_max), 
  m_x_medium_min(x_medium_min), m_x_medium_max(x_medium_max), 
  m_x_slow_min(x_slow_min), m_x_slow_max(x_slow_max) {

  ValidateInput(m_n_fast, m_n_medium, m_n_slow, 
		m_x_fast_min, m_x_fast_max,
		m_x_medium_min, m_x_medium_max,
		m_x_slow_min, m_x_slow_max);
		
}
