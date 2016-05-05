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

// Dummy value for padded areas of storage (for debug).
const RealT PADDING_FILL = - 3.1415926;
// Dummy value for areas of storage (for detecting non-properly initialized data).
const RealT INIT_FILL = 3.1415926;


#include <cassert>
#include <memory>

#include "multidimensional_storage.hpp"

#include "common.hpp"

static void ValidateInput(int n1, int n2, int n3, int padding) {

  assert(0 < n1);
  assert(0 < n2);
  assert(0 < n3);
  assert(0 <= padding);

}

MultiDimensionalStorage3D::MultiDimensionalStorage3D(int n1, int n2, int n3):
  m_n_fast(n1), m_n2(n2), m_n_slow(n3), m_padding(0), m_data(NULL) {
  
  ValidateInput(m_n_fast, m_n2, m_n_slow, m_padding);
  
}

MultiDimensionalStorage3D::MultiDimensionalStorage3D(int n1, int n2, int n3, int padding):
  m_n_fast(n1), m_n2(n2), m_n_slow(n3), m_padding(padding), m_data(NULL) {
  
  ValidateInput(m_n_fast, m_n2, m_n_slow, m_padding);
  
}

void MultiDimensionalStorage3D::Allocate() {

  LOG_DEBUG << "Allocating 3D storage with "
	    << "[n_fast(+padding), n2, n_slow] = "
	    << "[" << m_n_fast << "(+" << m_padding << "), "
	    << m_n2 << ", "
	    << m_n_slow << "]...";

  const size_t nb_elements = (m_n_fast + m_padding) * m_n2 * m_n_slow;
  const size_t size_in_bytes = nb_elements * sizeof(RealT);

  m_data = (RealT*)malloc(size_in_bytes);

  LOG_DEBUG << "Allocation done, memory used = "
	    << static_cast<double>(size_in_bytes / (1024.0 * 1024.0)) 
	    << " MBytes.";

  LOG_DEBUG << "Setting storage...";

  for (size_t i = 0; i < nb_elements; ++i)
    m_data[i] = INIT_FILL;

  // Now, Fill only the padded parts of the data. They should be left
  // untouched by the code.
  for (int islow = 0; islow < m_n_slow; ++islow) {
    for (int i2 = 0; i2 < m_n2; ++i2) {
      for (int ifast = m_n_fast; ifast < m_n_fast + m_padding; ++ifast) {

	const int n_fast_pad = m_n_fast + m_padding;
	const int n_medium = m_n2;

	const size_t index = n_medium * n_fast_pad * islow + n_fast_pad * i2 + ifast;
	m_data[index] = PADDING_FILL;

      }
    }
  }
  
  LOG_DEBUG << "Setting storage done.";
  
}

void MultiDimensionalStorage3D::DeAllocate() {

  LOG_DEBUG << "Deallocating 3D storage...";

  if (m_data != NULL)
    free(m_data);

  LOG_DEBUG << "Deallocation done.";
 
}

void MultiDimensionalStorage3D::Validate() {

  LOG_DEBUG << "Validating 3D storage...";

  size_t error_count = 0;

  // Padded area should be left untouched by the code.
  for (int islow = 0; islow < m_n_slow; ++islow) {
    for (int i2 = 0; i2 < m_n2; ++i2) {
      for (int ifast = m_n_fast; ifast < m_n_fast + m_padding; ++ifast) {
	
	const int n_fast_pad = m_n_fast + m_padding;
	const int n_medium = m_n2;

	const size_t index = n_medium * n_fast_pad * islow + n_fast_pad + ifast;
	if (m_data[index] != PADDING_FILL)
	  error_count += 1;

      }
    }
  }

  if (error_count != 0) {

    LOG_ERROR << "in MultiDimensional3D instance, "
	      << error_count << " elements in padded area have been touched.";

    std::abort();

  }

  LOG_DEBUG << "Validating 3D storage done.";

}
