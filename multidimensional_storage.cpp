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

#include <cassert>
#include <memory>

// Dummy value for padded areas of storage (for debug).
const RealT PADDING_FILL = - 3.1415926;
// Dummy value for areas of storage (for detecting non-properly initialized data).
const RealT INIT_FILL = 3.1415926;

#include "common.hpp"
#include "multidimensional_storage.hpp"

static void ValidateInput(int n1, int n2, int n3, int n4, int padding) {

  assert(0 < n1);
  assert(0 < n2);
  assert(0 < n3);
  assert(0 < n4);
  assert(0 <= padding);

}

MultiDimensionalStorage4D::MultiDimensionalStorage4D(int n1, int n2, int n3, int n4):
  m_n_fast(n1), m_n2(n2), m_n3(n3), m_n_slow(n4), m_padding(0), m_data(NULL) {
  
  ValidateInput(m_n_fast, m_n2, m_n3, m_n_slow, m_padding);
  
}

  MultiDimensionalStorage4D::MultiDimensionalStorage4D(int n1, int n2, int n3, int n4, int padding):
    m_n_fast(n1), m_n2(n2), m_n3(n3), m_n_slow(n4), m_padding(padding), m_data(NULL) {
  
    ValidateInput(m_n_fast, m_n2, m_n3, m_n_slow, m_padding);
    
}

void MultiDimensionalStorage4D::Allocate() {

  LOG_DEBUG << "Allocating 4D storage with "
	    << "[n_fast(+padding), n2, n3, n_slow] = "
	    << "[" << m_n_fast << "(+" << m_padding << "), "
    	    << m_n2 << ", "
	    << m_n3 << ", "
	    << m_n_slow << "]...";

  const size_t nb_elements = (m_n_fast + m_padding) * m_n2 * m_n3 * m_n_slow;
  const size_t size_in_bytes = nb_elements * sizeof(RealT);

  m_data = (RealT*)malloc(size_in_bytes);

  LOG_DEBUG << "Allocation done, memory used = "
	    << static_cast<double>(size_in_bytes / (1024.0 * 1024.0)) 
	    << " MBytes."
	    << "\n";

  LOG_DEBUG << "Setting storage...";

  for (size_t i = 0; i < nb_elements; ++i)
    m_data[i] = INIT_FILL;

  // Now, Fill only the padded parts of the data. They should be left
  // untouched by the code.
  const int n_fast = m_n_fast;
  const int n_fast_pad = n_fast + m_padding;
  const int n2 = m_n2;
  const int n3 = m_n3;
  const int n_slow = m_n_slow;
  
  for (int islow = 0; islow < n_slow; ++islow) {
    for (int i3 = 0; i3 < n3; ++i3) {
      for (int i2 = 0; i2 < n2; ++i2) {
	for (int ifast = n_fast; ifast < n_fast_pad; ++ifast) {

	  const size_t index = 
	    n3 * n2 * n_fast_pad * islow + n2 * n_fast_pad * i3 + n_fast_pad * i2 + ifast;

	  m_data[index] = PADDING_FILL;

	}
      }
    }
  }
  
  LOG_DEBUG << "Setting storage done."
	    << "\n";
  
}

void MultiDimensionalStorage4D::DeAllocate() {

  LOG_DEBUG << "Deallocating 4D storage...";

  if (m_data != NULL)
    free(m_data);

  LOG_DEBUG << "Deallocation done."
	    << "\n";
 
}

void MultiDimensionalStorage4D::Validate() {

  LOG_DEBUG << "Validating 4D storage...";

  size_t error_count = 0;
  
  const int n_fast = m_n_fast;
  const int n_fast_pad = n_fast + m_padding;
  const int n2 = m_n2;
  const int n3 = m_n3;
  const int n_slow = m_n_slow;

  // Padded area should be left untouched by the code.
  for (int islow = 0; islow < n_slow; ++islow) {
    for (int i3 = 0; i3 < n3; ++i3) {
      for (int i2 = 0; i2 < n2; ++i2) {
	for (int ifast = n_fast; ifast < n_fast_pad; ++ifast) {

	const size_t index = 
	  n3 * n2 * n_fast_pad * islow + n2 * n_fast_pad * i3 + n_fast_pad * i2 + ifast;

	if (m_data[index] != PADDING_FILL)
	  error_count += 1;

	}
      }
    }
  }

  if (error_count != 0) {

    LOG_ERROR << "in MultiDimensional4D instance, "
	      << error_count << " elements in padded area have been touched.";

    std::abort();

  }

  LOG_DEBUG << "Validating 4D storage done."
	    << "\n";

}

RealT* MultiDimensionalStorage4D::RawDataSlowDimension(int i) {

  assert((0 <= i) && (i < m_n_slow));

  const size_t offset = m_n3 * m_n2 * (m_n_fast + m_padding) * i;

  RealT* result = &(m_data[offset]);
  assert(result != NULL);

  return result;

}
