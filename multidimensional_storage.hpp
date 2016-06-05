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

#ifndef MULTIDIMENSIONAL_STORAGE_HPP
#define MULTIDIMENSIONAL_STORAGE_HPP

#include <array>
#include <cmath>

#include "common.hpp"
#include "variable_definitions.hpp"

// Dummy value for padded areas of storage (for debug).
const RealT PADDING_FILL = - 3.1415926;
// Dummy value for areas of storage (for detecting non-properly initialized data).
const RealT INIT_FILL = 0.0;

template <int nb_dimensions>
static void ValidateStorageDimensions(const std::array<size_t, nb_dimensions>& dimensions, int padding_fast) {
  
  for (auto dim: dimensions)
    assert(0 < dim);

  assert(0 <= padding_fast);
  
}

template <typename T, int nb_dimensions>
class MultiDimensionalStorage {
public:

  MultiDimensionalStorage():
    m_padding_fast(0) { 
    
    for (size_t& dim: m_dimensions)
      dim = 0;
  
  };

  MultiDimensionalStorage(std::array<size_t, nb_dimensions> dimensions):
    m_padding_fast(0), m_dimensions(dimensions) {

    ValidateStorageDimensions<nb_dimensions>(m_dimensions, m_padding_fast);
  
  };

  MultiDimensionalStorage(std::array<size_t, nb_dimensions> dimensions, size_t padding_fast):
    m_padding_fast(padding_fast), m_dimensions(dimensions) {

    ValidateStorageDimensions<nb_dimensions>(m_dimensions, m_padding_fast);
    
  };

  size_t dim_fast() { return m_dimensions.at(0); }
  size_t dim_fast() const { return m_dimensions.at(0); }

  size_t dim_slow() { return m_dimensions.at(nb_dimensions - 1); }
  size_t dim_slow() const { return m_dimensions.at(nb_dimensions - 1); }

  size_t padding_fast() { return m_padding_fast; }
  size_t padding_fast() const { return m_padding_fast; }

  std::array<size_t, nb_dimensions>& dimensions() {return m_dimensions;}
  const std::array<size_t, nb_dimensions>& dimensions() const {return m_dimensions;}

  void Allocate() {

    static_assert(0 < nb_dimensions, "Number of dimensions should be positive.");

    std::stringstream sstr;
    sstr << "Allocated multidimensional storage with (";

    size_t nb_elements = m_dimensions.at(0) + m_padding_fast;
    sstr << m_dimensions.at(0) << "+" << m_padding_fast;
  
    for (size_t i = 1; i < nb_dimensions; ++i) {
      
      nb_elements *= m_dimensions.at(i);
      sstr << "," << m_dimensions.at(i);

    }
        
    m_data.assign(nb_elements, INIT_FILL);

    sstr << ") elements, memory used: ";

    const double memory_used_in_mb = 
      static_cast<double>(nb_elements * sizeof(T)) / (1024.0 * 1024.0); 

    sstr << memory_used_in_mb << " MBytes.";

    LOG_DEBUG << sstr.str();

    // Now, fill the padding part of the data with specific values, do
    // detect modification later on.
    size_t n_slow = m_dimensions.at(nb_dimensions - 1);

    size_t n_medium = 1;

    for (size_t i = 1; i < nb_dimensions - 1; ++i)
      n_medium *= m_dimensions.at(i);
    
    const size_t n_fast = m_dimensions.at(0);
    const size_t n_fast_pad = n_fast + m_padding_fast;
    
    for (size_t islow = 0; islow < n_slow; ++islow) {

      T* data = RawDataSlowDimension(islow);

      for (size_t imedium = 0; imedium < n_medium; ++imedium) {
      
        const size_t index_base = imedium * n_fast_pad;
      
        for (size_t ifast = n_fast; ifast < n_fast_pad; ++ifast) {
        
          data[index_base + ifast] = PADDING_FILL;
        
        }
      }
    }
  }

  void DeAllocate() {

    // Do nothing, since memory is handled through a std::vector
    // container.

  }

  void ComputeRanges() const {

    // Compute min/max of storage values along the slowest dimension.
    const size_t n_slow = m_dimensions.at(nb_dimensions - 1);
    
    const T max_T = std::numeric_limits<T>::max();
    const T min_T = std::numeric_limits<T>::min();

    std::vector<T> minimums(n_slow, max_T);
    std::vector<T> maximums(n_slow, min_T);

    size_t n_medium = 1;

    for (size_t i = 1; i < nb_dimensions - 1; ++i)
      n_medium *= m_dimensions.at(i);

    const size_t n_fast = m_dimensions.at(0);
    const size_t n_fast_pad = n_fast + m_padding_fast;

    for (size_t islow = 0; islow < n_slow; ++islow) {

      const T* data = RawDataSlowDimension(islow);

      for (size_t imedium = 0; imedium < n_medium; ++imedium) {
      
        const size_t index_base = imedium * n_fast_pad;
      
        for (size_t ifast = 0; ifast < n_fast; ++ifast) {
        
          minimums.at(islow) = std::min(minimums.at(islow), data[index_base + ifast]);
          maximums.at(islow) = std::max(maximums.at(islow), data[index_base + ifast]);

        }
      }
    }
    
    std::stringstream msg;
    msg << "Variables range:\n\n";

    for (int islow = 0; islow < n_slow; ++islow) {

      msg << "____ " << variable::VARIABLE_NAMES[islow] << " in " 
          << "[" << minimums.at(islow) << ", " << maximums.at(islow) << "]\n";
    
    }

    LOG_DEBUG << msg.str();

  }

  void ComputeRanges() {

    return static_cast<const MultiDimensionalStorage* >(this)->ComputeRanges();

  }

  void Validate() {

    return static_cast<const MultiDimensionalStorage* >(this)->Validate();

  }

  void Validate() const {

    LOG_INFO << "Validating storage...";

    this->ComputeRanges();
  
    // Padded area should be left untouched by any client code.
    size_t pad_error_count = 0;

    const size_t n_slow = m_dimensions.at(nb_dimensions - 1);

    size_t n_medium = 1;

    for (size_t i = 1; i < nb_dimensions - 1; ++i)
      n_medium *= m_dimensions.at(i);

    const size_t n_fast = m_dimensions.at(0);
    const size_t n_fast_pad = n_fast + m_padding_fast;

    for (size_t islow = 0; islow < n_slow; ++islow) {

      const T* data = RawDataSlowDimension(islow);

      for (size_t imedium = 0; imedium < n_medium; ++imedium) {
      
        const size_t index_base = imedium * n_fast_pad;
      
        for (size_t ifast = n_fast; ifast < n_fast_pad; ++ifast) {
        
          if (data[index_base + ifast] != PADDING_FILL)
            pad_error_count += 1;
         
        }
      }
    }

    if (pad_error_count != 0) {

      LOG_ERROR << "in MultiDimensionalStorage instance, "
                << pad_error_count << " elements in padded area have been touched.";

      std::abort();

    }

    // Check for NaNs.
    size_t nan_error_count = 0;

    for (size_t islow = 0; islow < n_slow; ++islow) {

      const T* data = RawDataSlowDimension(islow);

      for (size_t imedium = 0; imedium < n_medium; ++imedium) {
      
        const size_t index_base = imedium * n_fast_pad;
      
        for (size_t ifast = 0; ifast < n_fast; ++ifast) {

          if (std::isnan(data[index_base + ifast]))
            nan_error_count += 1;
         
        }
      }
    }

    if (nan_error_count != 0) {

      LOG_ERROR << "in MultiDimensionalStorage instance, "
                << nan_error_count << " elements have NaN value.";

      std::abort();

    }

    LOG_INFO << "Validating 4D storage done."
             << "\n";

  }

  T* RawDataSlowDimension(size_t islow) {

    return const_cast<T*>(static_cast<const MultiDimensionalStorage*>(this)->RawDataSlowDimension(islow));
    
  }

  const T* RawDataSlowDimension(size_t islow) const {

    size_t offset = (m_dimensions.at(0) + m_padding_fast) * islow;
      
    for (size_t i = 1; i < nb_dimensions - 1; ++i)
      offset *= m_dimensions.at(i);
    
    if (offset >= m_data.size()) {

      LOG_ERROR << "Invalid offset (" << offset 
                << ") in accessing data with "
                << m_data.size() << " elements";

      std::abort();

    }

    return &(m_data.at(offset));

  }

private:

  size_t m_padding_fast;
  std::array<size_t, nb_dimensions> m_dimensions;
  std::vector<T> m_data;

};


typedef MultiDimensionalStorage<RealT, 4> MultiDimensionalStorage4D;

#endif // MULTIDIMENSIONAL_STORAGE_HPP
