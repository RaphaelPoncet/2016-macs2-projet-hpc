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

#include <sstream>
#include <string>

#include "common.hpp"
#include "array_vtk_io.hpp"

void WriteVTKXmlVariableHeader(const std::string& variable_name,
			       VTKDataFormat format, int nb_components,
			       size_t data_size_in_bytes,
			       size_t* data_offset_in_bytes_ptr,
			       std::ostream* os_ptr) {
  
  size_t data_offset_in_bytes = *data_offset_in_bytes_ptr;

  const std::string data_type_str = (sizeof(RealT) == 4 ? "Float32": "Float64");
  // VTK XML has two binary formats: 'binary' aka base64, and
  // 'appended'. We choose the second one.
  const std::string format_str = (format == ASCII ? "ascii" : "appended");
  
  std::stringstream stringstream;
  stringstream << nb_components;
  const std::string nb_components_str = stringstream.str();

  *os_ptr << "type=\"" << data_type_str << "\" "
	  << "Name=\"" << variable_name << "\" "
	  << "NumberOfComponents=\"" << nb_components_str << "\" "
	  << "format=\"" << format_str << "\"";

  if (format == BINARY) {

    *os_ptr << " offset=\"" << data_offset_in_bytes << "\"";

    // 4 is for an extra integer at the beginning of each array,
    // holding the length of it.
    data_offset_in_bytes += data_size_in_bytes + 4;
    *data_offset_in_bytes_ptr = data_offset_in_bytes;

  }
}

void WriteVTKXmlVariable(VTKDataFormat format, int nb_components,
			 int n_fast, int n_fast_padding, 
			 int n_medium, int n_slow, 
			 const RealT* data,
			 size_t* data_offset_in_bytes_ptr,
			 std::ostream* os_ptr) {

  if (nb_components == 1) {

    if (format == ASCII) {

      for (int islow = 0; islow < n_slow; ++islow) {
        for (int imedium = 0; imedium < n_medium; ++imedium) {
          for (int ifast = 0; ifast < n_fast; ++ifast) {

            const size_t index =
              n_medium * (n_fast + n_fast_padding) * islow + (n_fast + n_fast_padding) * imedium + ifast;
            *os_ptr << data[index] << "\n";

          }
        }
      }

    } else if (format == BINARY) {
      
      const size_t data_size_in_bytes = n_slow * n_medium * n_fast * sizeof(RealT);

      // VTK needs a 4 byte data size number written.
      assert(sizeof(int) == 4);
      const int data_size = (int)data_size_in_bytes;

      os_ptr->write(reinterpret_cast<const char*>(&data_size), 4);

      for (int islow = 0; islow < n_slow; ++islow) {
      	for (int imedium = 0; imedium < n_medium; ++imedium) {

	    const size_t index =
	      n_medium * (n_fast + n_fast_padding) * islow + (n_fast + n_fast_padding) * imedium;

      os_ptr->write(reinterpret_cast<const char*>(&(data[index])), n_fast * sizeof(RealT));

      	}
      }

    } else {

      assert(0);

    }
    
  } else if ((nb_components == 2) || (nb_components == 3)) {

    LOG_ERROR << "Multi component (e.g. vector) VTK output not implemented yet.";

    std::abort();

  } else {

    LOG_ERROR << "Invalid number of components: " << nb_components
	      << " (should be 1, 2 or 3).";
    std::abort();
    
  }

}
