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

static void WriteVTKXmlVariableHeader(const std::string& variable_name,
				      DataFormat format, int nb_components,
				      std::ostream* os_ptr) {
  
  const std::string data_type_str = (sizeof(RealT) == 4 ? "Float32": "Float64");
  const std::string format_str = (format == ASCII ? "ascii" : "binary");
  
  std::stringstream stringstream;
  stringstream << nb_components;
  const std::string nb_components_str = stringstream.str();

  *os_ptr << "type=\"" << data_type_str << "\" "
	  << "Name=\"" << variable_name << "\" "
	  << "NumberOfComponents=\"" << nb_components_str << "\" "
	  << "format=\"ascii\"";

}

void WriteVTKXmlVariable(const std::string& variable_name,
			 DataFormat format, int nb_components,
			 int n_fast, int n_fast_padding, 
			 int n_medium, int n_slow, 
			 RealT* data,
			 std::ostream* os_ptr) {

  *os_ptr << "<DataArray ";
  WriteVTKXmlVariableHeader(variable_name, format, nb_components, os_ptr);
  *os_ptr << ">\n";

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

      LOG_ERROR << "VTK binary output not implemented yet.";
      std::abort();

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

  *os_ptr << "</DataArray>\n";

}
