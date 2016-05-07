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

#include "array_binary_io.hpp"
#include "io_utils.hpp"

#include <fstream>
#include <iostream>

static void ValidatePositiveInteger(int parsed_integer, int* output_integer_ptr) {

  if (0 < parsed_integer) {
    
    *output_integer_ptr = parsed_integer;

  } else {
    
    LOG_ERROR << "Invalid value for n_fast: expecting a positive integer, got "
	      << "\'" << parsed_integer << "\'";
    std::abort();
    
  }
}

// Read an input file stream looking for a header of the form:
//
// # variable_name datatype nx ny nz nb_components
//
void ReadBinaryVariableHeader(std::istream& input_stream, 
			      std::string* name_ptr, DataType* datatype_ptr, 
			      int* n_fast_ptr, int* n_medium_ptr, int* n_slow_ptr, 
			      int* nb_components_ptr) {

  // Detect a header as a line starting with the following string.
  const std::string header_detect = "#";

  if (!input_stream) {

    LOG_ERROR << "Invalid input stream";
    std::abort();

  }

  LOG_INFO << "Parsing binary variable header...";

  // Will hold the current line.
  std::string line = "";

  int header_found_count = 0;

  while (std::getline(input_stream, line)) {

    if (line.substr(0, header_detect.size()) == header_detect) {

      LOG_VERBOSE << "Header line detected: \"" 
		  << line << "\""; 
      header_found_count += 1;

      std::string parsed_name = "";
      std::string parsed_datatype = "";
      int parsed_n_fast = -1;
      int parsed_n_medium = -1;
      int parsed_n_slow = -1;
      int parsed_nb_components = -1;
      
      std::istringstream isstream(line);
      std::string token;

      // Header magic string. Do nothing with it.
      std::getline(isstream, token, ' ');
      assert(token == header_detect);

      // Variable name.
      std::getline(isstream, token, ' ');
      if (!ParseToken<std::string>(token, &parsed_name)) {

	LOG_ERROR << "Bad input for variable name (expecting a string): "
		  << "\'" << token << "\'";
	std::abort();

      }

      // Variable datatype.
      std::getline(isstream, token, ' ');
      if (!ParseToken<std::string>(token, &parsed_datatype)) {

	LOG_ERROR << "Bad input for variable datatype (expecting a string): "
		  << "\'" << token << "\'";
	std::abort();

      }
      
      // n_fast.
      std::getline(isstream, token, ' ');
      if (!ParseToken<int>(token, &parsed_n_fast)) {

	LOG_ERROR << "Bad input for variable n_fast (expecting an integer): "
		  << "\'" << token << "\'";
	std::abort();

      }

      // n_medium.
      std::getline(isstream, token, ' ');
      if (!ParseToken<int>(token, &parsed_n_medium)) {

	LOG_ERROR << "Bad input for variable n_medium (expecting an integer): "
		  << "\'" << token << "\'";
	std::abort();

      }

      // n_slow.
      std::getline(isstream, token, ' ');
      if (!ParseToken<int>(token, &parsed_n_slow)) {

	LOG_ERROR << "Bad input for variable n_slow (expecting an integer): "
		  << "\'" << token << "\'";
	std::abort();

      }

      // nb_components.
      std::getline(isstream, token, ' ');
      if (!ParseToken<int>(token, &parsed_nb_components)) {

	LOG_ERROR << "Bad input for variable number of components (expecting an integer): "
		  << "\'" << token << "\'";
	std::abort();

      }

      // Input validation.

      // Name.
      *name_ptr = parsed_name;

      // Datatype.
      if (parsed_datatype == std::string("Float32")) {

	*datatype_ptr = FLOAT32;

      } else if (parsed_datatype == std::string("Float64")) {

	*datatype_ptr = FLOAT64;

      } else {

	LOG_ERROR << "Invalid value for datatype: expecting 'Float32' or 'Float64', got "
		  << "\'" << parsed_datatype << "\'";

	std::abort();

      }

      ValidatePositiveInteger(parsed_n_fast, n_fast_ptr);
      ValidatePositiveInteger(parsed_n_medium, n_medium_ptr);
      ValidatePositiveInteger(parsed_n_slow, n_slow_ptr);
      ValidatePositiveInteger(parsed_nb_components, nb_components_ptr);

      LOG_VERBOSE << "Parsed the following informations from header:";

      LOG_VERBOSE << "---{name:" << "\'" << *name_ptr << "\',"
		  << "datatype:" << "\'" << parsed_datatype << "\',"
		  << "n_fast:" << "\'" << *n_fast_ptr << "\',"
		  << "n_medium:" << "\'" << *n_medium_ptr << "\',"
		  << "n_slow:" << "\'" << *n_slow_ptr << "\',"
		  << "nb_components:" << "\'" << *nb_components_ptr << "\'"
		  << "}";

    }

  }

  if (header_found_count == 0) {

    LOG_ERROR << "Header line (e.g. a line starting with \"" 
	      << header_detect << "\") not found.";
    std::abort();

  } else if (header_found_count > 1) {

    LOG_ERROR << "Found " << header_found_count 
	      << " header lines, expected 1";

    std::abort();

  }

  LOG_INFO << "Parsing done.\n";

}

void ReadBinearyVariable() {

}
