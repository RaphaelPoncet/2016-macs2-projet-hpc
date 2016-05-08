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

// Read a simple binary input file stream looking for a header of the form:
//
// # variable_name datatype nb_components nx ny nz
//
void ReadBinaryVariable(std::istream& input_stream, 
			int n_fast, int n_fast_padding, int n_medium, 
			int n_slow, int nb_components, RealT* data) {

  assert(0 < n_fast);
  assert(0 < n_medium);
  assert(0 < n_slow);
  assert(0 <= n_fast_padding);
  assert(data != NULL);

  std::string variable_name;
  DataType variable_datatype;
  int variable_n_fast;
  int variable_n_medium;
  int variable_n_slow;
  int variable_nb_components;

  // Detect a header as a line starting with the following string.
  const std::string header_detect = "#";

  if (!input_stream) {

    LOG_ERROR << "Invalid input stream";
    std::abort();

  }

  LOG_DEBUG << "Parsing binary variable...";

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

      // nb_components.
      std::getline(isstream, token, ' ');
      if (!ParseToken<int>(token, &parsed_nb_components)) {
        
        LOG_ERROR << "Bad input for variable number of components (expecting an integer): "
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

      // Input validation.

      // Name.
      variable_name = parsed_name;

      // Datatype.
      if (parsed_datatype == std::string("Float32")) {

        variable_datatype = FLOAT32;

      } else if (parsed_datatype == std::string("Float64")) {

        variable_datatype = FLOAT64;

      } else {

        LOG_ERROR << "Invalid value for datatype: expecting 'Float32' or 'Float64', got "
                  << "\'" << parsed_datatype << "\'";

        std::abort();

      }

      ValidatePositiveInteger(parsed_n_fast, &variable_n_fast);
      ValidatePositiveInteger(parsed_n_medium, &variable_n_medium);
      ValidatePositiveInteger(parsed_n_slow, &variable_n_slow);
      ValidatePositiveInteger(parsed_nb_components, &variable_nb_components);

      if ((variable_nb_components != 1) && 
          (variable_nb_components != 2) && 
          (variable_nb_components != 3)) {

        LOG_ERROR << "Invalid value for nb_components: expecting '1', '2' or '3'"
                  << ", got \'" << variable_nb_components << "\'"; 

        std::abort();

      }

      LOG_VERBOSE << "Parsed the following informations from header:";

      LOG_VERBOSE << "---{name:" << "\'" << variable_name << "\',"
		  << "datatype:" << "\'" << parsed_datatype << "\',"
		  << "nb_components:" << "\'" << variable_nb_components << "\',"
		  << "n_fast:" << "\'" << variable_n_fast << "\',"
		  << "n_medium:" << "\'" << variable_n_medium << "\',"
		  << "n_slow:" << "\'" << variable_n_slow << "\'"
		  << "}";

      // Validate that the parsed variable header is compatible with
      // the parameters passed to this function.
      if (sizeof(RealT) != DATATYPE_SIZE[variable_datatype]) {

        LOG_ERROR << "Incompatible parsed datatype size: "
                  << DATATYPE_SIZE[variable_datatype] 
                  << " vs "
                  << sizeof(RealT) << ".";
        
        std::abort();

      }

      if (nb_components != variable_nb_components) {

        LOG_ERROR << "Incompatible parsed number of components: "
                  << variable_nb_components
                  << " vs "
                  << nb_components << ".";
        
        std::abort();

      }

      if ((n_fast != variable_n_fast) ||
          (n_medium != variable_n_medium) ||
          (n_slow != variable_n_slow)) {

        LOG_ERROR << "Incompatible parsed grid size: ("
                  << variable_n_fast << ", " 
                  << variable_n_medium << ", " 
                  << variable_n_slow << ") vs (" 
                  << n_fast << ", " 
                  << n_medium << ", " 
                  << n_slow << ").";

        std::abort();

      }

      // We only support scalar variables for now.
      assert(nb_components == 1);

      int index = 0;

      for (int islow = 0; islow < n_slow; ++islow) {
        for (int imedium = 0; imedium < n_medium; ++imedium) {

          input_stream.read(reinterpret_cast<char*>(&data[index]), 8 * n_fast);

          // we must account for padding in the fast direction.
          index += n_fast + n_fast_padding;
          
        }
      }
      
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

  LOG_DEBUG << "Parsing done.";

}
