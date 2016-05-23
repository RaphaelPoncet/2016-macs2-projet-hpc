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

#include "parameter_parser.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <streambuf>
#include <utility>
#include <vector>

#include "picojson.h"

#include "array_binary_io.hpp"
#include "grid_output.hpp"
#include "location_output.hpp"
#include "mathematical_parser.hpp"
#include "multidimensional_storage.hpp"
#include "rectilinear_grid.hpp"
#include "timeloop_manager.hpp"
#include "variable_definitions.hpp"

const std::string TMP_BINARY_FILENAME = "./tmp.dat";
const std::string MAGIC_ITER_STRING = "\%i";

OutputEvent::OutputEvent():
  m_type(""), m_format(""), m_stream_name(""), m_rhythm(0), 
  m_output_stream_ptr(NULL), m_file_ptr(NULL), m_index_slow(-1) {};

void OutputEvent::ParseFromJSON(const picojson::value& v) {

  // check if the type of the value is "object"
  if (! v.is<picojson::object>()) {

    LOG_ERROR << "Malformed JSON string. Check yout parameter file.";
    std::abort();

  }

  picojson::object o = v.get<picojson::object>();

  if (o["type"].is<std::string>())
    m_type = o["type"].get<std::string>();
  
  if (o["format"].is<std::string>())
    m_format = o["format"].get<std::string>();

  if (o["file"].is<std::string>())
    m_stream_name = o["file"].get<std::string>();

  if (o["rhythm"].is<double>())
    m_rhythm = o["rhythm"].get<double>();

  if (o["iz"].is<double>())
      m_index_slow = o["iz"].get<double>();

  if (o["name"].is<std::string>())
    m_var_name = o["name"].get<std::string>();

  if (o["formula"].is<std::string>()) {

    std::stringstream formula;
    formula << m_var_name << " = " << o["formula"].get<std::string>();
    m_formula = formula.str();

  }

}

void OutputEvent::Validate() {

  // Type is in a white list.
  if (m_type == "OutputVariables") {

    if (m_format == "binary" || m_format == "VTK") {

    } else {

      LOG_ERROR << "Invalid value for output event type: expected "
                << "\'VTK\'"
                << "or \'binary\'"
                << ", got \'" << m_format << "\'";

      std::abort();

    }

  } else if (m_type == "OutputReceivers") {

    if (m_format == "binary" || m_format == "VTK") {

    } else {

      LOG_ERROR << "Invalid value for output event type: expected "
                << "\'VTK\'"
                << "or \'binary\'"
                << ", got \'" << m_format << "\'";

      std::abort();

    }

  } else if (m_type == "CheckVariables") {

  } else if (m_type == "EvalVariable") {

  } else {

    LOG_ERROR << "Invalid value for output event type: expected "
              << "\'OutputVariables\'"
              << "or \'OutputReceivers\'"
              << ", got \'" << m_type << "\'";

    std::abort();

  }
  
}

void OutputEvent::Create(int nb_iter,
                         const MultiDimensionalStorage4D& storage,
                         const RectilinearGrid3D& grid) {

  if (this->type() == "OutputVariables") {

  } else if (this->type() == "OutputReceivers") {

    m_file_ptr = new std::ofstream;
    
    // Whatever the format (VTK or gridded binary), we write the
    // receivers data in a binary file, and convert it at the end if
    // needed.

    const std::string filename = 
      (m_format != "binary" ? TMP_BINARY_FILENAME : m_stream_name);

    m_file_ptr->open(filename, std::ios::out | std::ios::binary);

    m_location_output_ptr = new LocationOutput(m_index_slow);
    assert(0 < m_rhythm);
    const int nb_outputs = std::ceil((RealT)nb_iter / (RealT)m_rhythm);
    m_location_output_ptr->WriteHeader(grid, nb_outputs, m_file_ptr);

  } else if (this->type() == "CheckVariables") {

  } else if (this->type() == "EvalVariable") {

  } else {

    LOG_ERROR << "Unrecognized event "
              << "\'" << this->type() << "\'";
    std::abort();

  }

  UNUSED(storage);

}

void OutputEvent::Destroy() {

  if (this->type() == "OutputVariables") {

  } else if (this->type() == "OutputReceivers") {

    m_file_ptr->close();
    delete m_file_ptr;
    delete m_location_output_ptr;

    // Convert the receiver file, if needed.
    
    assert(m_format == "binary" || m_format == "VTK");

    if (m_format == "VTK") {

      std::string name = "";
      DataType datatype = NONE;
      int nb_components = - 1;
      int n_fast = - 1;
      int n_medium = - 1;
      int n_slow = - 1;
      RealT x_fast_min = 0.0;
      RealT x_fast_max = 0.0;
      RealT x_medium_min = 0.0;
      RealT x_medium_max = 0.0;
      RealT x_slow_min = 0.0;
      RealT x_slow_max = 0.0;

      LOG_VERBOSE << "Opening temporary file "
                  << "\'" << TMP_BINARY_FILENAME << "\'"
                  << " for VTK conversion...";

      std::ifstream tmp_binary_file(TMP_BINARY_FILENAME, std::ifstream::in | std::ifstream::binary);

      if (!tmp_binary_file) {
    
        LOG_ERROR << "Invalid input stream";
        std::abort();

      }

      ReadBinaryVariableHeader(&tmp_binary_file,
                               &name, 
                               &datatype,
                               &nb_components,
                               &n_fast, 
                               &n_medium,
                               &n_slow,
                               &x_fast_min,
                               &x_fast_max,
                               &x_medium_min,
                               &x_medium_max,
                               &x_slow_min,
                               &x_slow_max);

      const RectilinearGrid3D receiver_grid = 
        RectilinearGrid3D(x_fast_min, x_fast_max, n_fast, 
                          x_medium_min, x_medium_max, n_medium, 
                          x_slow_min, x_slow_max, n_slow);

      const int n_fast_padding = 0;

      MultiDimensionalStorage4D variable_storage_receivers = 
        MultiDimensionalStorage4D(n_fast, n_medium, n_slow, variable::NB_VARIABLES, n_fast_padding);
      
      variable_storage_receivers.Allocate();

      RealT* data = variable_storage_receivers.RawDataSlowDimension(variable::PRESSURE_0);

      tmp_binary_file.clear();
      tmp_binary_file.seekg(0);

      ReadBinaryVariable(tmp_binary_file, n_fast, n_fast_padding, n_medium, n_slow, nb_components, data);

      tmp_binary_file.close();

      LOG_VERBOSE << "Reading temporary file done.";

      LOG_VERBOSE << "Destroying temporary file...";

      std::remove(TMP_BINARY_FILENAME.c_str());

      LOG_VERBOSE << "Destroying temporary file done.";
  
      variable_storage_receivers.Validate();
      
      const std::string output_filename = m_stream_name;

      LOG_VERBOSE << "Writing VTK receiver file \"" << output_filename << "\"...";
  
      std::ofstream receivers_vtk_file(output_filename.c_str(), std::ofstream::binary);
  
      OutputGridAndData("vtk", receiver_grid, variable_storage_receivers, &receivers_vtk_file);
  
      receivers_vtk_file.close();
  
      LOG_VERBOSE << "Writing VTK receiver file done.\n";
      
      variable_storage_receivers.DeAllocate();
      
    }

  } else if (this->type() == "CheckVariables") {

  } else if (this->type() == "EvalVariable") {

  } else {

    assert(0);

  }

}

void OutputEvent::Execute(int iter,
                          const RectilinearGrid3D& grid,
                          MathematicalParser* parser_ptr,
                          MultiDimensionalStorage4D* variable_storage_ptr) {

  if (this->type() == std::string("OutputVariables")) {

    m_file_ptr = new std::ofstream;

    std::string base_name = m_stream_name;

    // Look for special patterns in the filename and replace them by
    // the iteration number.
    std::string::size_type position = 0;

    std::stringstream iteration_stream;
    iteration_stream << std::setfill('0') << std::setw(5) << iter; 
    std::string iter_string = iteration_stream.str();

    while ((position = base_name.find(MAGIC_ITER_STRING, position)) != std::string::npos) {

      base_name.replace(position, MAGIC_ITER_STRING.size(), iter_string);
      ++position;
      
    }
    
    const std::string output_filename = base_name;
      
    LOG_INFO << "Writing output file \"" << output_filename << "\"...";

    m_file_ptr->open(output_filename.c_str(), std::ofstream::binary);
    OutputGridAndData(m_format, grid, *variable_storage_ptr, m_file_ptr);

    LOG_INFO << "Writing output file done.\n";

    m_file_ptr->flush();
    m_file_ptr->close();
    delete m_file_ptr;

  } else if (this->type() == std::string("OutputReceivers")) {

    if (m_location_output_ptr == NULL) {

      LOG_ERROR << "in output event "
                << "\'" << m_type << "\'"
                << " (stream name \'" << m_stream_name << "\')"
                << ", attempting to write a NULL file object. "
                << "Did you initialize the event with Create() ?";

      std::abort();

    }
    
    m_location_output_ptr->Write(grid, *variable_storage_ptr, m_file_ptr);

  } else if (this->type() == "CheckVariables") {

    variable_storage_ptr->Validate();

  } else if (this->type() == "EvalVariable") {

    parser_ptr->EvaluateExpression(this->formula(), grid, variable_storage_ptr);

  } else {

    assert(0);

  }

}

static void ParseFormulaFromJSON(picojson::object& o, 
                                 const std::string& var_name,
                                 const RectilinearGrid3D& grid,
                                 MathematicalParser* math_parser_ptr,
                                 MultiDimensionalStorage4D* storage_ptr) {

  const std::string formula = o["formula"].get<std::string>();
            
  LOG_VERBOSE << "Evaluating variable "
              << "\'" << var_name << "\'"
              << " with formula "
              << "\'" << formula << "\'";

  std::stringstream parser_formula;
  parser_formula << var_name << "=" << formula;
  
  math_parser_ptr->EvaluateExpression(parser_formula.str(), grid, storage_ptr);

}

static void ParseGridFromJSON(picojson::object& v_grid, 
                              RectilinearGrid3D* grid_ptr) {

  UNUSED(v_grid);
  UNUSED(grid_ptr);

  const bool has_file = v_grid["file"].is<std::string>();
  assert(has_file);

  if (has_file) {
    
    const std::string grid_filename = v_grid["file"].get<std::string>();
    
    std::ifstream tmp_binary_file(grid_filename, std::ifstream::in | std::ifstream::binary);

    if (!tmp_binary_file) {
      
      LOG_ERROR << "Invalid input stream";
      std::abort();
      
    }

    std::string dummy_name = "";
    DataType dummy_datatype = NONE;
    int dummy_nb_components = - 1;
    int n_fast = - 1;
    int n_medium = - 1;
    int n_slow = - 1;
    RealT x_fast_min = 0.0;
    RealT x_fast_max = 0.0;
    RealT x_medium_min = 0.0;
    RealT x_medium_max = 0.0;
    RealT x_slow_min = 0.0;
    RealT x_slow_max = 0.0;

    ReadBinaryVariableHeader(&tmp_binary_file,
                             &dummy_name, 
                             &dummy_datatype,
                             &dummy_nb_components,
                             &n_fast, 
                             &n_medium,
                             &n_slow,
                             &x_fast_min,
                             &x_fast_max,
                             &x_medium_min,
                             &x_medium_max,
                             &x_slow_min,
                             &x_slow_max);

    // Initialize grid from geometry informations. 
    
    // Depending on variable support (cell or node), grid extent
    // and size will change a little.
    //
    // It is merely a matter of convenience for visualization: by
    // default, Paraview outputs cell variables with no interpolation,
    // and node variables with interpolation.

    const VariableSupport variable_support = variable::VARIABLE_SUPPORT;

    // The number of nodes is the number of cells + 1, except if the
    // number of nodes is 1, when the number of cells is also 1.
    const int n_fast_grid = 
      (variable_support == CELL ? (n_fast == 1 ? n_fast : n_fast + 1) : n_fast);
    const int n_medium_grid = 
      (variable_support == CELL ? (n_medium == 1 ? n_medium : n_medium + 1) : n_medium);
    const int n_slow_grid = 
      (variable_support == CELL ? (n_slow == 1 ? n_slow : n_slow + 1) : n_slow);
    
    assert(0 < n_fast);
    const RealT dx_fast = (x_fast_max - x_fast_min) / (RealT)n_fast;

    assert(0 < n_medium);
    const RealT dx_medium = (x_medium_max - x_medium_min) / (RealT)n_medium;

    assert(0 < n_slow);
    const RealT dx_slow = (x_slow_max - x_slow_min) / (RealT)n_slow;

    const RealT x_fast_min_grid = 
      (variable_support == CELL ? x_fast_min - 0.5 * dx_fast: x_fast_min);
    const RealT x_fast_max_grid = 
      (variable_support == CELL ? x_fast_max + 0.5 * dx_fast: x_fast_max);
    const RealT x_medium_min_grid = 
      (variable_support == CELL ? x_medium_min - 0.5 * dx_medium: x_medium_min);
    const RealT x_medium_max_grid = 
      (variable_support == CELL ? x_medium_max + 0.5 * dx_medium: x_medium_max);
    const RealT x_slow_min_grid = 
      (variable_support == CELL ? x_slow_min - 0.5 * dx_slow: x_slow_min);
    const RealT x_slow_max_grid = 
      (variable_support == CELL ? x_slow_max + 0.5 * dx_slow: x_slow_max);

    *grid_ptr = RectilinearGrid3D(x_fast_min_grid, x_fast_max_grid, n_fast_grid, 
                                  x_medium_min_grid, x_medium_max_grid, n_medium_grid, 
                                  x_slow_min_grid, x_slow_max_grid, n_slow_grid);
    
  }

}

void ParseParameterFile(int n_fast_padding,
                        std::istream* input_stream_ptr,
                        RectilinearGrid3D* grid_ptr,
                        MultiDimensionalStorage4D* storage_ptr,
                        MathematicalParser* math_parser_ptr,
                        std::vector<OutputEvent>* output_events_ptr,
                        TimeloopManager* timeloop_manager_ptr) {

  std::string json_string;

  input_stream_ptr->seekg(0, std::ios::end);   

  json_string.reserve(input_stream_ptr->tellg());

  input_stream_ptr->seekg(0, std::ios::beg);
  
  json_string.assign((std::istreambuf_iterator<char>(*input_stream_ptr)),
                     std::istreambuf_iterator<char>());

  std::stringstream sstr(json_string);

  picojson::value v;
  sstr >> v;

  picojson::object o = v.get<picojson::object>();

  const bool has_grid = o["grid"].is<picojson::object>();

  if (has_grid) {

    picojson::object v_grid = o["grid"].get<picojson::object>();

    ParseGridFromJSON(v_grid, grid_ptr);
    
  }

  math_parser_ptr->SetGridConstants(*grid_ptr);

  const bool has_parameters = o["parameters"].is<picojson::array>();

  if (has_parameters) {

    picojson::array parameters = o["parameters"].get<picojson::array>();

    for (auto i = parameters.begin(); i != parameters.end(); ++i) {

      const std::string formula = i->to_str();

      std::stringstream sstream(formula);
      std::string constant_name = "";
      std::getline(sstream, constant_name, '=');
      std::string constant_formula = "";
      std::getline(sstream, constant_formula, '=');

      math_parser_ptr->SetExpression(constant_formula);
      const double constant_value = math_parser_ptr->SafeEval();
      math_parser_ptr->AddConstant(constant_name, constant_value);
      
    }    
  }

  *storage_ptr = MultiDimensionalStorage4D(grid_ptr->n_fast(),
                                           grid_ptr->n_medium(),
                                           grid_ptr->n_slow(),
                                           variable::NB_VARIABLES,
                                           n_fast_padding);

  storage_ptr->Allocate();

  const bool has_init = o["init"].is<picojson::object>();

  if (has_init) {
    
    std::map<std::string, int> variable_database;
    
    for (int i = 0; i < variable::NB_VARIABLES; ++i) {

      auto pair = std::make_pair(std::string(variable::VARIABLE_NAMES[i]), i);
      variable_database.insert(pair);
      
    }
    
    picojson::object variable_names_json = o["init"].get<picojson::object>();

    for (auto i = variable_names_json.begin(); i != variable_names_json.end(); ++i) {

      const std::string var_name = i->first;
      int index_var = - 1;

      const bool var_found_in_database = 
        (variable_database.find(var_name) != variable_database.end());

      if (var_found_in_database) {

        index_var = variable_database[var_name];
        RealT* data = storage_ptr->RawDataSlowDimension(index_var);

        LOG_VERBOSE << "Variable "
                    << "\'" << var_name << "\'"
                    << " found in database with index "
                    << index_var;

        if (!i->second.is<picojson::object>()) {

          LOG_ERROR << "Malformed variable object";
          std::abort();

        }

        picojson::object o = i->second.get<picojson::object>();

        const bool has_filename = 
          o["file"].is<std::string>();

        const bool has_formula = 
          o["formula"].is<std::string>();

        if (!(has_formula || has_filename)) {

          LOG_ERROR << "Invalid initialization for variable "
                    << "\'" << var_name << "\'"
                    << ": expected "
                    << "\'formula\' or \'file\' field";

          std::abort();

        }

        if (has_formula) {

          ParseFormulaFromJSON(o, var_name, *grid_ptr,
                               math_parser_ptr, storage_ptr);

        }

        if (has_filename) {

          const std::string filename = o["file"].get<std::string>();
            
          LOG_VERBOSE << "Reading variable "
                      << "\'" << var_name << "\'"
                      << " from file "
                      << "\'" << filename << "\'";

          std::ifstream variable_file(filename, std::ifstream::in | std::ifstream::binary);
          
          if (!variable_file) {
      
            LOG_ERROR << "Invalid input stream";
            std::abort();
      
          }

          const int nb_components = 1;

          ReadBinaryVariable(variable_file, 
                             storage_ptr->n_fast(),
                             storage_ptr->n_fast_padding(),
                             storage_ptr->n2(),
                             storage_ptr->n3(),
                             nb_components, data);

          variable_file.close();

        }
          
      } else {

        LOG_ERROR << "Variable "
                  << "\'" << var_name << "\'"
                  << " not found in database";
        
        std::abort();
        
      }

    }

    const bool has_timeloop = o["timeloop"].is<picojson::object>();
    
    if (has_timeloop) {

      picojson::object v = o["timeloop"].get<picojson::object>();

      const bool has_niter = v["niter"].is<double>();
      const bool has_tfinal = v["tfinal"].is<double>();

      // Exactly 1 of the above options must be input.
      const bool correct_input =
        ((has_niter || has_tfinal) && !(has_niter && has_tfinal));

      if (!correct_input) {

        LOG_ERROR << "In \'timeloop\', one and only one option among "
                  << "\'niter\' and \'tfinal\' must be set";
        std::abort();

      }

      const bool has_dt = 
        (v["dt"].is<double>() || v["dt"].is<std::string>());

      if (!has_dt) {

        LOG_ERROR << "In \'timeloop\', missing mandatory "
                  << "\'dt\' (timestep) value";
        std::abort();

      }

      // CFL condition estimation should be put out of this function.
      RealT maximum_velocity = std::numeric_limits<RealT>::min();
      
      const RealT* velocity = storage_ptr->RawDataSlowDimension(variable::VELOCITY);

      const int n3 = storage_ptr->n3();
      const int n2 = storage_ptr->n2();
      const int n_fast = storage_ptr->n_fast();
      const int n_fast_pad = n_fast + storage_ptr->n_fast_padding();

      for (int i3 = 0; i3 < n3; ++i3) {
        for (int i2 = 0; i2 < n2; ++i2) {
          for (int ifast = 0; ifast < n_fast; ++ifast) {

            const size_t index = 
              n2 * n_fast_pad * i3 + n_fast_pad * i2 + ifast;

            maximum_velocity = fmax(maximum_velocity, velocity[index]);

          }
        }
      }

      RealT min_grid_size = std::numeric_limits<RealT>::max();

      if (grid_ptr->n_fast() > 1)
        min_grid_size = fmin(min_grid_size, grid_ptr->dx_fast());

      if (grid_ptr->n_medium() > 1)
        min_grid_size = fmin(min_grid_size, grid_ptr->dx_medium());

      if (grid_ptr->n_slow() > 1)
        min_grid_size = fmin(min_grid_size, grid_ptr->dx_slow());

      assert(0.0 < min_grid_size);
      assert(0.0 < maximum_velocity);
      
      const RealT cfl_stability = 0.5 * sqrt(2.0) * min_grid_size / maximum_velocity;

      LOG_INFO << "Maximum velocity:" << maximum_velocity 
               << ", minimum grid size: " << min_grid_size
               << ", CFL:" << cfl_stability;
      
      math_parser_ptr->AddConstant("CFL", cfl_stability);

      RealT dt = 0.0;

      if (v["dt"].is<double>()) {

        dt = v["dt"].get<double>();

      } else if (v["dt"].is<std::string>()) {

        math_parser_ptr->SetExpression(v["dt"].get<std::string>());
        dt = math_parser_ptr->SafeEval();
      
      }

      RealT final_time = 0.0;
      int nb_iter = 0;

      if (has_niter) {

        nb_iter = int(v["niter"].get<double>());
        final_time = nb_iter * dt;

      }

      if (has_tfinal) {

        final_time = v["tfinal"].get<double>();
        assert(0.0 < dt);
        nb_iter = ceil(final_time / dt);

      }

      *timeloop_manager_ptr = TimeloopManager(dt, final_time, nb_iter);
      
      LOG_VERBOSE << "Timeloop will be executed with "
                  << "nb_iter=" << timeloop_manager_ptr->nb_iter() << ", "
                  << "final_time=" << timeloop_manager_ptr->final_time() << ", "
                  << "dt=" << timeloop_manager_ptr->dt();

    }

    storage_ptr->Validate();
  }

  const bool has_output = o["output"].is<picojson::array>();

  if (has_output) {

    picojson::array output_events_json = o["output"].get<picojson::array>();

    for (auto it = output_events_json.begin(); it != output_events_json.end(); ++it) {

      OutputEvent event = OutputEvent();
      event.ParseFromJSON(*it);
      event.Validate();
      output_events_ptr->push_back(event);

      LOG_INFO << "Found output event: "
               << "type=\'" << event.type() << "\', "
               << "format=\'" << event.format() << "\', "
               << "filename=\'" << event.stream_name() << "\', "
               << "rhythm=\'" << event.rhythm() << "\', "
               << "var_name=\'" << event.var_name() << "\', "
               << "formula=\'" << event.formula() << "\'";

    }
  }

}
