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

#include "mathematical_parser.hpp"

#include <cassert>

#include "muParser.h"

#include "common.hpp"
#include "multidimensional_storage.hpp"
#include "rectilinear_grid.hpp"
#include "timer.hpp"
#include "variable_definitions.hpp"

MathematicalParser::MathematicalParser():
  m_parser() {}

void MathematicalParser::AddConstant(const std::string& constant_name,
                                     RealT constant_value) {

  m_parser.DefineConst(constant_name, constant_value);
  
}

void MathematicalParser::SetGridConstants(const RectilinearGrid3D& grid) {

  AddConstant("dx", grid.dx_fast());
  AddConstant("dy", grid.dx_medium());
  AddConstant("dz", grid.dx_slow());

  AddConstant("nx", grid.n_fast());
  AddConstant("ny", grid.n_medium());
  AddConstant("nz", grid.n_slow());

  AddConstant("xmin", grid.x_fast_min());
  AddConstant("xmax", grid.x_fast_max());
  AddConstant("ymin", grid.x_medium_min());
  AddConstant("ymax", grid.x_medium_max());
  AddConstant("zmin", grid.x_slow_min());
  AddConstant("zmax", grid.x_slow_max());

  // Nyquist spatial frequency.
  if (0.0 < grid.dx_fast())
    AddConstant("nyquist_x", 1.0 / (2.0 * grid.dx_fast()));

  if (0.0 < grid.dx_medium())
    AddConstant("nyquist_y", 1.0 / (2.0 * grid.dx_medium()));

  if (0.0 < grid.dx_slow())
    AddConstant("nyquist_z", 1.0 / (2.0 * grid.dx_slow()));

}

void MathematicalParser::PrintConstants() {

  auto constant_database = m_parser.GetConst();

  if (!constant_database.empty())
    for (auto item = constant_database.begin(); item != constant_database.end(); ++item)
      LOG_INFO << "  " << item->first << " =  " << item->second;

}

void MathematicalParser::EvaluateExpression(const std::string& expression,
                                            const RectilinearGrid3D& grid,
                                            MultiDimensionalStorage4D* variable_storage_ptr) {

  const int n_fast = variable_storage_ptr->n_fast();
  const int n_fast_padding = variable_storage_ptr->n_fast_padding();
  const int n_fast_padded = n_fast + n_fast_padding;
  const int n_medium = variable_storage_ptr->n2();
  const int n_slow = variable_storage_ptr->n3();

  // Temporary buffer to hold parser result.
  RealT* buffer = (RealT*) malloc(n_fast * sizeof(RealT));

  // Temporary buffer to hold grid coordinates in the fast
  // direction. We have to copy them because we do not want the parser
  // to modify them by user input.
  // RealT* x_grid = (RealT*) malloc(n_fast * sizeof(RealT));

  for (int islow = 0; islow < n_slow; ++islow) {
    for (int imedium = 0; imedium < n_medium; ++imedium) {

      std::vector<RealT> x_grid = grid.fast_coordinates();

      const size_t index_base = 
        n_medium * n_fast_padded * islow + n_fast_padded * imedium;

      // for (int ifast = 0; ifast < n_fast; ++ifast) 
      //   x_grid[ifast] = grid.fast_coordinates()[ifast];

      // Add pointers to variables to the parser.
      for (int ivar = 0; ivar < variable::NB_VARIABLES; ++ ivar) {

        RealT* data = variable_storage_ptr->RawDataSlowDimension(ivar);
        m_parser.DefineVar(variable::VARIABLE_NAMES[ivar], &(data[index_base]));

      }

      // Add grid coordinates to the parser.
      m_parser.DefineVar("x", VectorRawData<RealT>(x_grid));
      AddConstant("y", grid.medium_coordinates()[imedium]);
      AddConstant("z", grid.slow_coordinates()[islow]);

      m_parser.SetExpr(expression);

      try {

        m_parser.Eval(buffer, n_fast);
        
      } catch(mu::Parser::exception_type &e) {
    
        LOG_ERROR << "In formula "
                  << "\'" << e.GetExpr() << "\'"
                  << ", " << e.GetMsg();
  
        std::abort();
      }
    }
  }

  free(buffer);

}

std::vector<RealT> MathematicalParser::ReduceExpression(const std::string& expression,
                                                        const RectilinearGrid3D& grid,
                                                        MultiDimensionalStorage4D* variable_storage_ptr) {

  std::vector<RealT> result;
  // L1, L2 and Linfinite norms.
  result.assign(3, 0.0);

  RealT l1 = 0.0;
  RealT l2 = 0.0;
  RealT linf = 0.0;

  const int n_fast = variable_storage_ptr->n_fast();
  const int n_fast_padding = variable_storage_ptr->n_fast_padding();
  const int n_fast_padded = n_fast + n_fast_padding;
  const int n_medium = variable_storage_ptr->n2();
  const int n_slow = variable_storage_ptr->n3();

  // Temporary buffers to hold parser result.
  RealT* buffer = (RealT*) malloc(n_fast * sizeof(RealT));

  // Temporary buffer to hold grid coordinates in the fast
  // direction. We have to copy them because we do not want the parser
  // to modify them by user input.
  // RealT* x_grid = (RealT*) malloc(n_fast * sizeof(RealT));

  for (int islow = 0; islow < n_slow; ++islow) {
    for (int imedium = 0; imedium < n_medium; ++imedium) {

      std::vector<RealT> x_grid = grid.fast_coordinates();

      const size_t index_base = 
        n_medium * n_fast_padded * islow + n_fast_padded * imedium;

      // Add pointers to variables to the parser.
      for (int ivar = 0; ivar < variable::NB_VARIABLES; ++ ivar) {

        RealT* data = variable_storage_ptr->RawDataSlowDimension(ivar);
        m_parser.DefineVar(variable::VARIABLE_NAMES[ivar], &(data[index_base]));

      }

      // Add grid coordinates to the parser.
      m_parser.DefineVar("x", VectorRawData<RealT>(x_grid));
      AddConstant("y", grid.medium_coordinates()[imedium]);
      AddConstant("z", grid.slow_coordinates()[islow]);

      m_parser.SetExpr(expression);

      try {

        m_parser.Eval(buffer, n_fast);
        
      } catch(mu::Parser::exception_type &e) {
    
        LOG_ERROR << "In formula "
                  << "\'" << e.GetExpr() << "\'"
                  << ", " << e.GetMsg();
  
        std::abort();
      }

      for (int i = 0; i < n_fast; ++i) {

        l1 += fabsf(buffer[i]);
        l2 += buffer[i] * buffer[i];
        linf = fmax(linf, buffer[i]);

      }

    }
  }

  if (grid.n_fast() > 1) {

    l1 *= grid.dx_fast();
    l2 *= grid.dx_fast() * grid.dx_fast();

  }

  if (grid.n_medium() > 1) {

    l1 *= grid.dx_medium();
    l2 *= grid.dx_medium() * grid.dx_medium();

  }

  if (grid.n_slow() > 1) {

    l1 *= grid.dx_slow();
    l2 *= grid.dx_slow() * grid.dx_slow();

  }

  free(buffer);

  LOG_INFO << "Norm of "
           << "\'" << expression << "\' (L1, L2, Linf):"
           << "(" << l1 << ", " << l2 << ", " << linf << ")";

  result.at(0) = l1;
  result.at(1) = l2;
  result.at(2) = linf;

  return result;

}

void MathematicalParser::SetExpression(const std::string& expression) {

  m_parser.SetExpr(expression);

}
 
double MathematicalParser::SafeEval() {

  double value = 0.0;

  try {

    value = m_parser.Eval();
        
  } catch(mu::Parser::exception_type &e) {
    
    LOG_ERROR << "In formula "
              << "\'" << e.GetExpr() << "\'"
              << ", " << e.GetMsg();
  
    std::abort();
  }
  
  return value;
}
