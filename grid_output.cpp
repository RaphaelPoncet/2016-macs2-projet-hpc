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

#include "grid_output.hpp"

#include <cassert>
#include <fstream>
#include <sstream>

#include "array_vtk_io.hpp"
#include "common.hpp"
#include "multidimensional_storage.hpp"
#include "rectilinear_grid.hpp"
#include "variable_definitions.hpp"

void OutputGridAndData(const RectilinearGrid3D& grid, 
                       const MultiDimensionalStorage4D& variable_storage,
                       std::ofstream* os_ptr) {

  const VariableSupport variable_support = variable::VARIABLE_SUPPORT;
  const VTKDataFormat format = BINARY;

  grid.WriteHeaderVTKXml(os_ptr);
  grid.WriteVTKXmlAscii(os_ptr);

  if (format == ASCII) {

    *os_ptr << "<CellData>\n";
  
    if (variable_support == CELL) {

      const int n_slow = (grid.n_slow() == 1 ? 1 : grid.n_slow() - 1);
      assert(0 < n_slow);
      const int n_medium = (grid.n_medium() == 1 ? 1 : grid.n_medium() - 1);
      assert(0 < n_medium);
      const int n_fast = (grid.n_fast() == 1 ? 1 : grid.n_fast() - 1);
      assert(0 < n_fast);

      const size_t nb_grid_elements = n_slow * n_medium * n_fast;

      const size_t nb_elements_storage_3d = 
        variable_storage.n3() * variable_storage.n2() * variable_storage.n_fast();

      assert(nb_elements_storage_3d == nb_grid_elements);

      for (int ivar = 0; ivar < variable::NB_VARIABLES; ++ivar) {

        const std::string variable_name = variable::VARIABLE_NAMES[ivar];
        const int nb_components = 1;

        if (variable::VARIABLE_FLAGS[ivar] & WRITTEN) {

          size_t grid_data_offset_in_bytes = 0;
          const size_t grid_data_size_in_bytes = nb_grid_elements * sizeof(RealT);
          const RealT* data = variable_storage.RawDataSlowDimension(ivar);

          *os_ptr << "<DataArray ";
          WriteVTKXmlVariableHeader(variable_name, format, nb_components, grid_data_size_in_bytes, 
                                    &grid_data_offset_in_bytes, os_ptr);
          *os_ptr << ">\n";

          WriteVTKXmlVariable(format, nb_components, 
                              variable_storage.n_fast(), variable_storage.n_fast_padding(),
                              variable_storage.n2(), variable_storage.n3(), 
                              data, &grid_data_offset_in_bytes, os_ptr);

          *os_ptr << "</DataArray>\n";
    
        }
      }
    }

    *os_ptr << "</CellData>\n";

    *os_ptr << "<PointData>\n";

    if (variable_support == NODE) {

      const size_t nb_grid_elements = 
        grid.n_slow() * grid.n_medium() * grid.n_fast();

      const size_t nb_elements_storage_3d = 
        variable_storage.n_fast() * variable_storage.n2() * variable_storage.n3();

      assert(nb_elements_storage_3d == nb_grid_elements);

      for (int ivar = 0; ivar < variable::NB_VARIABLES; ++ivar) {

        const std::string variable_name = variable::VARIABLE_NAMES[ivar];
        const int nb_components = 1;
        UNUSED(nb_components);

        if (variable::VARIABLE_FLAGS[ivar] & WRITTEN) {

          const size_t grid_data_size_in_bytes = nb_grid_elements * sizeof(RealT);
          size_t grid_data_offset_in_bytes = 0;     
          const RealT* data = variable_storage.RawDataSlowDimension(ivar);
          const int nb_components = 1;

          *os_ptr << "<DataArray ";
          WriteVTKXmlVariableHeader(variable_name, format, nb_components, grid_data_size_in_bytes, 
                                    &grid_data_offset_in_bytes, os_ptr);
          *os_ptr << ">\n";

          WriteVTKXmlVariable(format, nb_components, 
                              variable_storage.n_fast(), variable_storage.n_fast_padding(),
                              variable_storage.n2(), variable_storage.n3(), 
                              data, &grid_data_offset_in_bytes, os_ptr);

          *os_ptr << "</DataArray>\n";

        }
      }
    }

    *os_ptr << "</PointData>\n";

  } else if (format == BINARY) {

    std::stringstream cell_data_arrays_header;
    std::stringstream point_data_arrays_header;

    cell_data_arrays_header << "<CellData>\n";
    point_data_arrays_header << "<PointData>\n";

    size_t grid_data_offset_in_bytes = 0;

    for (int ivar = 0; ivar < variable::NB_VARIABLES; ++ivar) {

      const std::string variable_name = variable::VARIABLE_NAMES[ivar];
      const int nb_components = 1;

      if (variable::VARIABLE_FLAGS[ivar] & WRITTEN) {

        if (variable_support == CELL) {

          const int n_slow = (grid.n_slow() == 1 ? 1 : grid.n_slow() - 1);
          assert(0 < n_slow);
          const int n_medium = (grid.n_medium() == 1 ? 1 : grid.n_medium() - 1);
          assert(0 < n_medium);
          const int n_fast = (grid.n_fast() == 1 ? 1 : grid.n_fast() - 1);
          assert(0 < n_fast);

          const size_t nb_grid_elements = n_slow * n_medium * n_fast;

          const size_t nb_elements_storage_3d = 
            variable_storage.n3() * variable_storage.n2() * variable_storage.n_fast();

          assert(nb_elements_storage_3d == nb_grid_elements);

          const size_t grid_data_size_in_bytes = nb_grid_elements * sizeof(RealT);

          cell_data_arrays_header << "<DataArray ";
          WriteVTKXmlVariableHeader(variable_name, format, nb_components, grid_data_size_in_bytes, 
                                    &grid_data_offset_in_bytes, &cell_data_arrays_header);
          cell_data_arrays_header << ">\n";
          cell_data_arrays_header << "</DataArray>\n";

        } else if (variable_support == NODE) {

          const size_t nb_grid_elements = 
            grid.n_slow() * grid.n_medium() * grid.n_fast();

          const size_t nb_elements_storage_3d = 
            variable_storage.n_fast() * variable_storage.n2() * variable_storage.n3();

          if (nb_elements_storage_3d != nb_grid_elements) {

            LOG_ERROR << "Incompatible grid and storage sizes: "
                      << "grid size = (" 
                      << grid.n_fast() << ", " << grid.n_medium() << ", " << grid.n_slow() << ")  vs "
                      << "storage size = (" << variable_storage.n_fast() << ", " 
                      << variable_storage.n2() << ", " 
                      << variable_storage.n3() << ")";
            std::abort();

          }

          const size_t grid_data_size_in_bytes = nb_grid_elements * sizeof(RealT);

          point_data_arrays_header << "<DataArray ";
          WriteVTKXmlVariableHeader(variable_name, format, nb_components, grid_data_size_in_bytes, 
                                    &grid_data_offset_in_bytes, &point_data_arrays_header);
          point_data_arrays_header << ">\n";
          point_data_arrays_header << "</DataArray>\n";

        }
      }
    }

    cell_data_arrays_header << "</CellData>\n";
    point_data_arrays_header << "</PointData>\n";

    *os_ptr << cell_data_arrays_header.str();
    *os_ptr << point_data_arrays_header.str();
    
    *os_ptr << "</Piece>\n"
                << "</RectilinearGrid>\n";

    // The underscore is intentional.
    *os_ptr << "<AppendedData encoding=\"raw\">\n_";
  
    if (variable_support == CELL) {

      const int n_slow = (grid.n_slow() == 1 ? 1 : grid.n_slow() - 1);
      assert(0 < n_slow);
      const int n_medium = (grid.n_medium() == 1 ? 1 : grid.n_medium() - 1);
      assert(0 < n_medium);
      const int n_fast = (grid.n_fast() == 1 ? 1 : grid.n_fast() - 1);
      assert(0 < n_fast);

      for (int ivar = 0; ivar < variable::NB_VARIABLES; ++ivar) {

        const int nb_components = 1;
	  
        size_t grid_data_offset_in_bytes = 0;
        const RealT* data = variable_storage.RawDataSlowDimension(ivar);

        if (variable::VARIABLE_FLAGS[ivar] & WRITTEN) {

          WriteVTKXmlVariable(format, nb_components, 
                              variable_storage.n_fast(), variable_storage.n_fast_padding(),
                              variable_storage.n2(), variable_storage.n3(), 
                              data, &grid_data_offset_in_bytes, os_ptr);

        }
      }

    } else if (variable_support == NODE) {

      const size_t nb_grid_elements = 
        grid.n_slow() * grid.n_medium() * grid.n_fast();

      const size_t nb_elements_storage_3d = 
        variable_storage.n_fast() * variable_storage.n2() * variable_storage.n3();

      assert(nb_elements_storage_3d == nb_grid_elements);

      for (int ivar = 0; ivar < variable::NB_VARIABLES; ++ivar) {

        const int nb_components = 1;
	  
        size_t grid_data_offset_in_bytes = 0;
        const RealT* data = variable_storage.RawDataSlowDimension(ivar);

        if (variable::VARIABLE_FLAGS[ivar] & WRITTEN) {

          WriteVTKXmlVariable(format, nb_components, 
                              variable_storage.n_fast(), variable_storage.n_fast_padding(),
                              variable_storage.n2(), variable_storage.n3(), 
                              data, &grid_data_offset_in_bytes, os_ptr);

        }
      }
    }

    *os_ptr << "\n</AppendedData>\n";

  }

  grid.WriteFooterVTKXml(os_ptr);

}
