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

#include "common.hpp"
#include "rectilinear_grid.hpp"

#include <fstream>

static void ValidateInput(int n1, int n2, int n3,
			  RealT x1_min, RealT x1_max,
			  RealT x2_min, RealT x2_max,
			  RealT x3_min, RealT x3_max) { 

  assert(0 < n1);
  assert(0 < n2);
  assert(0 < n3);

  assert(x1_min < x1_max);
  assert(x2_min < x2_max);
  assert(x3_min < x3_max);

}

static void LinearArray(RealT value_min, RealT value_max, int nb_elements, RealT* array) {

  assert(value_min <= value_max);
  assert(1 < nb_elements);
  assert(array != NULL);

  const RealT sampling = (value_max - value_min) / (RealT)(nb_elements - 1);
  assert(0 <= sampling);

  for (int i = 0; i < nb_elements; ++i)
    array[i] = value_min + (RealT)i * sampling;
  
}

RectilinearGrid3D::RectilinearGrid3D(RealT x_fast_min, RealT x_fast_max, int n_fast,
				     RealT x_medium_min, RealT x_medium_max, int n_medium,
				     RealT x_slow_min, RealT x_slow_max, int n_slow):
  m_n_fast(n_fast), m_n_medium(n_medium), m_n_slow(n_slow), 
  m_x_fast_min(x_fast_min), m_x_fast_max(x_fast_max), 
  m_x_medium_min(x_medium_min), m_x_medium_max(x_medium_max), 
  m_x_slow_min(x_slow_min), m_x_slow_max(x_slow_max) {

  LOG_DEBUG << "Initializing 3D rectilinear grid with (n1, n2, n3) = "
	    << "(" << m_n_fast << ", "
	    << m_n_medium << ", "
	    << m_n_slow << ")"
	    << " nodes...";

  ValidateInput(m_n_fast, m_n_medium, m_n_slow, 
		m_x_fast_min, m_x_fast_max,
		m_x_medium_min, m_x_medium_max,
		m_x_slow_min, m_x_slow_max);

  const RealT dummy_fill = - 3.1415926;

  m_fast_coordinates.assign(m_n_fast, dummy_fill);
  LinearArray(m_x_fast_min, m_x_fast_max, m_n_fast, VectorRawData<RealT>(m_fast_coordinates));

  m_medium_coordinates.assign(m_n_medium, dummy_fill);
  LinearArray(m_x_medium_min, m_x_medium_max, m_n_medium, VectorRawData<RealT>(m_medium_coordinates));
  m_slow_coordinates.assign(m_n_slow, dummy_fill);
  LinearArray(m_x_slow_min, m_x_slow_max, m_n_slow, VectorRawData<RealT>(m_slow_coordinates));
		
  LOG_DEBUG << "Initialization done.\n";

}

void RectilinearGrid3D::WriteHeaderVTKXml(std::ofstream* os_ptr) {

  *os_ptr << "<?xml version=\"1.0\"?>\n"
	  << "<VTKFile type=\"RectilinearGrid\" "
	  << "version=\"0.1\" "
	  << "byte_order=\"LittleEndian\">\n";
  
  *os_ptr << "<RectilinearGrid WholeExtent=\""
	  << 0 << " " << m_n_fast - 1 << " "
	  << 0 << " " << m_n_medium - 1<< " "
	  << 0 << " " << m_n_slow - 1
	  << "\">\n";

  *os_ptr << "<Piece Extent=\""
	  << 0 << " " << m_n_fast - 1 << " "
	  << 0 << " " << m_n_medium - 1 << " "
	  << 0 << " " << m_n_slow - 1 << "\""
	  << " Dimensions=\"" 
	  << m_n_fast << " " 
	  << m_n_medium << " " 
	  << m_n_slow << "\""
	  << ">\n";

}

void RectilinearGrid3D::WriteFooterVTKXml(std::ofstream* os_ptr) {

  *os_ptr << "</Piece>\n"
	  << "</RectilinearGrid>\n";
  
  *os_ptr << "</VTKFile>\n";

}

void RectilinearGrid3D::WriteVTKXmlAscii(std::ofstream* os_ptr) {

  *os_ptr << "<Coordinates>\n";

  *os_ptr << "<DataArray "
	  << "type=\"Float32\" "
	  << "format=\"ascii\">\n";

  const int nb_points_fast = m_fast_coordinates.size();

  for (int i = 0; i < nb_points_fast; ++i)
    *os_ptr << m_fast_coordinates.at(i) << "\n";

  *os_ptr << "</DataArray>\n";

  *os_ptr << "<DataArray "
	  << "type=\"Float32\" "
	  << "format=\"ascii\">\n";

  const int nb_points_medium = m_medium_coordinates.size();

  for (int i = 0; i < nb_points_medium; ++i)
    *os_ptr << m_medium_coordinates.at(i) << "\n";

  *os_ptr << "</DataArray>\n";

  *os_ptr << "<DataArray "
	  << "type=\"Float32\" "
	  << "format=\"ascii\">\n";

  const int nb_points_slow = m_slow_coordinates.size();
  
  for (int i = 0; i < nb_points_slow; ++i)
    *os_ptr << m_slow_coordinates.at(i) << "\n";
  
  *os_ptr << "</DataArray>\n";

  *os_ptr << "</Coordinates>\n";

}
