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

#ifndef ARRAY_VTK_IO_HPP
#define ARRAY_VTK_IO_HPP

typedef enum _DataFormat {ASCII, BINARY} DataFormat;

void WriteVTKXmlVariable(const std::string& variable_name,
			 DataFormat format, int nb_components,
			 int n_fast, int n_fast_padding, 
			 int n_medium, int n_slow, 
			 RealT* data,
			 std::ostream* os_ptr);

#endif // ARRAY_VTK_IO_HPP
