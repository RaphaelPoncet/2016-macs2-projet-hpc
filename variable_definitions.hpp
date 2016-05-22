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

#ifndef VARIABLE_DEFINITIONS_HPP
#define VARIABLE_DEFINITIONS_HPP

typedef enum _VariableSupport {NODE, CELL} VariableSupport;

enum VariableFlags {
  INITIALIZABLE = 1, // << 0
  WRITTEN = 1 << 1,
  DUMMY1 = 1 << 2,
  DUMMY2 = 1 << 3
};

namespace variable {

  const VariableSupport VARIABLE_SUPPORT = NODE;
  
  enum {VELOCITY,
        PRESSURE_0,
        PRESSURE_1,
        PRESSURE_REF,
        LAPLACE_PRESSURE,
        NB_VARIABLES};

  static const char* VARIABLE_NAMES[NB_VARIABLES] = 
    {"velocity", "pressure_0", "pressure_1", "pressure_ref", "laplace_pressure"};
  
  static const unsigned int VARIABLE_FLAGS[NB_VARIABLES] = 
    {WRITTEN | INITIALIZABLE,
     WRITTEN | INITIALIZABLE,
     WRITTEN | INITIALIZABLE,
     WRITTEN | INITIALIZABLE,
     WRITTEN | INITIALIZABLE};
    
}

#endif // VARIABLE_DEFINITIONS_HPP
