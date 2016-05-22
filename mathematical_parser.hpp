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

#ifndef MATHEMATICAL_PARSER_HPP
#define MATHEMATICAL_PARSER_HPP

#include <iosfwd>

#include "muParser.h"

class RectilinearGrid3D;
class MultiDimensionalStorage4D;

class MathematicalParser {
public:
  MathematicalParser();
  void SetGridConstants(const RectilinearGrid3D&);
  void AddConstant(const std::string& constant_name, RealT constant_value);
  void PrintConstants();
  void EvaluateExpression(const std::string& expression,
                          const RectilinearGrid3D& grid,
                          MultiDimensionalStorage4D* variable_storage_ptr);
  void SetExpression(const std::string& expression);
  double SafeEval();
private:
  mu::Parser m_parser;
};

#endif // MATHEMATICAL_PARSER_HPP
