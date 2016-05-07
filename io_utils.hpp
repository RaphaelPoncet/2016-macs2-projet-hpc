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

#ifndef IO_UTILS_HPP
#define IO_UTILS_HPP

template <typename T>
T ByteSwap(T t) {
  
  const int number_of_bytes = sizeof(T);

  union {
    T t;
    unsigned char bytes[number_of_bytes];
  } dat1, dat2;

  dat1.t = t;

  for (int id_byte = 0; id_byte < number_of_bytes / 2; ++id_byte)
    dat2.bytes[id_byte] = dat1.bytes[number_of_bytes - 1 - id_byte];

  return dat2.t;

}

// Convert a string into any type, returning false for invalid input
template <class T>
bool ParseToken(const std::string& s, T* t_ptr) {
  
  std::string garbage;
  std::istringstream iss(s);
  bool begins_good = !(iss >> *t_ptr).fail();
  iss.clear();
  iss.seekg(0, std::ios::beg);
  iss >> *t_ptr >> garbage;

  return begins_good && garbage.empty();

}

#endif // IO_UTILS_HPP
