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

#ifndef TIMELOOP_MANAGER_HPP
#define TIMELOOP_MANAGER_HPP

#include <fstream>
#include <iosfwd>
#include <string>

class TimeloopManager {
public:
  TimeloopManager();
  TimeloopManager(RealT dt, RealT tfinal, int niter);
  int nb_iter() { return m_nb_iter; }
  int nb_iter() const { return m_nb_iter; }
  RealT dt() { return m_dt; }
  RealT dt() const { return m_dt; }
  RealT final_time() { return m_final_time; }
  RealT final_time() const { return m_final_time; }
  RealT current_time() { return m_current_time; }
  RealT current_time() const { return m_current_time; }
  int current_iter() { return m_current_iter; }
  int current_iter() const { return m_current_iter; }
  void Update();
private:
  int m_nb_iter;
  RealT m_dt;
  RealT m_final_time;
  RealT m_current_time;
  RealT m_current_iter;
};

#endif // TIMELOOP_MANAGER_HPP
