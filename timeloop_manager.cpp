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

#include "timeloop_manager.hpp"

TimeloopManager::TimeloopManager():
  m_nb_iter(0), m_dt(0.0), m_final_time(0.0), m_current_time(0.0), m_current_iter(0) {};

TimeloopManager::TimeloopManager(RealT dt, RealT tfinal, int niter):
  m_nb_iter(niter), m_dt(dt), m_final_time(tfinal), m_current_time(0.0), m_current_iter(0) {};

void TimeloopManager::Update() {

  m_current_iter += 1;
  m_current_time += m_dt;

}
