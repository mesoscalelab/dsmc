// Copyright (c) 2019 Santosh Ansumali @ JNCASR
// See LICENSE

#pragma once

#include "particle.h"

template<int DIM,typename real,typename dof>
class cellList
{
public:
  cellList() {}

  ~cellList() {}

  std::vector<int> cellParticleList;
};
