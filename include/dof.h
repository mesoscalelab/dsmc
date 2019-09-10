// Copyright (c) 2019 Santosh Ansumali @ JNCASR
// See LICENSE

#pragma once

#include <cmath>
#include <vector>

template <typename real>
class dofParticle
{
public:
  real velocity() const { return vel; }
  real position() const { return pos; }

  real& velocity() { return vel; }
  real& position() { return pos; }

  void
  setPositionPeriodic(real boundBox)
  {
    if(pos > (boundBox)) {
      pos = std::fmod(pos,(boundBox));
    }
    if(pos < 0) {
      pos = (boundBox) - std::fmod(std::fabs(pos), boundBox);
    }
  }

private:
  real vel;
  real pos;
};
