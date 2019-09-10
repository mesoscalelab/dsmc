// Copyright (c) 2019 Santosh Ansumali @ JNCASR
// See LICENSE

#pragma once

#include "dof.h"

#define noCells 50	//only Z direction cells
#define nParticles noCells*1000

template<int DIM,typename real,typename dof>
class particle
{
public:
  void initialiseParticle()
  {
      numParticles = nParticles;
      numElements = DIM*numParticles;
      data = new dof[(int)numElements];
      cellCoord = new int[(int)numElements];
  }

  dof& getData(int particlenum, int dim) {return data[(int)numParticles*(dim - 1) + particlenum];}

  int& getcellCoord(int particlenum, int dim) {return cellCoord[(int)numParticles*(dim - 1) + particlenum];}

  real getnumParticles() {return numParticles;}

  int getDim() {return DIM;}

  real getnumElements() {return numElements;}

private:
  real numParticles;
  real numElements;
  dof* data;
  int* cellCoord;
};
