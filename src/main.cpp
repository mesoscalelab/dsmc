// Copyright (c) 2019 Santosh Ansumali @ JNCASR
// See LICENSE

#include <hardSphere.h>

int
main()
{   
  srand48(5);
  hardSphere<3,double,dofParticle<double> > p;
  p.hardSphereMain();
  return 0;
}
