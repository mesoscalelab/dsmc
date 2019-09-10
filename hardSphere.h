/*********************************************************************************************
 *   Copyright (c) <2018>, <Santosh Ansumali@JNCASR>                                         *
 *   All rights reserved.                                                                    *
 *   Redistribution and use in source and binary forms, with or without modification, are    *
 *   permitted provided that the following conditions are met:                               *
 *                                                                                           *
 *    1. Redistributions of source code must retain the above copyright notice, this list of *
 *       conditions and the following disclaimer.                                            *
 *    2. Redistributions in binary form must reproduce the above copyright notice, this list *
 *       of conditions and the following disclaimer in the documentation and/or other        *
 *       materials provided with the distribution.                                           *
 *    3. Neither the name of the <JNCASR> nor the names of its contributors may be used to   *
 *       endorse or promote products derived from this software without specific prior       *
 *       written permission.                                                                 *
 *                                                                                           *
 *       THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND     *
 *       ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED       *
 *       WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  *
 *       IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,    *
 *       INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,      *
 *       BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,       *
 *       DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF     *
 *       LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE     *
 *       OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED   *
 *       OF THE POSSIBILITY OF SUCH DAMAGE.                                                  *
 *                                                                                           *
 *       Suggestions:          ansumali@jncasr.ac.in                                         *
 *       Bugs:                 ansumali@jncasr.ac.in                                         *
 *                                                                                           *
 *********************************************************************************************/


/*********************************************************************************************

 *  @Author: Akshay Chandran and Santosh Ansumali                                            *

 *********************************************************************************************/


#ifndef HARDSPHERE_H_INCLUDED
#define HARDSPHERE_H_INCLUDED

#include "cellList.h"

#define gy .00001

#define dia 1.
#define M 1.
#define kT 1.
#define density .001
#define mfp 225.

#define boundBox 10000.	//X and Y direction bounds for periodicity
#define cSize 1
#define deltaZed .2*mfp
#define heightZ noCells*deltaZed
#define deltatime .1
#define totVol heightZ*boundBox*boundBox	//total volume

#define pi 3.14

template<int DIM,typename real,typename dof>
class hardSphere
{
    public:
        hardSphere()
        {
            Oparticle = new particle<DIM,real,dof>;
            Oparticle->initialiseParticle();
            numCells = noCells;

            cList = new cellList<DIM,real,dof>[numCells];
            uWall = 0.;
            mass = M;
            volume = totVol;
            deltat = deltatime;
            height = heightZ;
            deltaZ = deltaZed;
            diameter = dia;
	    extent[0] = boundBox;
	    extent[1] = boundBox;
	    extent[2] = height;
        }
        ~hardSphere()
        {
            //delete Oparticle;
            //delete[] cList;
        }

        real totKinEnergy();

        void evolveSystem();

        void updateVelocityBottomWall(int particleIndex, int dim);

        void updateVelocityTopWall(int particleIndex);

	void applyBoundaryConditions(int,int);

	void wallBoundary(int,int);

        void periodicBoundary(int,int);

        void clearCellParticles();

        void hardSphereMain();

        void initPositions();

        void initVelocity();

        void updateCellCoord();

        real gaussianRandom();

        real relSpeed(int particleOne, int particleTwo);

        void Collision();

        void processCollision(int particleOne, int particleTwo);

        void postProcess(int it);

        real getDensity(int cellNumber);

    private:
        int numCells;
        cellList<DIM,real,dof>* cList;
	particle<DIM,real,dof>* Oparticle;

        real uWall;

        real mass;
        real volume;
        real deltat;
        real height;
        real deltaZ;
        real diameter;
	real extent[DIM];
};

template<int DIM,typename real,typename dof>
real hardSphere<DIM,real,dof>::totKinEnergy()
{
    real KE = 0.;
    for(int dim = 1; dim <= DIM; dim++)
    {
        for(int particleIndex = 0; particleIndex < nParticles; particleIndex++)
        {
            KE += (Oparticle->getData(particleIndex,dim).velocity())*(Oparticle->getData(particleIndex,dim).velocity());
        }
    }

    return KE;
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::applyBoundaryConditions(int particleIndex, int dim)
{
	switch(dim)
	{
		case 1: periodicBoundary(particleIndex,1); break;
		case 2: periodicBoundary(particleIndex,2); break;
		case 3: wallBoundary(particleIndex,3); break;
	}
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::wallBoundary(int particleIndex,int dim)
{
	real deltatTemp;
	real posNew;
	real posOld = (Oparticle->getData(particleIndex,dim).position());
        real velOld = Oparticle->getData(particleIndex,dim).velocity();

	if(velOld > 0.)
	{
		deltatTemp = fabs((height - posOld)/(velOld));
		posNew = posOld + velOld*deltatTemp;
		posNew = posNew - velOld*(deltat - deltatTemp);
	}
	else
	{
		deltatTemp = fabs((posOld)/(velOld));
		posNew = posOld + velOld*deltatTemp;
		posNew = posNew - velOld*(deltat - deltatTemp);
	}

	Oparticle->getData(particleIndex,dim).position() = posNew;
	Oparticle->getData(particleIndex,dim).velocity() = -velOld;
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::evolveSystem()
{
    real posNew;
    real posOld;
    real velOld;
    real velNew;
    real deltatwBottom;
    real deltatwTop;

    for(int dim = 1; dim <= DIM; dim++)
    {
        for(int particleIndex = 0; particleIndex < Oparticle->getnumParticles(); particleIndex++)
        {
            posOld = fabs(Oparticle->getData(particleIndex,dim).position());
            velOld = Oparticle->getData(particleIndex,dim).velocity();

	    posNew = posOld + velOld*deltat;

	    if(dim == 2)
	    {
		posNew = posOld + velOld*deltat + .5*gy*deltat*deltat;
		velNew = velOld + gy*deltat;
		Oparticle->getData(particleIndex,dim).velocity() = velNew;
	    }

	    if(((posNew > extent[dim - 1]) || (posNew < 0.)))
	    {
		applyBoundaryConditions(particleIndex,dim);
	    }
	    else
	    {
		Oparticle->getData(particleIndex,dim).position() = posNew;
	    }
        }
    }
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::updateVelocityBottomWall(int particleIndex, int dim)
{
    real kTbyM = sqrt(kT/mass);
    switch(dim)
    {
        case 1: Oparticle->getData(particleIndex,1).velocity() = kTbyM*gaussianRandom();
        case 2: Oparticle->getData(particleIndex,2).velocity() = kTbyM*gaussianRandom();
        case 3: Oparticle->getData(particleIndex,3).velocity() = kTbyM*pow(-2.*log(drand48()),.5);
    }
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::updateVelocityTopWall(int particleIndex)
{
    //X and Y velocity components remain unchanged
    Oparticle->getData(particleIndex,3).velocity() = -Oparticle->getData(particleIndex,3).velocity();
    //reversal of normal component after specular collision with wall
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::periodicBoundary(int particleIndex,int dim)
{
    //periodic only in the X and Y directions
    {
        Oparticle->getData(particleIndex,dim).setPositionPeriodic(boundBox);
    }
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::clearCellParticles()
{
    for(int currCell = 0; currCell < numCells; currCell++)
    {
        cList[currCell].cellParticleList.clear();
    }
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::initPositions()
{
    for(int particleIndex = 0; particleIndex < nParticles; particleIndex++)
    {
        Oparticle->getData(particleIndex,1).position() = boundBox*drand48();
    }
    for(int particleIndex = 0; particleIndex < nParticles; particleIndex++)
    {
        Oparticle->getData(particleIndex,2).position() = boundBox*drand48();
    }
    for(int particleIndex = 0; particleIndex < nParticles; particleIndex++)
    {
        Oparticle->getData(particleIndex,3).position() = height*drand48();
    }
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::initVelocity()
{
    real scaleFactor = 0.;  //shifting mean to scalefactor value
    for(int particleIndex = 0; particleIndex < nParticles; particleIndex++)
    {
        Oparticle->getData(particleIndex,1).velocity() = gaussianRandom() + scaleFactor;
    }
    for(int particleIndex = 0; particleIndex < nParticles; particleIndex++)
    {
        Oparticle->getData(particleIndex,2).velocity() = gaussianRandom() + scaleFactor;
    }
    for(int particleIndex = 0; particleIndex < nParticles; particleIndex++)
    {
        Oparticle->getData(particleIndex,3).velocity() = gaussianRandom() + scaleFactor;
    }
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::updateCellCoord()
{
    int cellNumber;
    real delZInv = 1./deltaZ;
    for(int particleIndex = 0; particleIndex < nParticles; particleIndex++)
    {
        cellNumber = (int)(Oparticle->getData(particleIndex,3).position()*delZInv);
        cList[cellNumber].cellParticleList.push_back(particleIndex);
    }
}

template<int DIM,typename real,typename dof>
real hardSphere<DIM,real,dof>::gaussianRandom() //Box-Mueller Transform
{
	real uniRand = drand48();
	real uniRand2 = drand48();

	return pow(-2.*log(uniRand),.5)*sin(2*pi*uniRand2);
}

template<int DIM,typename real,typename dof>
real hardSphere<DIM,real,dof>::relSpeed(int particleOne, int particleTwo)
{
    real relSp = 0.;
    real p1Vel;
    real p2Vel;
    for(int dim = 1; dim <= DIM; dim++)
    {
        p1Vel = Oparticle->getData(particleOne,dim).velocity();
        p2Vel = Oparticle->getData(particleTwo,dim).velocity();
        relSp += (p1Vel - p2Vel)*(p1Vel - p2Vel);
    }

    return pow(relSp,.5);
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::processCollision(int particleOne, int particleTwo)
{
    real p1Vel;
    real p2Vel;
    real p1VelNew;
    real p2VelNew;
    real vCOM[DIM];
    real vRel[DIM];
    real relSp = relSpeed(particleOne,particleTwo);

    for(int dim = 1; dim <= DIM; dim++)
    {
        p1Vel = Oparticle->getData(particleOne,dim).velocity();
        p2Vel = Oparticle->getData(particleTwo,dim).velocity();
        vCOM[dim - 1] = .5*(p1Vel + p2Vel);
    }

    real phiRand = drand48();
    real phi = 2.*pi*phiRand;

    real thetaRand = 2.*drand48() - 1.;
    real cosTheta = thetaRand;
    real sinTheta = pow(1. - cosTheta*cosTheta,.5);

    vRel[0] = relSp*sinTheta*cos(phi);
    vRel[1] = relSp*sinTheta*sin(phi);
    vRel[2] = relSp*cosTheta;

    for(int dim = 1; dim <= DIM; dim++)
    {
        p1VelNew = vCOM[dim - 1] + .5*vRel[dim -1];
        p2VelNew = vCOM[dim - 1] - .5*vRel[dim -1];

        Oparticle->getData(particleOne,dim).velocity() = p1VelNew;
        Oparticle->getData(particleTwo,dim).velocity() = p2VelNew;
    }
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::Collision()
{
    for(int currCell = 0; currCell < numCells; currCell++)
    {
        real vrMax = 2.;
        real selectCriteria;
        real relSp;
        int cellPartSize = cList[currCell].cellParticleList.size();
        real numCollisions = cellPartSize*(cellPartSize)*deltat*pi*diameter*diameter*vrMax*.5*numCells*100./volume;
        //100 effective particles per particle
        real uRandPairSelect[2];

        for(int collNo = 1; collNo <= (int)numCollisions; collNo++)
        {
            uRandPairSelect[0] = (int)(cellPartSize*drand48());
            uRandPairSelect[1] = (int)(cellPartSize*drand48());
            int particleOne = cList[currCell].cellParticleList[uRandPairSelect[0]];
            int particleTwo = cList[currCell].cellParticleList[uRandPairSelect[1]];
            relSp = relSpeed(particleOne,particleTwo);
            selectCriteria = drand48();

            if((relSp/vrMax) > selectCriteria)
            {
                processCollision(particleOne,particleTwo);
            }
        }
    }
}

double vCells[3][noCells];

template<int DIM,typename real, typename dof>
void hardSphere<DIM,real,dof>::postProcess(int it)
{
    real momentum = 0.;
    for(int currCell = 0; currCell < numCells; currCell++)
    {
        for(int dim = 1; dim <= DIM; dim++)
        {
            for(int particleIndex = 0; particleIndex < cList[currCell].cellParticleList.size(); particleIndex++)
            {
                int pIndex = cList[currCell].cellParticleList[particleIndex];
                vCells[dim - 1][currCell] += Oparticle->getData(pIndex,dim).velocity();
            }
        }
    }
}

template<int DIM,typename real, typename dof>
real hardSphere<DIM,real,dof>::getDensity(int cellNumber)
{
    return mass*(cList[cellNumber].cellParticleList.size())*numCells/volume;
}

int iterations;

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::hardSphereMain()
{
    real transKE;
    int count = 0;
    iterations = 0;

	while(iterations < 1)
	{

    		initPositions();
    		initVelocity();

		real initKE = totKinEnergy();
    		updateCellCoord();
		real totTime = 0.;
	    while(totTime <= 100000.)
	    {
		evolveSystem();

		Collision();

		transKE = totKinEnergy();   //Kinetic Energy at any time t

		cout<<"Ensemble: 5"<<endl;
		cout<<"iteration: "<<iterations;
		cout<<"\ntime: "<<totTime<<endl;
		cout<<"KE ratio: "<<transKE/initKE<<endl;

		if(count%10 == 0)   //updating cells every 10 steps
		{
		    clearCellParticles();
		    updateCellCoord();
		}

		totTime += deltat;
			count++;

		ofstream resultAv("./results/t100k/output_100k_n10k_bb10k_ens5.csv");

		postProcess(iterations);
		for(int currCell = 0; currCell < numCells; currCell++)
		{
			for(int dim = 1; dim <= DIM; dim++)
			{
				vCells[dim][currCell] = vCells[dim][currCell]/1.;
			}

			resultAv<<currCell<<"\t"<<vCells[1][currCell]<<endl;

			vCells[0][currCell] = 0.;
			vCells[1][currCell] = 0.;
			vCells[2][currCell] = 0.;		
		}
	    }
	    iterations++;
	}
}

#endif // HARDSPHERE_H_INCLUDED
