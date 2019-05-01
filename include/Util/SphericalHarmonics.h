/*************************************************
*	SphericalHarmonics.h
*
*	Release: August 2012
*	Update: February 2015
*
*	University of North Carolina at Chapel Hill
*	Department of Computer Science
*	
*	Ilwoo Lyu, ilwoolyu@cs.unc.edu
*************************************************/

#pragma once
#include "Geom.h"

class SphericalHarmonics
{
public:
	static void basis(int degree, float *p, float *Y);
};

