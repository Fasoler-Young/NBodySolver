#pragma once

#include "NBodyData.h"
#include <omp.h>

class NBodySolver
{
	NBodyData* data;
public:
	NBodySolver(NBodyData* data);
	NBodyData* get_data();
	virtual void step(value_type dt) = 0;
};


