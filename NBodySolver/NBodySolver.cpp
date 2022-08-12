#include "NBodySolver.h"

NBodySolver::NBodySolver(NBodyData* data)
	:data(data)
{
}

NBodyData* NBodySolver::get_data() 
{
	return data;
}
