// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "select_functions.h"

#include <stdlib.h>
#include <stdio.h>

#include "Parameters.h"
#include "Macros.h"

#include "bases.h"
#include "cubature.h"

/*
 *	Purpose:
 *		Select basis, grad_basis and cubature functions based on input type.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void select_functions(basis_tdef *basis, grad_basis_tdef *grad_basis, cubature_tdef *cubature, const unsigned int type)
{
	switch(type) {
	case LINE:
	case QUAD:
	case HEX:
		*basis      = basis_TP;
		*grad_basis = grad_basis_TP;
		*cubature   = cubature_TP;
		break;
	case TRI:
		*basis      = basis_SI;
		*grad_basis = grad_basis_SI;
		*cubature   = cubature_TRI;
		break;
	case TET:
		*basis      = basis_SI;
		*grad_basis = grad_basis_SI;
		*cubature   = cubature_TET;
		break;
	case WEDGE:
		printf("Error: WEDGE elements use a combination of TRI and LINE basis functions/nodes.\n"), EXIT_MSG;
		break;
	case PYR:
		*basis      = basis_PYR;
		*grad_basis = grad_basis_PYR;
		*cubature   = cubature_PYR;
		break;
	default:
		printf("Error: Unsupported type = %d.\n",type), EXIT_MSG;
		break;
	}
}

void select_functions_basis(basis_tdef *basis, const unsigned int type)
{
	switch(type) {
	case LINE:
	case QUAD:
	case HEX:
		*basis      = basis_TP;
		break;
	case TRI:
	case TET:
		*basis      = basis_SI;
		break;
	case WEDGE:
		printf("Error: WEDGE elements use a combination of TRI and LINE basis functions/nodes.\n"), EXIT_MSG;
		break;
	case PYR:
		*basis      = basis_PYR;
		break;
	default:
		printf("Error: Unsupported type = %d.\n",type), EXIT_MSG;
		break;
	}
}

void select_functions_grad_basis(grad_basis_tdef *grad_basis, const unsigned int type)
{
	switch(type) {
	case LINE:
	case QUAD:
	case HEX:
		*grad_basis = grad_basis_TP;
		break;
	case TRI:
	case TET:
		*grad_basis = grad_basis_SI;
		break;
	case WEDGE:
		printf("Error: WEDGE elements use a combination of TRI and LINE basis functions/nodes.\n"), EXIT_MSG;
		break;
	case PYR:
		*grad_basis = grad_basis_PYR;
		break;
	default:
		printf("Error: Unsupported type = %d.\n",type), EXIT_MSG;
		break;
	}
}

void select_functions_cubature(cubature_tdef *cubature, const unsigned int type)
{
	switch(type) {
	case LINE:
	case QUAD:
	case HEX:
		*cubature   = cubature_TP;
		break;
	case TRI:
		*cubature   = cubature_TRI;
		break;
	case TET:
		*cubature   = cubature_TET;
		break;
	case WEDGE:
		printf("Error: WEDGE elements use a combination of TRI and LINE basis functions/nodes.\n"), EXIT_MSG;
		break;
	case PYR:
		*cubature   = cubature_PYR;
		break;
	default:
		printf("Error: Unsupported type = %d.\n",type), EXIT_MSG;
		break;
	}
}
