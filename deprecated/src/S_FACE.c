// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "S_FACE.h"
#include "S_VOLUME.h"

#include <stddef.h>
#include <limits.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"

#include "element_functions.h"

/*
 *	Purpose:
 *		Provide functions associated with S_FACE struct.
 */

static inline unsigned int       compute_Indmf   (unsigned int const Indfh);
static inline unsigned int       compute_Indlfh  (unsigned int const Indfh);
static        unsigned int       compute_Indsfh  (struct S_SIDE_DATA const *const data, unsigned int const type);
static        void               set_type_FACE   (struct S_FACE *const FACE);


void update_data_FACE (struct S_FACE *const FACE)
{
	/*
	 *	Purpose:
	 *		Transfer data to the S_SIDE_DATA structs.
	 *
	 *	Comments:
	 *		This function should no longer be required when the new notation is adopted. (ToBeDeleted)
	 */

	size_t const Nsides = (FACE->Boundary ? 1 : 2 );
	for (size_t side = 0; side < Nsides; side++) {
		struct S_SIDE_DATA *const data = &FACE->data[side];

		if (side == 0) {
			data->Indfh  = FACE->VfL;
			data->VOLUME = FACE->VL;
		} else if (side == 1) {
			data->Indfh  = FACE->VfR;
			data->VOLUME = FACE->VR;
		}

		data->Indmf  = compute_Indmf(data->Indfh);
		data->Indlfh = compute_Indlfh(data->Indfh);
		if (side == 0)
			set_type_FACE(FACE);

		data->Indsfh = compute_Indsfh(data,FACE->type);
	}
}


// static functions

static inline unsigned int compute_Indmf (unsigned int const Indfh)
{
	return Indfh / NFREFMAX;
}

static inline unsigned int compute_Indlfh (unsigned int const Indfh)
{
	return Indfh % NFREFMAX;
}

static unsigned int compute_Indsfh (struct S_SIDE_DATA const *const data, unsigned int const type)
{
	/*
	 *	Notation:
	 *		Indlsfh : (Ind)ex of (l)ocal (s)ub-(f)ace potentially (h)-refined.
	 */

	unsigned int const Indlfh = data->Indlfh;

	unsigned int Indlsfh = UINT_MAX;
	if (Indlfh == 0) {
		Indlsfh = 0;
	} else {
		switch (type) {
		case POINT:
			// Do nothing
			break;
			Indlsfh = 0;
		case LINE:
			Indlsfh = Indlfh-1;
			break;
		case TRI:
			if (Indlfh <= 4) // Isotropic refinement (1 -> 4)
				Indlsfh = Indlfh-1;
			else
				EXIT_UNSUPPORTED;
			break;
		case QUAD:
			if (Indlfh <= 4) // Isotropic refinement (1 -> 4)
				Indlsfh = Indlfh-1;
			else if (Indlfh <= 6) // Horizontal refinement (1 -> 2)
				Indlsfh = Indlfh-5;
			else if (Indlfh <= 8) // Vertical refinement (1 -> 2)
				Indlsfh = Indlfh-7;
			else
				EXIT_UNSUPPORTED;
			break;
		default:
			EXIT_UNSUPPORTED;
			break;
		}
	}
	return (data->Indmf)*NSUBFMAX+Indlsfh;
}

static void set_type_FACE (struct S_FACE *const FACE)
{
	/*
	 *	Purpose:
	 *		Set the type of the FACE based on the neighbouring VOLUME.
	 */

	if (DB.d == 1) {
		FACE->type = POINT;
	} else if (DB.d == 2) {
		FACE->type = LINE;
	} else if (DB.d == 3) {
		struct S_ELEMENT const *const ELEMENT = get_ELEMENT_F_type(FACE->data[0].VOLUME->type,FACE->data[0].Indmf);
		FACE->type = ELEMENT->type;
	} else {
		EXIT_UNSUPPORTED;
	}
}
