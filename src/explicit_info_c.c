// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "explicit_info_c.h"
#include "Parameters.h"
#include "S_VOLUME.h"
#include "S_FACE.h"


#include "Macros.h"

#include "S_ELEMENT.h"

#include "element_functions.h"

/*
 *	Purpose:
 *		Provide support functions for explicit_VOLUME/FACE_info_c functions.
 *
 *	Comments:
 *		The pointer to the perturbed VOLUME and a boolean indicating whether all RHS_c contributions should be computed
 *		are passed to the complex functions. When applicable, this eliminates redundant computations of the RHS in
 *		VOLUMEs which are not neighbours of the perturbed VOLUME.
 *
 *		May need modification for HDG.
 */

struct S_LOCAL_MESH_ELEMENTS compute_local_ELEMENT_list (struct S_VOLUME const *const VOLUME, char const object)
{
	/*
	 *	Purpose:
	 *		Return the list of local VOLUMEs (i.e. those adjacent to and including the current VOLUME).
	 */

	struct S_ELEMENT const *const ELEMENT = get_ELEMENT_type(VOLUME->type);

	struct S_LOCAL_MESH_ELEMENTS local_ELEMENTs;
	unsigned int size = 0;

	switch (object) {
	case 'V':
		local_ELEMENTs.Vlist[size++] = (struct S_VOLUME *const) VOLUME;
		for (size_t mf = 0; mf < ELEMENT->Nf; mf++) {
			for (size_t sf = 0; sf < NSUBFMAX; sf++) {
				struct S_FACE const *const FACE = VOLUME->FACE[mf*NSUBFMAX+sf];

				if (!FACE || FACE->Boundary)
					break;

				struct S_VOLUME *VNeigh = (FACE->data[0].VOLUME == VOLUME ? FACE->data[1].VOLUME : FACE->data[0].VOLUME);

				local_ELEMENTs.Vlist[size++] = VNeigh;
			}
		}
		local_ELEMENTs.NVl = size;
		break;
	case 'F':
		for (size_t mf = 0; mf < ELEMENT->Nf; mf++) {
			for (size_t sf = 0; sf < NSUBFMAX; sf++) {
				struct S_FACE *const FACE = VOLUME->FACE[mf*NSUBFMAX+sf];

				if (!FACE)
					break;

				local_ELEMENTs.Flist[size++] = FACE;
			}
		}
		local_ELEMENTs.NFl = size;
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}

	return local_ELEMENTs;
}

bool is_VOLUME_in_local_list (struct S_VOLUME const *const VOLUME,
                              struct S_LOCAL_MESH_ELEMENTS const *const local_ELEMENTs)
{
	for (size_t i = 0; i < local_ELEMENTs->NVl; i++) {
		if (VOLUME == local_ELEMENTs->Vlist[i])
			return 1;
	}
	return 0;
}

bool is_FACE_in_local_list (struct S_FACE const *const FACE,
                            struct S_LOCAL_MESH_ELEMENTS const *const local_ELEMENTs)
{
	for (size_t i = 0; i < local_ELEMENTs->NFl; i++) {
		if (FACE == local_ELEMENTs->Flist[i])
			return 1;
	}
	return 0;
}
