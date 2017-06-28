// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "explicit_info_c.h"
#include "Parameters.h"
#include "S_VOLUME.h"

#include "S_ELEMENT.h"
#include "S_FACE.h"

#include "element_functions.h"

/*
 *	Purpose:
 *		Provide support functions for explicit_VOLUME/FACE_info_c functions.
 *
 *	Comments:
 *		VOLUME is passed to the functions such that only the necessary RHS terms are evaluated when applicable. When
 *		compute_all is passed as TRUE, the RHS of every VOLUME is computed. May need modification for HDG.
 */

struct S_LOCAL_VOLUMES compute_local_VOLUME_list (struct S_VOLUME const *const VOLUME)
{
	/*
	 *	Purpose:
	 *		Return the list of local VOLUMEs (i.e. those adjacent to and including the current VOLUME).
	 */

	struct S_ELEMENT const *const ELEMENT = get_ELEMENT_type(VOLUME->type);

	struct S_LOCAL_VOLUMES local_VOLUMEs;

	unsigned int size = 0;

	local_VOLUMEs.Vlist[size++] = (struct S_VOLUME *const) VOLUME;
	for (size_t mf = 0; mf < ELEMENT->Nf; mf++) {
		for (size_t sf = 0; sf < NSUBFMAX; sf++) {
			struct S_FACE const *const FACE = VOLUME->FACE[mf*NSUBFMAX+sf];

			if (!FACE || FACE->Boundary)
				break;

			struct S_VOLUME *VNeigh = ( FACE->data[0].VOLUME == VOLUME ? FACE->data[1].VOLUME : FACE->data[0].VOLUME );

			local_VOLUMEs.Vlist[size++] = VNeigh;
		}
	}
	local_VOLUMEs.size = size;

	return local_VOLUMEs;
}

bool is_VOLUME_in_local_list (struct S_VOLUME const *const VOLUME, struct S_LOCAL_VOLUMES const *const local_VOLUMEs)
{
	for (size_t i = 0; i < local_VOLUMEs->size; i++) {
		if (VOLUME == local_VOLUMEs->Vlist[i])
			return 1;
	}
	return 0;
}
