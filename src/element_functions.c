#include <stdlib.h>
#include <stdio.h>

#include "database.h"
#include "parameters.h"

/*
 *	Purpose:
 *		Simple element-related functions:
 *			int        is_ELEMENT_present(const int type);
 *			*S_ELEMENT get_ELEMENT_type(const int type);
 *			*S_ELEMENT get_ELEMENT_Eclass(const int Eclass, const int Esubclass);
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

int is_ELEMENT_present(const int type)
{
	struct S_ELEMENT *ELEMENT = DB.ELEMENT;

	while (ELEMENT != NULL) {
		if (type == ELEMENT->type) {
			if (ELEMENT->present)
				return 1;
			else
				return 0;
		}
		ELEMENT = ELEMENT->next;
	}
	printf("Error: Element type not found.\n"), exit(1);
}


struct S_ELEMENT *get_ELEMENT_type(const int type)
{
	struct S_ELEMENT *ELEMENT = DB.ELEMENT;

	while (ELEMENT != NULL) {
		if (type == ELEMENT->type)
			return ELEMENT;

		ELEMENT = ELEMENT->next;
	}

	printf("Error: Element type not found.\n"), exit(1);
}

struct S_ELEMENT *get_ELEMENT_Eclass(const int Eclass, const int Esubclass)
{
	struct S_ELEMENT *ELEMENT = DB.ELEMENT;

	if (Eclass == C_TP || (Eclass == C_WEDGE && Esubclass == C_TP)) {
		while (ELEMENT != NULL) {
			if (ELEMENT->type == LINE)
				return ELEMENT;

			ELEMENT = ELEMENT->next;
		}
	} else if (Eclass == C_SI || (Eclass == C_WEDGE && Esubclass == C_SI)) {
		while (ELEMENT != NULL) {
			if (ELEMENT->type == TRI)
				return ELEMENT;

			ELEMENT = ELEMENT->next;
		}
	} else if (Eclass == C_WEDGE) {
		printf("Error: Operators for WEDGE elements must be taken from TP or SI elements.\n"), exit(1);
	} else if (Eclass == C_PYR) {
		while (ELEMENT != NULL) {
			if (ELEMENT->type == PYR)
				return ELEMENT;

			ELEMENT = ELEMENT->next;
		}
	}
	printf("Error: Element class not found.\n"), exit(1);
}
