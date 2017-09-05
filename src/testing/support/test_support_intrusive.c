// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 */

#include "test_support_intrusive.h"

#include "test_support_volume.h"
//#include "test_support_face.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Macros.h"

#include "volume.h"
#include "face.h"
#include "file_processing.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

struct Intrusive_List* constructor_file_name_IL
	(const char*const list_name, const char*const file_name, const struct const_Intrusive_List*const elements)
{
	FILE* file = fopen_checked(file_name); // closed

	struct Intrusive_List* intrusive_list = constructor_empty_IL();

	bool found_var = false;
	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),file)) {
		if (strstr(line,list_name)) {
			found_var = true;
			if (strstr(list_name,"Volume"))
				push_back_IL(intrusive_list,(struct Intrusive_Link*) constructor_Volume(file,line,elements));
//			else if (strstr(list_name,"Faces"))
//				push_back_IL(intrusive_list,(struct Intrusive_Link*) constructor_Face(file,line));
			else
				EXIT_UNSUPPORTED;
		}
	}

	fclose(file);

	if (!found_var)
		EXIT_ERROR("Did not find a '%s' member in the file: %s.",list_name,file_name);

/*for (const struct Intrusive_Link* curr = intrusive_list->first; curr; curr = curr->next) {
	struct Volume* volume = (struct Volume*) curr;
	printf("%d %d %d %d\n",volume->index,volume->boundary,volume->curved,volume->element->type);
}
EXIT_UNSUPPORTED;*/


	return intrusive_list;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
