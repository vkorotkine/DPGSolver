// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 */

#include "test_support_intrusive.h"

#include "test_support_volume.h"
#include "test_support_face.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "macros.h"
#include "constants_alloc.h"
#include "constants_intrusive.h"

#include "volume.h"
#include "face.h"
#include "file_processing.h"

// Static function declarations ************************************************************************************* //

/// Container for intrusive list related information.
struct IL_Info {
	int list_name; ///< \ref Intrusive_List::name.
};

/**	\brief Set the \ref IL_Info.
 *	\return A copy of the constructed \ref IL_Info. */
struct IL_Info set_IL_Info
	(const char*const list_name ///< The `char*` list name.
	);

// Interface functions ********************************************************************************************** //

struct Intrusive_List* constructor_file_name_IL
	(const char*const list_name, const char*const file_name, const struct const_Intrusive_List*const elements,
	 const struct Intrusive_List*const volumes)
{
	FILE* file = fopen_checked(file_name); // closed

	struct IL_Info il_info = set_IL_Info(list_name);

	struct Intrusive_List* intrusive_list = constructor_empty_IL(il_info.list_name);

	bool found_var = false;
	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),file)) {
		if (strstr(line,list_name)) {
			found_var = true;
			if (il_info.list_name == IL_VOLUME)
				push_back_IL(intrusive_list,(struct Intrusive_Link*) constructor_Volume(file,line,elements));
			else if (il_info.list_name == IL_FACE)
				push_back_IL(intrusive_list,(struct Intrusive_Link*) constructor_Face(file,line,elements,volumes));
			else
				EXIT_UNSUPPORTED;
		}
	}

	fclose(file);

	if (!found_var)
		EXIT_ERROR("Did not find a '%s' member in the file: %s.",list_name,file_name);

	return intrusive_list;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

struct IL_Info set_IL_Info (const char*const list_name)
{
	struct IL_Info il_info;

	if (strstr(list_name,"Volume")) {
		il_info.list_name = IL_VOLUME;
	} else if (strstr(list_name,"Face")) {
		il_info.list_name = IL_FACE;
	} else {
		EXIT_UNSUPPORTED;
	}

	return il_info;
}
