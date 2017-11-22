set	(SOURCE
	 adaptation.c
	 face_solver_adaptive.c
	 element_adaptation.c
	 volume_solver_adaptive.c
	)

set	(LIBS_DEPEND
	)

set (LIB Adaptation)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} PRIVATE ${LIBS_DEPEND})
