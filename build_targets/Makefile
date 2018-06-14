# Makefile used to add targets for helm in emacs.

default: all_debug

# Debug
.PHONY: all_debug fast_debug_1D fast_debug_2D fast_debug_3D debug_1D debug_2D debug_3D
all_debug:
	@echo "\nRunning make meshes && make for all debug build directories.\n"
	cd .. && ./make_debug_builds.sh

fast_debug_1D:
	@echo "\nRunning parallel make meshes && make for the 1D debug build directory.\n"
	cd ../build_debug_1D && make meshes && make -j

fast_debug_2D:
	@echo "\nRunning parallel make meshes && make for the 2D debug build directory.\n"
	cd ../build_debug_2D && make meshes && make -j

fast_debug_3D:
	@echo "\nRunning parallel make meshes && make for the 3D debug build directory.\n"
	cd ../build_debug_3D && make meshes && make -j

debug_1D:
	@echo "\nRunning make meshes && make for the 1D debug build directory.\n"
	cd ../build_debug_1D && make meshes && make -j1

debug_2D:
	@echo "\nRunning make meshes && make for the 2D debug build directory.\n"
	cd ../build_debug_2D && make meshes && make -j1

debug_3D:
	@echo "\nRunning make meshes && make for the 3D debug build directory.\n"
	cd ../build_debug_3D && make meshes && make -j1

# Release
.PHONY: all_release fast_release_1D fast_release_2D fast_release_3D release_1D release_2D release_3D
all_release:
	@echo "\nRunning make meshes && make for all release build directories.\n"
	cd .. && ./make_release_builds.sh

fast_release_1D:
	@echo "\nRunning parallel make meshes && make for the 1D release build directory.\n"
	cd ../build_release_1D && make meshes && make -j

fast_release_2D:
	@echo "\nRunning parallel make meshes && make for the 2D release build directory.\n"
	cd ../build_release_2D && make meshes && make -j

fast_release_3D:
	@echo "\nRunning parallel make meshes && make for the 3D release build directory.\n"
	cd ../build_release_3D && make meshes && make -j

release_1D:
	@echo "\nRunning make meshes && make for the 1D release build directory.\n"
	cd ../build_release_1D && make meshes && make -j1

release_2D:
	@echo "\nRunning make meshes && make for the 2D release build directory.\n"
	cd ../build_release_2D && make meshes && make -j1

release_3D:
	@echo "\nRunning make meshes && make for the 3D release build directory.\n"
	cd ../build_release_3D && make meshes && make -j1

# Documentation
.PHONY: doc
doc:
	@echo "\nRunning Doxygen in build_debug_3D.\n"
	cd ../build_debug_3D && make -j doc

# Execution
.PHONY: quick_debug_1D quick_debug_2D quick_debug_3D quick_1D quick_2D quick_3D memcheck_1D memcheck_2D memcheck_3D

quick_debug_1D: fast_debug_1D
	cd ../build_debug_1D/script_files && ./quick.sh

quick_debug_2D: fast_debug_2D
	cd ../build_debug_2D/script_files && ./quick.sh

quick_debug_3D: fast_debug_3D
	cd ../build_debug_3D/script_files && ./quick.sh

memcheck_debug_1D: fast_debug_1D
	cd ../build_debug_1D/script_files && ./memcheck.sh

memcheck_debug_2D: fast_debug_2D
	cd ../build_debug_2D/script_files && ./memcheck.sh

memcheck_debug_3D: fast_debug_3D
	cd ../build_debug_3D/script_files && ./memcheck.sh

quick_1D: fast_debug_1D
	cd ../build_1D/script_files && ./quick.sh

quick_2D: fast_debug_2D
	cd ../build_2D/script_files && ./quick.sh

quick_3D: fast_debug_3D
	cd ../build_3D/script_files && ./quick.sh
