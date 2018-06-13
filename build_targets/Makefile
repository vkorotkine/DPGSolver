# Makefile used to add targets for helm in emacs.

default: all_debug

# Debug
.PHONY: all_debug fast_debug_1d fast_debug_2d fast_debug_3d debug_1d debug_2d debug_3d
all_debug:
	@echo "\nRunning make meshes && make for all debug build directories.\n"
	cd .. && ./make_debug_builds.sh

fast_debug_1d:
	@echo "\nRunning parallel make meshes && make for the 1d debug build directory.\n"
	cd ../build_debug_1d && make meshes && make -j

fast_debug_2d:
	@echo "\nRunning parallel make meshes && make for the 2d debug build directory.\n"
	cd ../build_debug_2d && make meshes && make -j

fast_debug_3d:
	@echo "\nRunning parallel make meshes && make for the 3d debug build directory.\n"
	cd ../build_debug_3d && make meshes && make -j

debug_1d:
	@echo "\nRunning make meshes && make for the 1d debug build directory.\n"
	cd ../build_debug_1d && make meshes && make -j1

debug_2d:
	@echo "\nRunning make meshes && make for the 2d debug build directory.\n"
	cd ../build_debug_2d && make meshes && make -j1

debug_3d:
	@echo "\nRunning make meshes && make for the 3d debug build directory.\n"
	cd ../build_debug_3d && make meshes && make -j1

# Release
.PHONY: all_release fast_release_1d fast_release_2d fast_release_3d release_1d release_2d release_3d
all_release:
	@echo "\nRunning make meshes && make for all release build directories.\n"
	cd .. && ./make_release_builds.sh

fast_release_1d:
	@echo "\nRunning parallel make meshes && make for the 1d release build directory.\n"
	cd ../build_release_1d && make meshes && make -j

fast_release_2d:
	@echo "\nRunning parallel make meshes && make for the 2d release build directory.\n"
	cd ../build_release_2d && make meshes && make -j

fast_release_3d:
	@echo "\nRunning parallel make meshes && make for the 3d release build directory.\n"
	cd ../build_release_3d && make meshes && make -j

release_1d:
	@echo "\nRunning make meshes && make for the 1d release build directory.\n"
	cd ../build_release_1d && make meshes && make -j1

release_2d:
	@echo "\nRunning make meshes && make for the 2d release build directory.\n"
	cd ../build_release_2d && make meshes && make -j1

release_3d:
	@echo "\nRunning make meshes && make for the 3d release build directory.\n"
	cd ../build_release_3d && make meshes && make -j1

# Documentation
.PHONY: doc
doc:
	@echo "\nRunning Doxygen in build_debug_3d.\n"
	cd ../build_debug_3d && make -j doc

# Execution
.PHONY: quick_debug_1d quick_debug_2d quick_debug_3d quick_1d quick_2d quick_3d memcheck_1d memcheck_2d memcheck_3d

quick_debug_1d: fast_debug_1d
	cd ../build_debug_1d/script_files && ./quick.sh

quick_debug_2d: fast_debug_2d
	cd ../build_debug_2d/script_files && ./quick.sh

quick_debug_3d: fast_debug_3d
	cd ../build_debug_3d/script_files && ./quick.sh

memcheck_debug_1d: fast_debug_1d
	cd ../build_debug_1d/script_files && ./memcheck.sh

memcheck_debug_2d: fast_debug_2d
	cd ../build_debug_2d/script_files && ./memcheck.sh

memcheck_debug_3d: fast_debug_3d
	cd ../build_debug_3d/script_files && ./memcheck.sh

quick_1d: fast_debug_1d
	cd ../build_1d/script_files && ./quick.sh

quick_2d: fast_debug_2d
	cd ../build_2d/script_files && ./quick.sh

quick_3d: fast_debug_3d
	cd ../build_3d/script_files && ./quick.sh
