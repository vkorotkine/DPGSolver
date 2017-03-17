# Makefile for mesh generation

PYTHONC := python3

# run 'python3 -c "import generate_meshes.py"' in terminal and check filename in __pycache__ to check PYTHON_BINARY
PYTHON_BINARY := cpython-35
PYCACHE       := __pycache__

SOURCE      := generate_meshes.py
SOURCE_NAME := $(SOURCE:%.py=%)
COMPILED := $(PYCACHE)/$(SOURCE_NAME).$(PYTHON_BINARY).pyc


SHELL := /bin/bash

# Default goal (Generates compiled python code)
$(COMPILED) : $(SOURCE)
	@echo
	@echo Creating compiled code for: $^
	@$(PYTHONC) -c $(^:%.py="import %")
	@echo

.PHONY : clean
clean:
	rm $(COMPILED)

.PHONY : meshes_all


# .geo files
DCUBE_GEO := $(wildcard dCube/*.geo)


meshes_all:
	@echo Testing
	@echo $(DCUBE_GEO)

UPDATE_H := dCube2D_TRI0x.msh dCube2D_QUAD0x.msh
#                  dCube3D_TET0x.msh
UPDATE_H := $(addprefix dCube/Poisson/Dirichlet/Test/,$(UPDATE_H))

.PHONY : meshes_update_h
meshes_update_h : $(UPDATE_H)
$(UPDATE_H) : $(DCUBE_GEO)
	@$(PYTHONC) ./$(COMPILED) update_h
	

