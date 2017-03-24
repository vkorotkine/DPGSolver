"""
Purpose:
	Generate information needed for h-refinement in the code.

Comments:
	The values selected in set_VToVe are arbitrary in the sense that they correspond to specific figures for the
	h-refinement ELEMENTs which are shown in documentation/h_refinement_info_figures.pdf. This ordering does not
	correpond to that of the k = 2 geometry nodes used in the code.

Notation:
	Nv    : (N)umber of (v)olumes
	VToVe : (V)OLUME (To) (Ve)rtex array
"""

import sys
import numpy as np
import collections

### Classes ###
class Refinement_class:
	"""Class holding information relating to refinement options of each element type."""

	def __init__(self,name):
		self.name = name

		Nv = 0
		if (name.find('TRI') != -1):
			if (name == 'TRI4'):
				Nv = 4
				VToVe  = [[0,3,4],[3,1,5],[4,5,2],[5,4,3]]
				VTypes = ['TRI' for i in range(Nv)]
		elif (name.find('QUAD') != -1):
			if (name == 'QUAD4'):
				Nv = 4
				VToVe = [[0,4,5,6],[4,1,6,7],[5,6,2,8],[6,7,8,3]]
				VTypes = ['QUAD' for i in range(Nv)]
			elif (name == 'QUAD2r'):
				Nv = 2
				VToVe = [[0,4,2,8],[4,1,8,3]]
				VTypes = ['QUAD' for i in range(Nv)]
			elif (name == 'QUAD2s'):
				Nv = 2
				VToVe = [[0,1,5,7],[5,7,2,3]]
				VTypes = ['QUAD' for i in range(Nv)]
		elif (name.find('TET') != -1):
			if (name == 'TET8'):
				Nv = 8
				VToVe = [[0,4,5,7],[4,1,6,8],[5,6,2,9],[7,8,9,3],[4,9,8,6],[9,4,7,5],[8,7,9,4],[6,5,4,9]]
				VTypes = ['TET' for i in range(Nv)]
			elif (name == 'TET12'):
				Nv = 12
				VToVe = [[0,4,5,7],[4,1,6,8],[5,6,2,9],[7,8,9,3],[10,9,8,6],[9,10,7,5],\
				         [8,7,10,4],[6,5,4,10],[10,4,7,5],[4,10,8,6],[6,5,10,9],[8,7,9,10]]
				VTypes = ['TET' for i in range(Nv)]
			elif (name == 'TET6'):
				Nv = 6
				VToVe = [[0,4,5,7],[4,1,6,8],[5,6,2,9],[7,8,9,3],[8,7,6,5,4],[7,8,5,6,9]]
				VTypes = ['TET' for i in range(4)]
				for i in range(2): VTypes.append('PYR')
		elif (name.find('PYR') != -1):
			if (name == 'PYR10'):
				Nv = 10
				VToVe = [[0,5,6,7,10],[5,1,7,8,11],[6,7,2,9,12],[7,8,9,3,13],\
				         [6,7,12,10],[7,8,13,11],[11,10,7,5],[13,12,9,7],[11,10,13,12,7],[10,11,12,13,4]]
				VTypes = ['PYR' for i in range(4)]
				for i in range(4): VTypes.append('TET')
				for i in range(2): VTypes.append('PYR')
		else:
			print(name)
			EXIT_TRACEBACK()

		if (Nv == 0):
			print("Did not find a refinement type for the given input.\n")
			EXIT_TRACEBACK()

		self.Nv     = Nv
		self.VToVe  = VToVe
		self.VTypes = VTypes


	def compute_info(self,ELEMENTs):
		self.compute_FToVe(ELEMENTs)
		self.compute_EToE_EToF(ELEMENTs)
		self.compute_InOut(ELEMENTs)

	def compute_FToVe(self,ELEMENTs):
		FToVe = []
		for v in range(0,self.Nv):
			FToVe.append([])
			ELEMENT = get_ELEMENT_type(self.VTypes[v],ELEMENTs)
			for f in range(0,ELEMENT.Nf):
				FToVe[v].append([])
				for ve in range(0,ELEMENT.Nfve[f]):
					FToVe[v][f].append(self.VToVe[v][ELEMENT.FToVe[f][0][ve]])

		self.FToVe = FToVe

	def compute_EToE_EToF(self,ELEMENTs):
		FToVe = self.FToVe

		EToE = []
		EToF = []

		Found = []
		for v in range(0,self.Nv):
			Found.append([])
			ELEMENT = get_ELEMENT_type(self.VTypes[v],ELEMENTs)
			for f in range(0,ELEMENT.Nf):
				Found[v].append(0)

		for v in range(0,self.Nv):
			EToE.append([])
			EToF.append([])
			ELEMENT = get_ELEMENT_type(self.VTypes[v],ELEMENTs)
			for f in range(0,ELEMENT.Nf):
				EToE[v].append('x')
				EToF[v].append('x')
				for v2 in range(0,self.Nv):
					ELEMENT2 = get_ELEMENT_type(self.VTypes[v2],ELEMENTs)
					for f2 in range(0,ELEMENT2.Nf):
						if (not ([v,f] == [v2,f2]) and not Found[v][f] and  \
						    sorted(FToVe[v][f]) == sorted(FToVe[v2][f2])):
							EToE[v][f] = v2
							EToF[v][f] = f2

							Found[v][f]   = 1
							Found[v2][f2] = 1

		for v in range(0,self.Nv):
			ELEMENT = get_ELEMENT_type(self.VTypes[v],ELEMENTs)
			for f in range(0,ELEMENT.Nf):
				if (EToE[v][f] != 'x' and EToE[v][f] != 'xx'):
					v2 = EToE[v][f]
					f2 = EToF[v][f]

					EToE[v2][f2] = 'xx'
					EToF[v2][f2] = 'xx' 


		self.EToE = EToE
		self.EToF = EToF

	def compute_InOut(self,ELEMENTs):
		def get_IndOrd_Indices(Nve,vals):
			"""Return IndOrd indices for RL and LR. The convention used here corresponds to the cases in d = 2, 3
			   expected outputs of DPG_ROOT/src/test_unit_get_face_ordering.c (P = 1). The notation may be different
			   here as compared to the code (L == In, R == Out, e.g. LR == InOut).
			   
			   The concept behind the convention is:
			   	Given an index for each configuration of vals, find the index of the configuration as seen by the
				neighboring element.

			   Examples:
				A TRI FACE on a 3D TET is oriented such that its vertices have indices [0,2,1] (RL = 3 below) with
				respect to the reference TET FACE (i.e. it has the same position of vertex 0, but is oriented in the
				opposite direction).

				From the point of view of the neighboring TET, to get from [0,2,1] to the reference [0,1,2] requires an
				ordering of [0,2,1] (LR = 3)

				Had the ordering of the FACE with respect to the reference been [1,2,0] (RL = 1), to get back to the
				reference [0,1,2] would required an ordering of [2,0,1] (LR = 2).


			   Why this index is necessary is explained below.

			   In a 2D code, elements from the mesh file are generally (this is necessary in the code) defined such that
			   their orientations are all in the same direction (e.g. counter-clockwise for TRI elements). VOLUMEs with
			   incorrect orientation are "flipped-over" and will be seen in the code as having negative volume. If this
			   convention is maintained throughout the code, it would be known a priori that the FACEs of adjacent
			   ELEMENTs are ALWAYS oriented in the opposite direction to each other (i.e. RL == LR == 1) and this
			   information would be unnecessary. HOWEVER, the convention in the code changes this FACE orientation in
			   the reference TRI and QUAD elements so that this is no longer true in the present case.

			   Further, in 3D, special cases arise where the ordering cannot be simply defined as reversed (TRI: RL =
			   1,2; QUAD: RL = 5,6) and these indices become necessary to know the correct orientation of FACE data with
			   respect to each of the neighboring VOLUMEs.

			   Thus, in the interest of code generality, this information is used in 2D as well and we do not impose a
			   restriction on the orientation of the FACEs of the reference ELEMENTs.
			"""

			data = collections.namedtuple('data','RL LR')
			if (Nve == 2): # LINE (d == 2)
				if   (vals == [0,1]): IndOrdInfo = data(RL = 0, LR = 0)
				elif (vals == [1,0]): IndOrdInfo = data(RL = 1, LR = 1)
				else                        : EXIT_TRACEBACK()
			elif (Nve == 3): # TRI (d == 3)
				if   (vals == [0,1,2]): IndOrdInfo = data(RL = 0, LR = 0)
				elif (vals == [1,2,0]): IndOrdInfo = data(RL = 1, LR = 2)
				elif (vals == [2,0,1]): IndOrdInfo = data(RL = 2, LR = 1)
				elif (vals == [0,2,1]): IndOrdInfo = data(RL = 3, LR = 3)
				elif (vals == [2,1,0]): IndOrdInfo = data(RL = 4, LR = 4)
				elif (vals == [1,0,2]): IndOrdInfo = data(RL = 5, LR = 5)
			elif (Nve == 4): # QUAD (d == 3)
				if   (vals == [0,1,2,3]): IndOrdInfo = data(RL = 0, LR = 0)
				elif (vals == [1,0,3,2]): IndOrdInfo = data(RL = 1, LR = 1)
				elif (vals == [2,3,0,1]): IndOrdInfo = data(RL = 2, LR = 2)
				elif (vals == [3,2,1,0]): IndOrdInfo = data(RL = 3, LR = 3)
				elif (vals == [0,2,1,3]): IndOrdInfo = data(RL = 4, LR = 4)
				elif (vals == [2,0,3,1]): IndOrdInfo = data(RL = 5, LR = 6)
				elif (vals == [1,3,0,2]): IndOrdInfo = data(RL = 6, LR = 5)
				elif (vals == [3,1,2,0]): IndOrdInfo = data(RL = 7, LR = 7)
			else:
				EXIT_TRACEBACK()

			return IndOrdInfo

		PrintEnabled = 0

		FToVe = self.FToVe
		EToE  = self.EToE
		EToF  = self.EToF

		FToVeOut = []
		FToVeRef = []
		IndOrdRL = []
		IndOrdLR = []
		for v in range(0,self.Nv):
			FToVeOut.append([])
			FToVeRef.append([])
			IndOrdRL.append([])
			IndOrdLR.append([])

			ELEMENT = get_ELEMENT_type(self.VTypes[v],ELEMENTs)
			for f in range(0,ELEMENT.Nf):
				FToVeOut[v].append([])
				FToVeRef[v].append([])
				IndOrdRL[v].append('x')
				IndOrdLR[v].append('x')
				if (EToE[v][f] == 'x'):
					FToVeOut[v][f].append('x')
					FToVeRef[v][f].append('x')
				elif (EToE[v][f] == 'xx'):
					FToVeOut[v][f].append('xx')
					FToVeRef[v][f].append('xx')
					IndOrdRL[v][f] = 'xx'
					IndOrdLR[v][f] = 'xx'
				else:
					# Store FToVeOut for printing below if desired
					v2 = EToE[v][f]
					f2 = EToF[v][f]
					FToVeOut[v][f].append(FToVe[v2][f2])

					# Find indices relating to IndOrd
					len_IndOrd = len(FToVe[v][f])

					val_IndOrd = []
					for i in range(len_IndOrd):
						for j in range(len_IndOrd):
							if (FToVe[v][f][i] == FToVe[v2][f2][j]):
								val_IndOrd.append(j)

					# Store FToVeRef for printing below if desired
					FToVeRef[v][f].append(val_IndOrd)

					IndOrdReturn = get_IndOrd_Indices(len_IndOrd,val_IndOrd)

					IndOrdRL[v][f] = IndOrdReturn.RL
					IndOrdLR[v][f] = IndOrdReturn.LR

		if (PrintEnabled):
			list_print(FToVeOut,"v",3)
			list_print(FToVeRef,"v",3)

		self.IndOrdRL = IndOrdRL
		self.IndOrdLR = IndOrdLR


class ELEMENT_class:
	"""Class holding information relating to each element type."""
	def __init__(self,EType):
		self.EType = EType

		if (EType == 'TRI'):
			Nf        = 3
			Nfve      = [2 for i in range(0,Nf)]

			VToVe    = [[0,1,2]]
			RefTypes = ['TRI4']
		elif (EType == 'QUAD'):
			Nf        = 4
			Nfve      = [2 for i in range(0,Nf)]

			VToVe    = [[0,1,2,3]]
			RefTypes = ['QUAD4','QUAD2r','QUAD2s']
		elif (EType == 'TET'):
			Nf        = 4
			Nfve      = [3 for i in range(0,Nf)]

			VToVe    = [[0,1,2,3]]
			RefTypes = ['TET8','TET12','TET6']
		elif (EType == 'PYR'):
			Nf        = 5
			Nfve      = [3 for i in range(0,Nf-1)]
			for i in range(1): Nfve.append(4)

			VToVe    = [[0,1,2,3,4]]
			RefTypes = ['PYR10']
		else:
			EXIT_TRACEBACK()

		self.Nf    = Nf
		self.Nfve  = Nfve
		self.VToVe = VToVe

		self.set_FToVe()
		self.initialize_RefTypes(RefTypes)

	def set_FToVe(self):
		Nf    = self.Nf
		Nfve  = self.Nfve

		if (self.EType == 'TRI'):
			FToVe = [[[1,2]],[[0,2]],[[0,1]]]
		elif(self.EType == 'QUAD'):
			FToVe = [[[0,2]],[[1,3]],[[0,1]],[[2,3]]]
		elif(self.EType == 'TET'):
			FToVe = [[[2,3,1]],[[2,3,0]],[[0,1,3]],[[0,1,2]]]
		elif(self.EType == 'PYR'):
			FToVe = [[[0,2,4]],[[1,3,4]],[[0,1,4]],[[2,3,4]],[[0,1,2,3]]]
		else:
			EXIT_BACKTRACE()

		self.FToVe = FToVe

	def initialize_RefTypes(self,RefTypes):
		self.RefTypes = []
		for i in range(len(RefTypes)):
			RefTypeCurrent = Refinement_class(RefTypes[i])
			self.RefTypes.append(RefTypeCurrent)

	def set_RefTypes(self,ELEMENTs):
		for i in range(len(self.RefTypes)):
			self.RefTypes[i].compute_info(ELEMENTs)



### Functions ###
PYTHON_ROOT = '../'
sys.path.insert(0,PYTHON_ROOT)
from support_functions import EXIT_TRACEBACK

def np_array_print(np_array,name):
	ndim = np_array.ndim

	if (ndim == 1 or ndim == 2):
		if (ndim == 1):
			s = [[str(e) for e in row] for row in [np_array]]
		else:
			s = [[str(e) for e in row] for row in np_array]
		lens  = [max(map(len, col)) for col in zip(*s)]
		fmt   = '\t'.join('{{:{}}}'.format(x) for x in lens)
		table = [fmt.format(*row) for row in s]

		print('\n'.join(table),'\n')
	elif (ndim == 3):
		for i in range(np_array.shape[ndim-1]):
			print(name,":",i)
			np_array_print(np_array[:,:,i],"")
	else:
		print("Add support.\n")
		EXIT_TRACEBACK()

def list_print(list_array,name,dim):
	"""
	if (dim == 2):
		for item in list_array:
			print('\t'.join(map(str,item)))
		print("\n")
	elif (dim == 3):
		print(list_array)
		print(len(list_array))
		for i in range(len(list_array)):
			print(name,":",i)
			list_print(list_array[i],"",2)
	"""
	if (dim == 2):
		for item in list_array:
			print('\t'.join(map(str,item)))
		print("\n")
	elif (dim == 3):
		for i in range(len(list_array)):
			print(name,":",i)
			list_print(list_array[i],"",2)
	else:
		EXIT_TRACEBACK()

def get_ELEMENT_type(EType,ELEMENTs):
	for ELEMENT in ELEMENTs:
		if (ELEMENT.EType == EType):
			return ELEMENT

	print("Did not find the ELEMENT of type:",EType) 
	EXIT_TRACEBACK()


if __name__ == '__main__':
	"""Generate h refinement related information for each of the supported ELEMENT types in the code."""

	# Initialize all ELEMENT types (necessary as some refinements depend on several ELEMENTs)
	ELEMENTs = []
	EType = 'TRI';   ELEMENTs.append(ELEMENT_class(EType))
	EType = 'QUAD';  ELEMENTs.append(ELEMENT_class(EType))
	EType = 'TET';   ELEMENTs.append(ELEMENT_class(EType))
	EType = 'PYR';   ELEMENTs.append(ELEMENT_class(EType))

#	EType = 'TRI'
#	EType = 'QUAD'
	EType = 'TET'

	PrintAll = 1
	PrintInd = 0
	for ELEMENT in ELEMENTs:
		if (ELEMENT.EType != EType):
			continue
		ELEMENT.set_RefTypes(ELEMENTs)

		for i in range(len(ELEMENT.RefTypes)):
			if (not PrintAll and i != PrintInd):
				continue

			RefType = ELEMENT.RefTypes[i]
			print("\n\nh refinement information for RefType:",RefType.name,"\n\n")

			print("VToVe:")
#			list_print(ELEMENT.VToVe,"",2)
			list_print(RefType.VToVe,"",2)

			print("FToVe:")
#			list_print(ELEMENT.FToVe,"v",3)
			list_print(RefType.FToVe,"v",3)

			print("EToE:")
			list_print(RefType.EToE,"",2)

			print("EToF:")
			list_print(RefType.EToF,"",2)

			print("IndOrdRL :")
			list_print(RefType.IndOrdRL,"",2)

			print("IndOrdLR :")
			list_print(RefType.IndOrdLR,"",2)

			print("\n\n\n***Update the Notation.***\n\n\n")
