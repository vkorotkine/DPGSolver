"""
Purpose:
	Generate information needed for h-refinement in the code.

Comments:
	The values selected for VToVe (Parent and Refined), FVertices, and EVertices are arbitrary in the sense that they
	correspond to specific figures for the h-refined versions of the elements (shown in
	DPG_ROOT/documentation/hrefinement_info_figures.pdf) This ordering does not correspond to that of the k = 2 geometry
	nodes used in the code (Note for the TET12 refinement that an additional node is included at the centroid as
	compared to the set of equi-spaced k = 2 geometry nodes).

Notation:
	d     : (d)imension
	Nv    : (N)umber of (v)olumes
	Nf    : (N)umber of (f)aces
	Ne    : (N)umber of (e)dges
	Nfve  : (N)umber of (f)ace (ve)rtices
	Neve  : (N)umber of (e)dge (ve)rtices

	ETypes    : (E)LEMENT      (Types): TRI, QUAD, TET, HEX, WEDGE, PYR
	VTypes    : (V)OLUME       (Types)
	RefTypes  : h-(Ref)inement (Types) (i.e. types of how the reference element will be split during h-refinement)
	FVertices : (F)ACE (Vertices) (i.e. sets of vertices located on each FACE)
	EVertices : (E)DGE (Vertices)

	VToVe : (V)OLUME (To) (Ve)rtex connectivity
	FToVe : (F)ACE   (To) (Ve)rtex connectivity
	EToVe : (E)DGE   (To) (Ve)rtex connectivity
	VToV  : (V)OLUME (To) (V)OLUME connectivity
	VToF  : (V)OLUME (To) (F)ACE   connectivity
	VToE  : (V)OLUME (To) (E)DGE   connectivity. Note that only the 'x' entries are used (i.e. the external EDGEs) as
	        the EDGE connectivity is not required in the code.
	FToP  : (F)ACE   (To) (P)ARENT connectivity
	EToP  : (E)DGE   (To) (P)ARENT connectivity

	IndOrdRL : (Ind)ex (Ord)ering (R)ight to (L)eft. See comments in compute_InOut.
	IndOrdLR : (Ind)ex (Ord)ering (L)eft to (R)ight. See comments in compute_InOut.
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
		if (name.find('LINE') != -1):
			if (name == 'LINE2'):
				Nv = 2
				VToVe  = [[0,1],[1,2]]
				VTypes = ['LINE' for i in range(Nv)]
		elif (name.find('TRI') != -1):
			if (name == 'TRI4'):
				Nv = 4
				VToVe  = [[0,1,3],[1,2,4],[3,4,5],[4,3,1]]
				VTypes = ['TRI' for i in range(Nv)]
		elif (name.find('QUAD') != -1):
			if (name == 'QUAD4'):
				Nv = 4
				VToVe = [[0,1,3,4],[1,2,4,5],[3,4,6,7],[4,5,7,8]]
				VTypes = ['QUAD' for i in range(Nv)]
			elif (name == 'QUAD2r'):
				Nv = 2
				VToVe = [[0,1,6,7],[1,2,7,8]]
				VTypes = ['QUAD' for i in range(Nv)]
			elif (name == 'QUAD2s'):
				Nv = 2
				VToVe = [[0,2,3,5],[3,5,6,8]]
				VTypes = ['QUAD' for i in range(Nv)]
		elif (name.find('TET') != -1):
			if (name == 'TET8'):
				Nv = 8
				VToVe = [[0,1,3,6],[1,2,4,7],[3,4,5,8],[6,7,8,9],[1,8,7,4],[8,1,6,3],[7,6,8,1],[4,3,1,8]]
				VTypes = ['TET' for i in range(Nv)]
			elif (name == 'TET12'):
				Nv = 12
				VToVe = [[0,1,3,6],[1,2,4,7],[3,4,5,8],[6,7,8,9],[10,8,7,4],[8,10,6,3],\
				         [7,6,10,1],[4,3,1,10],[10,1,6,3],[1,10,7,4],[4,3,10,8],[7,6,8,10]]
				VTypes = ['TET' for i in range(Nv)]
			elif (name == 'TET6'):
				Nv = 6
				VToVe = [[0,1,3,6],[1,2,4,7],[3,4,5,8],[6,7,8,9],[7,6,4,3,1],[6,7,3,4,8]]
				VTypes = ['TET' for i in range(4)]
				for i in range(2): VTypes.append('PYR')
		elif (name.find('HEX') != -1):
			if (name == 'HEX8'):
				Nv = 8
				VToVe = [[0,1, 3, 4, 9, 10,12,13],[1, 2, 4, 5, 10,11,13,14],[3, 4, 6, 7, 12,13,15,16],[4, 5, 7, 8, 13,14,16,17],\
				         [9,10,12,13,18,19,21,22],[10,11,13,14,19,20,22,23],[12,13,15,16,21,22,24,25],[13,14,16,17,22,23,25,26]]
				VTypes = ['HEX' for i in range(Nv)]
		elif (name.find('PYR') != -1):
			if (name == 'PYR10'):
				Nv = 10
				VToVe = [[0,1,3,4,9],[1,2,4,5,10],[3,4,6,7,11],[4,5,7,8,12],\
				         [3,4,11,9],[4,5,12,10],[10,9,4,1],[12,11,7,4],[10,9,12,11,4],[9,10,11,12,13]]
				VTypes = ['PYR' for i in range(4)]
				for i in range(4): VTypes.append('TET')
				for i in range(2): VTypes.append('PYR')
		elif (name.find('WEDGE') != -1):
			if (name == 'WEDGE8'):
				Nv = 8
				VToVe = [[0,1,3,6,7,9],[1,2,4,7,8,10],[3,4,5,9,10,11],[4,3,1,10,9,7],\
				         [6,7,9,12,13,15],[7,8,10,13,14,16],[9,10,11,15,16,17],[10,9,7,16,15,13]]
				VTypes = ['WEDGE' for i in range(Nv)]
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
		self.compute_EToVe(ELEMENTs)
		self.compute_VToX(ELEMENTs)
		self.compute_InOut(ELEMENTs)
		self.compute_FE_correspondence(ELEMENTs)

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

	def compute_EToVe(self,ELEMENTs):
		EToVe = []

		ELEMENT = get_ELEMENT_type(self.name,ELEMENTs)
		if (ELEMENT.d <= 2):
			self.EToVe = EToVe
			return

		for v in range(0,self.Nv):
			EToVe.append([])
			ELEMENT = get_ELEMENT_type(self.VTypes[v],ELEMENTs)
			for e in range(0,ELEMENT.Ne):
				EToVe[v].append([])
				for ve in range(0,ELEMENT.Neve):
					EToVe[v][e].append(self.VToVe[v][ELEMENT.EToVe[e][0][ve]])

		self.EToVe = EToVe

	def compute_VToX(self,ELEMENTs):
		FToVe = self.FToVe
		EToVe = self.EToVe

		VToV = []
		VToF = []

		# (V)OLUME (To) (V)OLUME (E)DGE is not unique here as EDGEs may be shared by multiple VOLUMEs, but this is not
		# important for its usage here.
		VToVE = []
		VToE  = []

		FoundE = []
		for v in range(0,self.Nv):
			ELEMENT = get_ELEMENT_type(self.VTypes[v],ELEMENTs)

			FoundE.append([])
			for e in range(0,ELEMENT.Ne):
				FoundE[v].append(0)

		for v in range(0,self.Nv):
			ELEMENT = get_ELEMENT_type(self.VTypes[v],ELEMENTs)

			VToV.append([])
			VToF.append([])
			for f in range(0,ELEMENT.Nf):
				VToV[v].append('x')
				VToF[v].append('x')
				for v2 in range(0,self.Nv):
					ELEMENT2 = get_ELEMENT_type(self.VTypes[v2],ELEMENTs)
					for f2 in range(0,ELEMENT2.Nf):
						if (not ([v,f] == [v2,f2]) and sorted(FToVe[v][f]) == sorted(FToVe[v2][f2])):
							VToV[v][f] = v2
							VToF[v][f] = f2

			VToVE.append([])
			VToE.append([])
			for e in range(0,ELEMENT.Ne):
				VToVE[v].append('x')
				VToE[v].append('x')
				for v2 in range(0,self.Nv):
					ELEMENT2 = get_ELEMENT_type(self.VTypes[v2],ELEMENTs)
					for e2 in range(0,ELEMENT2.Ne):
						if (not ([v,e] == [v2,e2]) and not FoundE[v][e] and \
						    sorted(EToVe[v][e]) == sorted(EToVe[v2][e2])):
							VToVE[v][e] = v2
							VToE[v][e]  = e2

							FoundE[v][e]   = 1
							FoundE[v2][e2] = 1

		EDGE_set = set()
		for v in range(0,self.Nv):
			ELEMENT = get_ELEMENT_type(self.VTypes[v],ELEMENTs)
			for f in range(0,ELEMENT.Nf):
				if (VToV[v][f] != 'x' and VToV[v][f] != 'xx'):
					v2 = VToV[v][f]
					f2 = VToF[v][f]

					VToV[v2][f2] = 'xx'
					VToF[v2][f2] = 'xx'

			# Note: Ne == 0 for (d <= 2)
			for e in range(0,ELEMENT.Ne):
				if (tuple(sorted(EToVe[v][e])) in EDGE_set):
					VToE[v][e]  = 'xx'
				else:
					EDGE_set.add(tuple(sorted(EToVe[v][e])))

		self.VToV = VToV
		self.VToF = VToF
		self.VToE = VToE

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
			if (Nve == 1): # POINT (d = 1)
				IndOrdInfo = data(RL = 0, LR = 0)
			elif (Nve == 2): # LINE (d == 2)
				if   (vals == [0,1]): IndOrdInfo = data(RL = 0, LR = 0)
				elif (vals == [1,0]): IndOrdInfo = data(RL = 1, LR = 1)
				else                : EXIT_TRACEBACK()
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
		VToV  = self.VToV
		VToF  = self.VToF

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
				if (VToV[v][f] == 'x'):
					FToVeOut[v][f].append('x')
					FToVeRef[v][f].append('x')
				elif (VToV[v][f] == 'xx'):
					FToVeOut[v][f].append('xx')
					FToVeRef[v][f].append('xx')
					IndOrdRL[v][f] = 'xx'
					IndOrdLR[v][f] = 'xx'
				else:
					# Store FToVeOut for printing below if desired
					v2 = VToV[v][f]
					f2 = VToF[v][f]
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

	def compute_FE_correspondence(self,ELEMENTs):

		# Parent ELEMENT
		ELEMENT = get_ELEMENT_type(self.name,ELEMENTs)

		Nf        = ELEMENT.Nf
		Ne        = ELEMENT.Ne
		FVertices = ELEMENT.FVertices
		EVertices = ELEMENT.EVertices


		VToF  = self.VToF
		VToE  = self.VToE
		FToVe = self.FToVe
		EToVe = self.EToVe

		# FACE correspondence
		FToP = []
		for v in range(self.Nv):
			FToP.append([])
			for f in range(len(VToF[v])):
				FToP[v].append([])
				if (VToF[v][f] != 'x'): # Internal FACEs
					FToP[v][f].append('x')
					continue

				# Find the parent FACE on which the child FACE lies
				for fP in range(Nf):
					Found = 1
					for ve in range(len(FToVe[v][f])):
						if (not FToVe[v][f][ve] in FVertices[fP]):
							Found = 0
							break;

					if (Found):
						break;

				if (not Found):
					list_print(self.VToV,"",2)
					list_print(VToF,"",2)
					print(FVertices)
					print(FToVe[v][f])
					print("Error: Did not find the associated parent FACE.\n")
					EXIT_TRACEBACK()

				FToP[v][f].append(fP)
				FToP[v][f].append('F')

		self.FToP = FToP

		if (ELEMENT.d == 2):
			# EDGEs == FACEs
			return

		# EDGE correspondence
		EToP = []
		for v in range(self.Nv):
			EToP.append([])
			for e in range(len(VToE[v])):
				EToP[v].append([])

				# If the EDGE lies on a parent EDGE, find the parent EDGE
				for eP in range(Ne):
					Found = 1
					for ve in range(len(EToVe[v][e])):
						if (not EToVe[v][e][ve] in EVertices[eP]):
							Found = 0
							break;

					if (Found):
						EToP[v][e].append(eP)
						EToP[v][e].append('E')
						break;

				# If the EDGE lies on a parent FACE, find the parent FACE
				if (not Found):
					for fP in range(Nf):
						Found = 1
						for ve in range(len(EToVe[v][e])):
							if (not EToVe[v][e][ve] in FVertices[fP]):
								Found = 0
								break;

						if (Found):
							EToP[v][e].append(fP)
							EToP[v][e].append('F')
							break;

				if (not Found):
					EToP[v][e].append('x')


		self.EToP = EToP


class ELEMENT_class:
	"""Class holding information relating to each element type."""
	def __init__(self,EType):
		self.EType = EType

		Neve = 2
		Nfve = []
		if (EType == 'LINE'):
			d         = 1
			Nf        = 2
			Ne        = 0
			Nfve      = [1 for i in range(0,Nf)]

			VToVe    = [[0,2]]
			RefTypes = ['LINE2']

			FVertices = [set([0]),set([2])]
			EVertices = []
		elif (EType == 'TRI'):
			d         = 2
			Nf        = 3
			Ne        = 0
			Nfve      = [2 for i in range(0,Nf)]

			VToVe    = [[0,2,5]]
			RefTypes = ['TRI4']

			FVertices = [set([2,4,5]),set([0,3,5]),set([0,1,2])]
			EVertices = []
		elif (EType == 'QUAD'):
			d         = 2
			Nf        = 4
			Ne        = 0
			Nfve      = [2 for i in range(0,Nf)]

			VToVe    = [[0,2,6,8]]
			RefTypes = ['QUAD4','QUAD2r','QUAD2s']

			FVertices = [set([0,3,6]),set([2,5,8]),set([0,1,2]),set([6,7,8])]
			EVertices = []
		elif (EType == 'TET'):
			d         = 3
			Nf        = 4
			Ne        = 6
			Nfve      = [3 for i in range(0,Nf)]

			VToVe    = [[0,2,5,9]]
			RefTypes = ['TET8','TET12','TET6']

			FVertices = [set([2,4,5,7,8,9]),set([0,3,5,6,8,9]),set([0,1,2,6,7,9]),set([0,1,2,3,4,5])]
			EVertices = [set([2,4,5]),set([0,3,5]),set([0,1,2]),set([0,6,9]),set([2,7,9]),set([5,8,9])]
		elif (EType == 'HEX'):
			d         = 3
			Nf        = 6
			Ne        = 12
			Nfve      = [4 for i in range(0,Nf)]

			VToVe    = [[0,2,6,8,18,20,24,26]]
			RefTypes = ['HEX8']

			FVertices = [set([0,3,6,9,12,15,18,21,24]), set([2,5,8,11,14,17,20,23,26]),set([0,1,2,9,10,11,18,19,20]),\
			             set([6,7,8,15,16,17,24,25,26]),set([0,1,2,3,4,5,6,7,8]),      set([18,19,20,21,22,23,24,25,26])]
			EVertices = [set([0,1,2]),   set([6,7,8]),   set([18,19,20]),set([24,25,26]),set([0,3,6]),  set([2,5,8]),\
			             set([18,21,24]),set([20,23,26]),set([0,9,18]),  set([2,11,20]), set([6,15,24]),set([8,17,26])]
		elif (EType == 'PYR'):
			d         = 3
			Nf        = 5
			Ne        = 8
			Nfve      = [3 for i in range(0,Nf-1)]
			for i in range(1): Nfve.append(4)

			VToVe    = [[0,2,6,8,13]]
			RefTypes = ['PYR10']

			FVertices = [set([0,3,6,9,11,13]),set([2,5,8,10,12,13]),set([0,1,2,9,10,13]),set([6,7,8,11,12,13]),\
			             set([0,1,2,3,4,5,6,7,8])]
			EVertices = [set([0,3,6]),set([2,5,8]),set([0,1,2]),set([6,7,8]),\
			             set([0,9,13]),set([2,10,13]),set([6,11,13]),set([8,12,13])]
		elif (EType == 'WEDGE'):
			d  = 3
			Nf = 5
			Ne = 9
			Nfve.extend((4 for i in range(0,3)))
			Nfve.extend((3 for i in range(0,2)))

			VToVe    = [[0,2,5,12,14,17]]
			RefTypes = ['WEDGE8']

			FVertices = [set([2,4,5,8,10,11,14,16,17]),set([0,3,5,6,9,11,12,15,17]),set([0,1,2,6,7,8,12,13,14]),\
			             set([0,1,2,3,4,5]),set([12,13,14,15,16,17])]
			EVertices = [set([2,4,5]),set([0,3,5]),set([0,1,2]),set([14,16,17]),set([12,15,17]),set([12,13,14]),\
			             set([0,6,12]),set([2,8,14]),set([5,11,17])]
		else:
			EXIT_TRACEBACK()

		self.d         = d
		self.Nf        = Nf
		self.Ne        = Ne
		self.Nfve      = Nfve
		self.Neve      = Neve
		self.VToVe     = VToVe
		self.FVertices = FVertices
		self.EVertices = EVertices

		self.set_FToVe()
		self.set_EToVe()
		self.initialize_RefTypes(RefTypes)

	def set_FToVe(self):
		Nf    = self.Nf
		Nfve  = self.Nfve

		if (self.EType == 'LINE'):
			FToVe = [[[0]],[[1]]]
		elif (self.EType == 'TRI'):
			FToVe = [[[1,2]],[[0,2]],[[0,1]]]
		elif(self.EType == 'QUAD'):
			FToVe = [[[0,2]],[[1,3]],[[0,1]],[[2,3]]]
		elif(self.EType == 'TET'):
			FToVe = [[[2,3,1]],[[2,3,0]],[[0,1,3]],[[0,1,2]]]
		elif(self.EType == 'HEX'):
			FToVe = [[[0,2,4,6]],[[1,3,5,7]],[[0,1,4,5]],[[2,3,6,7]],[[0,1,2,3]],[[4,5,6,7]]]
		elif(self.EType == 'PYR'):
			FToVe = [[[0,2,4]],[[1,3,4]],[[0,1,4]],[[2,3,4]],[[0,1,2,3]]]
		elif(self.EType == 'WEDGE'):
			FToVe = [[[1,2,4,5]],[[0,2,3,5]],[[0,1,3,4]],[[0,1,2]],[[3,4,5]]]
		else:
			EXIT_BACKTRACE()

		self.FToVe = FToVe

	def set_EToVe(self):
		Ne   = self.Ne
		Neve = self.Neve

		if (self.d <= 2):
			EToVe = []
		elif(self.EType == 'TET'):
			EToVe = [[[1,2]],[[0,2]],[[0,1]],[[0,3]],[[1,3]],[[2,3]]]
		elif(self.EType == 'HEX'):
			EToVe = [[[0,1]],[[2,3]],[[4,5]],[[6,7]],[[0,2]],[[1,3]],[[4,6]],[[5,7]],[[0,4]],[[1,5]],[[2,6]],[[3,7]]]
		elif(self.EType == 'PYR'):
			EToVe = [[[0,2]],[[1,3]],[[0,1]],[[2,3]],[[0,4]],[[1,4]],[[2,4]],[[3,4]]]
		elif(self.EType == 'WEDGE'):
			EToVe = [[[1,2]],[[0,2]],[[0,1]],[[4,5]],[[3,5]],[[3,4]],[[0,3]],[[1,4]],[[2,5]]]
		else:
			EXIT_BACKTRACE()

		self.EToVe = EToVe

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

""" Not used but potentially useful for another function in the future. (ToBeModified)
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
"""

def list_print(list_array,name,dim):
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
		if (EType.find(ELEMENT.EType) != -1):
			return ELEMENT

	print("Did not find the ELEMENT of type:",EType)
	EXIT_TRACEBACK()

def list_print_BC(list_array,names,dim):
	"""Print list_array in a format which can be copied directly into the c-code."""
	if (dim == 2):
		if (names[1] == "FToP"):
			Pnames = ['IndF','IndFP','IndBC']
		elif (names[1] == "EToP"):
			Pnames = ['IndE','IndEP','IndBC']
		else:
			EXIT_TRACEBACK()

		NEntries = 3
		Plists = ['' for i in range(NEntries)]

		count = 0
		item_count = 0
		for item in list_array:
			if item[0] != 'x':
				Plists[0] += Pnames[0]+'['+str(count)+']  = '+str(item_count)+'; '
				Plists[1] += Pnames[1]+'['+str(count)+'] = '+str(item[0])+'; '
				Plists[2] += Pnames[2]+'['+str(count)+'] = '+str(item[1])+'; '
				count += 1
			item_count += 1

		for i in range(NEntries):
			print(Plists[i])
		print("\n")
	elif (dim == 3):
		for i in range(len(list_array)):
			if (names[0] == "vh"):
				print(names[0],":",i+1)
			else:
				EXIT_TRACEBACK()
			list_print_BC(list_array[i],names,2)
	else:
		EXIT_TRACEBACK()


if __name__ == '__main__':
	"""Generate h refinement related information for each of the supported ELEMENT types in the code."""

	# Initialize all ELEMENT types (necessary as some refinements depend on several ELEMENTs)
	ELEMENTs = []
	EType = 'LINE';  ELEMENTs.append(ELEMENT_class(EType))
	EType = 'TRI';   ELEMENTs.append(ELEMENT_class(EType))
	EType = 'QUAD';  ELEMENTs.append(ELEMENT_class(EType))
	EType = 'TET';   ELEMENTs.append(ELEMENT_class(EType))
	EType = 'HEX';   ELEMENTs.append(ELEMENT_class(EType))
	EType = 'PYR';   ELEMENTs.append(ELEMENT_class(EType))
	EType = 'WEDGE'; ELEMENTs.append(ELEMENT_class(EType))

#	EType = 'LINE'
	EType = 'TRI'
#	EType = 'QUAD'
#	EType = 'TET'
#	EType = 'HEX'
#	EType = 'PYR'
#	EType = 'WEDGE'

	PrintAll = 0
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
			list_print(RefType.VToVe,"",2)

			print("FToVe:")
			list_print(RefType.FToVe,"v",3)

			print("VToV:")
			list_print(RefType.VToV,"",2)

			print("VToF:")
			list_print(RefType.VToF,"",2)

#			print("VToE:")
#			list_print(RefType.VToE,"",2)

			print("IndOrdRL :")
			list_print(RefType.IndOrdRL,"",2)

			print("IndOrdLR :")
			list_print(RefType.IndOrdLR,"",2)

			print("FToP :")
#			list_print(RefType.FToP,"v",3)
			list_print_BC(RefType.FToP,["vh","FToP"],3)

			if (ELEMENT.d == 3):
				print("EToP :")
#				list_print(RefType.EToP,"v",3)
				list_print_BC(RefType.EToP,["vh","EToP"],3)
