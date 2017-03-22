"""
Purpose:
	Generate information needed for h-refinement in the code.

Comments:
	The values selected in set_VToVe are arbitrary in the sense that they correspond to specific figures for the
	h-refinement ELEMENTs which are shown in documentation/h_refinement_info_figures.pdf. This ordering does not
	correpond to that of the k = 2 geometry nodes used in the code.

Notation:
	NV    : (N)umber of (V)OLUMEs
	VToVe : (V)OLUME (To) (Ve)rtex array
"""

import sys
import numpy as np


### Classes ###
class Refinement_class:
	"""Class holding information relating to refinement options of each element type."""
	def __init__(self,name,ELEMENT):
		self.name    = name
		self.ELEMENT = ELEMENT

		self.FToVe = []

		NV = 0
		if (name.find('TRI') != -1):
			if (name == 'TRI4'):
				NV = 4
				VToVe = np.array([[0,3,4],[3,1,5],[4,5,2],[5,4,3]],dtype=np.int)
		elif (name.find('QUAD') != -1):
			if (name == 'QUAD4'):
				NV = 4
				VToVe = np.array([[0,4,5,6],[4,1,6,7],[5,6,2,8],[6,7,8,3]],dtype=np.int)
			elif (name == 'QUAD2r'):
				NV = 2
				VToVe = np.array([[0,4,2,8],[4,1,8,3]],dtype=np.int)
			elif (name == 'QUAD2s'):
				NV = 2
				VToVe = np.array([[0,1,5,7],[5,7,2,3]],dtype=np.int)
		else:
			EXIT_TRACEBACK()

		if (NV == 0):
			print("Did not find a refinement type for the given input.\n")
			EXIT_TRACEBACK()

		self.NV = NV
		self.VToVe = VToVe

		self.compute_FToVe()
		self.compute_EToE_EToF()

	def compute_FToVe(self):
		ELEMENT = self.ELEMENT

		FToVe = np.zeros((self.NV,max(ELEMENT.Nfve),ELEMENT.Nf),dtype=np.int)
		for vh in range(0,self.NV):
			for f in range(0,ELEMENT.Nf):
				FToVe[vh,:,f] = self.VToVe[vh,ELEMENT.FToVe[0,:,f]]

		self.FToVe = FToVe

	def compute_EToE_EToF(self):
		ELEMENT = self.ELEMENT

		FToVe = self.FToVe

		EToE = np.zeros((self.NV,ELEMENT.Nf),dtype=np.int)
		EToF = np.zeros((self.NV,ELEMENT.Nf),dtype=np.int)

		# Make list of faces with indices sorted from lowest to highest
		# Loop through the list searching for duplicates
		# If no duplicate is found, place a blank entry in EToE and EToF. If found, place correct index for first
		# occurence.
		FToVe_list = []
		for vh in range(0,self.NV):
			for f in range(0,ELEMENT.Nf):
				FToVe_list.append(FToVe[vh,:,f])
				print("vhf",vh,f)
				print(FToVe[vh,:,f])

		EXIT_TRACEBACK()
		



class ELEMENT_class:
	"""Class holding information relating to each element type."""
	def __init__(self,EType):
		self.EType = EType

		if (EType == 'TRI'):
			Nf        = 3
			Nfve      = [2 for i in range(0,Nf)]

			VToVe    = np.array([0,1,2],dtype=np.int)
			RefTypes = ['TRI4']
		elif (EType == 'QUAD'):
			Nf        = 4
			Nfve      = [2 for i in range(0,Nf)]

			VToVe    = np.array([0,1,2,3],dtype=np.int)
			RefTypes = ['QUAD4','QUAD2r','QUAD2s']
#			RefTypes = ['QUAD4']
		else:
			EXIT_TRACEBACK()

		self.Nf        = Nf
		self.Nfve      = Nfve

		self.VToVe     = VToVe

		self.set_FToVe()

		self.RefTypes = []
		for i in range(len(RefTypes)):
			RefTypeCurrent = Refinement_class(RefTypes[i],self)
			self.RefTypes.append(RefTypeCurrent)

	def set_FToVe(self):
		Nf    = self.Nf
		Nfve  = self.Nfve

		FToVe = np.zeros((1,max(Nfve),Nf),dtype=np.int)
		if (self.EType == 'TRI'):
			FToVe[0,:,0] = [1,2]
			FToVe[0,:,1] = [0,2]
			FToVe[0,:,2] = [0,1]
		elif(self.EType == 'QUAD'):
			FToVe[0,:,0] = [0,2]
			FToVe[0,:,1] = [1,3]
			FToVe[0,:,2] = [0,1]
			FToVe[0,:,3] = [2,3]
		else:
			EXIT_BACKTRACE()

		self.FToVe = FToVe



### Functions ###
PYTHON_ROOT = '../'
sys.path.insert(0,PYTHON_ROOT)
from support_functions import EXIT_TRACEBACK

def np_array_print(np_array):
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
	else:
		print("Add support.\n")
		EXIT_TRACEBACK()


if __name__ == '__main__':
	"""Generate h refinement related information for each of the supported ELEMENT types in the code."""

#	EType = 'TRI'
	EType = 'QUAD'

	ELEMENT = ELEMENT_class(EType)

	for RefType in range(len(ELEMENT.RefTypes)):
		print("VToVe:")
		np_array_print(ELEMENT.VToVe)
		np_array_print(ELEMENT.RefTypes[RefType].VToVe)

		print("FToVe:")
		for f in range(ELEMENT.Nf):
			print("f: ",f)
			np_array_print(ELEMENT.FToVe[:,:,f])
			np_array_print(ELEMENT.RefTypes[RefType].FToVe[:,:,f])
