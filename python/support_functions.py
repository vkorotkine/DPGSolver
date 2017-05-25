import sys
import traceback

### Special ###
def EXIT_TRACEBACK():
	print('\nError: Unsupported')
	traceback.print_stack()
	sys.exit()

### File Input/Output Related ###
def f_write(f,Ntabs,Nnewlines,string):
	"""Write string to file with specified number of tabs and newlines."""

	output = ''
	for n in range(0,Ntabs):
		output += '\t'
	output += string
	for n in range(0,Nnewlines):
		output += '\n'
	f.write(output)

### Printing ###
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
