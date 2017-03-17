import sys 
import traceback

### Special Functions ###
def EXIT_TRACEBACK():
	print('\nError: Unsupported')
	traceback.print_stack()
	sys.exit()
