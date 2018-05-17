
def skip_lines(N,f):
	""" Simple. """
	for n in range(N):
		next(f)

def f_write(f,n_tabs,string):
	"""Write string to file with specified number of tabs and newline"""
	for n in range(0,n_tabs):
		f.write('\t')
	f.write(string + '\n')

