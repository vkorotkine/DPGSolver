clear

cd /Users/philipzwanenburg/Desktop/Research_Codes/DPGC/cases

#valgrind --leak-check=yes /Users/philipzwanenburg/Desktop/Research_Codes/DPGC/bin/DPGSolver.exe PeriodicVortex
valgrind --track-origins=yes --leak-check=yes /Users/philipzwanenburg/Desktop/Research_Codes/DPGC/bin/DPGSolver.exe PeriodicVortex

#valgrind --leak-check=yes /Users/philipzwanenburg/Desktop/Research_Codes/DPGC/bin/DPGSolver.exe SupersonicVortex
