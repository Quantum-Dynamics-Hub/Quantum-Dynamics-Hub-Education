# Contents of PROJECT folder:
___
Nonadiabatic MD of vinylchloride
1. (NewtonX) DFTB+ and Newton X inputs for geometry optimization, single point calculation, initial condition generation. Useful scripts (should be run in TRAJECTORIES/ folder):
	1. read_TRAJ.sh <start> <end> - grep and print last STEP for sequence of trajectories 
	2. start.sh <start> <end> - submit trajectories from TRAJ<start> to TRAJ<end>
	3. clean_TRAJ.sh <start> <end> - clean all trajectories from <start> to <end>
2. (NEXMD) NEXMD outputs and input
3. (figures) Contains images that were obtained by analyzing the results
- LaTeX and PDF versions of Report
