C:\Users\narma\OneDrive - University of North Carolina at Chapel Hill\PhD_research\Swine_flu\GISAID_work\sequences\clade_subsets

Total number of sequences prior to cleaning: 
	# 4551

# Generated clade subsets by:
	1) Using mcc tree to assign 'unassigned' sequences based on shared parent nodes with other sequences that were pre-assigned clades from BV-BRC
	2) re-do clade groupings, grouping together:
		- 1990.4-like into 1990.4,
  		- 2010.1-like into 2010.1,
  		- Other-Human-2000-like into Other-Human-2000,
 		- Other-Human-2010-like into Other-Human-2010,
  		- Other-Human-1970-like into Other-Human-1970

	2) Doing a simple downsampling and equal proportion downsampling based on unique combination of state, date, shared parent node and assigned clade 

		a) Simple downsampling: just keeping unique sequences for each state, date, and shared parent node (n=4403)

		b) Equal downsampling: first only keeping one sequence for each unique location, date, clade, and shared parent grouping, then for each location, date, and clade grouping making sure there is one sequence representative, randomly sampling from these. This was repeated 3 times to make sure the same number of samples were being retained/kept in each iteration. 
			n = 3823