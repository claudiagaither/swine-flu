C:\Users\narma\OneDrive - University of North Carolina at Chapel Hill\PhD_research\Swine_flu\GISAID_work\sequences

2/24/2026

1) gisaid_epiflu_sequence.fasta- Narmada and Claudia obtained
	
	- Contains all H3N2 viruses and sub clades from swine in the USA collected and sampled as of 2/24/2026
	- Sequence headers have virus name, collection date, clade, segment, and segment number
2) align_flu_longleaf.sh - Claudia
	- created shell script to align sequences in LongLeaf


2/25/2026

1) H3N2_aligned_mafft.fasta- Narmada
	- Aligned input gisaid_epiflu_sequence.fasta using align_flu_longleaf.sh and Longleaf
	- MAFFT alignment with first deduplicating sequence names before aligning
2) .gitattributes- Narmada
	- I had to convert the file endings to Unix in Longleaf which Git (and Windows) will try to convert back to Windows line endings
	- This file will specify to keep these endings if we need to put these sequences back through Longleaf