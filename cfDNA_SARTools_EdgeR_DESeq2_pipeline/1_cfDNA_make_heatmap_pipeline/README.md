README — run_make_heatmap.sh
=============================

How to Run:
-----------
To run the script in the background using `nohup`:

---------------------------------------------------------------------------------------------
nohup bash run_make_heatmap.sh <input_directory> <reference_file> <output_prefix> > heatmap_run.out 2>&1 &
---------------------------------------------------------------------------------------------

- `<input_directory>`: Path to the folder containing your `.singleFrag.bed` files  
- `<reference_file>`: Path to the reference regions file (tab-separated, no header)  
- `<output_prefix>`: A label added to all output files  

Example:
---------------------------------------------------------------------------------------------
nohup bash run_make_heatmap.sh ./cfDNA_samples reference_regions.txt ATAC_10kbp > heatmap_run.out 2>&1 &
---------------------------------------------------------------------------------------------

You can monitor progress using:
---------------------------------------------------------------------------------------------
tail -f heatmap_run.out
---------------------------------------------------------------------------------------------

Purpose:
--------
This script automates the generation of heatmaps/metaplots using the `make_heatmap` tool 
to visualize cfDNA fragment enrichment across a set of genomic regions. 
It processes multiple `.singleFrag.bed` files in parallel, using a provided reference file, 
and records which samples were processed successfully or failed.

Input Files:
------------
1. **cfDNA Fragment Files (`*.singleFrag.bed`):**
   - Each file should contain cfDNA fragment data aligned to the genome.
   - These files will be looped over and processed individually.

2. **Reference File (e.g., `reference_regions.txt`):**
   - A tab-separated file with the following columns (no header line):
     Region_ID    Midpoint    Chromosome    Start    End

   Example:
   chr1:10000-10500    10250    chr1    10000    10500
   chr1:20000-20800    20400    chr1    20000    20800

Output Files:
-------------
For each input sample:
- Heatmap/metaplot file:
  SAMPLE_NAME_OUTPUTPREFIX.txt
- Command log for each sample:
  SAMPLE_NAME_OUTPUTPREFIX.log

Additionally:
- successful_samples.log → List of samples that completed successfully
- failed_samples.log → List of samples that failed to process

Script Logic:
-------------
- Loops through all `.singleFrag.bed` files in the given `<input_directory>`.
- Extracts the sample name from each filename.
- Constructs a custom output file name for each sample.
- Runs `make_heatmap` on each sample in the background.
- Allows up to 5 jobs to run concurrently (customizable).
- Logs output of each sample to its own log file.
- Updates success/failure logs based on exit status.

Important Notes:
----------------
- The reference file must NOT contain a header line.
- Ensure `make_heatmap` is executable:
  chmod +x make_heatmap
- Adjust MAX_CONCURRENT_JOBS and THREADS in the script if needed.
- Default resolution is set by: -5000 20 500
  → 5 kb upstream and downstream, 20 bins of 500 bp each

make_heatmap Options
-------------
	-t =  specify number of threads to use
	-b = specify method of defining relative bin locations:
		c = command line, fixed bin size
	-h = specify hit file type:
		b = basic bed
	-l = specify location within hit to match to gene list:
		c = center
	-a = specify location within each gene to anchor relative bin locations:
		u = user-specified, contained in second column of gene list file
	-v = specify method of determining value for each bin
		t = compute total of hit values
	-d = specify method of determining location of each bin
		p = interpret relative bin locations as physical distance
	-s = specify handling of strand identifiers
		b = utilize hits from both strands, also for use with data that lacks strand information
	-- End of command options

Author:
-------
Sakuntha Devaka Gunarathna  
Published: July 31, 2025
