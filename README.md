# POWEREDBiSeq
Scripts for POWEREDBiSeq, a tool which estimates power in bisulfite sequencing (BS) studies


To use POWEREDBiSeq, you must first download the script POWEREDBiSeqFunctions.R
In R, set the working directory to the folder that contains it and run <i>source("POWEREDBiSeqFunctions.R")</i>

You bisulfite sequencing data must be in the format: 
- Column 1: Site IDs (e.g. "chr1_100000478")
- Columns 2 to n+1 (where n is number of samples): DNA methylation values for each sample, scaled between 0-100
- Columns n+2 to 2n +1: Read depth value for each sample

If your data is in .cov files, function RRBSMatrixMaker can be used:

<i>RRBSMatrixMaker(filePath, coordinatesOfUniqueCharactersInSampleName = "all")</i>

Where <b>filePath</b> is the folder path to the coverage files (they must be in the same folder) and <b>coordinatesOfUniqueCharactersInSampleName</b> if the column names are unique sample ids. 


To apply POWEREDBiSeq, now that your BS data is prepared, you will need to input 
- <b>rd</b>: the read depth threshold you plan to use
- <b>meanDiff</b>: the mean difference in DNAm expected between your groups which must be between 0 - 1 and represents a proportion, not a percentage
- <b>nSamplesNeeded</b>: the minimum number of samples per group
- <b>pheno</b>: if you want to, you can include a factor vector of group assignments, otherwise samples will be randomly assigned

The output contains 
- the read depth
- power
- min and max power of 5 estimations
- the Bonferroni p-value used in the power calculation
- the proportion and number of sites that would be tested with given the filtering applied
