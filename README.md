# Cactus_Alignments_Tools
The toolbox used to deal with cactus whole genome alignments. Because the efficiency of perl is limitted, I am happy to work with any person how would like to improve this pipeline. 

BTW, one anthoer tool, [halSynteny](https://github.com/ComparativeGenomicsToolkit/hal.git) would also be a good way to extracted the reasonable single copy alignments by a pairwise way. Here my pipeline support a chance to deal with all synteny chain based on one uniq reference species.

Requirement:
https://github.com/ComparativeGenomicsToolkit/hal.git
And add the bin to your PATH

How to use:
You could use the extract_shell.pl to create the work shell. 


# Ho to work
To obtain single-copy alignments from Cactus alignments, a pipeline to get the best synteny blocks was designed as follows:
![image](https://user-images.githubusercontent.com/9262911/139777521-502855a7-4b0d-46e6-a2d2-9a18132f111c.png)

1. Determination of the reference species. The reference species would provide the positional coordinates for synteny identification. 

2. For each query species, aligned blocks from the same chromosome and the same strand were regarded as a synteny chain.

3. For the synteny chain with insertions, the original synteny chain was subdivided into several synteny chains with locally colinear blocks to make sure each synteny chain had a continuous order within a single strand.

4. For each query species, all synteny chains were ordered by length without gaps and clustered to the different types at 100-bp intervals.

5. For each type, all synteny chains were further ordered by aligned scores (S). MA is the number of matched bases, and MI is the number of aligned but unmatched bases. Insertions and deletions were also calculated. The number of break points in the query is denoted by B, and the number of gaps in the query is denoted by G. So, S was calculated as

![image](https://user-images.githubusercontent.com/9262911/139777602-83cdbc40-1347-478b-b7a4-728fb48afde2.png)

6. For each query species, the base from the best synteny chain (from higher length and score to lower) was selected to constitute the final best synteny blocks.

7. To prevent high chimeric synteny blocks, all insertions shorter than 10 bp were eliminated.

Finally, single-copy synteny blocks were obtained from each query species. To optimize this process, the whole-genome alignment was split into multiple 100-kb (for penguins) or 10-kb (for all birds) windows based on reference and were filtered in parallel.
