--------------------------------------------------------------------------- 
 uKIN - Using Knowledge In Networks
---------------------------------------------------------------------------- 

 I. Input 

 There are three required input files:
 1) a network file 
 2) a prior knowledge file containing a list of nodes (genes) known to be
 disease associated, possibly with weights on them
 3) a file of newly implicated genes, each with a weight

 !!! Plus a required path to a Matlab executable since uKIN needs Matlab to do
 some heavy lifting of big matrices. The path should be specified as:
 matlab=/path/to/matlab after the three input files.
 
 Additionally, the user may provide:
 5) a custom value of the restart parameter alpha controlling how much 
 influence each type of information (prior and new) has.
 Default is alpha = 0.5
 6) a custom value of the diffusion parameter gamma controlling how far the prior
 knowledge spreads. Default is gamma = 1.
 7) output prefix which is used in the beginning of the name of 
 the output file

 II. Output

 output_prefix_results.txt is written in the uKIN directory. The file
 contains a list of candidate genes ranked by how frequently they are
 visited as the guided walks reach the stationary distribution.

 III. How to run

 Simply issue:
 ruby uKIN.rb network_file.txt prior_knowledge.txt new_genes.txt matlab=/path/to/matlab {alpha=0.4 gamma=0.8 output_prefix=my_output}

 
 Note that uKIN is implemented in Ruby but requires Matlab for the heavy matrix 
 operations. 
 You can simply install Ruby via: sudo apt-get install ruby-full
 Matlab is typically provided for free to academic institutions. You can also
 download a 30 day free trial from their website.

 IV. Input File Formats

 1. Network file: each line specifies an edge, white space delimited:
 GENE_ID GENE_ID

 2. Prior knowledge: list of genes
 GENE_ID
 GENE_ID

 3. Newly implicated genes: each line is white space delimited
 GENE_ID WEIGHT
 GENE_ID WEIGHT