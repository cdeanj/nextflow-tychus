Output Files
============
Output files from each module (assembly and alignment) are described below. The ``directory`` subsections below indicate which directory results are stored in once the workflow completes. File names will be prefixed with the name of the FASTQ file used as input.

Assembly Module
---------------
Output from the assembly module includes the assemblies output from each assembler, as well as the hybrid assembly produced from the contig integration step from CISA. Also included are the genome annotations from Prokka.

Directory: Contigs
~~~~~~~~~~~~~~~~~~
============================ ===============================================================
File Name                    Description
============================ ===============================================================
abyss-contigs.fa             
velvet-contigs.fa	     
spades-contigs.fa            
idba-contigs.fa              
master_integrated_contigs.fa 
============================ ===============================================================

Directory: Annotations
~~~~~~~~~~~~~~~~~~~~~~
========= ===============================================================
File Name                    Description
========= ===============================================================
*.faa
*.ffn
*.fna
*.fsa
*.gbk
*.gff
*.sqn
*.tbl
*.txt
========= ===============================================================

Alignment Module
----------------
Output from the alignment module includes the results for alignment against the user inputted resistance, plasmid, and virulence databases. Results from the variant calling and phylogenetic analyses are also recorded.

Directory: Alignment
~~~~~~~~~~~~~~~~~~~~
============================ ===============================================================
File Name                    Description
============================ ===============================================================
*
============================ ===============================================================

Directory: Variants
~~~~~~~~~~~~~~~~~~~
============================ ===============================================================
File Name                    Description
============================ ===============================================================
*
============================ ===============================================================

Directory: Phylogeny
~~~~~~~~~~~~~~~~~~~~
============================ ===============================================================
File Name                    Description
============================ ===============================================================
*
============================ ===============================================================
