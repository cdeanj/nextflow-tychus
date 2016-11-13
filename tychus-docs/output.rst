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
abyss-contigs.fa             Abyss assembly contigs
velvet-contigs.fa	     Velvet assembly contigs
spades-contigs.fa            SPades assembly contigs
idba-contigs.fa              IDBA assembly contigs
master_integrated_contigs.fa CISA integrated contigs
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

============================ ===============================================================
File Name                    Description
============================ ===============================================================





============================ ===============================================================
