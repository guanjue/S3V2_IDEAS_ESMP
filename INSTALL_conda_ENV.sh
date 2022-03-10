
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels mvdbeek

conda update conda

conda create -n s3v2_ideas gsl=2.2.1 python=3 r=3.6 bedtools ucsc_tools numpy scipy multiprocess r-doParallel r-foreach r-data.table

conda activate s3v2_ideas

