
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels mvdbeek

conda update conda

conda create -n s3v2_ideas
conda install -n s3v2_ideas gsl=2.2.1 r python=3 ucsc_tools
conda install -n s3v2_ideas bedtools numpy scipy multiprocess
conda install -n s3v2_ideas r-doParallel r-foreach r-data.table

conda activate s3v2_ideas
