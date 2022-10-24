# Download from ENCODE data portal
# https://www.encodeproject.org/files/ENCFF356LFX/
wget https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz
gunzip ENCFF356LFX.bed.gz
mv ENCFF356LFX.bed hg38-blacklist.ENCODE.2020-05-05.ENCFF356LFX.bed

# Download from the following link:
# https://github.com/Boyle-Lab/Blacklist/tree/master/lists
wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz
wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg19-blacklist.v2.bed.gz
wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
gunzip *

