bw_file=$1
bed_file_200bp=$2
output=$3
script_folder=$4

echo get bed
cat $bed_file_200bp | awk -F '\t' -v OFS='\t' '{if (int(($2+$3)/2)-2500 > 0) print $1, int(($2+$3)/2)-2500, int(($2+$3)/2)+2500, $1"_"$2"_"$3; else print $1, 0, int(($2+$3)/2)+2500, $1"_"$2"_"$3 }' > bed_file_5kb
cat $bed_file_200bp | awk -F '\t' -v OFS='\t' '{if (int(($2+$3)/2)-5000 > 0) print $1, int(($2+$3)/2)-5000, int(($2+$3)/2)+5000, $1"_"$2"_"$3; else print $1, 0, int(($2+$3)/2)+5000, $1"_"$2"_"$3 }' > bed_file_10kb

echo get average signal
time $script_folder/bigWigAverageOverBed $bw_file $bed_file_200bp chip.mean.200bp.tab
time $script_folder/bigWigAverageOverBed $bw_file bed_file_5kb chip.mean.5kb.tab
time $script_folder/bigWigAverageOverBed $bw_file bed_file_10kb chip.mean.10kb.tab
echo sort by peak id
sort -k1,1 chip.mean.200bp.tab > chip.mean.200bp.sort.tab
sort -k1,1 chip.mean.5kb.tab > chip.mean.5kb.sort.tab
sort -k1,1 chip.mean.10kb.tab > chip.mean.10kb.sort.tab

echo get input mat
cut -f1 chip.mean.200bp.sort.tab | awk -F '_' -v OFS='\t' '{print $1, $2, $3}' > chip.meanrc.bg.txt
time cut -f6 chip.mean.200bp.sort.tab > chip.mean.200bp.sort.tab.txt
time cut -f6 chip.mean.5kb.sort.tab > chip.mean.5kb.sort.tab.txt
time cut -f6 chip.mean.10kb.sort.tab > chip.mean.10kb.sort.tab.txt
time paste chip.meanrc.bg.txt chip.mean.5kb.sort.tab.txt chip.mean.10kb.sort.tab.txt | sort -k1,1 -k2,2n > chip.meanrc.bg.5_10kb.tab.txt

echo get output signal
time Rscript $script_folder/get_local_bg_rc.atac.R chip.mean.200bp.sort.tab.txt chip.meanrc.bg.5_10kb.tab.txt $output

echo clean folder
rm bed_file_5kb bed_file_10kb
rm chip.mean.200bp.tab chip.mean.5kb.tab chip.mean.10kb.tab
rm chip.mean.200bp.sort.tab chip.mean.5kb.sort.tab chip.mean.10kb.sort.tab

rm chip.mean.200bp.sort.tab.txt chip.mean.5kb.sort.tab.txt chip.mean.10kb.sort.tab.txt
rm chip.meanrc.bg.txt chip.meanrc.bg.5_10kb.tab.txt

#time bash ../get_macs_bg_script/get_macs_bg.sh T_CD8_SPL.atac.meanrc.raw.bw test.1.bed test.bedgraph /storage/home/gzx103/scratch/vision/all_final_data/rc_norm/get_macs_bg_script


