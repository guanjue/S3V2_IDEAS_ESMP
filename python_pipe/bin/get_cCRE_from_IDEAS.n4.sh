id_name=$1
script_dir=$2
working_dir=$3
### go to IDEAS epigenetic state folder
#cd $id_name'_IDEAS_output'
#id_name='run_IDEAS_S3V2_mean'
#script_dir='/storage/home/gzx103/scratch/S3norm_IDEAS_cCRE/cCRE_callbin/'
#time bash get_cCRE_from_IDEAS.sh 'run_IDEAS_S3V2_mean' '/storage/home/gzx103/scratch/S3norm_IDEAS_cCRE/cCRE_callbin/'

cd $working_dir
### get unmerge bin bed
time Rscript $script_dir/get_bins_from_IDEAS.n4.R $id_name

### define state combinations
declare -a state_combination=("3" "32" "321")

### get AVERAGE merge peak
for st in "${state_combination[@]}"
do
	echo $st
	bedtools merge -i $id_name'.AVERAGE.'$st'.bed' > $id_name'.AVERAGE.'$st'.M.bed'
done

### get ct peaks
for st in "${state_combination[@]}"
do
	echo $st
	cp $id_name'.AVERAGE.'$st'.M.bed' $id_name'.AVERAGE_allct.'$st'.M.bed'
done

### get ct average intersect & uniq peaks
for ct in $(tail -n+2 $id_name'.ct.list.txt')
do
	echo $ct
	for st in "${state_combination[@]}"
	do
		echo $st
		### mergine ct bins
		bedtools merge -i $id_name'.'$ct'.'$st'.bed' > $id_name'.'$ct'.'$st'.M.bed'
		### intersect
		bedtools intersect -a $id_name'.'$ct'.'$st'.M.bed' -b $id_name'.AVERAGE_allct.'$st'.M.bed' > tmp1.bed
		### ct uniq
		bedtools intersect -a $id_name'.'$ct'.'$st'.M.bed' -b $id_name'.AVERAGE_allct.'$st'.M.bed' -v > tmp2.bed
		### ave uniq
		bedtools intersect -b $id_name'.'$ct'.'$st'.M.bed' -a $id_name'.AVERAGE_allct.'$st'.M.bed' -v > tmp3.bed
		### pool all
		cat tmp1.bed tmp2.bed tmp3.bed | sort -k1,1 -k2,2n > tmp4.bed && mv tmp4.bed $id_name'.AVERAGE_allct.'$st'.M.bed'
	done
done

#cp $id_name'.AVERAGE.321.M.bed' $id_name'.AVERAGE_allct.321.M.bed'

### get ct binary mat
for st in "${state_combination[@]}"
do
	echo $st
	cp $id_name'.AVERAGE_allct.'$st'.M.bed' $id_name'.AVERAGE_allct.'$st'.M.mat.txt'
done


### get columns
for ct in $(tail -n+2 $id_name'.ct.list.txt')
do
	echo $ct
	for st in "${state_combination[@]}"
	do
		echo $st
		bedtools intersect -a $id_name'.AVERAGE_allct.'$st'.M.bed' -b $id_name'.'$ct'.321.bed' -c > tmp1.txt
		#bedtools intersect -a $id_name'.AVERAGE_allct.'$st'.M.bed' -b $id_name'.'$ct'.32.bed' -c > tmp1.txt
		cut -f4 tmp1.txt > tmp2.txt
		paste $id_name'.AVERAGE_allct.'$st'.M.mat.txt' tmp2.txt > $id_name'.AVERAGE_allct.'$st'.M.mat.txt.tmp' \
		&& mv $id_name'.AVERAGE_allct.'$st'.M.mat.txt.tmp' $id_name'.AVERAGE_allct.'$st'.M.mat.txt'
	done
done

### remove unreproducible ones
for st in "${state_combination[@]}"
do
	echo $st
	time Rscript $script_dir/get_reproducible_pk.R $id_name'.ct.list.txt' $id_name'.AVERAGE_allct.'$st'.M.mat.txt' $id_name'.AVERAGE_allct.'$st'.M.rep0.bed'
	bedtools merge -i $id_name'.AVERAGE_allct.'$st'.M.rep0.bed' > $id_name'.AVERAGE_allct.'$st'.M.rep1.bed'
	bedtools intersect -a $id_name'.AVERAGE_allct.'$st'.M.rep1.bed' -b $id_name'.AVERAGE.'$st'.M.bed' -v > tmp5.bed
	cat $id_name'.AVERAGE.'$st'.M.bed' tmp5.bed | sort -k1,1 -k2,2n > $id_name'.AVERAGE.'$st'.M.rep.bed'
done

### get final cCRE list
bedtools intersect -a $id_name'.AVERAGE_allct.32.M.rep1.bed' -b $id_name'.AVERAGE_allct.3.M.rep1.bed' -v > $id_name'.AVERAGE_allct.32.NOT.3.M.rep.bed'
bedtools intersect -a $id_name'.AVERAGE_allct.321.M.rep1.bed' -b $id_name'.AVERAGE_allct.32.M.rep1.bed' -v > $id_name'.AVERAGE_allct.321.NOT.32.M.rep.bed'
cat $id_name'.AVERAGE.3.M.rep.bed' $id_name'.AVERAGE_allct.32.NOT.3.M.rep.bed' $id_name'.AVERAGE_allct.321.NOT.32.M.rep.bed' | sort -k1,1 -k2,2n > $id_name'.cCRE.bed'
#cat $id_name'.AVERAGE.3.M.rep.bed' $id_name'.AVERAGE_allct.32.NOT.3.M.rep.bed' | sort -k1,1 -k2,2n > $id_name'.cCRE.bed'

bedtools merge -d 200 -i $id_name'.cCRE.bed' > $id_name'.cCRE.M.bed'

### get ct cCREs
for ct in $(tail -n+2 $id_name'.ct.list.txt')
do
	echo $ct
	bedtools intersect -a $id_name'.'$ct'.32.M.bed' -b $id_name'.'$ct'.3.M.bed' -v > $id_name'.'$ct'.32.NOT.3.M.bed'
	bedtools intersect -a $id_name'.'$ct'.321.M.bed' -b $id_name'.'$ct'.32.M.bed' -v > $id_name'.'$ct'.321.NOT.32.M.bed'
	cat $id_name'.'$ct'.3.M.bed' $id_name'.'$ct'.32.NOT.3.M.bed' $id_name'.'$ct'.321.NOT.32.M.bed' | sort -k1,1 -k2,2n > $id_name'.'$ct'.cCRE.bed'
#	cat $id_name'.'$ct'.3.M.bed' $id_name'.'$ct'.32.NOT.3.M.bed' | sort -k1,1 -k2,2n > $id_name'.'$ct'.cCRE.bed'
done

### clean folder
for ct in $(cat $id_name'.ct.list.txt')
do
	echo $ct
	rm $id_name'.'$ct*3*
done


