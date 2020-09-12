#!/bin/bash

## print a message and exit
function die () {
        echo $1 >> $log
        exit 100
}

function binSignal () {
   #get average signal for bins for experiment and control if available,
   # if no control use all one's
   # output is bedgraph files (2)

   #expects bed3&4 files of windows named:
   w3="windowsNoBlack.noid.bed"
   w4="windowsNoBlack.withid.bed"
   #tmpControl.bed
   #writes bedgraphs named:
   #<celltype>.<mark>.<id>.ip.idsort.bedgraph
   #<celltype>.<mark>.<id>.ctrl.idsort.bedgraph

   line=$1
   outdir=$2
   log=$4
   script_dir=$3
   read -a f <<< "$line"
   #format: cell<tab>mark<tab>id<tab>file.bw<tab>input.bw
   
   cd $outdir
   test $? -eq 0 || die "Failed cd $outdir"

   command=$script_dir"bigWigAverageOverBed ${f[3]} ${w4} ${f[0]}.${f[1]}.${f[2]}.tmp"
   echo $command
   $command
   test $? -eq 0 || die "Failed $command"
   
   sort -k1,1n ${f[0]}.${f[1]}.${f[2]}.tmp | cut -f 5 | paste $w3 - > ${f[0]}.${f[1]}.${f[2]}.ip.idsort.bedgraph
   test $? -eq 0 || die "Failed to sort, cut, and paste .tmp to .bedgraph"

   if [ ${f[4]+1} ]
      then
         command=$script_dir"bigWigAverageOverBed ${f[4]} $w4 ${f[0]}.${f[1]}.${f[2]}.control.tmp"
         $command
         test $? -eq 0 || die "Failed $command"

         sort -k1,1n ${f[0]}.${f[1]}.${f[2]}.control.tmp | cut -f 5 | paste $w3 - > ${f[0]}.${f[1]}.${f[2]}.ctrl.idsort.bedgraph
         test $? -eq 0 || die "Failed to sort, cut, paste from .tmp to .bedgraph for control"
   else #copy all 1s file from outside loop
         command="cp tmpControl.bed ${f[0]}.${f[1]}.${f[2]}.ctrl.idsort.bedgraph"
         $command
         test $? -eq 0 || die "Failed $command"
   fi

#do cleanup outside of threads
   return 0
}

OUTDIR=$1
GENOMESIZES=$2
BLACK=$3
METADATA=$4
WORKFLOWDIR=$5
nthreads=$6
logfile=$7

## create the output directory if necessary
if [ ! -d "$OUTDIR" ]; then
        mkdir -p "$OUTDIR"
        test $? -eq 0 || die "Insufficient privileges to create $OUTDIR"
fi

export PATH=/gpfs/group/rch8/legacy/group/bin:$PATH


if [ ! -f "$OUTDIR/windowsNoBlack.noid.bed" ]; then
   echo "makewindow noid"
   command="bedtools makewindows -g $GENOMESIZES -w 200"
   $command > $OUTDIR/windows.bed
   test $? -eq 0 || die "failed: $command"
   echo $command

   command="bedtools subtract -a $OUTDIR/windows.bed -b $BLACK"
   $command > $OUTDIR/windowsNoBlack.noid.bed
   test $? -eq 0 || die "failed: $command"
   echo $command
fi

if [ ! -f "$OUTDIR/windowsNoBlack.withid.bed" ]; then
   echo "makewindow withid"
   i=1
   while read -r line; do
      echo -e "$line\t${i}"
      let i++
   done < $OUTDIR/windowsNoBlack.noid.bed > $OUTDIR/windowsNoBlack.withid.bed
   test $? -eq 0 || die "failed to add ids to windows"
fi
echo "Windows done"

#print a 1 foreach window, default control if none given
#create file outside of loop, copy as needed
while read -r line; do printf "$line\t1\n"; done < $OUTDIR/windowsNoBlack.noid.bed > $OUTDIR/tmpControl.bed
test $? -eq 0 || die "failed: to create default control"

#convert bw to bedgraphs
#read metadata file for bigwigs and control bigwigs
i=0
while read -r line 
   do
   binSignal "$line" $OUTDIR $WORKFLOWDIR &
   ((i++))
   if [ $i -ge $nthreads ]; then
      wait
      i=0
   fi
done < $METADATA

wait
#when all threads finished clean up
rm tmpControl.bed
rm *.tmp

exit 0
