#!/usr/bin/env bash
#set -e
#set -o pipefail

## print a message and exit
function die () {
        echo $1 >> $logfile
        exit 100
}

## function to convert one file
function get_bw {
    GENOMESIZES=$1
    bedgraph_file=$2
    OUTDIR=$3
    SCRIPTDIR=$4

    sort -k1,1 -k2,2n ${bedgraph_file} > ${bedgraph_file}.pksort.bedgraph
    test $? -eq 0 || die "Failed to sort ${bedgraph_file}"
    #$SCRIPTDIR'bedGraphToBigWig' ${bedgraph_file}.pksort.bedgraph ${GENOMESIZES} ${bedgraph_file}.bw
    bedGraphToBigWig ${bedgraph_file}.pksort.bedgraph ${GENOMESIZES} ${bedgraph_file}.bw
    test $? -eq 0 || die "Failed to convert ${bedgraph_file} to bigwig"
    mv ${bedgraph_file}.bw $OUTDIR
    test $? -eq 0 || die "Failed to move ${bedgraph_file}"
    rm ${bedgraph_file}.pksort.bedgraph
    return 0
}

SCRIPTDIR=$1
OUTPREFIX=$2
nthreads=$3
GENOMESIZES=$4
logfile=$5

## create the output directory if necessary
if [ ! -d "${OUTPREFIX}_NBP" ]; then
        mkdir -p "${OUTPREFIX}_NBP"
        test $? -eq 0 || die "Insufficient privileges to create ${OUTPREFIX}_NBP"
fi
if [ ! -d "${OUTPREFIX}_RC" ]; then
        mkdir -p "${OUTPREFIX}_RC"
        test $? -eq 0 || die "Insufficient privileges to create ${OUTPREFIX}_RC"
fi

export PATH=$SCRIPTDIR:$PATH

#assume working directory is the output directory
echo "log errors to $logfile"
echo "get bw......"
echo "for S3V2 RC"
i=0
shopt -s nullglob
for f in *.S3V2.bedgraph
do
   get_bw $GENOMESIZES $f "${OUTPREFIX}_RC/" $SCRIPTDIR &
   ((i++))
   if [ $i -ge $nthreads ]; then
      wait
      i=0
   fi
done

#get average
echo "for average signal"
for f in *.average_sig.bedgraph.S3V2.ave.bedgraph
do
   get_bw $GENOMESIZES $f "${OUTPREFIX}_RC/" $SCRIPTDIR &
   ((i++))
   if [ $i -ge $nthreads ]; then
      wait
      i=0
   fi
done

echo "for S3V2 -log10 p-value"
for f in *.S3V2.bedgraph.NBP.bedgraph
do
   get_bw $GENOMESIZES $f "${OUTPREFIX}_NBP" $SCRIPTDIR &
   ((i++))
   if [ $i -ge $nthreads ]; then
      wait
      i=0
   fi
done

echo "for average signal"
for f in *.average_sig.bedgraph.S3V2.ave.bedgraph.NBP.bedgraph
do
   get_bw $GENOMESIZES $f "${OUTPREFIX}_NBP" $SCRIPTDIR &
   ((i++))
   if [ $i -ge $nthreads ]; then
      wait
      i=0
   fi
done

wait
echo "get bw......Done"

exit 0

