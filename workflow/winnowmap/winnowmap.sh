#!/bin/bash
# Usage: ./winnowmap.sh <working-dir> <in-read> <in_asm> <n-threads>

WDIR=$1
IN_READ=$2
IN_ASM=$3
OUT_PREFIX=$4
N_THREADS=$5

_IN_READ=$(basename ${IN_READ})
_IN_ASM=$(basename ${IN_ASM})
OUT_BAM=${OUT_PREFIX}.bam
OUT_MERYL=${OUT_PREFIX}.meryl
OUT_REP=${OUT_PREFIX}.rep

mkdir -p ${WDIR}
for FILE in ${IN_READ} ${IN_ASM}; do
    ln -sf $(readlink -f ${FILE}) ${WDIR}/$(basename ${FILE})
done

CDIR=$(pwd)
cd ${WDIR}
meryl count k=15 output ${OUT_MERYL} ${_IN_ASM}
meryl print greater-than distinct=0.9998 ${OUT_MERYL} >${OUT_REP}
winnowmap -W ${OUT_REP} -t${N_THREADS} -ax "map-pb" --eqx --secondary=no -Y ${_IN_ASM} ${_IN_READ} |
    samtools sort -@${N_THREADS} -o ${OUT_BAM}
samtools index -@${N_THREADS} ${OUT_BAM}
cd ${CDIR}
