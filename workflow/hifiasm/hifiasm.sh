#!/bin/bash
# Usage: ./hifiasm.sh <working-dir> <in-read> <n-threads> [<hifiasm-options>]

WDIR=$1
IN_READ=$2
N_THREADS=$3
OPTIONS=${@:4}

_IN_READ=$(basename ${IN_READ})
OUT_PREFIX=${_IN_READ%.gz}
OUT_PREFIX=${OUT_PREFIX%.*}.hifiasm

mkdir -p ${WDIR}
ln -sf $(readlink -f ${IN_READ}) ${WDIR}/${_IN_READ}

CDIR=$(pwd)
cd ${WDIR}
hifiasm -o ${OUT_PREFIX} -t ${N_THREADS} ${OPTIONS} ${_IN_READ}
for DATA in *tg.gfa; do
    gfatools gfa2fa ${DATA} >${DATA%.gfa}.fasta
done
cd ${CDIR}
