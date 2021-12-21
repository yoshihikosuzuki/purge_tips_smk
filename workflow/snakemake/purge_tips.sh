#!/bin/bash
shopt -s expand_aliases && source ~/.bashrc && set -e || exit 1
alias error='printf "\033[0;35m[ERROR]\033[0m ";'
alias warn='printf "\033[0;33m [WARN]\033[0m ";'
alias info='printf "\033[0;32m [INFO]\033[0m ";'
_exit_w_msg() {
    MSG=$1
    error echo -e $MSG
    kill -INT $$
}

USAGE="Usage: purge_tips [-T<n_threads(64)>] [-e<e_thres(5)>] <hifi.fastq>"
N_THREADS=64 E_THRES=5
while getopts :T:e:h OPT; do
    case $OPT in
    T) N_THREADS=$OPTARG ;;
    e) E_THRES=$OPTARG ;;
    h) _exit_w_msg "$USAGE" ;;
    \?) _exit_w_msg "Illegal option: -$OPTARG" ;;
    esac
done
READS=${@:$OPTIND:1}

run_hifiasm() {
    # Dependency
    ml hifiasm gfatools
    # Input data/parameter
    READS_FASTQ=$1
    N_THREADS=$2
    # Temporary variable
    OUT_PREFIX=$(basename ${READS_FASTQ} .gz)
    OUT_PREFIX=${OUT_PREFIX%.*}
    # Output
    OUT_BAM=${REF_PREFIX}.${READS_PREFIX}.sorted.bam
    # Command
    hifiasm -o ${OUT_PREFIX} -t ${N_THREADS} ${READS_FASTQ}
    for DATA in *tg.gfa; do
        gfatools gfa2fa ${DATA} > ${DATA%.gfa}.fasta
    done
    # Return value
    ASM_GFA=${OUT_PREFIX}.bp.p_utg.noseq.gfa
    ASM_FASTA=${OUT_PREFIX}.bp.p_utg.fasta
}

run_winnowmap() {
    # Dependency
    ml samtools winnowmap
    # Input data/parameter
    REF_FASTA=$1
    READS_FASTQ=$2
    N_THREADS=$3
    TYPE="map-pb"
    # Temporary variable
    REF_PREFIX=${REF_FASTA%.*}
    READS_PREFIX=$(basename ${READS_FASTQ} .gz)
    READS_PREFIX=${READS_PREFIX%.*}
    # Intermediate output
    REF_MERYL=${REF_PREFIX}.meryl
    REP=${REF_PREFIX}.rep
    # Output
    OUT_BAM=${REF_PREFIX}.${READS_PREFIX}.sorted.bam
    # Command
    meryl count k=15 output ${REF_MERYL} ${REF_FASTA}
    meryl print greater-than distinct=0.9998 ${REF_MERYL} >${REP}
    winnowmap -W ${REP} -t${N_THREADS} -ax ${TYPE} --eqx --secondary=no -Y ${REF_FASTA} ${READS_FASTQ} |
        samtools sort -@${N_THREADS} -o ${OUT_BAM}
    samtools index -@${N_THREADS} ${OUT_BAM}
    # Return value
    READS_BAM=${OUT_BAM}
}

purge_tips() {
    # Dependency
    ml seqkit samtools
    # Input data/parameter
    NOSEQ_GFA=$1
    READS_FASTQ=$2
    READS_BAM=$3
    E_THRES=$4
    # Intermediate output
    NG_CNAMES=${NOSEQ_GFA}.e${E_THRES}.ng_cnames
    NG_RNAMES=${NOSEQ_GFA}.e${E_THRES}.ng_rnames
    # Output
    NG_FASTQ=${READS_FASTQ/.fastq/.ng.fastq}
    OK_FASTQ=${READS_FASTQ/.fastq/.ok.fastq}
    # Command
    awk -v thres_cov=${E_THRES} '/^S/ && int(substr($5,6)) <= thres_cov {print $2}' ${NOSEQ_GFA} >${NG_CNAMES}
    if [ ! -s ${NG_CNAMES} ]; then
        info echo "No low-cov contigs found. Finish purge_tips."
        exit 1
    fi
    samtools view ${READS_BAM} $(cat ${NG_CNAMES}) |
        cut -f1 |
        sort -u >${NG_RNAMES}
    seqkit grep -v -n -f ${NG_RNAMES} ${READS_FASTQ} > ${OK_FASTQ}
    seqkit grep -n -f ${NG_RNAMES} ${READS_FASTQ} > ${NG_FASTQ}
    # Return value
    FILTERED_FASTQ="filtered.fastq"
    ln -sf ${OK_FASTQ} ${FILTERED_FASTQ}
}

# run_hifiasm ${READS} ${N_THREADS}
# run_winnowmap ${ASM_FASTA} ${READS} ${N_THREADS}
# purge_tips ${ASM_GFA} ${READS} ${READS_BAM} ${E_THRES}
echo "hoge"
