#!/bin/bash

# FN OSTRAVA
# CREATED BY PETR BROŽ, Bioxsys s.r.o.

exec &>> rna_fuze.log

T="$(date +%s)"

mkdir -p tmp
rm -f FINISH.TXT

####################VARIABLES
#source /var/anaconda3/etc/profile.d/conda.sh

#CONDA_BASE=$(conda info --base)
cpu_total=$(nproc)
cpu=$(nproc)
full_path=$(readlink -e fusion.sh)


####################VARIABLES

fastp_env=/home/finalistdx/anaconda3/envs/fastp_env/bin/fastp
arriba=/home/finalistdx/anaconda3/envs/arriba_env/bin/arriba
star=/home/finalistdx/anaconda3/envs/arriba_env/bin/STAR
plot=/home/finalistdx/anaconda3/envs/arriba_env/bin/draw_fusions.R


dict=/media/finalistdx/Analysis/Databases/iGenomes/hg38/genome.dict
STAR_INDEX=/home/finalistdx/Databases/iGenomes/
REFERENCE=/media/finalistdx/Analysis/Databases/iGenomes/hg38/genome.fa
GENCODE=/media/finalistdx/Analysis/fusions/databases_fusion/gencode.v28.annotation.gtf
BLACK_LIST=/media/finalistdx/Analysis/fusions/databases_fusion/blacklist_hg38_GRCh38_2018-11-04.tsv.gz
COSMIC_FUSION=/media/finalistdx/Analysis/fusions/databases_fusion/CosmicFusionExport.tsv
LIB_DIR=/home/finalistdx/Databases/iGenomes/hg38/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir/
CYTO=/media/finalistdx/Analysis/fusions/databases_fusion/cytobands_hg38_GRCh38_2018-02-23.tsv
DOM=/media/finalistdx/Analysis/fusions/databases_fusion/protein_domains_hg38_GRCh38_2018-03-06.gff3

#star fusion
resource_lib=/media/finalistdx/Analysis/fusions/databases_fusion/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play/ctat_genome_lib_build_dir

#fusion catcher db
fusioncatcher_db=/home/finalistdx/anaconda3/envs/fusioncatcher_env/share/fusioncatcher-1.33/db/human_v102


function merge_multilanes() {

    for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F | rev | cut -c 22- | rev ; done | sort | uniq)

        do

        echo "INFO: Processing: $i"

    cat "$i"_L00*_R1_001.fastq.gz > "$i"_ME_L001_R1_001.fastq.gz


    cat "$i"_L00*_R2_001.fastq.gz > "$i"_ME_L001_R2_001.fastq.gz

    done;

    find . -type f \( -iname "*.gz" ! -iname "*_ME_L00*" \) -exec rm -rf {} +


    }


    function trimming() {


    for i in $(ls *.fastq.gz | rev | cut -c 22- | rev | sort | uniq)

    do

    echo "INFO: Analyzing file: $i"

    $fastp_env -i ${i}_L001_R1_001.fastq.gz -I ${i}_L001_R2_001.fastq.gz -o ${i}_trim_L001_R1_001.fastq.gz -O ${i}_trim_L001_R2_001.fastq.gz -j ${i}.json -h ${i}.html -w $cpu_total

    done


    mkdir -p original_fastq

    find -maxdepth 1 -iname "*fastq.gz" ! -iname "*trim*fastq.gz" -exec mv {} original_fastq/ \;


}


function star_fusion() {

    for i in $(ls *trim*.fastq.gz | rev | cut -c 22- | rev | sort | uniq)

    do

    #conda activate arriba_env

    echo "INFO: RUN STAR ALIGNER"

    echo "INFO: Analyzing file: $i"

    $star \
        --runThreadN "$cpu" \
        --genomeDir "$STAR_INDEX" \
        --genomeLoad NoSharedMemory \
        --readFilesIn ${i}_L001_R1_001.fastq.gz ${i}_L001_R2_001.fastq.gz \
        --readFilesCommand "gunzip -c" \
        --twopassMode Basic \
        --outStd BAM_Unsorted --outSAMtype BAM Unsorted \
        --outSAMunmapped Within \
        --outBAMcompression 0 \
        --outFilterMultimapNmax 1 --outFilterMismatchNmax 3 \
        --outFileNamePrefix ${i} \
        --chimSegmentMin 10 \
        --alignMatesGapMax 100000 \
        --alignIntronMax 100000 \
        --chimOutType WithinBAM SoftClip \
        --chimJunctionOverhangMin 10 \
        --chimScoreMin 1 \
        --chimOutJunctionFormat 1 \
        --chimScoreDropMax 30 \
        --chimScoreJunctionNonGTAG 0 \
        --chimScoreSeparation 1 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --chimNonchimScoreDropMin 10 \
        --peOverlapNbasesMin 12 \
        --peOverlapMMp 0.1 \
        --chimSegmentReadGapMax 3 | samtools sort -@ "$cpu" -T tmp -O BAM -o ${i}_out_sorted.bam -
    rm -f ${i}_out.bam

    done;

    #conda deactivate

    T="$(($(date +%s)-T))"

    echo "INFO: Time of STAR alignment in seconds: ${T} s"

}


function arriba() {

    #conda activate arriba_env

    for i in $(ls *sorted.bam)

    do

    echo "INFO: RUN FUSION DETECTION"

    echo "INFO: Analyzing file: $i"

    $arriba \
        -x $i \
        -o ${i%.bam}_fusions.tsv \
        -O ${i%.bam}_fusions.discarded.tsv \
        -a "$REFERENCE" \
        -g "$GENCODE" \
        -b "$BLACK_LIST" \
        -k $COSMIC_FUSION # see section "Complete Fusion Export" at http://cancer.sanger.ac.uk/cosmic/download
    #   -d structural_variants_from_WGS.tsv \

    done;

    #conda deactivate

    # concatanate all results with header to one final final report (include sample name)

    touch tmp_result.xls

    for i in *_fusions.tsv;

    do name=${i%_out*}

    awk -v name=$name '{OFS="\t"}NR>1{print name,$0}' $i >> tmp_result.xls

    done

    file=$(ls -1 *_fusions.tsv | awk 'NR==1{print $0}')

    head -1 $file | awk '{OFS="\t"}{print "NAME",$0}' | tr -d "#" | cat - tmp_result.xls > final_fusions.xls

    rm -f tmp_result.xls

    T="$(($(date +%s)-T))"
    echo "INFO: Time of fusions in seconds: ${T} s"

}


function plot_arriba() {

    #conda activate arriba_env

    for i in *fusions.tsv

    do

    echo "INFO: Plotting fusion: $i"

    $plot \
        --fusions=$i \
        --output=${i%.tsv}.pdf \
        --annotation=$GENCODE \
        --cytobands=$CYTO \
        --proteinDomains=$DOM \
    #    --alignments=${i%_fusions.tsv}.bam \

    done

    T="$(($(date +%s)-T))"

    echo "INFO: Time of fusion pipeline in seconds: ${T} s"

    touch FINISH_ALGORITMUS_1.TXT
    #conda deactivate

}


function star_fusion_all() {

    source /home/finalistdx/anaconda3/bin/activate star-fusion_env

    echo "INFO: Starting new star fusion analysis"

    # Loop through trimmed fastq files
    for sample_id in $(find . -maxdepth 1 \( -type f -o -type l \) -name "*trim*.fastq.gz" | while read F; do
        basename "$F" | rev | cut -c 22- | rev
    done | sort | uniq); do

    # Construct paired-end read paths
    left_fq="${sample_id}_L001_R1_001.fastq.gz"
    right_fq="${sample_id}_L001_R2_001.fastq.gz"

    # Run first part of STAR-Fusion with error checking
    if STAR-Fusion --genome_lib_dir "$resource_lib" \
                    --left_fq "$left_fq" \
                    --right_fq "$right_fq" \
                    --CPU $cpu \
                    --tmpdir tmp \
                    --examine_coding_effect \
                    --output_dir "${sample_id}_out_star_fusion" &> /dev/null; then
        echo "STAR-Fusion completed first step for sample: $sample_id"

        # Run second part of STAR-Fusion (assuming Chimeric.out.junction exists in current directory)
        STAR-Fusion --genome_lib_dir "$resource_lib" \
                -J Chimeric.out.junction \
                --output_dir "${sample_id}_out_star_fusion" &> /dev/null

        mv "${sample_id}_out_star_fusion"/star-fusion.fusion_predictions.abridged.tsv "${sample_id}_star_fusion.tsv"

        if [ $? -eq 0 ]; then
        echo "STAR-Fusion completed for sample: $sample_id"

        else
        echo "Error running second part of STAR-Fusion for sample: $sample_id"
        fi
    else
        echo "Error running first part of STAR-Fusion for sample: $sample_id"
    fi

    done

    echo "STAR-Fusion analysis completed for all samples."
    conda deactivate

    for i in *_star_fusion.tsv;

    do name=${i%_trim*}

    awk -v name=$name '{OFS="\t"}NR>1{print name,$0}' $i >> tmp_result.xls

    done

    file=$(ls -1 *_star_fusion.tsv | awk 'NR==1{print $0}')

    head -1 $file | awk '{OFS="\t"}{print "NAME",$0}' | tr -d "#" | cat - tmp_result.xls > final_star_fusions.xls

    rm -f tmp_result.xls

    touch FINISH_ALGORITMUS_2.TXT

}


function fusioncatcher_all() {

    source /home/finalistdx/anaconda3/bin/activate fusioncatcher_env

    for i in $(find . -maxdepth 1  \( -type f -o -type l \) -name "*trim*.fastq.gz" | while read F; do basename $F | rev | cut -c 22- | rev ; done | sort | uniq); do
    echo "INFO: Processing FUSIONCATCHER sample: $i"

    fusioncatcher -d $fusioncatcher_db \
    -i "${i}"_L001_R1_001.fastq.gz,"${i}"_L001_R2_001.fastq.gz \
    -o "${i}"_fusion_catcher \
    --skip-conversion-grch37 \
    --skip-fastqtk \
    -p $cpu_total

    mv "${i}_fusion_catcher"/final-list_candidate-fusion-genes.txt "${i}_fusion_catcher.tsv"

    done

    conda deactivate

    touch FINISH_ALGORITMUS_3.TXT

}

##########
#merge_multilanes
# trimming  #
# star_fusion #
# arriba  #
# plot_arriba  #
star_fusion_all
# fusioncatcher_all
##########


T="$(($(date +%s)-T))"
echo "INFO: Time of fusion pipeline in seconds: ${T} s"

