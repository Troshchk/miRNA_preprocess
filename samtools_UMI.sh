#!/bin/bash


files=" ./46635Q_mapped/46635Q_hairpinh_bowtie.bam ./46635Q_mapped/46635Q_mature_non_h_bowtie.bam ./46635Q_mapped/46635Q_hairpin_non_h_bowtie.bam ./47356Q_mapped/47356Q_matureh_bowtie.bam ./47356Q_mapped/47356Q_hairpinh_bowtie.bam ./47356Q_mapped/47356Q_mature_non_h_bowtie.bam ./47356Q_mapped/47356Q_hairpin_non_h_bowtie.bam ./47370Q_mapped/47370Q_matureh_bowtie.bam ./47370Q_mapped/47370Q_hairpinh_bowtie.bam ./47370Q_mapped/47370Q_mature_non_h_bowtie.bam ./47370Q_mapped/47370Q_hairpin_non_h_bowtie.bam ./47371Q_mapped/47371Q_matureh_bowtie.bam ./47371Q_mapped/47371Q_hairpinh_bowtie.bam ./47371Q_mapped/47371Q_mature_non_h_bowtie.bam ./47371Q_mapped/47371Q_hairpin_non_h_bowtie.bam ./47442Q_mapped/47442Q_matureh_bowtie.bam ./47442Q_mapped/47442Q_hairpinh_bowtie.bam ./47442Q_mapped/47442Q_mature_non_h_bowtie.bam ./47442Q_mapped/47442Q_hairpin_non_h_bowtie.bam"


for file in $files  ; do samtools sort $file -o ${file/.bam/_sorted.bam}; samtools index ${file/.bam/_sorted.bam}; umi_tools dedup -I ${file/.bam/_sorted.bam} --output-stats=deduplicated -S ${file/.bam/_sorted_deduplicated.bam}; done
