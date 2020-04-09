#!/usr/bin/env nextflow
//params.sample_csv = "HCC-1187_manifest.csv"
sample_attributes = []
file(params.sample_csv).readLines().each { it -> sample_attributes << it.tokenize(',') }

fastq_channel = Channel.fromFilePairs(params.input_fastqs)
fastq_channel1 = Channel.fromFilePairs(params.input_fastqs)

attr_channel = Channel.from(sample_attributes)
attr_channel1 = Channel.from(sample_attributes)



fastq_channel1.join(attr_channel1).subscribe {println it}
params.inputs = fastq_channel.join(attr_channel)

// process fake_align {

//   input:
//   tuple val(sample_name), file(fastqs),  val(tn) from params.input_channels

//   output:
//   tuple val(sample_name), val(tn), file("unaligned.bam") into ubam_ch1

//   """
//   echo ${sample_name}, ${tn}
//   head -n 2 ${fastqs[0]} > "unaligned.bam"
//   head -n 2 ${fastqs[1]} >> "unaligned.bam"
//   """
// }

// process print_align {
//   input:
//   tuple val(sample_name), val(tn), file(ubam) from ubam_ch1

//   """
//   echo ${sample_name}, ${tn}
//   cat ${ubam} | echo
//   """
// }

ref_fasta = Channel.value(params.genome_fasta)
ref_alt = Channel.value(params.genome_alt)
wes_coverage_interval_list = Channel.value(params.agilent_coverage_interval_list)
known_indels_sites_VCFs = Channel.value(params.genome_known_indels_sites_VCFs)
group_interval = Channel.value(params.sequence_group_interval)


process FastqToSam {

    input:
    tuple val(prefix), file(fastqs),  val(sample_name), val(tn), val(analyte) from params.inputs

    when:
    analyte == 'dna'

    output:
    tuple val(sample_name), val(tn), file("${sample_name}_unaligned.bam") into ubam_ch1, ubam_ch2

    cpus 2

    memory "150 GB"

    container "157538628385.dkr.ecr.us-west-2.amazonaws.com/bfx_docker:latest"

    errorStrategy 'retry'
    """
    java -XX:ParallelGCThreads=1 \
    -Xmx120G -jar /usr/picard/picard.jar \
    FastqToSam \
      FASTQ=${fastqs[0]} \
      FASTQ2=${fastqs[1]} \
      OUTPUT=${sample_name}_unaligned.bam \
      READ_GROUP_NAME=${sample_name} \
      SAMPLE_NAME=${sample_name} \
      LIBRARY_NAME=${sample_name} \
      PLATFORM=illumina \
      SEQUENCING_CENTER=PSNL
    """
}

process CollectQualityYieldMetrics {
  input:
  tuple val(sample_name), val(tn), file(ubam) from ubam_ch2

  output:
  tuple val(sample_name), val(tn), file("${sample_name}.unmapped.quality_yield_metrics") into ubam_metrics_ch

  cpus 1

  memory "12 GB"

  container "157538628385.dkr.ecr.us-west-2.amazonaws.com/bfx_docker:latest"

  errorStrategy 'retry'

  """
  java \
  -Xmx9G -jar /usr/picard/picard.jar \
  CollectQualityYieldMetrics \
    INPUT=${ubam} \
    OQ=true \
    OUTPUT=${sample_name}.unmapped.quality_yield_metrics
  """
}

process SamtoFastqAndBwaMemAndMba {


  cpus 32

  memory "12 GB"

  container "157538628385.dkr.ecr.us-west-2.amazonaws.com/bfx_docker:latest"

  errorStrategy 'retry'

  input:
  tuple val(sample_name), val(tn), file(ubam) from ubam_ch1
  val ref_fasta
  val ref_alt

  output:
  tuple val(sample_name), val(tn), file("${sample_name}.aligned.unsorted.bam") into unsort_bam_ch

  """
    set -o pipefail
    set -e
    export bwa_version=`bwa 2>&1 |  grep -e '^Version' |  sed 's/Version: //'`
    export bwa_commandline="bwa mem -K 100000000 -p -v 3 -t 16 -Y"

    # set the bash variable needed for the command-line
    # if ref_alt has data in it,
    if [ -s ${ref_alt} ]; then
      java \
      -Xmx5G -jar /usr/picard/picard.jar \
        SamToFastq \
        INPUT=${ubam} \
        FASTQ=/dev/stdout \
        INTERLEAVE=true \
        NON_PF=true | \
      /usr/bin/bwa mem -K 10000000 -p -v 3 -t 32 -Y ${ref_fasta} /dev/stdin - 2> >(tee ${sample_name}.aligned.unsorted.bwa.stderr.log >&2) | \
      java -Dsamjdk.compression_level=${params.compression_level} -Xmx3G \
      -jar /usr/picard/picard.jar \
      MergeBamAlignment \
        VALIDATION_STRINGENCY=SILENT \
        EXPECTED_ORIENTATIONS=FR \
        ATTRIBUTES_TO_RETAIN=X0 \
        ATTRIBUTES_TO_REMOVE=NM \
        ATTRIBUTES_TO_REMOVE=MD \
        ALIGNED_BAM=/dev/stdin \
        UNMAPPED_BAM=${ubam} \
        OUTPUT=${sample_name}.aligned.unsorted.bam \
        REFERENCE_SEQUENCE=${ref_fasta} \
        PAIRED_RUN=true \
        SORT_ORDER="unsorted" \
        IS_BISULFITE_SEQUENCE=false \
        ALIGNED_READS_ONLY=false \
        CLIP_ADAPTERS=false \
        MAX_RECORDS_IN_RAM=2000000 \
        ADD_MATE_CIGAR=true \
        MAX_INSERTIONS_OR_DELETIONS=-1 \
        PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
        PROGRAM_RECORD_ID="bwamem" \
        PROGRAM_GROUP_VERSION="v0.7.15" \
        PROGRAM_GROUP_COMMAND_LINE="bwa mem -K 100000000 -p -v 3 -t 32 -Y" \
        PROGRAM_GROUP_NAME="bwamem" \
        UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
        ALIGNER_PROPER_PAIR_FLAGS=true \
        UNMAP_CONTAMINANT_READS=true \
        ADD_PG_TAG_TO_READS=false

      grep -m1 "read .* ALT contigs" ${sample_name}.aligned.unsorted.bwa.stderr.log | \
      grep -v "read 0 ALT contigs"

    # else ref_alt is empty or could not be found
    else
      echo "oh no"
      exit 1;
    fi
  """
}

process MarkDuplicates {

  cpus 4

  memory "16 GB"

  container "157538628385.dkr.ecr.us-west-2.amazonaws.com/bfx_docker:latest"

  errorStrategy 'retry'

  input:
  tuple val(sample_name), val(tn), file(unsort_bam) from unsort_bam_ch

  output:
  tuple val(sample_name), val(tn), file("${sample_name}.aligned.unsorted.duplicates_marked.bam") into dupmark_bam_ch

  """
      java -Dsamjdk.compression_level=${params.compression_level} \
      -Xmx14G -jar /usr/picard/picard.jar \
      MarkDuplicates \
        INPUT=${unsort_bam} \
        OUTPUT=${sample_name}.aligned.unsorted.duplicates_marked.bam \
        METRICS_FILE=${sample_name}.dupmetrics.txt \
        VALIDATION_STRINGENCY=SILENT \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
        ASSUME_SORT_ORDER="queryname" \
        CLEAR_DT="false" \
        ADD_PG_TAG_TO_READS=false
  """
}

process SortSam {
  cpus 1

  memory "12 GB"

  container "157538628385.dkr.ecr.us-west-2.amazonaws.com/bfx_docker:latest"

  errorStrategy 'retry'

  input:
  tuple val(sample_name), val(tn), file(mrkdup_bam) from dupmark_bam_ch

  output:
  tuple val(sample_name), val(tn), file("${sample_name}.aligned.duplicate_marked.sorted.bam") into dupmark_sorted_bam_ch

  """
  java -Dsamjdk.compression_level=${params.compression_level} \
  -Xmx8G -jar /usr/picard/picard.jar \
  SortSam \
    INPUT=${mrkdup_bam} \
    OUTPUT=${sample_name}.aligned.duplicate_marked.sorted.bam \
    SORT_ORDER="coordinate" \
    MAX_RECORDS_IN_RAM=300000
  """
}

process FixTags {
  cpus 1

  memory "14 GB"

  container "157538628385.dkr.ecr.us-west-2.amazonaws.com/bfx_docker:latest"

  errorStrategy 'retry'

  input:
  tuple val(sample_name), val(tn), file(dupmark_sorted_bam) from dupmark_sorted_bam_ch
  val ref_fasta

  output:
  tuple val(sample_name), val(tn), file("${sample_name}.aligned.duplicate_marked.sorted.fixed.bam") into dup_sort_fix_bam_ch1
  tuple val(sample_name), val(tn), file("${sample_name}.aligned.duplicate_marked.sorted.fixed.bai") into dup_sort_fix_bam_idx_ch1


  """
  java -XX:ParallelGCThreads=1 \
  -Xmx11G -jar /usr/picard/picard.jar \
  SetNmMdAndUqTags \
    INPUT=${dupmark_sorted_bam} \
    OUTPUT=${sample_name}.aligned.duplicate_marked.sorted.fixed.bam \
    REFERENCE_SEQUENCE=${ref_fasta} \
    CREATE_INDEX=true
  """

}

process BaseQualRecalWES {

  cpus 1

  memory "15 GB"

  container "broadinstitute/gatk:4.0.8.1"

  errorStrategy 'retry'


  input:
  tuple val(sample_name), val(tn), file(dup_sort_fix_bam) from dup_sort_fix_bam_ch1
  tuple val(sample_name), val(tn), file(dup_sort_fix_bam_idx) from dup_sort_fix_bam_idx_ch1
  val ref_fasta
  val known_indels_sites_VCFs
  val group_interval

  output:
  tuple val(sample_name), val(tn), file("${sample_name}.recal_report.csv") into bqsr_report_ch
  tuple val(sample_name), val(tn), file(dup_sort_fix_bam) into dup_sort_fix_bam_ch2
  """
  gatk --java-options -Xmx12G \
  BaseRecalibrator \
    -R ${ref_fasta} \
    -I ${dup_sort_fix_bam} \
    -O ${sample_name}.recal_report.csv \
    --use-original-qualities \
    --known-sites ${known_indels_sites_VCFs[0]} \
    --known-sites ${known_indels_sites_VCFs[1]} \
    --known-sites ${known_indels_sites_VCFs[2]} \
    -L ${group_interval}
  """

}

process ApplyBQSRWES {

  cpus 1

  memory "10 GB"

  container "broadinstitute/gatk:4.0.8.1"

  errorStrategy 'retry'

  input:
  tuple val(sample_name), val(tn), file(dup_sort_fix_bam) from dup_sort_fix_bam_ch2
  tuple val(sample_name), val(tn), file(bqsr_report) from bqsr_report_ch
  val ref_fasta

  output:
  tuple val(sample_name), val(tn), file("${sample_name}.dedup.recal.bam") into recal_bam_ch1, recal_bam_ch2
  tuple val(sample_name), val(tn), file("${sample_name}.dedup.recal.bai") into recal_bam_idx_ch

  """
  gatk --java-options -Xmx7G \
  ApplyBQSR \
    -R ${ref_fasta} \
    -I ${dup_sort_fix_bam} \
    -O ${sample_name}.dedup.recal.bam \
    -bqsr ${bqsr_report} \
    --static-quantized-quals 10 \
    --static-quantized-quals 20 \
    --static-quantized-quals 30 \
    --create-output-bam-index \
    --add-output-sam-program-record
  """
}

process ValidateBam {

  cpus 1

  memory "10 GB"

  container "157538628385.dkr.ecr.us-west-2.amazonaws.com/bfx_docker:latest"

  errorStrategy 'retry'

  input:
  tuple val(sample_name), val(tn), file(recal_bam) from recal_bam_ch1
  val ref_fasta

  output:
  tuple val(sample_name), val(tn), file("${recal_bam}"), file(recal_bam) into recal_bam_ch3
  tuple val(sample_name), val(tn), file("${sample_name}.dedup.recal_bam_validation_report.txt") into bam_valid_report_ch

  """
  java -XX:ParallelGCThreads=1 \
  -Xmx7G -jar /usr/picard/picard.jar \
  ValidateSamFile \
    INPUT=${recal_bam} \
    OUTPUT=${sample_name}.dedup.recal_bam_validation_report.txt \
    REFERENCE_SEQUENCE=${ref_fasta} \
    MODE=VERBOSE \
    IGNORE="MISSING_TAG_NM" \
    MAX_OUTPUT=1000000000 \
    IS_BISULFITE_SEQUENCED=false
  """

}

process CollectHSMetrics {

  cpus 1

  memory "10 GB"

  container "157538628385.dkr.ecr.us-west-2.amazonaws.com/bfx_docker:latest"

  errorStrategy 'retry'

  input:
  tuple val(sample_name), val(tn), file(recal_bam) from recal_bam_ch2
  val ref_fasta
  val wes_coverage_interval_list

  output:
  tuple val(sample_name), val(tn), file("${sample_name}.hsmetrics.txt") into hsmetrics_ch
  tuple val(sample_name), val(tn), file("${sample_name}.pertarget_metrics.txt") into per_target_metrics_ch

  """
  java \
  -jar -Xmx6G /usr/picard/picard.jar \
  CollectHsMetrics \
    INPUT=${recal_bam} \
    OUTPUT=${sample_name}.hsmetrics.txt \
    REFERENCE_SEQUENCE=${ref_fasta} \
    BAIT_INTERVALS=${wes_coverage_interval_list} \
    TARGET_INTERVALS=${wes_coverage_interval_list} \
    PER_TARGET_COVERAGE=${sample_name}.pertarget_metrics.txt \
    VALIDATION_STRINGENCY=SILENT
  """

}


// splitting the processed bams into tumor and normal channels
recal_bam_ch3.branch {
  tumor_bam_ch: it[1] =~ /^tumor/
  normal_bam_ch: it[1] =~ /^normal/
}.set {tn_bam_channel}

// Emits as a value, useful to run parallel Mutect 2 on split intervals
tumor_ch = tn_bam_channel.tumor_bam_ch.first()
normal_ch = tn_bam_channel.normal_bam_ch.first()

// channel from the split intervals, used to run mutect2 on parallel interval chunks
split_intervals = Channel.fromPath(params.split_intervals).splitText(by:20000, file:"scattered.intervals")

gnomad = Channel.value(params.gnomad)
gnomad_idx = Channel.value(params.gnomad_idx)

// replaced the splitintervals step with in built nextflow splitText Channel above

// process SplitIntervals {

//   cpus 1

//   memory "4 GB"

//   disk "40 GB"

//   container "broadinstitute/gatk:4.0.8.1"


//   input:
//   val split_intervals
//   val ref_fasta

//   output:
//   file("interval-files/*.intervals") into intervals_ch

//   """
//   set -e
//   mkdir interval-files
//   gatk --java-options -Xmx3G SplitIntervals \
//       -R ${ref_fasta} \
//       -L ${split_intervals} \
//       -scatter ${params.scatter_count} \
//       -O interval-files
//   """

// }

process Mutect2 {

  cpus 1

  memory "4 GB"

  container "broadinstitute/gatk:4.0.8.1"

  errorStrategy 'retry'

  input:
  val gnomad
  val ref_fasta
  file interval from split_intervals

  tuple val(tumor_name), val(tumor_status), tbam, tbam2 from tumor_ch
  tuple val(normal_name), val(normal_status), nbam, nbam2 from normal_ch

  output:
  val(tumor_name)
  file("${interval}-output.vcf") into dummy_vcf_ch
  file("${interval}-output.vcf") into unfiltered_vcf_ch
  file("${interval}-output.vcf.idx") into unfiltered_vcf_idx_ch

  """
  set -e
  echo ${tbam}
  echo ${tbam2}
  echo ${nbam}
  echo ${nbam2}
  ls
  gatk --java-options "-Xmx3G" Mutect2 \
      -R ${ref_fasta} \
      -I ${tbam} \
      -tumor ${tumor_name} \
      -I ${nbam} \
      -normal ${normal_name} \
      --germline-resource ${gnomad} \
      -L ${interval} \
      -O ${interval}-output.vcf
  """

}

process MergeVCFs {

  cpus 1

  memory "4 GB"

  container "broadinstitute/gatk:4.0.8.1"

  errorStrategy 'retry'

  input:
  val(tumor_name)
  // Not using this properly, using here to block process till all upstream Mutect2 are done
  file(dummy_merge) from dummy_vcf_ch.collectFile(name:"dummy.txt", skip:100000)

  file(unfiltered_vcfs) from unfiltered_vcf_ch.collect()
  file(unfiltered_vcf_idx) from unfiltered_vcf_idx_ch.collect()

  output:
  file "${tumor_name}.dedup.recal.unfiltered.vcf" into unfiltered_m2_vcf_ch

  """
  set -e
  ls *vcf > vcfs_to_merge.txt
  gatk --java-options "-Xmx3G" \
  MergeVcfs \
  -I vcfs_to_merge.txt \
  -O ${tumor_name}.dedup.recal.unfiltered.vcf
  """

}
