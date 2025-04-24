#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: Workflow
label: SVision-Pro workflow for structural variant detection in cancer samples
doc: |
  A workflow that uses SVision-Pro to identify structural variants in tumor/normal paired samples
  from Oxford Nanopore Technology (ONT) long-read sequencing data.

requirements:
  - class: ScatterFeatureRequirement
  - class: SubworkflowFeatureRequirement
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement
  - class: MultipleInputFeatureRequirement

inputs:
  reference_genome:
    type: File
    format: edam:format_1929  # FASTA format
    doc: Reference genome in FASTA format
  
  reference_genome_index:
    type: File?
    doc: Index file for reference genome (.fai)
  
  tumor_reads:
    type: File
    format: edam:format_3179  # FASTQ or BAM
    doc: Tumor sample ONT reads in FASTQ or BAM format
  
  normal_reads:
    type: File
    format: edam:format_3179  # FASTQ or BAM
    doc: Normal/control sample ONT reads in FASTQ or BAM format
  
  output_prefix:
    type: string
    doc: Prefix for output files
  
  threads:
    type: int?
    default: 16
    doc: Number of threads to use for alignment and SV calling

  min_sv_size:
    type: int?
    default: 50
    doc: Minimum structural variant size to detect
  
  min_support_reads:
    type: int?
    default: 3
    doc: Minimum number of supporting reads required for SV calling

outputs:
  tumor_bam:
    type: File
    outputSource: align_tumor/aligned_bam
  
  normal_bam:
    type: File
    outputSource: align_normal/aligned_bam
  
  sv_vcf:
    type: File
    outputSource: run_svisionpro_somatic/sv_vcf
  
  sv_report:
    type: File
    outputSource: run_svisionpro_somatic/summary_report

steps:
  # Step 1: Index reference genome if not provided
  index_reference:
    run: samtools-faidx.cwl
    in:
      reference_genome: reference_genome
      fai_exists: reference_genome_index
    out: [reference_index]
    when: $(inputs.fai_exists === null)

  # Step 2: Align tumor sample reads to reference
  align_tumor:
    run: minimap2-align.cwl
    in:
      reference: reference_genome
      reads: tumor_reads
      threads: threads
      output_prefix:
        valueFrom: $(inputs.output_prefix + ".tumor")
    out: [aligned_bam, aligned_bam_index]

  # Step 3: Align normal sample reads to reference
  align_normal:
    run: minimap2-align.cwl
    in:
      reference: reference_genome
      reads: normal_reads
      threads: threads
      output_prefix:
        valueFrom: $(inputs.output_prefix + ".normal")
    out: [aligned_bam, aligned_bam_index]

  # Step 4: Run SVision-Pro for somatic variant detection
  run_svisionpro_somatic:
    run: svisionpro-somatic.cwl
    in:
      reference_genome: reference_genome
      reference_index: 
        source: [index_reference/reference_index, reference_genome_index]
        pickValue: first_non_null
      tumor_bam: align_tumor/aligned_bam
      tumor_bam_index: align_tumor/aligned_bam_index
      normal_bam: align_normal/aligned_bam
      normal_bam_index: align_normal/aligned_bam_index
      min_sv_size: min_sv_size
      min_support_reads: min_support_reads
      threads: threads
      output_prefix: output_prefix
    out: [sv_vcf, summary_report]

# Below are the tool definitions for each step in the workflow