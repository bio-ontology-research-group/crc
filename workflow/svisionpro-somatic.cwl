#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool
label: SVision-Pro somatic structural variant caller
doc: |
  SVision-Pro is a structural variant (SV) caller designed for long-read sequencing data,
  with special optimization for Oxford Nanopore Technology reads.
  This tool definition is for somatic variant calling using tumor-normal paired samples.

baseCommand: [SVision-Pro]

requirements:
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: genomicslab/svision-pro:latest

inputs:
  reference_genome:
    type: File
    format: edam:format_1929  # FASTA format
    inputBinding:
      prefix: --ref
    doc: Reference genome in FASTA format
  
  reference_index:
    type: File
    doc: Index file for reference genome (.fai)
  
  tumor_bam:
    type: File
    format: edam:format_2572  # BAM format
    inputBinding:
      prefix: --tumor-bam
    doc: Tumor sample BAM file
    secondaryFiles:
      - .bai
  
  tumor_bam_index:
    type: File
    doc: Index file for tumor BAM
  
  normal_bam:
    type: File
    format: edam:format_2572  # BAM format
    inputBinding:
      prefix: --normal-bam
    doc: Normal/control sample BAM file
    secondaryFiles:
      - .bai
  
  normal_bam_index:
    type: File
    doc: Index file for normal BAM
  
  output_prefix:
    type: string
    inputBinding:
      prefix: --output-prefix
    doc: Prefix for output files
  
  min_sv_size:
    type: int?
    default: 50
    inputBinding:
      prefix: --min-sv-size
    doc: Minimum structural variant size to detect
  
  min_support_reads:
    type: int?
    default: 3
    inputBinding:
      prefix: --min-support-reads
    doc: Minimum number of supporting reads required for SV calling
  
  threads:
    type: int?
    default: 16
    inputBinding:
      prefix: --threads
    doc: Number of threads to use

arguments:
  - prefix: --mode
    valueFrom: "somatic"
  - prefix: --tech
    valueFrom: "ont"  # Specify Oxford Nanopore Technology
  - prefix: --output-dir
    valueFrom: "./output"

outputs:
  sv_vcf:
    type: File
    outputBinding:
      glob: output/$(inputs.output_prefix).somatic.vcf
    doc: VCF file containing somatic structural variants
  
  summary_report:
    type: File
    outputBinding:
      glob: output/$(inputs.output_prefix).somatic.summary.html
    doc: HTML summary report of somatic structural variants

stderr: svision_pro_somatic.log
stdout: svision_pro_somatic.out

# Note: You might need to adjust parameters based on the specific version of SVision-Pro being used
# Always refer to the SVision-Pro documentation for the most up-to-date parameter information