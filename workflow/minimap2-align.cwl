#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool
label: Minimap2 long read alignment tool
doc: |
  Minimap2 is a versatile sequence alignment program that aligns DNA or mRNA sequences against
  a large reference database. It is particularly optimized for Oxford Nanopore and PacBio reads.

baseCommand: [minimap2]

requirements:
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/minimap2:2.24--h7132678_1

inputs:
  reference:
    type: File
    format: edam:format_1929  # FASTA format
    inputBinding:
      position: 1
    doc: Reference genome in FASTA format
  
  reads:
    type: File
    format: edam:format_3179  # FASTQ or BAM format
    inputBinding:
      position: 2
    doc: ONT reads in FASTQ or BAM format
  
  threads:
    type: int?
    default: 16
    inputBinding:
      prefix: -t
    doc: Number of threads to use
  
  output_prefix:
    type: string
    doc: Prefix for output files

arguments:
  - position: 0
    prefix: -ax
    valueFrom: "map-ont"  # Preset for Oxford Nanopore reads
  - position: 0
    prefix: -L
    valueFrom: "--secondary=no"  # Skip secondary alignments
  - position: 0
    prefix: --MD
    valueFrom: "y"  # Output the MD tag
  - position: 3
    shellQuote: false
    valueFrom: |
      | samtools sort -@ $(inputs.threads) -o $(inputs.output_prefix).bam -

outputs:
  aligned_bam:
    type: File
    format: edam:format_2572  # BAM format
    outputBinding:
      glob: $(inputs.output_prefix).bam
    secondaryFiles:
      - .bai
  
  aligned_bam_index:
    type: File
    outputBinding:
      glob: $(inputs.output_prefix).bam.bai

stdout: minimap2_alignment.log

# After alignment, we need to index the BAM file
# Note that samtools is required in the same environment
# or can be called through a separate tool definition