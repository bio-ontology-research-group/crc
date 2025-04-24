#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool
label: Samtools faidx index generator
doc: |
  Tool to generate FASTA index files using samtools faidx

baseCommand: [samtools, faidx]

requirements:
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/samtools:1.15--h1170115_1

inputs:
  reference_genome:
    type: File
    format: edam:format_1929  # FASTA format
    inputBinding:
      position: 1
    doc: Reference genome in FASTA format
  
  fai_exists:
    type: File?
    doc: Optional existing index file to avoid re-indexing if already available

outputs:
  reference_index:
    type: File
    outputBinding:
      glob: $(inputs.reference_genome.basename).fai
    doc: Index file for the reference genome

# This tool will need to be run conditionally in the workflow if an index is not provided