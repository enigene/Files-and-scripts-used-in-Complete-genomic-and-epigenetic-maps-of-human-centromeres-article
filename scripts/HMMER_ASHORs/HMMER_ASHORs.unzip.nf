#!/usr/bin/env nextflow

// Alpha satellite annotation of long DNA FASTA sequences (version for uncompressed input)

params.input_fasta='./fasta/*.fa'
params.HMM="./HMMs/AS-HORs-hmmer3.0-290621.hmm"
params.output="./HMMER_output"
hmmertblout2bed="./scripts/hmmertblout2bed.awk"
add_length_gap="./scripts/add_length_gap.R"
params.score_to_length_threshold=0.7

Channel.fromPath(params.input_fasta)
  .map{file->tuple(file.baseName,
                   file)}
  .set{input_fasta_ch}

hmm = file("${params.HMM}").getBaseName()

process hmmer {
//  errorStrategy 'ignore'
  input:
    tuple id, file('input.fasta') from input_fasta_ch

  output:
    set id, file("${id}-vs-${hmm}.tbl.out") into hmmout2bed

  publishDir "${params.output}", mode: 'link', pattern: '*.{gz}'
  module 'HMMER/3.1b2-foss-2018a'
  scratch true
  cpus 32
  script:
  """
    nhmmer --notextw --noali --tblout /dev/stdout -o /dev/null ${params.HMM} input.fasta > ${id}-vs-${hmm}.tbl.out
  """
}

process bedops {
  errorStrategy 'ignore'
  input:
    tuple id,
          file("${id}-vs-${hmm}.tbl.out") from hmmout2bed

  output:
    set file("${id}-vs-${hmm}.bed"), // into publishDir,
        file("${id}-vs-${hmm}.tbl.lengap.tsv") into publishDir

  publishDir "${params.output}/bed", mode: 'link', pattern: '*.{bed,tsv}'
  module 'BEDOPS'
  module 'R/3.5.3'
  scratch true
  shell:
  '''
    awk -v th=!{params.score_to_length_threshold} -f !{hmmertblout2bed} "!{id}-vs-!{hmm}.tbl.out" | sed "s/ /_SPACE_/g" | sort-bed - | bedmap --max-element --fraction-both 0.3 - | awk '!a[\$0]++' | bedmap --max-element --fraction-both 0.3 - | awk '!a[\$0]++' | sed "s/_SPACE_/ /g" > "!{id}-vs-!{hmm}.tbl.out.bed"
    awk -v FS="\t" -v OFS="\t" '{$5=sprintf("%d",$5);print}' "!{id}-vs-!{hmm}.tbl.out.bed" > "!{id}-vs-!{hmm}.bed"
    Rscript !{add_length_gap} "!{id}-vs-!{hmm}.tbl.out.bed" .
  '''
}
