process filter_sequences {
  publishDir "${params.runinfo}/", mode: 'copy', pattern: "filter_sequences.log"

  input:
  tuple val(chunk_id), path(multifasta), path(selected_ids)

  output:
  tuple val(chunk_id), path("filtered_${multifasta}"), emit: filtered_fasta
  path ".command.log"

  script:
  """
  #!/bin/bash
  echo "--------------------\nFilter sequences from multifasta\n--------------------"

  seqtk subseq $multifasta $selected_ids > filtered_$multifasta
  """

}
