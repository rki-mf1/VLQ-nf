process filter_sequences {
  //publishDir "${params.runinfo}/", mode: 'copy', pattern: "filter_sequences.log"

  input:
  tuple val(chunk_id), path(multifasta), path(selected_ids)

  output:
  tuple val(chunk_id), path("filtered_${multifasta}"), emit: filtered_fasta
  path ".command.log", emit:log

  script:
  """
  #!/bin/bash
  echo "_______________________________________________________________"
  echo "Filter sequences from multifasta"

  seqtk subseq $multifasta $selected_ids > filtered_$multifasta

  echo "_______________________________________________________________"
  """

}
