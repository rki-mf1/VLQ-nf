process filter_sequences {


  input:
  tuple val(chunk_id), path(multifasta), path(selected_ids)

  output:
  tuple val(chunk_id), path("filtered_${multifasta}"), emit: filtered_fasta
  //path ".command.log", emit:log
  //path "sample_counter.txt", emit: sample_counter

  script:
  """
  #!/bin/bash
  echo "_______________________________________________________________"
  echo "Filter sequences from multifasta"

  seqtk subseq $multifasta $selected_ids > filtered_$multifasta

  #counter=\$(grep '>' filtered_$multifasta | wc -l)
  #echo \$counter >> sample_counter.txt
  #grep '>' filtered_$multifasta >> sample_counter.txt

  echo "_______________________________________________________________"
  """

}
