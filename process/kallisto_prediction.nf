process kallisto_prediction {
  maxForks 1
  //publishDir "${params.runinfo}/", mode: 'copy', pattern: ".command.log", saveAs: {filename -> "kallisto_prediction.log"}
  publishDir "${params.output}/${fastq.simpleName}/", mode: 'copy', pattern: "kallisto_out/*.{tsv,json,h5}"

  input:
  tuple path(fastq), path(ref_index), path(final_selection)

  output:
  path "kallisto_out/predictions.tsv", emit: prediction_ch
  path "kallisto_out/*"
  path ".command.log", emit: log

  script:
  """
  #!/bin/bash
  echo "_______________________________________________________________"
  echo "Run kallisto"
  kallisto quant -i $ref_index -o kallisto_out/ --single -l $params.fragment_length -s $params.fragment_length_sd -t $params.kallisto_threads $fastq

  output_predictions.py  kallisto_out/abundance.tsv --metadata $final_selection -m $params.min_ab --voc $params.vocs -o kallisto_out/predictions.tsv
  """

}
