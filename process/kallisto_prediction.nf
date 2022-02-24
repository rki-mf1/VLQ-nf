process kallisto_prediction {
  publishDir "${params.output}/${sample.baseName}", mode: 'copy', pattern: "kallisto_out/*.{tsv,json,h5}"


  input:
  tuple path(sample), path(ref_index), path(final_selection)

  output:
  path "kallisto_out/predictions.tsv", emit: prediction_ch
  path "kallisto_out/*"
  path ".command.log", emit: log

  script:
  """
  #!/bin/bash
  echo "_______________________________________________________________"
  echo "Run kallisto"
  kallisto quant -i $ref_index -o kallisto_out/ --single -l $params.fragment_length -s $params.fragment_length_sd -t $params.kallisto_threads $sample

  output_predictions.py  kallisto_out/abundance.tsv --metadata $final_selection -m $params.min_ab -o kallisto_out/predictions.tsv
  """

}
