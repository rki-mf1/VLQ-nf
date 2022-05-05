process kallisto_prediction_single {
  publishDir "${params.output}/${sample.baseName}", mode: 'copy', pattern: "kallisto_out/*.{tsv,json,h5}"

  maxForks 1

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
  echo $sample
  kallisto quant -t $params.kallisto_threads -i $ref_index -b $params.bootstrap -o kallisto_out/ --single -l $params.fragment_length -s $params.fragment_length_sd $sample

  output_predictions.py  kallisto_out/abundance.tsv --metadata $final_selection -m $params.min_ab -o kallisto_out/predictions.tsv
  """
}
