process kallisto_prediction_paired {
  publishDir "${params.output}/${sample_name}/", mode: 'copy', pattern: "kallisto_out/*.{tsv,json,h5}"

  maxForks 1

  input:
  tuple path(sample_1), path(sample_2), path(ref_index), path(final_selection)
  val sample_name

  output:
  path "kallisto_out/predictions.tsv", emit: prediction_ch
  path "kallisto_out/*"
  path ".command.log", emit: log

  script:

  """
  #!/bin/bash
  echo "_______________________________________________________________"
  echo "Run kallisto"
  kallisto quant -t $params.kallisto_threads -i $ref_index -b $params.bootstrap -o kallisto_out/ $sample_1 $sample_2

  output_predictions.py  kallisto_out/abundance.tsv --metadata $final_selection -m $params.min_ab -o kallisto_out/predictions.tsv
  """
}
