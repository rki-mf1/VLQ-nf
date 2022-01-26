process kallisto_prediction {
  publishDir "${params.runinfo}/", mode: 'copy', pattern: ".command.log", saveAs: {filename -> "kallisto_prediction.log"}
  publishDir "${params.output}/${fastq.simpleName}/", mode: 'copy', pattern: "kallisto_out/*"

  input:
  path fastq
  path ref_index
  path final_selection

  output:
  path "kallisto_out/predictions.tsv", emit: prediction_ch
  path "kallisto_out/*"
  path ".command.log"

  script:
  """
  #!/bin/bash

  kallisto quant -i $ref_index -o kallisto_out/ --single -l 200 -s 20 -t 20 $fastq

  output_predictions.py  kallisto_out/abundance.tsv --metadata $final_selection -m $params.min_ab --voc $params.vocs -o kallisto_out/predictions.tsv 


  """

}
