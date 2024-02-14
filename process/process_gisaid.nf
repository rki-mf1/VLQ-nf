process process_gisaid {
  publishDir "${params.runinfo}/", mode: 'copy', pattern: ".command.log", saveAs: {filename -> "process_gisaid.log"}

  input:
  path meta

  output:
  path "processed_gisaid_metadata.tsv", emit: gisaid_processed
  path ".command.log"


  script:
  """
  #!/bin/bash
  echo "_______________________________________________________________"
  echo "--- Process GISAID data"

  process_gisaid.py -meta $meta

  """

}
