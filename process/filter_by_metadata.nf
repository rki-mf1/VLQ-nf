process filter_by_metadata {
  publishDir "${params.runinfo}/", mode: 'copy', pattern: ".command.log", saveAs: {filename -> "filter_by_metadata.log"}

  
  input:
  path gisaid_df

  output:
  path "selection_by_metadata.csv", emit: selection_df
  path ".command.log"

  script:
  """
  #!/bin/bash
  echo "_______________________________________________________________"
  echo "--- Filter samples by metadata"

  filter_by_metadata.py -gisaid $gisaid_df --continent $params.continent --country $params.country --startdate $params.startdate --enddate $params.enddate --min_len $params.min_len -k $params.k --seed $params.seed
  """
}
