process filter_by_metadata {
  publishDir "${params.runinfo}/", mode: 'copy', pattern: '.command.log', saveAs: {filename -> "filter_by_metadata.log"}

  input:
  path desh_df
  path gisaid_df

  output:
  path "selection_by_metadata.csv", emit: selection_df
  path ".command.log"

  script:
  if ("${desh_df.simpleName}" == 'DUMMY')
    """
    filter_references.py -gisaid $gisaid_df --country $params.country --startdate $params.startdate --enddate ${params.enddate} --min_len $params.min_len -k $params.k --seed $params.seed
    """

  else if ("${gisaid_df.simpleName}" == 'DUMMY')
    """
    filter_references.py -desh $desh_df --country $params.country --startdate $params.startdate --enddate $params.enddate --min_len $params.min_len -k $params.k --seed $params.seed
    """

  else
    """
    filter_references.py -gisaid ${gisaid_df} -desh ${desh_df} --country ${params.country} --startdate ${params.startdate} --enddate ${params.enddate} --min_len ${params.min_len} -k ${params.k} --seed {$params.seed}
    """

}
