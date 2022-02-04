process process_gisaid {
  publishDir "${params.runinfo}/", mode: 'copy', pattern: '.command.log', saveAs: {filename -> "process_gisaid.log"}

  input:
  path meta
  path map_tsv

  output:
  path "processed_gisaid_metadata.tsv", emit: gisaid_processed
  path ".command.log"


  script:
  if ("${map_tsv.simpleName}" != "DUMMY")
    """
    #!/bin/bash
    echo "-------------------\nProcess GISAID data\n--------------------"

    tar --exclude='readme.txt' -xhf ${meta}

    process_gisaid.py -meta metadata.tsv -epi ${map_tsv}

    rm metadata.tsv
    """

  else
    """
    #!/bin/bash
    echo "-------------------\nProcess GISAID data\n--------------------"

    tar --exclude='readme.txt' -xhf ${meta}

    process_gisaid.py -meta metadata.tsv

    rm metadata.tsv
    """

}
