process process_desh {
  publishDir "${params.runinfo}/", mode: 'copy', pattern: ".command.log", saveAs: {filename -> "process_desh.log"}

  input:
  tuple path(meta), path(lineage)
  path fasta

  output:
  path "processed_desh_metadata.tsv", emit: desh_processed
  path ".command.log"

  script:
  """
  #!/bin/bash
  echo "_______________________________________________________________"
  echo "--- Process DESH data"

  cp ${meta} cp_${meta}
  cp ${lineage} cp_${lineage}
  unxz -d cp_${meta}
  unxz -d cp_${lineage}

  process_desh.py -meta cp_${meta.baseName} -lineage cp_${lineage.baseName} -fasta $fasta

  rm cp_${lineage.baseName} cp_${meta.baseName}
  """

}
