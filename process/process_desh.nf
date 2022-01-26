process process_desh {
  publishDir "${params.runinfo}/", mode: 'copy', pattern: '.command.log', saveAs: {filename -> "process_desh.log"}

  input:
  tuple path(meta), path(lineage)
  path fasta

  output:
  path "processed_desh_metadata.tsv", emit: desh_processed
  path ".command.log"

  script:
  """
  cp ${meta} cp_${meta}
  cp ${lineage} cp_${lineage}
  cp ${fasta} cp_${fasta}
  unxz -d cp_${meta}
  unxz -d cp_${lineage}
  gzip -d cp_${fasta}

  process_desh.py -meta cp_${meta.baseName} -lineage cp_${lineage.baseName} -fasta cp_${fasta.baseName}

  rm cp_${fasta.baseName} cp_${lineage.baseName} cp_${meta.baseName} 
  """

}
