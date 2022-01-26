/*
in order to decompress .xz and .tar.xz, xz-utils must be installed
decompress .fasta.xz and .tar.xz to .fasta.gz
*/
process unxz {
  publishDir "${params.runinfo}/", mode: 'copy', pattern: '.command.log', saveAs: {in_filename -> "unxz.log"}

  input:
  path in_file

  output:
  path "cp_${in_file.baseName}.gz", emit: unxz_file
  path ".command.log"

  script:
    """
    cp $in_file cp_${in_file}
    unxz -d cp_${in_file}
    gzip cp_${in_file.baseName}
    # rm cp_${in_file.baseName}
    """

}
