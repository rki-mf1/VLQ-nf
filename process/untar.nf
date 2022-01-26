/*
in order to decompress .xz and .tar.xz, xz-utils must be installed
decompress .fasta.xz and .tar.xz to .fasta.gz
*/
process untar {
  publishDir "${params.runinfo}/", mode: 'copy', pattern: '.command.log', saveAs: {in_filename -> "untar.log"}

  input:
  path in_file

  output:
  path "*.gz", emit: untar_file
  path ".command.log"

  script:
    """
    #cp $in_file cp_${in_file}
    # mkdir untar_out
    tar --exclude='readme.txt' -xhf ${in_file}
    gzip *.fasta
    #rm cp_${in_file}
    """

}
