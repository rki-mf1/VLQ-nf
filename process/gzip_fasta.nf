/*
* In order to decompress .xz and .tar.xz, xz-utils must be installed
* Idea: run process on condition that infile != *.fasta.gz, else give infile as output
*/
process gzip_fasta {
  storeDir "${params.databases}/INPUT/"

  input:
  path in_file

  output:
  path "*.fasta.gz", emit: gz_fasta
  path ".command.log", emit: log


  script:
    """
    #!/bin/bash

    if [[ ${in_file} == *.tar* ]]
    then
      tar --exclude='readme.txt' -xhf ${in_file}
      gzip *.fasta

    elif [[ ${in_file} == *.xz ]]
    then
      cp $in_file cp_${in_file}
      unxz -d cp_${in_file}
      gzip cp_${in_file.baseName}

    elif [[ ${in_file} == *.fasta ]]
    then
      cp $in_file cp_${in_file}
      gzip cp_${in_file}

    fi

    """

}
