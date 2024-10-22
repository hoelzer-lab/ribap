/*Comment section: */

process generate_html {
  label 'python3'
  publishDir "${params.output}/", mode: 'copy', pattern: "web" 

  input:
    tuple val(ident), file(holy_table)
    file(roary)
    file(individual_annotations)
    file(tree_svg) 

  output:
    file("web")

  script:
    """
    # Try to copy the file
    if cp "$baseDir/data/web.tar.gz" .; then
      echo "File copied successfully."
    else
      # the command might fail when singularity is used, see: https://github.com/hoelzer-lab/ribap/issues/67
      echo "File copy failed. Attempting to download the file..."
      if wget --no-check-certificate https://osf.io/fcyjn/download -O web.tar.gz; then
        echo "File downloaded successfully."
      else
        echo "File download failed."
      fi
    fi
    tar zxvf web.tar.gz
    
    mkdir tree
    cp *.svg tree/

    # clarify github issue 11 if this is correct and intended behavior
    #cat ribap_individual_annotation*.csv > ribap_individual_annotation.csv
    cp ribap_individual_annotation95.csv ribap_individual_annotation.csv

    generate_html.py . > web/ribap.html
    mv tree web/

    """
}

