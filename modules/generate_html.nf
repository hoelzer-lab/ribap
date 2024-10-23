/*Comment section: */

process generate_html {
  label 'python3'
  publishDir "${params.output}/", mode: 'copy', pattern: "web" 

  input:
    tuple val(ident), file(holy_table)
    file(roary)
    file(individual_annotations)
    file(tree_svg) 
    file(web_dir)

  output:
    file("web")

  script:
    """
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

