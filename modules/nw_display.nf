/*Comment section: */

process nw_display {
  label 'newick_utils'
  publishDir "${params.output}/fasttree", mode: 'copy', pattern: "*.svg" 

  input: 
    file(nwk)

  output:
    file("*.svg")

  script:
    """
    for NWK in *.nwk; do
      nw_display -v 25 -i "font-size:6" -l "font-size:12;font-family:helvetica;font-style:italic" -Il -w 750 -b "opacity:0" -s \${NWK} > \$(basename \${NWK} .nwk).svg
    done
    """
}

