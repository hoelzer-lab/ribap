/*Comment section: */

process combine_roary_ilp {
  label 'python3'
  publishDir "${params.output}/05-combine", mode: 'copy', pattern: "*.csv" 
  publishDir "${params.output}/05-combine", mode: 'copy', pattern: "*.txt" 

  input: 
    tuple val(ident), file(roary), file(strain_ids), file(prokka_gff)
    file(solved_ilps)

  output:
    tuple val(ident), file("holy*.csv")
    tuple val(ident), file("ribap*.csv")
    tuple val(ident), file("*.txt")

  script:
    """
    mkdir solved
    cp *.simple solved/

    mkdir prokka
    cp *.gff prokka/

    combine_roary_ilp.py ${strain_ids} ${ident}/gene_presence_absence.csv solved/ holy_python_ribap_"${ident}".csv ${ident} > ribap_roary"${ident}"_summary.txt
    """
}

