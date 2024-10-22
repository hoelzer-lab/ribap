/*Comment section: */

process combine_roary_ilp {
  label 'python3'
  publishDir "${params.output}/05-combine", mode: 'copy', pattern: "*.csv" 
  publishDir "${params.output}/05-combine", mode: 'copy', pattern: "*.txt" 

  input: 
    tuple val(ident), file(roary), file(strain_ids), file(prokka_gff)
    file(solved_ilps)
    file(script)

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

    # copy the py script to execute it here with the extended command. 
    # However this might fail if baseDir can not be reached from the container.
    # If so, download the script from GitHub as a fallback.
    #if cp "$baseDir/bin/combine_roary_ilp.py" .; then
    #  echo "Script copied successfully."
    #else
    #  echo "Script could not be copied. Download from GitHub."
    #  git clone https://github.com/hoelzer-lab/ribap.git ribap-tmp
    #  cp ribap-tmp/bin/combine_roary_ilp.py .
    #  rm -r ribap-tmp
    #fi

    # setrecursionlimit see: https://github.com/hoelzer-lab/ribap/issues/66
    python -c "import sys;sys.setrecursionlimit(${params.set_recursion_limit});exec(open('combine_roary_ilp.py').read())" ${strain_ids} ${ident}/gene_presence_absence.csv solved/ holy_python_ribap_"${ident}".csv ${ident} > ribap_roary"${ident}"_summary.txt
    """
}

