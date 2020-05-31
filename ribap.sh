#!/bin/bash
set -e 

trap "clean_tmp" SIGINT SIGHUP SIGTERM

function clean_tmp {
    echo -e 'Received an interrupt!'
    if [ $KEEP -eq 0 ]; then
        echo -e 'Cleaning temporary files.'
        rm -r "$OUTDIR"/tsv
        rm -r "$OUTDIR"/ilp
        rm -r "$OUTDIR"/mmseq2
        rm -r "$OUTDIR"/fasta
        rm -r "$OUTDIR"/roary
        echo -e 'Done'
    else
        echo -e 'Temporary files are kept, since you set the -k flag!'
    fi
    echo -e 'RIBAP is exiting now.'
}


function check_dependencies {
    [[ $(python3 --version &>/dev/null; echo $?) -gt 0 ]] && echo -e "\e[1;91mERROR: \e[0mMissing dependency - \e[38;5;32mpython3" && return 1
    [[ $(prokka -v &>/dev/null; echo $?) -gt 0 ]] && echo -e "\e[1;91mERROR: \e[0mMissing dependency - \e[38;5;32mprokka" && return 1
    [[ $(roary --version &>/dev/null; echo $?) -gt 0 ]] && echo -e "\e[1;91mERROR: \e[0mMissing dependency - \e[38;5;32mroary" && return 1
    [[ $(which cd-hit-est &>/dev/null; echo $?) -gt 0 ]] && echo -e "\e[1;91mERROR: \e[0mMissing dependency - \e[38;5;32mcd-hit-est" && return 1
    [[ $(mcl -h &>/dev/null; echo $?) -gt 0 ]] && echo -e "\e[1;91mERROR: \e[0mMissing dependency - \e[38;5;32mmcl" && return 1
    [[ $(blastp -version &>/dev/null; echo $?) -gt 0 ]] && echo -e "\e[1;91mERROR: \e[0mMissing dependency - \e[38;5;32mblastp" && return 1
    [[ $(makeblastdb -version &>/dev/null; echo $?) -gt 0 ]] && echo -e "\e[1;91mERROR: \e[0mMissing dependency - \e[38;5;32mmakeblastdb" && return 1
    [[ $(parallel --version &>/dev/null; echo $?) -gt 0 ]] && echo -e "\e[1;91mERROR: \e[0mMissing dependency - \e[38;5;32mparallel" && return 1
    [[ $(raxmlHPC-PTHREADS-SSE3 -v &>/dev/null; echo $?) -gt 0 ]] && echo -e "\e[1;91mERROR: \e[0mMissing dependency - \e[38;5;32mraxmlHPC" && return 1
    [[ $(mafft --version &>/dev/null; echo $?) -gt 0 ]] && echo -e "\e[1;91mERROR: \e[0mMissing dependency - \e[38;5;32mmafft" && return 1
    [[ $(which fasttree &>/dev/null; echo $?) -gt 0 ]] && echo -e "\e[1;91mERROR: \e[0mMissing dependency - \e[38;5;32mFastTree" && return 1
    [[ $(nw_display -h &>/dev/null; echo $?) -gt 0 ]] && echo -e "\e[1;91mERROR: \e[0mMissing dependency - \e[38;5;32mnw_display // newick_utils" && return 1
    [[ $(glpsol --version &>/dev/null; echo $?) -gt 0 ]] && echo -e "\e[1;91mERROR: \e[0mMissing dependency - \e[38;5;32mglpk" && return 1

    return 0
}

function show_help {
    echo ""
    echo "RIBAP - Roary ILP Bacterial Annotation Pipeline"
    echo ""
    echo "Annotate your protein sequences with Prokka and determine a pan genome with Roary."
    echo "This genome is refined with the usage of ILPs that solve the best matching for each pairwise"
    echo "strain blastp comparison."
    echo ""
    echo "Usage:"
    echo "ribap [-h] [-v] [-k] [-t] [-p CPU] [-o OUTPUT] <fasta1> <fasta2> ..."
}

function log_message {
    >&2 echo -e "\e[38;5;29m[LOG] \e[0m$1"
}

function print_citation {
    echo -en "          "; for i in {28..33} {33..28} {28..33}; do echo -en "\e[38;5;${i}m#\e[0m" ; done ; echo
    echo -e "             \e[38;5;29mPlease cite :)\e[0m"
    echo -en "          "; for i in {28..33} {33..28} {28..33}; do echo -en "\e[38;5;${i}m#\e[0m" ; done ; echo
}

function rename_fasta_for_prokka {
    for FILE in $@; do
        
        if [ ! -f "$FILE" ]; then
            continue
        fi
        i=0
        IFS=$'\n'
        BN="${FILE%.*}"
        BN="${BN##*/}"
        #EXT="${FILE##*.}"
        EXT=fasta
        PREFIX="${BN:0:25}"
        RENAMED="$OUTDIR"/fasta/"$BN"_RENAMED."$EXT"
        RENAMED=$(echo "$RENAMED" | sed 's/,/_/g')
        cp "$FILE" "$OUTDIR"/fasta
        cp "$FILE" "$RENAMED"
        for entry in $(grep '>' $RENAMED); do
            i=$((i+1))
            entry=$(echo $entry | sed 's/|/\\|/g')
            #sed -r -i 's@'"$entry"'@>'"$PREFIX"'_'"$i"'@' "$RENAMED"
            #sed -r -i 's@'"$entry"'@>contig'"$i"'@' "$RENAMED"
            sed -r 's@'"$entry"'@>contig'"$i"'@' "$RENAMED" > "$RENAMED"_TMP
            mv "$RENAMED"_TMP "$RENAMED"
        done
    done
}

check_dependencies || (echo "Exiting..." && exit 1)

COMMAND="$0 $@"
DIR=$(dirname $(which "$0"))
OPTIND=1
CPUS=1
OUTDIR="$PWD"
KEEP=0

while getopts "hktvo:p:" opts; do
    case "$opts" in
        h|\?)
            show_help
            exit 0
            ;;
        k)  KEEP=1
            ;;
        t)  BUILDTREE=1
            ;;
        v)  VERBOSE=1
            ;;
        o)  OUTDIR="$(cd $(dirname $OPTARG); pwd)"/"$(basename $OPTARG)"
            ;;
        p)  CPUS="$OPTARG"
            ;;
    esac
done


shift $((OPTIND-1))

NUMSTRAINS="$#"

if [ $# == 0 ]; then
    show_help
    >&2 echo ""
    >&2 echo ""
    >&2 echo -e "\e[1;91mERROR: \e[0mIt seems that you forgot to enter input files. The core gene set of 0 strains is most likely empty."
    >&2 echo "Exiting"
    exit 1
fi

if [ "$VERBOSE" ]; then
    echo ''
    echo -en "          \e[38;5;28mW\e[38;5;29mE\e[38;5;30mL\e[38;5;31mC\e[38;5;32mO\e[38;5;33mM\e[38;5;32mE \e[38;5;31mT\e[38;5;30mO \e[38;5;29mR\e[38;5;28mI\e[38;5;29mB\e[38;5;30mA\e[38;5;32mP\e[38;5;33m!"; echo
    echo -en "     "; for i in {33..28} {28..33} {33..28} {28..33} {33..28}; do echo -en "\e[38;5;${i}m#\e[0m" ; done ; echo
    echo ''
    log_message 'This is your command:'
    echo ''
    log_message "$COMMAND"
    log_message "Output files will be written into $OUTDIR."
    echo ''
    echo '---------------------------------'
    echo ''
    log_message 'Annotating your genes based on genome fastas you provided.'
    echo ''
fi


## Prokka
## conda install -c conda-forge -c bioconda prokka
print_citation
echo ''
echo '---------------------------------'
echo 'RIBAP uses Prokka for gene annotation. Please cite this awesome work:'
echo 'Seemann T.; Prokka: rapid prokaryotic genome annotation; Bioinformatics 2014 Jul 15;30(14):2068-9'
echo 'PMID:24642063'
echo '---------------------------------'
echo ''


mkdir -p "$OUTDIR"
mkdir -p "$OUTDIR"/fasta
mkdir -p "$OUTDIR"/log

TMPTEST=$IFS
rename_fasta_for_prokka "$@"
IFS=$TMPTEST

for STRAIN in "$OUTDIR"/fasta/*RENAMED*; do
    BN=$(basename "$STRAIN" .fasta)
    #prokka --outdir "$OUTDIR"/prokka/"$BN" --prefix "$BN" --genus Chlamydia --usegenus --cpus "$CPUS" --force "$STRAIN" 2>"$OUTDIR"/log/ribap_"$BN"_prokka.log #--gram neg ??
    prokka --outdir "$OUTDIR"/prokka/"$BN" --prefix "$BN" --cpus "$CPUS" --force "$STRAIN" 2>"$OUTDIR"/log/ribap_"$BN"_prokka.log
done

if [ "$VERBOSE" ]; then
    log_message 'Gene annotation with Prokka is done.'
    log_message 'I will now continue with roary in order to'
    log_message 'calculate a first version of the pan genome.'
    echo ''
fi

for ANNO in "$OUTDIR"/prokka/*/*.gff; do
    BN=$(basename "$ANNO" .gff)
    echo $(grep -vm1 '^#' $ANNO | awk '{print $9}' | cut -d'=' -f2 | cut -d'_' -f1),"$BN"
done > "$OUTDIR"/strain_ids.txt

## Roary w/ identity 60, 70, 80, 90, 95
## installed via Miniconda

print_citation
echo ''
echo '---------------------------------'
echo 'RIBAP uses roary for a scaffold core gene set. If you use and like RIBAP, give roary a citation as well:'
echo 'Andrew J. Page, Carla A. Cummins et al., "Roary: Rapid large-scale prokaryote pan genome analysis", Bioinformatics, 2015'
echo 'doi:10.1093/bioinformatics/btv421'
echo '---------------------------------'
echo ''

mkdir -p "$OUTDIR"/roary && cd "$OUTDIR"/roary
for IDENT in 60 70 80 90 95; do
    OUT="$OUTDIR"/roary/"$IDENT"
    roary -e --mafft -p "$CPUS" -v -i "$IDENT" -r "$OUTDIR"/prokka/*/*.gff -f "$OUT" &> "$OUTDIR"/log/ribap_roary_"$IDENT".log
done
cd - >/dev/null

if [ "$VERBOSE" ]; then
    log_message 'Pan genome calculation by roary done.'
    log_message 'I will now perform a individual blastp analysis similar to roary, but without clustering sequences.'
    log_message 'This is needed for the ILPs and might take a while... stay tuned!'
    log_message ''
    log_message 'Using the following blastp filter:'
    log_message 'Sequence Similarity between hits: 60%'
    log_message 'Alignment hit length: 40% of query length'
    echo ''
fi

## ILPs (on perc=60 roary output)
#mkdir -p "$OUTDIR"/blast
#ALL="$OUTDIR"/blast/all.blast

mkdir -p "$OUTDIR"/mmseq2
ALL="$OUTDIR"/mmseq2/all_proteins.fa
cat "$OUTDIR"/prokka/*/*.faa > "$ALL"
MMSEQDB="$OUTDIR"/mmseq2/mmseq2.db

MMDIR=/home/co68mol/miniconda3/envs/mmseq
PATH=$PATH:/home/co68mol/miniconda3/envs/mmseq/bin/

#makeblastdb -in "$ALL" -dbtype prot -parse_seqids >/dev/null
#blastp -task blastp -num_threads "$CPUS" -query "$ALL" -db "$ALL" -evalue 1e-10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send evalue bitscore slen" | awk '{if($3>60 && $4>($9*0.4)){print $0}}' > "$OUTDIR"/blast/all_vs_all.blast 2>/dev/null
echo "Creating DB"
mmseqs createdb "$ALL" "$MMSEQDB"
echo "Creating Index"
mmseqs createindex "$MMSEQDB" "$OUTDIR"/mmseq2/tmp
echo "Starting MMSeqs2 Search"
mmseqs search "$MMSEQDB" "$MMSEQDB" "${MMSEQDB%.*}_result" "$OUTDIR"/mmseq2/tmp -a
echo "Converting results"
mmseqs convertalis "$MMSEQDB" "$MMSEQDB" "${MMSEQDB%.*}_result" "${MMSEQDB%.*}_result.csv"

#create single .tsv from blastout
MMSEQ="${MMSEQDB%.*}_result.csv"
#BLAST="$OUTDIR"/blast/all_vs_all.blast
STRAIN_IDS="$OUTDIR"/strain_ids.txt # this needs to be automatically generated after the prokka runs # Done ;)

if [ "$VERBOSE" ]; then
    log_message 'All done with the blastp command. Now let me prepare'
    log_message 'the pairwise ILPs in order to refine the cluster generated by roary.'
    echo ''
fi

mkdir -p "$OUTDIR"/tsv
ulimit -n 4096
"$DIR"/bin/mmseq2tsv.py "$MMSEQ" "$STRAIN_IDS" "$OUTDIR"/tsv
ulimit -n 1024
mkdir -p "$OUTDIR"/ilp
#pairwise ILP comparison
print_citation
echo ''
echo '---------------------------------'
echo 'The ILPs of RIBAP are based on the FFDCJ model by Martinez et al. and modified by us. To honour their nice craftmanship of theory and modelling, please cite:'
echo 'Martinez, Fábio V et al., 2015 “On the Family-Free DCJ Distance and Similarity.” Algorithms for Molecular Biology'
echo 'PMID: 25859276'
echo '---------------------------------'
echo ''

parallel -j "$CPUS" ''"$DIR"'/bin/ILP.py --max --indel {}' ::: "$OUTDIR"/tsv/*tsv #2>/dev/null

if [ "$VERBOSE" ]; then
    log_message 'ILPs are prepared. Solving all of them might'
    log_message 'take a while.'
    echo ''
fi

#print_citation
echo ''
echo 'ILPs are solved using the freely-available ILP solver GLPK.'
echo 'Copyright (C) 2000-2017 Andrew Makhorin, Department for Applied'
echo 'Informatics, Moscow Aviation Institute, Moscow, Russia. All rights'
echo 'reserved. E-mail: <mao@gnu.org>.'
echo ''
echo 'This program has ABSOLUTELY NO WARRANTY.'
echo ''
echo 'This program is free software; you may re-distribute it under the terms'
echo 'of the GNU General Public License version 3 or later.'
echo ''

parallel -j "$CPUS" 'glpsol --lp {} --mipgap 0.01 --memlim 16834 --tmlim 120 -o {.}.sol >/dev/null' ::: "$OUTDIR"/ilp/*ilp 2>/dev/null


# hardcoded stuff is hardcoded
#CPLEX="/home/kevin/program/CPLEX_Studio128/cplex/bin/x86-64_linux/cplex"
#for ILP in "$OUTDIR"/ilp/*ilp; do
#    "$CPLEX" -c "read $ILP" "lp" "set timelimit 3600" "set mip tolerances mipgap 0.01" "set threads $CPUS" "set workmem 16384" "opt" "set output writelevel 3" "write ${ILP%.*}.sol" "y" "quit" 2>/dev/null
#done

#for SOL in "$OUTDIR"/ilp/*sol; do
#    python3 "$DIR"/bin/cplexsol_to_simple.py "$SOL"
#done

#for SIMPLE in "$OUTDIR"/ilp/*simple.1; do
#    #grep -E '^x_' "$SIMPLE" > "${SIMPLE%.*}"
#    awk '/^x_/ {print $2,$1}' "$SIMPLE" > "${SIMPLE%.*}"
#done

#function awk_parallel_magic {
#    awk '$2 ~ /x_A.*_B/ && $4 == 1 {print}' $1 > ${1%.*}.simple
#}
#export -f awk_parallel_magic
#parallel -j "$CPUS" 'awk_parallel_magic {}' ::: "$OUTDIR"/ilp/*sol 2>/dev/null


#for SOL in for SOL in "$OUTDIR"/ilp/*sol; do
#    sed -E -i '/x_A[^[:space:]]+$/ N;s/\n//g' "$SOL"
#done

function sed_parallel_remove_newline_in_sol {
    sed -E -i '/x_A[^[:space:]]+$/ N;s/\n//g' "$1"
}

export -f sed_parallel_remove_newline_in_sol
parallel -j "$CPUS" 'sed_parallel_remove_newline_in_sol {}' ::: "$OUTDIR"/ilp/*sol

for tsv in "$OUTDIR"/tsv/*tsv; do
    BN=$(basename $tsv .tsv)
    for SOL in "$OUTDIR"/ilp/"$BN"*sol; do
        #sed -E -i '/x_A[^[:space:]]+$/ N;s/\n//g' "$SOL"
        awk '$2 ~ /x_A.*_/ && $4 == 1 {print}' "$SOL"
    done > "$OUTDIR"/ilp/"$BN".ilp.simple
done

if [ "$VERBOSE" ]; then
    log_message 'ILPs are solved.'
    log_message 'Now I will combine the roary and ILP results for you.'
    log_message 'Stay tuned!'
    echo ''
fi

for IDENT in 60 70 80 90 95; do
    python3 "$DIR"/bin/combine_roary_ilp.py "$OUTDIR"/strain_ids.txt "$OUTDIR"/roary/"$IDENT"/gene_presence_absence.csv "$OUTDIR"/ilp/ "$OUTDIR"/holy_python_ribap_"$IDENT".csv "$IDENT" > "$OUTDIR"/ribap_roary"$IDENT"_summary.txt
done

if [ "$VERBOSE" ]; then
    log_message 'All roary clusters have been refined with the ILPs'
    log_message 'Now I will do multiple sequence alignments and trees'
    log_message 'for all RIBAP groups. Almost there.'
    echo ''
fi

mkdir -p "$OUTDIR"/msa
python3 "$DIR"/bin/create_msa_tree.py "$OUTDIR" "$OUTDIR"/holy_python_ribap_95.csv

#The Newick Utilities: High-throughput Phylogenetic tree Processing in the UNIX Shell
#Thomas Junier and Evgeny M. Zdobnov
#Bioinformatics 2010 26:1669-1670
#http://bioinformatics.oxfordjournals.org/content/26/13/1669.full
#doi:10.1093/bioinformatics/btq243

parallel -j "$CPUS" 'mafft {} > {.}_mafft.aln 2>/dev/null' ::: "$OUTDIR"/msa/group*faa 2>/dev/null
parallel -j "$CPUS" 'fasttree {} > {.}_tree.nwk 2>/dev/null' ::: "$OUTDIR"/msa/group*aln 2>/dev/null
parallel -j "$CPUS" 'nw_display -v 25 -i "font-size:6" -l "font-size:12;font-family:helvetica;font-style:italic" -Il -w 750 -b "opacity:0" -s {} > {.}.svg 2>/dev/null' ::: "$OUTDIR"/msa/*nwk 2>/dev/null
#parallel -j "$CPUS" 'nw_display {} > {.}.ascii' ::: "$OUTDIR"/msa/*nwk

python3 "$DIR"/bin/concat_coreMSA.py "$OUTDIR"/msa "$NUMSTRAINS"

mkdir -p "$OUTDIR"/coreGenome/
mv "$OUTDIR"/msa/coreGenome_mafft.aln "$OUTDIR"/coreGenome/

if [ "$VERBOSE" ]; then
    log_message 'You fill find a nice multiple sequence alignment of all core genes here:'
    log_message ""$OUTDIR"/coreGenome/coreGenome_mafft.aln"
    echo ''
fi

if [ "$BUILDTREE" ]; then
    if [ "$VERBOSE" ]; then
        echo '---------------------------------'
        log_message 'I reckon you used the -t parameter - thus, you would like to have'
        log_message 'a tree based on the core genome. Sure thing, I will use RAxML for this.'
        log_message 'Be aware, this will take a lot of time.'
        echo '---------------------------------'
    fi
    raxmlHPC-PTHREADS-SSE3 -T "$CPUS" -f a -x 1234 -p 1234 -s "$OUTDIR"/coreGenome/coreGenome_mafft.aln -n aa -m PROTGAMMAWAG -N 100 -w "$OUTDIR"/coreGenome/
fi

rm "$OUTDIR"/msa/*faa
mkdir -p "$OUTDIR"/web
mkdir -p "$OUTDIR"/tree
wget https://www.rna.uni-jena.de/supplements/ribap/web.tar.gz
tar zxvf web.tar.gz
gunzip -r web
cp -r "$DIR"/web/* "$OUTDIR"/web
cp "$OUTDIR"/msa/*svg "$OUTDIR"/tree

if [ "$VERBOSE" ]; then
    log_message 'I am basically done. Let me prepare the final output files for you.'
    echo ''
fi

python3 "$DIR"/bin/generate_html.py "$OUTDIR" > "$OUTDIR"/web/ribap.html
cd "$OUTDIR"/
ln -s web/ribap.html .

if [ "$VERBOSE" ]; then
    log_message 'You will find the resulting .html file here:'
    log_message "$OUTDIR/ribap.html"
    echo ''
fi

mv "$OUTDIR"/holy_python_ribap_"$IDENT".csv "$OUTDIR"/ribap_gene_absence_presence.csv
mv "$OUTDIR"/ribap_roary"$IDENT"_summary.txt "$OUTDIR"/ribap_summary.txt
rm "$OUTDIR"/ribap_roary*.txt "$OUTDIR"/holy_*

if [ "$VERBOSE" ]; then
    log_message 'All done!'
    echo '---------------------------------'
fi

if [ "$KEEP" -eq 0 ]; then
    if [ "$VERBOSE" ]; then
        log_message 'Let me clean up all the mess I did to get here.'
    fi
    clean_tmp
fi

echo ''
echo ''
print_citation
echo ''
echo '---------------------------------'
echo 'Thank you for using RIBAP'
echo 'We hope we were able to help you with your study.'
echo 'If you are using RIBAP for results you would like to publish'
echo 'please be a kind soul and leave a small citation for us in your'
echo 'manuscript.'
echo ''
echo 'Lamkiewicz, Kevin et al., 2020, "RIBAP: An all-in-one annotation and core gene clustering pipeline for bacterial genomes", JOURNAL'
echo 'PMID: XXXXXXXXXX'
echo ''
echo 'This program is free software; you may re-distribute it under the terms'
echo 'of the GNU General Public License version 3 or later.'
echo '---------------------------------'
echo ''
echo 'Stay awesome, dude! :)'
echo ''

cat "$OUTDIR"/ribap_summary.txt
exit 0
