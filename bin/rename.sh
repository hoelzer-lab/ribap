#!/bin/bash
set -e 

function rename_fasta_for_prokka {
    for FILE in $@; do
        
        if [ ! -f "$FILE" ]; then
            continue
        fi
        i=0
        IFS=$'\n'
        BN="${FILE%.*}"
        BN="${BN##*/}"
        BN="${BN//./_}"
        #EXT="${FILE##*.}"
        EXT=fasta
        PREFIX="${BN:0:25}"
        RENAMED="$OUTDIR"/"$BN"_RENAMED."$EXT"
        RENAMED=$(echo "$RENAMED" | sed 's/,/_/g')
        #cp "$FILE" "$OUTDIR"/fasta
        cp "$FILE" "$RENAMED"
        for entry in $(grep '>' $RENAMED); do
            i=$((i+1))
            entry=$(echo $entry | sed 's/|/\\|/g')
            sed -r 's@'"$entry"'@>contig'"$i"'@' "$RENAMED" > "$RENAMED"_TMP
            mv "$RENAMED"_TMP "$RENAMED"
        done
    done
}

OUTDIR='./'
rename_fasta_for_prokka "$@"
