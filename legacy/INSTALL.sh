#!/usr/bin/env bash

trap "clean_up" SIGINT SIGTERM SIGHUP EXIT


#if [[ $OSTYPE =~ "darwin" ]]; then
# alias realpath='readlink -f'
#fi 

function clean_up {
    echo 'ERROR: I noticed an interruption signal.'
    echo 'Before exiting, let me clean up the mess I did.'
    #rm -r "$RIBAPDIR"
    #mv "$LOG" "$backup"
    echo 'Okay - everything I installed so far is removed.'
    echo "You will find a .log file in $HOME."
    echo 'Bye Bye!'
    echo "$PATH" >> $LOG
}

function log_block {
    echo '' >> "$LOG"
    echo '##################################' >> "$LOG"
    echo '' >> "$LOG"
}

function show_help {
    echo 'Best help ever! yaaay'
}

function test_installation {
    echo "Testing installed software -- stay tuned!"
    log_block
    echo 'TESTING SOFTWARE' >> "$LOG"
    if ! "$RIBAPDIR"/bin/prokka --version >> "$LOG" 2>&1; then
        echo 'Something went wrong with prokka.'
        echo 'Please see if your system meets all requirements and Perl libraries'
        echo 'https://github.com/tseemann/prokka'
        echo ''
        echo 'sudo apt-get install libdatetime-perl libxml-simple-perl libdigest-md5-perl git default-jre bioperl'
        echo 'sudo cpan Bio::Perl'
        echo ''
        echo "For more information, have a look at $LOG"
        #rm -r "$RIBAPDIR"
        exit 1
    fi
    log_block
    if ! "$RIBAPDIR"/bin/roary -h >> "$LOG" 2>&1; then
        echo 'I had an issue with roary.'
        echo 'Please see if your system meets all requirements and Perl libraries'
        echo 'https://github.com/sanger-pathogens/Roary/blob/master/README.md#required-dependencies'
        echo ''
        echo 'and check whether required Perl libraries are installed.'
        echo ''
        echo 'sudo cpanm  Array::Utils Bio::Perl Exception::Class File::Basename File::Copy File::Find::Rule File::Grep File::Path File::Slurper File::Spec File::Temp File::Which FindBin Getopt::Long Graph Graph::Writer::Dot List::Util Log::Log4perl Moose Moose::Role Text::CSV PerlIO::utf8_strict Devel::OverloadInfo Digest::MD5::File'
        echo ''
        echo "For more information, have a look at $LOG"
        #rm -r "$RIBAPDIR"
        exit 1
    fi
    log_block
    if ! "$RIBAPDIR"/bin/glpsol --help >> "$LOG" 2>&1; then
        echo ''
        echo ''
        echo 'Somehow the GLPK package is not working.'
        echo 'Please check whether all system requirements are met.'
        echo 'https://en.wikibooks.org/wiki/GLPK/Linux_OS#Installation_procedure'
        echo ''
        echo "For more information, have a look at $LOG"
        #rm -r "$RIBAPDIR"
        exit 1
    fi

    echo 'It seems like everything went well.'
    echo 'Please copy the following lines to your .bashrc / .profile / .bash_profile:'
    echo 'export PATH=$PATH:'"$RIBAPDIR"'/bin'
    echo 'export PERL5LIB=$PERL5LIB:'"$RIBAPDIR"'/lib'
    echo ''
    echo 'Thank you for using RIBAP! Stay awesome.'
    exit 0 

}

RIBAPDIR="$HOME"/ribap
LOG="$HOME"/ribap_installation.log
SCRIPTDIR=$(dirname $(realpath $0))

while getopts "htp:" opts; do
    case "$opts" in
        h|\?)
            show_help
            exit 0
            ;;
        p)  
            PREFIX="$OPTARG"
            RIBAPDIR="$(cd $(dirname $PREFIX); pwd)"/"$(basename $PREFIX)"/ribap
            ;;
        t)
            test_installation
            exit 0
            ;;

    esac
done

shift $((OPTIND-1))


echo "Welcome to the RIBAP Installation Pipeline."
echo "Please make sure you have a stable internet connection"
echo "and the GNU wget program installed."
echo ""
echo "RIBAP will install all things (except Python3 packages) into: "
echo "$RIBAPDIR"
echo "You will further find a log file of all installation processes"
echo "in your home directory $HOME ."

if [ ! -d "$RIBAPDIR" ]; then
    mkdir -p "$RIBAPDIR"
    export PATH="$RIBAPDIR"/bin:"$RIBAPDIR"/binaries/linux:"$PATH"
fi

echo 'RIBAP Installation Log' > $LOG
log_block

echo 'Installing Python packages for RIBAP, if needed.'
echo 'Make sure that your pip3 binary is connected to your /usr/bin/env python3 variable!'
echo 'RIBAP will install networkx and docopt for you.'
sleep 1

if [ ! $(which pip3) ]; then
    echo 'I was not able to find pip3! Please install pip3 yourself and link it to your python3 interpreter'
    exit 1
fi

if ! python3 -c "import networkx" &>/dev/null; then
    log_block
    pip3 install networkx >> "$LOG"
fi
echo 'networkx is installed'

if ! python3 -c "import docopt" &>/dev/null; then
    log_block
    pip3 install docopt >> "$LOG"
fi    
echo 'docopt is installed'

if ! python3 -c "import Bio" &>/dev/null; then
    log_block
    pip3 install BioPython >> "$LOG"
fi
echo 'BioPython is installed'

echo 'Python packages done.'
echo ''
echo 'Let me install Prokka, roary and GLPK for you.'

#################################################
RIBAPTMP="$RIBAPDIR"/tmp/
cd "$RIBAPDIR" && mkdir -p "$RIBAPTMP" && mkdir -p "$RIBAPDIR"/bin
cd "$RIBAPTMP"
#################################################

PATH="$RIBAPDIR"/bin:"$PATH"

#################################################

cd "$RIBAPTMP"
echo 'Downloading newick utilities'
sleep 1
echo ''
wget http://cegg.unige.ch/pub/newick-utils-1.6-Linux-x86_64-disabled-extra.tar.gz
tar -xf newick-utils-1.6-Linux-x86_64-disabled-extra.tar.gz
cp newick-utils-1.6/src/nw_* "$RIBAPDIR"/bin
echo ''
echo 'Done with newick utilities'
echo ''
sleep 1

#################################################

cd "$RIBAPTMP"
echo 'Downloading cd-hit'
sleep 1
echo ''
wget https://github.com/weizhongli/cdhit/releases/download/V4.8.1/cd-hit-v4.8.1-2019-0228.tar.gz
tar -xf cd-hit-v4.8.1-2019-0228.tar.gz
log_block
cd cd-hit-v4.8.1-2019-0228/ && make openmp=no >> $LOG
cd cd-hit-auxtools/ && make >> $LOG && cd ..
log_block
find . -executable -type f | xargs -I{} cp {} "$RIBAPDIR"/bin
echo ''
echo 'Done with cd-hit'
echo ''
sleep 1

#################################################

cd "$RIBAPTMP"
echo 'Downloading MCL'
sleep 1
echo ''
wget https://micans.org/mcl/src/mcl-14-137.tar.gz
tar -xf mcl-14-137.tar.gz
log_block
cd mcl-14-137 && ./configure --prefix="$RIBAPDIR" >> $LOG
make install >> $LOG
log_block
echo ''
echo 'Done with MCL'
echo ''
sleep 1

#################################################

cd "$RIBAPTMP"
echo 'Downloading BLAST+'
sleep 1
echo ''
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.9.0+-x64-linux.tar.gz
tar -xf ncbi-blast-2.9.0+-x64-linux.tar.gz
cp ncbi-blast-2.9.0+/bin/* "$RIBAPDIR"/bin
echo ''
echo 'Done with BLAST+'
echo ''
sleep 1

#################################################

cd "$RIBAPTMP"
echo 'Downloading GNU parallel'
sleep 1
echo ''
wget ftp://ftp.gnu.org/gnu/parallel/parallel-20190622.tar.bz2
tar -xf parallel-20190622.tar.bz2
cd parallel-20190622
log_block
./configure --prefix="$RIBAPDIR" >> $LOG && make >> $LOG && make install >> $LOG
log_block
echo ''
echo 'Done with GNU Parallel'
echo ''
sleep 1

#################################################

cd "$RIBAPTMP"
echo 'Downloading RAxML'
sleep 1
echo ''
wget https://github.com/stamatak/standard-RAxML/archive/master.zip
unzip master.zip && cd standard-RAxML-master
for FILE in *SSE*.gcc; do
    make -f "$FILE"
    rm *.o
done

# eeeewwww - windows! ;)
find . -executable -type f | grep -v "Windows" | xargs -I{} cp {} "$RIBAPDIR"/bin/

echo ''
echo 'Done with RAxML'
echo ''
sleep 1

#################################################

cd "$RIBAPTMP"
echo 'Downloading MAFFT'
sleep 1
echo ''
wget https://mafft.cbrc.jp/alignment/software/mafft-7.427-without-extensions-src.tgz
tar -xf mafft-7.427-without-extensions-src.tgz && cd mafft-7.427-without-extensions/core/
sed -i -e '/PREFIX = / s:PREFIX = .*:PREFIX = '"$RIBAPDIR"/bin':' Makefile
sed -i -e '/BINDIR = / s:BINDIR = .*:BINDIR = '"$RIBAPDIR"/bin':' Makefile
make clean && make && make install
echo ''
echo 'Done with MAFFT'
echo ''
sleep 1

#################################################

cd "$RIBAPTMP"
echo 'Downloading FastTree'
sleep 1
echo ''
wget http://www.microbesonline.org/fasttree/FastTree
chmod 775 FastTree && cp FastTree "$RIBAPDIR"/bin
ln -s "$RIBAPDIR"/bin/FastTree "$RIBAPDIR"/bin/fasttree
echo ''
echo 'Done with FastTree'
echo ''
sleep 1
#################################################

echo 'Downloading and Installing Prokka'
sleep 1
echo ''
git clone https://github.com/tseemann/prokka.git "$RIBAPTMP"/prokka
rsync -a "$RIBAPTMP"/prokka/* "$RIBAPDIR"
log_block
"$RIBAPDIR"/bin/prokka --setupdb >> "$LOG" 2>&1
if [ $? -eq 0 ]; then
    echo ''
    echo 'Prokka is ready to use.'
else
    echo ''
    echo 'Uh-oh. Something went wrong while the database setup of Prokka!'
    echo "Please have a look at $LOG"
    echo 'https://github.com/tseemann/prokka'
    #rm -rf "$RIBAPDIR"
    exit 1
fi
echo ''

#################################################

cd "$RIBAPTMP"

echo 'Downloading and Installing roary'
sleep 1
echo ''
wget  https://github.com/sanger-pathogens/Roary/tarball/master
tar -xf master
rsync -a sanger-pathogens-Roary-*/* "$RIBAPDIR"
#cd "$RIBAPDIR"/lib/Bio
#find . -iname "*pm" | xargs -n 1 sed -i 's/DIR => getcwd/DIR => $self->output_directory/'
echo ''
echo 'I am done with roary'
echo ''

#################################################

cd "$RIBAPTMP"

echo 'Downloading and Installing GLPK'
sleep 1
echo ''
wget ftp://ftp.gnu.org/gnu/glpk/glpk-4.65.tar.gz
tar -xf glpk-4.65.tar.gz
cd glpk-4.65
log_block
./configure --prefix="$RIBAPDIR" >> "$LOG"
log_block
make --jobs=4 >> "$LOG"
log_block
make install >> "$LOG"
log_block
make clean >> "$LOG"
echo ''
echo 'GLPK is good to go!'
echo ''
sleep 1

#   "CD-HIT: a fast program for clustering and comparing large sets of protein or nucleotide sequences", Weizhong Li & Adam Godzik. Bioinformatics, (2006) 22:1658-1659
#   "CD-HIT: accelerated for clustering the next generation sequencing data", Limin Fu, Beifang Niu, Zhengwei Zhu, Sitao Wu & Weizhong Li. Bioinformatics, (2012) 28:3150-3152

cd "$RIBAPTMP"
git clone https://github.com/soedinglab/MMseqs2.git
cd MMseqs2
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
make
make install 
cp bin/* "$RIBAPDIR"/bin/


cp -r "$SCRIPTDIR"/scripts "$SCRIPTDIR"/ribap.sh "$SCRIPTDIR"/web "$RIBAPDIR"/bin

echo 'Cleaning up...'
cd "$RIBAPDIR"
rm -rf "$RIBAPTMP"

test_installation

