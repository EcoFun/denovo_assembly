#!/usr/bin/env bash

__author__="Ludovic Duvaux"
__maintainer__="Ludovic Duvaux"
__license__="GPL_v3"

### usage:
usage="SYNOPSIS
$(basename "$0") [-cDemh] [-d data_dir] [-r suffix_res] [-t n_cores] -k 30,50,70 -p sample1,sample2 -s sample_suffix

FUNCTION
prepare merged and connected reads pairs for assembly using 'abyss-pe'.
The assembler can be different but work well with 'ABYSS'.

OPTIONS
    -c              pass Konnector analysis
    -d              raw data directory 
    -D              Debbug mode (i.e. dry run)
    -e              extend read using konnector
                    [BEWARE: not functional, miss arguments hard coded for bloom filter]
    -k              kmer list for konnector
                    can be specified as  '-k 30,50,70'
                    or '-k \"30 50 70\"'
                    or '-k \"$(echo seq 30 20 70)\"' and so on.
    -m              pass abyss-mergepair analysis
    -p              prefixes of the samples to be processed.
                    can be specified as '-p F,G,H')
                    or '-p \"F G H\"'
                    or '-p \"$(echo {F..H})\"' and so on.
    -r              single suffix for results files.
    -s              single suffix corresponding to the sample names to be
    -t              number of cores for fastqc and konnector [1]
    -h              show this help and exit"


### print out passed command line
echo "Submitted command line:"
echo `basename $0` $*
echo ""
echo ""

### command line options
# read the options
TEMP="$(getopt -o cd:eDk:mp:r:s:t:h -n 'abyss_prepare_data.sh' -- "$@")"

# if invalid option, print usage and exit
if [ $? -ne 0 ] ; then  # $? is the status of the last I don't remember what
    echo ""
    echo "$usage"
    exit
fi

# extract options and their arguments into variables.
eval set -- "$TEMP"
unset TEMP

# set initial values for some parameters and extract option values
pass_kon="false"
dir="./"
pass_merg="false"
DEB="false"
konextd="false"
while true ; do
    case "$1" in
        -h)
            echo "$usage"
            exit
            ;;
        -c)
            pass_kon="true"
            shift
            ;;
        -d)
            case "$2" in
                "") shift 2 ;;
                *) dir=$2 ; shift 2 ;;
            esac ;;
        -D)
            DEB="true"
            shift
            ;;
        -e)
            konextd="true"
            shift
            ;;
        -k)
            case "$2" in
                "") shift 2 ;;
                *) ks=$2 ; shift 2 ;;
            esac ;;
        -m)
            pass_merg="true"
            shift
            ;;
        -p)
            case "$2" in
                "") shift 2 ;;
                *) samp=$2 ; shift 2 ;;
            esac ;;
        -r)
            case "$2" in
                "") shift 2 ;;
                *) sufres=$2 ; shift 2 ;;
            esac ;;
        -s)
            case "$2" in
                "") shift 2 ;;
                *) suf=$2 ; shift 2 ;;
            esac ;;
        -t)
            case "$2" in
                "") ncpu="1" ; shift 2 ;;
                *) ncpu=$2 ; shift 2 ;;
            esac ;;
        --) shift ; break ;;
        *) echo "Internal error!" 
        exit 1 ;;
    esac
done

# test if required parameters are specified
if  [ -z "$ks" ] ; then
    echo "### ERROR: Please specify option '-k'"
    printf "\n\n"
    echo "$usage"
    exit
elif [ -z "$samp" ] ; then
    echo "### ERROR: Please specify option '-p'"
    printf "\n\n"
    echo "$usage"
    exit
elif [ -z "$suf" ] ; then
    echo "### ERROR: Please specify option '-s'"
    printf "\n\n"
    echo "$usage"
    exit
elif [ -z "$sufres" ] ; then
    echo "### ERROR: Please specify option '-r'"
    printf "\n\n"
    echo "$usage"
    exit
fi

# prepare options
ks=$(echo $ks | sed 's/,/ /g')
samp=$(echo $samp | sed 's/,/ /g')

if [ "$DEB" = "true" ]; then
    echo "option -c is" $pass_kon
    echo "option -d is" $dir
    echo "option -D is" $DEB
    echo "option -e is" $konextd
    echo "option -k is" $ks
    echo "option -m is" $pass_merg
    echo "option -p is" ${samp[@]}
    echo "option -r is" $sufres
    echo "option -s is" $suf
    echo "option -t is" $ncpu
    echo ""

    echo 'print sample'
    for i in $samp; do
        echo $i
    done
    echo ""
fi

####### script
for k in $ks; do
    for sa in ${samp[@]}; do
        echo "#### Process sample $sa for kmer=$k ####"
        dmerg=01_merged_fastq/${sa}
        mkdir -p $dmerg
        pref1=${sa}_${sufres}

        # 1) mergepairs
        if [ "$pass_merg" = "false" ]; then
            echo "# 1) mergepairs"
            date
            if [ ! -f "$dmerg/${pref1}_merged.fastq.gz" ] || [ ! -f "$dmerg/${pref1}_reads_1.fastq.gz" ] || [ ! -f "$dmerg/${pref1}_reads_2.fastq.gz" ] ; then  # if files does not exist

                echo "## Run abyss-mergepairs for $sa ##"
                echo abyss-mergepairs -o $dmerg/$pref1 $dir/${sa}$suf
                if [ "$DEB" = "false" ]; then
                    time (abyss-mergepairs -o $dmerg/$pref1 $dir/${sa}$suf > /dev/null) 2>&1
                fi
                
                # compress new files
                echo ""
                echo "find $dmerg -name \"*fastq\" -exec gzip {} \;"
                if [ "$DEB" = "false" ]; then
                    time find $dmerg -name "*fastq" -exec gzip {} \;
                fi
            else # pass
                echo "## Abyss-mergepairs already done for $sa ##"
            fi
            
            # 2) fastqc on mergepairs files
            echo ""
            echo "# 2) fastqc on mergepairs files"
            date
            dfastqc=$dmerg/fastqc_mergepairs
            mkdir -p $dfastqc
            if [ ! -f "$dfastqc/${pref1}_merged_fastqc.html" ] || [ ! -f "$dfastqc/${pref1}_reads_1_fastqc.html" ] || [ ! -f "$dfastqc/${pref1}_reads_2_fastqc.html" ] ; then

                echo "## Run fastqc on mergepairs files ##"
                echo fastqc --outdir=$dfastqc -t $ncpu $dmerg/$pref1*fastq.gz
                if [ "$DEB" = "false" ]; then
                    time (
                        fastqc --outdir=$dfastqc -t $ncpu $dmerg/$pref1*fastq.gz
                    ) 2>&1
                fi
            else
                echo "## fastqc on mergepairs already done for $sa ##"
            fi
        else
            echo "# 1) Pass mergepairs for $sa"
            echo "# 2) Pass fastqc on mergepair files for $sa"
        fi
        
        #~# 3) konnector
        dkonn=02_konnector_fastq/${sa}/$k
        mkdir -p $dkonn
        
        if [ "$pass_kon" = "false" ]; then
            echo ""
            echo "# 3) konnector"
            date
            pref2=${sa}_${sufres}_k${k}

            if [ ! -f "$dkonn/${pref2}_pseudoreads.fa.gz" ] || [ ! -f "$dkonn/${pref2}_reads_1.fq.gz" ] || [ ! -f "$dkonn/${pref2}_reads_1.fq.gz" ]; then 
                echo "## Run konnector (kmer=$k) ##"
                if [ "$konextd" = true ] ; then
                    echo konnector -j $ncpu -k $k -D 100 -E -o $dkonn/$pref2 $dir/${sa}$suf
                    if [ "$DEB" = "false" ]; then
                        time (
                            konnector -j $ncpu -k $k -D 100 -E -o $dkonn/$pref2 $dir/${sa}$suf
                        ) 2>&1
                    fi
                    echo ""
                    echo ""
                elif  [ "$konextd" = false ] ; then
                    echo konnector -j $ncpu -k $k -o $dkonn/$pref2 $dir/${sa}$suf
                    if [ "$DEB" = "false" ]; then
                        time (
                            konnector -j $ncpu -k $k -o $dkonn/$pref2 $dir/${sa}$suf
                        ) 2>&1
                    fi
                    echo ""
                else
                    echo "ERROR: bad value for konextd"
                    exit 
                fi

                # compress new files
                echo "find $dkonn -name \"*fastq\" -exec gzip {} \;"
                if [ "$konextd" = true ] ; then
                    if [ "$DEB" = "false" ]; then
                        time find $dkonn -name "*f[qa]" -exec gzip {} \;
                    fi
                fi
            else
                echo "## konnector for $sa already done with kmer=$k ##"
            fi
        else
            echo "# 3) Pass konnector for $sa"
        fi
        echo "#### Sample $sa processed for kmer=$k :) ####"
        echo ""
        echo ""
    done
done
