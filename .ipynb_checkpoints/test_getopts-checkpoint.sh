#!/bin/bash

echo
echo -e "\t#################################### "
echo -e "\t##                                ## "
echo -e "\t##    Imputation / QC Pipeline    ## "
echo -e "\t##          Torkamani Lab         ## "
echo -e "\t##                                ## "
echo -e "\t##         Author: Raquel Dias    ## "
echo -e "\t##                 Shaun Chen     ## "
echo -e "\t##  Last modified: 12/24/19       ## "
echo -e "\t##                                ## "
echo -e "\t#################################### "
echo 

### preallocate variables
required=
optional_hasarg=
got_opt_arg=false
optional_noarg=false
has_default="DEFAULT"

### initialize (clear) opt parsing variables
OPTIND=
OPTARG=
opt=
myinput=$1


while getopts "hv:b:cd:" opt; do

    case "$opt" in
        h)
            echo "\
            usage:
            ------
            getoptsdemo [ -h ] -a REQ [ -b OPTHAS ] [ -c ] [ -d OPTDEF ]
            description:
            ------------
            This demonstrates parsing arguments with getopts
            optional arguments:
            -------------------
            -h          Print this help message and exit.
            -run        CONFIRM        Optional opt that does not take an argument.
            -start
            -end
            -ref
            OPTDEF   If given, sets value of has_default to OPTDEF.
                        Otherwise, the value of has_default defaults to
                        DEFAULT.
            "
            exit 0
            ;;
        v)
            myinput=${OPTARG}
            ;;
        b)
            optional_hasarg=${OPTARG}
            got_opt_arg=true
            ;;
        c)
            optional_noarg=true
            ;;
        d)
            has_default=${OPTARG}
            ;;
        ?)
            echo "Error: did not recognize option, ${OPTARG}."
            echo "Please try -h for help."
            exit 1
            ;;
    esac
done


if [[ "$myinput" == "" ]]; then
    echo -e "ERROR: Required VCF_PATH not set."
    echo -e "Please try -h for help."
    echo
    exit 1
fi
echo "myinput: $myinput"


if [[ "$got_opt_arg" == true ]]; then
    if [[ "$optional_hasarg" == "" ]]; then
        echo "-b option set but no argument given"
        echo "Please try -h for help."
        exit 1
    else
        echo "optional_hasarg: $optional_hasarg"
    fi
else
    echo "optional_hasarg not set"
fi

if [[ "$optional_noarg" == true ]]; then
    echo "optional_noarg set"
else
    echo "optional_noarg not set"
fi

echo "has_default: $has_default"