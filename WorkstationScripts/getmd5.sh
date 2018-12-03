#!/bin/bash

# Copyright (c) 2018
# Amy L. Olex, Virginia Commonwealth University
# alolex at vcu.edu

## Usage Function
usage()
{
cat << EOF
usage: $0 options
This script is used to calculate the MD5SUM of all files in a given directory.

REQUIRED ARGUMENTS:
   -f      Input file 

-f INPUT DIRECTORY:
Path to directory where files are located.

EXAMPLE USAGE:
   >> ~/getmd5.sh -f /path/to/folder/

EOF
}

## Parsing input arguments

FILE=
while getopts “hf:” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         f)
             FILE=$OPTARG
             ;;
         ?)
             usage
             exit
             ;;
     esac
done

## Handle empty arguments

if [[ -z $FILE ]]; then usage; exit 1; fi



cd $FILE

for file in `ls`
do
   md5sum $file > $file.md5
done
