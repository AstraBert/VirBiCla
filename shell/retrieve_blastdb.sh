#!/bin/bash

# Usage function
usage() {
    echo "Usage: retrieve_blastdb.sh -db DBNAME -o FASTAFILE

    REQUIRED ARGUMENTS:
    -db | --database: Provide the name of the database to download
    -o | --output: Provide the path to the output fasta file
  Input retrieve_blastdb.sh -h,--help to show the help message"
    exit 1
}


db="None"
output="None"

while [ $# -gt 0 ]; do
    case "$1" in
    -h | --help)
        usage
        ;;
    -db | --database) ######## <- new samples argument
        db="$2"
        shift 2
        ;;
    -o | --output)
        output="$2"
        shift 2
        ;;
    *)
        echo "Unknown option: $1"
        exit 1
        ;;
    esac
done


if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
    usage
fi

if [ "${db}" = "None"  ] && [ ! "$1" = "-h" ] && [ ! "$1" = "--help" ]; then
    echo "Missing required argument: DBNAME with the name of the db to download"
    usage
fi

if [ "${output}" = "None"  ] && [ ! "$1" = "-h" ] && [ ! "$1" = "--help" ]; then
    echo "Missing required argument: OUTPUT with the path to the output fasta file"
    usage
fi

update_blastdb.pl --decompress $db
blastdbcmd -entry all -db $db -out $output
rm ${db}.n*
rm *.tar.gz.md5
