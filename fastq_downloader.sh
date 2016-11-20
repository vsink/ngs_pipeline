#!/usr/bin/env bash
counter=1

if [ -z "$1" ];
then
echo "Fastq accession file list is not found!"
acc_list=10
else
acc_list=1
acc_file=$1

fi

if [ -z "$2" ];
then
echo "By default files will be downloaded as SRA"
ftype=sra
else
if  [ $2 == "sra" ]; 
then
ftype=sra
else
ftype=fastq
fi
fi


function sra_download {

while read access_nbr

do
        three_letter=${access_nbr:0:3}
        six_letter=${access_nbr:0:6}
wget "ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/${three_letter}/${six_letter}/${access_nbr}/${access_nbr}.sra"

done < $acc_file
}

function fastq_download {

while read access_nbr

do
fastq-dump $access_nbr
printf "$counter\t$access_nbr\tOK\n"
((counter++))
done < $acc_file
}


if  [[ $acc_list == 1 && $ftype == "sra" ]];
then
sra_download
fi


if  [[ $acc_list == 1 && $ftype == "fastq" ]];
then
fastq_download
fi
 
