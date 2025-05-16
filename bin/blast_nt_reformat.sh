#!/bin/bash

max_num=`curl --list-only https://ftp.ncbi.nlm.nih.gov/blast/db/ | grep -Po '\"nt.[0-9]+?\.tar.gz' | grep -Po '[0-9]+' | tail -1`

seq -w 000 $max_num | parallel wget ${1}/nt.{}.tar.gz
ls | grep 'tar.gz' | parallel tar -xzf
ls | grep 'tar.gz' | parallel rm -f

mkdir blast_nt_db
mv nt* blast_nt_db
