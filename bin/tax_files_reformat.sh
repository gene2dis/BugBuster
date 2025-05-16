#!/bin/bash

wget $1
mkdir tax_files
tar -xzf *.tar.gz
mv nodes.dmp tax_files
mv names.dmp tax_files
rm *.dmp gc.prt readme.txt taxdump.tar.gz
