#!/bin/bash

wget -O checkm2_database.tar.gz ${1}

tar -xvzf checkm2_database.tar.gz
rm checkm2_database.tar.gz
rm CONTENTS.json
find . -type 'f' -name '*.dmnd' | xargs -I {} mv {} .
rm -rf CheckM2_database
