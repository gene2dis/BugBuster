#!/bin/sh

wget $1
mkdir k2_ref_db
tar -xzf *.tar.gz -C k2_ref_db && rm *.tar.gz
