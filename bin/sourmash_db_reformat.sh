#!/bin/bash

for f in *.csv; do
    gzip $f
done
