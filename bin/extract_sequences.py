#!/usr/bin/env python3
# extract_sequences.py <proteins.faa> <hit list>: writes the proteins listed in
# the hit list (tab-separated: read id, ARG, predicted class, best hit) to
# <bin>_ARGs.faa, the header carrying the ARG and its class.
# Copied unchanged (apart from this header) from /usr/bin/extract_sequences.py
# in quay.io/ffuentessantander/biopython_utils:1.0, ARG_FASTA_FORMATTER's former
# image, which was the only copy (roadmap U6, Phase 10).
from Bio import SeqIO
import pandas as pd
import sys

# Generar vector con los ID que hay que buscar

df = pd.read_csv(sys.argv[2], sep='\t', header=None)

seq_id = df[0].tolist()
ARG_gene = df[1].tolist()
predicted_ARG_class = df[2].tolist()
match = df[3].tolist()

# Leer fasta
with open(sys.argv[1], "r") as inF:
    for record in SeqIO.parse(inF, "fasta"):
        # Checkear si los ID estan en el fasta
        name = record.description
        id = str(name.split(" ")[0])
        file_name = str(sys.argv[1]).replace("_proteins.faa","")
        new_file_name = str(sys.argv[1]).replace("_proteins.faa","_ARGs.faa")
        if id in seq_id:
            # cosas varias para el header y que sera util mas adelante
            index_number = seq_id.index(id)
            ARG_gene_find = ARG_gene[index_number]
            predicted_ARG_class_find = predicted_ARG_class[index_number]
#            match_find = match[index_number]
            with open(new_file_name, 'a') as archivo:
                archivo.write(f">%s_from_%s_gene_%s_ARG_class_%s\n" % (id, file_name, ARG_gene_find, predicted_ARG_class_find))  # Escribir el header
                archivo.write(f"{record.seq}\n")  # Escribir la secuencia
