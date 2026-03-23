library("seqinr")
library("tidyverse")

path_list <- list.files(pattern="*mmseq", full.names = F)

for (path in path_list){
  assign(path, list.files(path, full.names = TRUE) %>%
    map(function(x) {
      print(x)
      read_tsv(x)
    }
    ) %>%
    reduce(rbind) %>%
    mutate(Seq_names = paste(sample_bin,target,sep = "__")) %>%
      group_by(Seq_names) %>%
      filter(evalue == min(evalue)) %>%
      filter(qcov == max(qcov)) %>%
      filter(pident == max(pident)) %>%
      filter(tcov == max(tcov)))
  
  assign(paste0("Final_",path), get(path)[!duplicated(get(path)[,c('tseq')]),] %>%
           mutate(header = paste(Seq_names, query, sep = "__")) %>%
           ungroup() %>%
           select(header, tseq))
  
  write.fasta(sequences = as.list(get(paste0("Final_",path))$tseq),
              names = get(paste0("Final_",path))$header,
              nbchar = 80,
              file.out = paste0(path,".faa"))
}
