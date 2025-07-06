annotate_exons_introns <- function(txt_file) {
  # Read the sequence
  seq <- readLines(txt_file, warn = FALSE)
  seq <- paste(seq, collapse = "")  # In case the sequence spans multiple lines
  
  # Use regex to identify regions of exons (uppercase) and introns (lowercase)
  matches <- stringr::str_match_all(seq, "([A-Z]+|[a-z]+)")[[1]][,1]
  
  pos <- 1
  exon_count <- 0
  intron_count <- 0
  codon_counter <- 1
  res <- data.frame()
  
  for (frag in matches) {
    len <- nchar(frag)
    start <- pos
    end <- pos + len - 1
    
    if (grepl("^[A-Z]+$", frag)) {
      exon_count <- exon_count + 1
      codons <- ceiling(len / 3)
      first_codon <- codon_counter
      last_codon <- codon_counter + codons - 1
      codon_counter <- last_codon + 1
      type <- paste0("exon", exon_count)
      res <- rbind(res, data.frame(
        name = type,
        start = start,
        end = end,
        first_codon = first_codon,
        last_codon = last_codon
      ))
    } else {
      intron_count <- intron_count + 1
      type <- paste0("intron", intron_count)
      res <- rbind(res, data.frame(
        name = type,
        start = start,
        end = end,
        first_codon = NA,
        last_codon = NA
      ))
    }
    
    pos <- end + 1
  }
  
  return(res)
}

# Human Gene APC (ENST00000257430.9) from GENCODE V48
# Description: APC regulator of WNT signaling pathway, transcript variant 31 (from RefSeq NM_001407467.1)
# Gencode Transcript: ENST00000257430.9
# Gencode Gene: ENSG00000134982.19
# Transcript (Including UTRs)
# Position: hg38 chr5:112,737,885-112,846,239 Size: 108,355 Total Exon Count: 16 Strand: +
#   Coding Region
# Position: hg38 chr5:112,754,891-112,844,126 Size: 89,236 Coding Exon Count: 15 
# Downloaded sequences for introns and CDR exons 
# link: https://genome.ucsc.edu/cgi-bin/hgc?hgsid=2706277374_d7hIO0Lu5U5hakfta2i6iCRVXPym&g=htcDnaNearGene&i=ENST00000257430.9&c=chr5&l=112737884&r=112846239&o=knownGene&boolshad.hgSeq.promoter=0&hgSeq.promoterSize=1000&boolshad.hgSeq.utrExon5=0&hgSeq.cdsExon=on&boolshad.hgSeq.cdsExon=0&boolshad.hgSeq.utrExon3=0&hgSeq.intron=on&boolshad.hgSeq.intron=0&boolshad.hgSeq.downstream=0&hgSeq.downstreamSize=1000&hgSeq.granularity=gene&hgSeq.padding5=0&hgSeq.padding3=0&boolshad.hgSeq.splitCDSUTR=0&hgSeq.casing=exon&boolshad.hgSeq.maskRepeats=0&hgSeq.repMasking=lower&submit=submit

annotated_regions <- annotate_exons_introns("APC_gene.txt")
write.csv(annotated_regions, "APC_gene_annotated.csv", row.names = FALSE)