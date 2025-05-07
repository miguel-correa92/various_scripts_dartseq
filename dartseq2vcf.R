# Turn dartseq report into vcf with per-allele depth included



setwd('~/bean_collection/DCob20-5617/OrderAppendix_1_250kcounts_DCob20-5617/')

require(data.table)
require(dtplyr)


rename_loci <- function(locus) {
  # locus is a vector of strings. Usually the AlleleID column from DArTseq reports
  a <- str_extract(string = locus, pattern = '^(\\d)+(|)') #extraer el numero del clon
  b <- str_extract(string = locus, pattern = '-[0-9]{0,3}:[A-Z]>[A-Z]') %>% 
    str_replace(':','-') %>% 
    str_replace('>','/')
  
  str_c(a,b)
}

recode.snps <- function(snp.file, AlleleID, regex.samples, chr.col, pos.col) {
  
  # leer sin encabezado
  temp <- fread(file = snp.file, header = FALSE, nrows = 10)
  nremove <- sum(temp[, 1] == "*") 
  # leer completo
  snps <- fread(file = snp.file, header = TRUE, skip = nremove)
  
  # recodificar los datos
  snps2 <- snps %>% 
    lazy_dt() %>% 
    # select(AlleleID, all_of(c(chr.col, pos.col)), matches(regex.samples)) %>%
    select(CHROM = all_of(chr.col), POS = all_of(pos.col), ID = AlleleID, matches(regex.samples)) %>%
    pivot_longer(cols = matches(regex.samples)) %>% 
    mutate(value = fifelse(value == 0, '0/0',
                           fifelse(value == 1, '1/1', 
                                   fifelse(value == 2, '0/1',
                                           './.')))) %>% 
    pivot_wider(names_from = name, values_from = value) %>% 
    mutate(CHROM = str_replace(CHROM, '^$', '.'),
           POS = str_replace(POS, '^0$', '.'),
           ID = rename_loci(ID), 
           REF = str_extract(ID, '[ACGT](?=/)'),
           ALT = str_extract(ID, '(?<=/)[ACGT]'),
           QUAL = '.', 
           FILTER = '.',
           .after = ID) %>%
    as.data.table()
  
  return(snps2)

}

gl.recode.snps <- function(gl, chr.col, pos.col){
  require(data.table)
  require(tidyverse)
  require(dartR)
  
  # Alignment 
  align <- gl$other$loc.metrics %>% 
    select(CHROM = all_of(chr.col), POS = all_of(pos.col), ID = AlleleID) %>% 
    lazy_dt() %>% 
    mutate(ID = rename_loci(ID), CHROM = as.character(CHROM)) %>% 
    arrange(CHROM, POS)
  
  gl$other$loc.metrics$AlleleID %>% length()
  
  align
  
  snps2 <- gl %>% 
    as.matrix() %>% 
    as.data.frame() %>% 
    rownames_to_column('name') %>% 
    lazy_dt() %>% 
    # select(AlleleID, all_of(c(chr.col, pos.col)), matches(regex.samples)) %>%
    # select(CHROM = all_of(chr.col), POS = all_of(pos.col), ID = AlleleID, matches(regex.samples)) %>%
    pivot_longer(cols = -name, names_to = 'ID') %>% 
    mutate(value = fifelse(is.na(value), './.',
                           fifelse(value == 0, '0/0',
                                   fifelse(value == 1, '0/1', 
                                           '1/1')))) %>% 
    pivot_wider(names_from = name, values_from = value) %>% 
    left_join(align, ., by = 'ID') %>% 
    # as.data.table()
    mutate(CHROM = str_replace(CHROM, '^$', '.'),
           POS = str_replace(POS, '^0$', '.'),
           #ID = rename_loci(ID), 
           REF = str_extract(ID, '[ACGT](?=/)'),
           ALT = str_extract(ID, '(?<=/)[ACGT]'),
           QUAL = '.', 
           FILTER = '.',
           .after = ID) %>%
    as.data.table()
  
  return(snps2)
  
   
}

consolidate.SEQ <- function(count.file, names.row = 4) {
  
  require(tidyverse)
  require(data.table)
  require(dtplyr)
  
  # leer sin encabezado
  temp <- fread(file = count.file, header = FALSE, nrows = 10)
  nremove <- sum(temp[, 1] == "*") 
  # old.names <- temp[4, ] %>% unlist(use.names = FALSE) %>% enframe() %>% count(value)
  new.names <- temp[c(names.row, nremove+1), ] %>% 
    data.table::transpose() %>% 
    filter(V1 != '*')
  
  # new.names %>% enframe() %>% view()
  # leer completo
  counts <- fread(file = count.file, header = TRUE, skip = nremove)
  
  # manener unicamente columnas de muestras y de nombre de alelo
  counts2 <- counts %>% 
    lazy_dt() %>% 
    # select(AlleleID, all_of(c(chr.col, pos.col)), matches(regex.samples)) %>%
    select(ID = AlleleID, all_of(new.names$V2)) 
  
  # renombrar columnas
  
  # x.names <- colnames(counts2) %>% unlist(use.names = FALSE) %>% enframe() %>% count(value)
  
  # 1. Consolidar replicas tecnicas ----------------------------------------------
  
  df.final <- counts2 %>% 
    pivot_longer(cols = -ID, names_to = 'V2', values_to = 'reads') %>% 
    left_join(new.names, by = 'V2') %>% 
    select(-V2) %>% 
    group_by(V1, ID) %>% 
    summarise(reads = sum(reads), .groups = 'drop') %>% 
    pivot_wider(names_from = V1, values_from = reads) %>% 
    as.data.table()
    
  
  return(df.final)

}

global.depth <- function(consolidated.seq, sample.list) {
  
  consolidated.seq <- consolidated.seq %>% 
    lazy_dt() %>% 
    select(ID, all_of(sample.list))
  
  DP <- consolidated.seq %>% 
    mutate(ID = rename_loci(ID)) %>% 
    pivot_longer(cols = -ID) %>% 
    group_by(ID) %>% 
    summarise(DP = sum(value), .groups = 'drop') %>% 
    mutate(DP = str_c('DP=', DP)) 
  
  AD <- consolidated.seq %>% 
    pivot_longer(cols = -ID) %>% 
    group_by(ID) %>% 
    summarise(value = sum(value), .groups = 'drop') %>% 
    arrange(ID) %>% 
    mutate(ID = rename_loci(ID)) %>%
    arrange(ID) %>% 
    mutate(allele = rep(c('REF', 'ALT'), nrow(.)/2)) %>% 
    pivot_wider(names_from = allele, values_from = value) %>% 
    unite('AD', c(REF, ALT), sep = ',') 
  
  df.global.depth <- left_join(DP, AD, by = 'ID') %>% 
    unite('INFO', DP, AD, sep = ';') 
  
  return(df.global.depth)

}

sample.depth <- function(consolidated.seq) {

  DP.sample <- consolidated.seq %>% 
    lazy_dt() %>% 
    mutate(ID = rename_loci(ID)) %>% 
    pivot_longer(cols = -ID) %>% 
    group_by(ID, name) %>% 
    summarise(DP = sum(value), .groups = 'drop')
  
  AD.sample <- consolidated.seq %>% 
    lazy_dt() %>% 
    pivot_longer(cols = -ID) %>% 
    mutate(allele = fifelse(str_detect(ID, '--'), 'REF', 'ALT')) %>% 
    mutate(ID = rename_loci(ID)) %>% 
    pivot_wider(names_from = allele, values_from = value) %>%
    unite('AD', c(REF, ALT), sep = ',')   
  
  df.sample.depth <- left_join(DP.sample, AD.sample, by = c('ID', 'name')) %>% 
    unite('INFO', DP, AD, sep = ':') %>% 
    as.data.table()
  
  return(df.sample.depth)
  
}

# join marker and depth data

join.GT.DP <- function(snp.df, depth.df) {
  
  depth.df <- lazy_dt(depth.df)
  #Split marker data
  meta.df <- snp.df %>% 
    select(CHROM:FILTER) %>% 
    lazy_dt()
  
  snp.df <- snp.df %>% 
    lazy_dt() %>% 
    select(-c(CHROM,POS, REF:FILTER)) %>% 
    pivot_longer(cols = -ID)
  
  # join marker and depth data
  
  all.df <- inner_join(snp.df, depth.df, by = c('ID', 'name')) %>% 
    unite('GT', c(value, INFO), sep = ':') %>% 
    pivot_wider(names_from = name, values_from = GT)
    
  # add meta
  
  all.df <- left_join(meta.df, all.df, by = 'ID')
  
  return(all.df)
  
}


DArTSeq2vcf <- function(snp.file, count.file, regex.samples, chr.col, pos.col, names.row = 4, out.name = 'temp.vcf') {
  # 1. recode snps
  # 2. consolidate SEQ files (merge technical replicates)
  # 3. global DP and AD
  # 4. per sample DP and AD
  # 4.1 check that SNP and SEQ files have the same number of markers
  # 5. join marker and depth data
  # 6. add info column with global DP and AD 
  # 6.1 Add FORMAT column
  # 7. prepare header of vcf
  # 8. Save file
  
  step1 <- recode.snps(snp.file = snp.file, regex.samples = regex.samples, 
                      chr.col = chr.col, pos.col = pos.col)
  
  step2 <- consolidate.SEQ(count.file = count.file, names.row = names.row)
  
  step3 <- global.depth(consolidated.seq = step2, sample.list = str_subset(colnames(step1), pattern = regex.samples))
  
  x1 <- length(unique(step1$ID))
  x2 <- length(unique(step2$ID))/2
  cat(x1, 'unique SNP IDs in SNP report\n')
  cat(x2, 'unique SNP IDs in Read Counts report \n')
  if(x1 != x2) cat('Different number of markers in SNP and SEQ files.\n')
  
  step4 <- sample.depth(consolidated.seq = step2)
  
  step5 <- join.GT.DP(snp.df = step1, depth.df = step4)
  
  step6 <- step5 %>% 
    inner_join(step3, by = 'ID') %>% 
    relocate(INFO, .after = FILTER) %>% 
    mutate(FORMAT = 'GT:DP:AD', .after = INFO) %>% 
    rename(`#CHROM` = CHROM) %>% 
    as.data.table()
  
  genome.size <- c("Chr01" = 51433939,"Chr02" = 49670989,"Chr03" = 53438756,"Chr04" = 48048378,"Chr05" = 40923498,"Chr06" = 31236378,
                   "Chr07" = 40041001,"Chr08" = 63048260,"Chr09" = 38250102,"Chr10" = 44302882,"Chr11" = 53580169) %>% 
    enframe('chromosome', 'length.chr') %>% 
    mutate(field = str_c('##contig=<ID=',chromosome, ',length=',length.chr,',URL=https://phytozome-next.jgi.doe.gov/info/Pvulgaris_v2_1'))
  
  '##contig=<ID=ctg1,length=81195210,URL=ftp://somewhere.example/assembly.fa,md5=f126cdf8a6e0c7f379d618ff66beb2da,...>'
  
  header <- c('##fileformat=VCFv4.5', 
              '##assembly=https://phytozome-next.jgi.doe.gov/info/Pvulgaris_v2_1', 
              '##source=DArTseq', 
              '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
              '##INFO=<ID=AD,Number=1,Type=Integer,Description="Allele Total Depth">',
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype"',
              '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
              '##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Allele Read Depth">', 
              "##contig=<ID=Chr01,length=51433939,URL=https://phytozome-next.jgi.doe.gov/info/Pvulgaris_v2_1", 
              "##contig=<ID=Chr02,length=49670989,URL=https://phytozome-next.jgi.doe.gov/info/Pvulgaris_v2_1", 
              "##contig=<ID=Chr03,length=53438756,URL=https://phytozome-next.jgi.doe.gov/info/Pvulgaris_v2_1", 
              "##contig=<ID=Chr04,length=48048378,URL=https://phytozome-next.jgi.doe.gov/info/Pvulgaris_v2_1", 
              "##contig=<ID=Chr05,length=40923498,URL=https://phytozome-next.jgi.doe.gov/info/Pvulgaris_v2_1", 
              "##contig=<ID=Chr06,length=31236378,URL=https://phytozome-next.jgi.doe.gov/info/Pvulgaris_v2_1", 
              "##contig=<ID=Chr07,length=40041001,URL=https://phytozome-next.jgi.doe.gov/info/Pvulgaris_v2_1", 
              "##contig=<ID=Chr08,length=63048260,URL=https://phytozome-next.jgi.doe.gov/info/Pvulgaris_v2_1", 
              "##contig=<ID=Chr09,length=38250102,URL=https://phytozome-next.jgi.doe.gov/info/Pvulgaris_v2_1", 
              "##contig=<ID=Chr10,length=44302882,URL=https://phytozome-next.jgi.doe.gov/info/Pvulgaris_v2_1", 
              "##contig=<ID=Chr11,length=53580169,URL=https://phytozome-next.jgi.doe.gov/info/Pvulgaris_v2_1"
              ) %>% 
    cbind()
    
  
  write(header, out.name, sep = '\t')
  fwrite(step6, file = out.name, append = TRUE, quote = FALSE, sep = '\t', col.names = TRUE)
  cat('VCF saved to ', out.name, '\n')

}

gl2vcf <- function(gl, count.file, chr.col, pos.col, names.row = 4, out.name = 'gl.temp.vcf') {
  # 1. recode snps
  # 2. consolidate SEQ files (merge technical replicates)
  # 3. global DP and AD
  # 4. per sample DP and AD
  # 4.1 check that SNP and SEQ files have the same number of markers
  # 5. join marker and depth data
  # 6. add info column with global DP and AD 
  # 6.1 Add FORMAT column
  # 7. prepare header of vcf
  # 8. Save file
  
  step1 <- gl.recode.snps(gl, chr.col = chr.col, pos.col = pos.col)
  
  step2 <- consolidate.SEQ(count.file = count.file, names.row = names.row)
  
  step3 <- global.depth(consolidated.seq = step2, sample.list = indNames(gl))
  
  x1 <- length(unique(step1$ID))
  x2 <- length(unique(step2$ID))/2
  cat(x1, 'unique SNP IDs in SNP genlight\n')
  cat(x2, 'unique SNP IDs in Read Counts report \n')
  if(x1 != x2) cat('Different number of markers in SNP and SEQ files.\n')
  
  step4 <- sample.depth(consolidated.seq = step2)
  
  step5 <- join.GT.DP(snp.df = step1, depth.df = step4)
  
  step6 <- step5 %>% 
    inner_join(step3, by = 'ID') %>% 
    relocate(INFO, .after = FILTER) %>% 
    mutate(FORMAT = 'GT:DP:AD', .after = INFO) %>% 
    rename(`#CHROM` = CHROM) %>% 
    as.data.table()
  
  genome.size <- c("Chr01" = 51433939,"Chr02" = 49670989,"Chr03" = 53438756,"Chr04" = 48048378,"Chr05" = 40923498,"Chr06" = 31236378,
                   "Chr07" = 40041001,"Chr08" = 63048260,"Chr09" = 38250102,"Chr10" = 44302882,"Chr11" = 53580169) %>% 
    enframe('chromosome', 'length.chr') %>% 
    mutate(field = str_c('##contig=<ID=',chromosome, ',length=',length.chr,',URL=https://phytozome-next.jgi.doe.gov/info/Pvulgaris_v2_1'))
  
  '##contig=<ID=ctg1,length=81195210,URL=ftp://somewhere.example/assembly.fa,md5=f126cdf8a6e0c7f379d618ff66beb2da,...>'
  
  header <- c('##fileformat=VCFv4.5', 
              '##assembly=https://phytozome-next.jgi.doe.gov/info/Pvulgaris_v2_1', 
              '##source=DArTseq', 
              '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
              '##INFO=<ID=AD,Number=1,Type=Integer,Description="Allele Total Depth">',
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype"',
              '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
              '##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Allele Read Depth">', 
              "##contig=<ID=Chr01,length=51433939,URL=https://phytozome-next.jgi.doe.gov/info/Pvulgaris_v2_1", 
              "##contig=<ID=Chr02,length=49670989,URL=https://phytozome-next.jgi.doe.gov/info/Pvulgaris_v2_1", 
              "##contig=<ID=Chr03,length=53438756,URL=https://phytozome-next.jgi.doe.gov/info/Pvulgaris_v2_1", 
              "##contig=<ID=Chr04,length=48048378,URL=https://phytozome-next.jgi.doe.gov/info/Pvulgaris_v2_1", 
              "##contig=<ID=Chr05,length=40923498,URL=https://phytozome-next.jgi.doe.gov/info/Pvulgaris_v2_1", 
              "##contig=<ID=Chr06,length=31236378,URL=https://phytozome-next.jgi.doe.gov/info/Pvulgaris_v2_1", 
              "##contig=<ID=Chr07,length=40041001,URL=https://phytozome-next.jgi.doe.gov/info/Pvulgaris_v2_1", 
              "##contig=<ID=Chr08,length=63048260,URL=https://phytozome-next.jgi.doe.gov/info/Pvulgaris_v2_1", 
              "##contig=<ID=Chr09,length=38250102,URL=https://phytozome-next.jgi.doe.gov/info/Pvulgaris_v2_1", 
              "##contig=<ID=Chr10,length=44302882,URL=https://phytozome-next.jgi.doe.gov/info/Pvulgaris_v2_1", 
              "##contig=<ID=Chr11,length=53580169,URL=https://phytozome-next.jgi.doe.gov/info/Pvulgaris_v2_1"
  ) %>% 
    cbind()
  
  
  write(header, out.name, sep = '\t')
  fwrite(step6, file = out.name, append = TRUE, quote = FALSE, sep = '\t', col.names = TRUE)
  cat('VCF saved to ', out.name, '\n')
  
}

























