# Read a list of dart reports and join them into a single data.table
## path_reports (STRING): directory where all files are stored
## marker (STRING): type of report being read. The options are
##   'snp_onerow' searches for the following pattern in the filename: "[Mm]apping"
##   'snp_tworow' searches for the following pattern in the filename: "SNP_[0-9]"
##   'silico'     searches for the following pattern in the filename: "[Ss]ilicoDArT"
## recursive: should directories within path_reports be searched as well? 
## pattern (STRING): use if the files have a pattern besides the ones metiones above
join.dart.raw <- function(path_reports, marker = NULL, pattern = NULL, recursive = FALSE){
  require(purrr)
  require(data.table)
  
  if (is.null(pattern)) {
    if (!marker %in% c('snp_onerow', 'snp_tworow', 'silico')) stop("Marker can only be one of 'snp_onerow', 'snp_tworow', 'silico'")
    
    if (marker == 'snp_onerow'){
      pattern <- '[Mm]apping'
    } else if(marker == 'snp_tworow'){
      pattern <- 'SNP_[0-9]'  
    } else {
      pattern <- '[Ss]ilicoDArT'
    }
    
  }
  
  list_files <- list.files(pattern = pattern, path = path_reports, full.names = TRUE, recursive = recursive)
  
  df1 <- fread(file = list_files[[1]], header = FALSE)
  nremove <- sum(df1[, 1] == "*") + 1
  list_files2 <- list_files[-1]
  dfs <- map(list_files2, ~fread(file = .x, header = FALSE, skip = nremove))
  
  list_all <- vector(mode = 'list', length = length(list_files))
  list_all[[1]] <- df1
  list_all[2:length(list_all)] <- dfs
  
  df_all <- rbindlist(list_all, use.names = FALSE, idcol = NULL)
  
  return(df_all)
  
} #end joind.dart.raw



# calculate the call rate of a genlight/DartR object and return it as a table with alleleID and call rate

callrate.loci.gl <- function(gl){
  require(dplyr)
  require(magrittr)
  require(tibble)
  cr <- gl %>%
    as.matrix %>% 
    is.na(.) %>% 
    not() %>% 
    colMeans() %>%  
    as_tibble(rownames = 'LocusID', .name_repair = 'unique') %>% 
    rename(callrate.loci = value)
  
  return(cr)
}

# calculate the call rate of a genlight/DartR object and return it as a table with:
## sample id
## call rate
## number of markers called 
## total number of markers across the dataset

callrate.sample.gl <- function(gl){
  require(dplyr)
  require(magrittr)
  require(tibble)
  cr <- gl %>%
    as.matrix %>% 
    is.na(.) %>% 
    not() %>% 
    rowMeans() %>%  
    as_tibble(rownames = 'sample') %>% 
    repair_names() %>% 
    rename(callrate.sample = value)
  
  n.markers <- gl %>%
    as.matrix %>% 
    is.na(.) %>% 
    not() %>% 
    rowSums() %>%  
    as_tibble(rownames = 'sample', .name_repair = 'unique') %>% 
    rename(n.called.markers = value)
  
  df <- left_join(cr, n.markers, by = 'sample') %>% 
    mutate(n.total.markers = nLoc(gl))
  
  return(df)
}

#Clean AlleleIDs from DArT to make make them equal to LocNames in dartR
rename_loci <- function(locus) {
  # locus is a vector of strings. Usually the AlleleID column from DArTseq reports
  a <- str_extract(string = locus, pattern = '^(\\d)+(|)') #extraer el numero del clon
  b <- str_extract(string = locus, pattern = '-[0-9]{0,3}:[A-Z]>[A-Z]') %>% 
    str_replace(':','-') %>% 
    str_replace('>','/')
  
  str_c(a,b)
}

# Combinar las replicas tecnicas de los reportes SEQ de DArTseq en una sola columna
SEQ.consolidate.tech.replicates <- function(x, 
                                            first.meta.col, 
                                            last.meta.col, 
                                            suffix, 
                                            col.AlleleID, regex = FALSE) {
  
  # Uso

  # x =  Reporte SEQ, sin encabezado
  # first.meta.col = Nombre de la primera columna del metadata 
  # last.meta.col  = Nombre de la ultima columna del metadata 
  # suffix  = Sufijo usado para diferenciar las muestras replicadas de igual nombre
  # col.AlleleID = Nombre de la columna con los AlleleIDs
  # regex = el patron (suffix) es una expresion regular o debe ser tomada literalmente

  
  # require(data.table)
  # require(dtplyr)
  require(stringr)
  require(tidyselect)
  require(dplyr)
  
  # 1. separar los metadatos de los marcadores de los conteos de lecturas, dejando los AlleleIDs
  x <- x #%>% lazy_dt()
  
  x.meta <- x %>% select(all_of(first.meta.col):all_of(last.meta.col)) %>% as.data.table()
  
  meta.colnames <- colnames(x.meta)
  meta.colnames <- meta.colnames[meta.colnames != col.AlleleID]
  
  x.mks <- x %>% select(!all_of(meta.colnames)) %>% as.data.table()
  
  if (regex) {
    # 2. Extraer la lista de muestras repetidas usando el sufijo
    rep.samples <- colnames(x.mks) %>% str_subset(pattern = regex(suffix))  
    # 3. limpiar nombres borrando el sufijo
    samples <- str_remove(rep.samples, pattern = regex(suffix))
  
  } else{
    rep.samples <- colnames(x.mks) %>% str_subset(pattern = coll(suffix))  
    
    samples <- str_remove(rep.samples, pattern = coll(suffix))
  }
  
  # 4.  Extraer los pares de muestras
  # 4.1 Hacer la suma de lecuras
  # 4.2 Renombrar la columna final con el nombre de la muestra
  # 4.3 guardar en un df
  
  # s <- 1
  
  list.pairs <- vector(mode = 'list', length = length(samples))
  
  x.new.cols <- x.mks %>% select(all_of(col.AlleleID))
  
  
  for (s in seq_along(list.pairs)) {
    cat(samples[s], '\n')
    
    sample.name <- samples[s]
    
    
    sample.name.replace <- ('sum.reads')
    names(sample.name.replace) <- sample.name
    
    temp <- x.mks %>% 
      select(AlleleID, matches(samples[s])) %>% 
      mutate(sum.reads = rowSums(across(-all_of(col.AlleleID))), .keep = 'unused') %>% 
      rename(all_of(sample.name.replace))
    
    x.new.cols <- left_join(x.new.cols, temp, by = col.AlleleID)
    
  }
  
  # 5. Borrar las columnas viejas
  
  x.mks.no.repeats <- x.mks %>% select(-any_of(c(rep.samples, samples)))
  
  # 6. agregar las columnas nuevas
  
  x.mks.new <- left_join(x.mks.no.repeats, x.new.cols, by = col.AlleleID)
  
  # 7. Agregar las columnas del metadata
  
  x.final <- left_join(x.meta, x.mks.new, by = col.AlleleID)
}# end SEQ.consolidate.tech.replicates

# Recalcular la profundidad promedio por lecturas de un subset de datos
recalc.rdepth <- function(x, AlleleID.col ='AlleleID') {
  # Uso
  ## x: df con los AlleleIDs originales y las muestras en formato de conteos de lecturas por alelo
  ## AlleleID.col: Nombre de la columna con los AlleleIDs originales
  ## Nota importante: Puede que haya problema con rename_loci. Lo mejor es cargarla por fuera de la funcion tambien
  require(dplyr)
  require(data.table)
  require(dtplyr)
  
  
  rename_loci <- function(locus) {
    a <- str_extract(string = locus, pattern = '^(\\d)+(|)') #extraer el numero del clon
    b <- str_extract(string = locus, pattern = '-[0-9]{0,3}:[A-Z]>[A-Z]') %>% 
      str_replace(':','-') %>% 
      str_replace('>','/')
    
    str_c(a,b)
  }
    
  x <- as.data.table(x)
  
  # 1. reemplazar 0 por NAs
  # https://stackoverflow.com/questions/11036989/replace-all-0-values-to-na
  for (j in names(x)) set(x,which(x[[j]] == 0),j,NA)
  
  
  # 2. calcular profundidad de secuenciacion por alelo
  
  rdepth.allele <- x %>% lazy_dt() %>% 
    mutate(r.depth.allele = rowMeans(across(-all_of(AlleleID.col)), na.rm = TRUE), .keep = 'unused') %>% 
    as_tibble()
  
  rdepth.allele <- replace_na(rdepth.allele, replace = list( r.depth.allele = 0))
  
  # 3. calcular profundidad de secuenciacion por marcador
  
  rdepth.mk <- x %>% 
    lazy_dt() %>% 
    mutate(AlleleID2 = rename_loci(.data[[AlleleID.col]]), .after = 1) %>% 
    pivot_longer(cols = -c(all_of(AlleleID.col), AlleleID2), names_to = 'sample', values_to = 'value') %>% 
    group_by(AlleleID2) %>% 
    summarise(rdepth.mk = mean(value, na.rm = TRUE)) %>% 
    as.data.table()
  
  # 4. poner todo junto
  
  rdepth.all <- rdepth.allele %>% lazy_dt() %>% 
    mutate(AlleleID2 = rename_loci(.data[[AlleleID.col]]), .after = 1) %>% 
    left_join(rdepth.mk, by = 'AlleleID2') %>% 
    as.data.table() %>% 
    select(-AlleleID2)
  
   return(rdepth.all)

}



# Calcular el minor allele count y agegarlo a un objeto genlight
## x es un objeto genlight/DartR
calc.mac.gl <- function(x) {
  
  if (!class(x) %in% c('genlight', 'dartR')) {stop('only objects of class genlight or dartR valid as input')}
  
  x <- gl.recalc.metrics(x, verbose = 1)
  
  # A = CallRate x N.individuos
  # B = FreqHomRef * A * 2 ## 2 es porque los homocigotos tienen dos copias de cada alelo
  # C = FreqHomSnp * A * 2 ## 2 es porque los homocigotos tienen dos copias de cada alelo
  # D = FreqHets   * A * 1 ## 1 es porque los heterocigotos tienen 1 copia de cada alelo
  # AlleleCountRef = B + D
  # AlleleCOuntSnp = C + D
  # AlleleCountRef + AlleleCountSnp = N.individuos * 2
  # MAC = min(AlleleCountRef, AlleleCOuntSnp)
  
  A <- round(x$other$loc.metrics$CallRate * nInd(x), digits = 0)
  B <- round(A * x$other$loc.metrics$FreqHomRef * 2, digits = 0)
  C <- round(A * x$other$loc.metrics$FreqHomSnp * 2, digits = 0)
  D <- round(A * x$other$loc.metrics$FreqHets   * 1, digits = 0)
  
  AlleleCountRef <- B + D
  AlleleCountSnp <- C + D
  
  # Check
  
  check.sum <- all(AlleleCountRef + AlleleCountSnp == A * 2)
  
  if (check.sum) {
    cat('Sums of allele counts match the expected value of (Number of individuals * 2 * Callrate)\n')
  } else {
    
    cat('The sums of allele counts of at least one marker does not meet the expected value of (Number of individuals * 2 * Callrate)')
  }
  
  # Hallar el conteo alelico menor por marcador
  MinorAlleleCount <- pmin(AlleleCountRef, AlleleCountSnp)
  
  x$other$loc.metrics$MinorAlleleCount <- MinorAlleleCount
  
  return(x)

} # fin de calc.mac.gl

# convertir un objeto genlight al formato geno que usa LEA (dartR tambien tiene esta funcion)
genlight2geno <- function(x, out.name){
  
  genotype <- as.matrix(x)
  genotype[is.na(genotype)] <- 9   
  
  # write geno
  
  require(data.table)
  require(readr)
  require(tibble)
  # Transponer la matriz de datos
  
  t.genotype <- transpose(as.data.table(genotype))
  # t.genotype <- transpose(as.data.table(genotype, keep.rownames = 'sample'), make.names = 'sample' )
  
  # fwrite(t.genotype, file = paste0(out.name, '.geno'), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '')
  write_delim(t.genotype, file = paste0(out.name, '.geno'), quote = 'none', col_names = FALSE, delim = '', eol = '\n')
  
  # all(indNames(x) == colnames(t.genotype))
  
  return(NULL)
  
}

# Extraer los valores de cross-entropy de un proyecto de snmf del paquete LEA en forma de data.frame

extract.cross.entropy <- function(snmf.project, K) {
  # snmf.project = proyect created by snmf (already read into R)
  # K =  numeric vector detailing range of K to extract the cross-entropy values
  
  
  all.ce <- map(K, function(x) cross.entropy(snmf.project, x))
  runs <- map_chr(all.ce, ~colnames(.x))
  names(all.ce) <- runs
  all.ce <- map(all.ce, function(x){
    colnames(x) <- 'ce'
    x <- as_tibble(x, rownames = 'rep.num')
    return(x)
  })
  
  all.ce <- bind_rows(all.ce, .id = 'run') %>% mutate(run = factor(run, levels = str_c("K", K, sep = ' = ')))
  
  return(all.ce)
}

# extraer las matrices q y los valores de cross-entropy de un proyecto de snmf
process.snmf <- function(x, y, k) {
  # x es el proyecto de snmf
  # y es el objeto genlight del que se obtuvo el archivo en formato geno 
  ##  (basicamente para extraer el nombre de las muestras)
  # k es un vector numerico indicando los valores de K a extraer
  sample.names <- indNames(y) 
  
  
  # Cross-entropy
  
  all.ce <- extract.cross.entropy(x, k)
  
  
  # mejor valor de cross entropy
  best.ce <- all.ce %>% 
    group_by(run) %>% 
    slice(which.min(ce)) %>% 
    mutate(best.run = as.numeric(str_remove(rep.num, 'run ')),
           numeric.k = as.numeric(str_remove(run, 'K = ')))
  
  
  # Extraer matrices Q
  q.matrices <- map2(best.ce$numeric.k, best.ce$best.run, function(.x, .y){
    Q(object = x, K = .x, run = .y)
  })
  
  names(q.matrices) <- best.ce$run
  
  # Reorganizar las matrices en una sola tabla
  q.matrices.df <- q.matrices %>% 
    map(function(x){
      x %>% 
        as_tibble() %>% 
        mutate(accession = sample.names, .before = 1) %>% 
        pivot_longer(cols = -accession, values_to = 'adm.coef', names_to = 'ancestral.pop')
      
    }) %>% 
    bind_rows(.id = 'run')
  
  
  final.result <- list(cross.entropy.df = all.ce, q.matrices = q.matrices.df)
  
  return(final.result)
} # final process.snmf

# convertir de genlight a hapmap
gl2hapmap <- function(gl, reference.name) {
    
  # gl <- test.gl.2
  # gl <- snp.data
  alleleID <- locNames(gl)
  sample.names <- indNames(gl)
  SNPs <- gl$other$loc.metrics$SNP %>% as.character() 
  
    # Detectar el formato del nombre del genoma de referencia.
  ## Algunos reportes dart incluyen la posicion del tag y del SNP. Anotados como ChromPosTag y ChromPosSnp, respectivamente
  
  check.reference.format <- str_detect(colnames(gl$other$loc.metrics), 'ChromPosSnp') %>% any()
  
  if (check.reference.format) {
    
    chrom.col <- str_c('Chrom_', reference.name)
    pos.col   <- str_c('ChromPosSnp_', reference.name)
    
  } else {
    
    chrom.col <- str_c('Chrom_', reference.name)
    pos.col   <- str_c('ChromPos_', reference.name)
    
  }
  
  ## Extraer informacion de cromosoma y position del SNP
  chromosome <- gl$other$loc.metrics[,chrom.col] %>% as.character()
  snp.position <- gl$other$loc.metrics[,pos.col]
  
  
  # Detectar la presencia de informacion de Strand
  ## Archivos recientes (2022- adelante) incluyen informacion de la hebra ('strand') a la que se alinea la lectura
  
  check.strand <- str_detect(colnames(gl$other$loc.metrics), '[Ss]trand') %>% any()
  
  ## extraer la columna de strand o dejarla como NA
  if (check.strand) {
    strand <- str_c('Strand_', reference.name)
    strand.col <- gl$other$loc.metrics[,strand] %>% as.character()
  } else {
    strand.col <- NA
  }
   
  # Girar la matriz y agregar nombres de muestras antes de cualquier cosa 
  # (la version 1.2.2 de dtplyr con la version 1.14.4 de data.table reorganizan las columnas en pivot_wide, 
  # por eso se usa 'transpose' de data.table)
  
  mk.matrix <- as.matrix(gl) %>% 
    as.data.table(keep.rownames = FALSE) %>% 
    data.table::transpose()
    
  
  colnames(mk.matrix) <- sample.names
  
  # El siguiente bloque:
  ## - Agrega ID del Alelo, 
  ## - Los SNP
  ## - Separa la columna SNP en Alelo ref y alt, en ese orden
  ## - Transforma la tabla al formato largo
  ## - Transforma los genotypos de conteos a letras: ej. 0 -> AA; 1 -> AT; 2 -> TT
  ## - Transforma la tabla al formato largo
  ## - Agrega las columnas adicionales del formato HapMap
  ## - Reemplaza los NA en los marcadores por NN
  
  hapmap.table <- mk.matrix %>% 
    # lazy_dt() %>% 
    mutate(alleleID = alleleID, 
           SNP.alleles = SNPs %>% str_extract('[ACGT]>[ACGT]'), 
           ref = str_sub(SNP.alleles, start = 1, end = 1), 
           alt = str_sub(SNP.alleles, start = 3, end = 3), 
           chrom = chromosome, 
           pos = snp.position,
           strand = strand.col, 
           assembly = reference.name, 
           .before = 1) %>% 
    # as_tibble()
    pivot_longer(cols = -c(alleleID:assembly), names_to = 'sample', values_to = 'genotype', ) %>% 
    mutate(genotype = if_else(genotype == 0, 
                              str_c(ref, ref, sep = ''),
                              if_else(genotype == 1, 
                                      str_c(ref, alt, sep = ''),
                                      str_c(alt, alt, sep = '')))) %>% 
    
    ## pivot wider tambien reorganiza las filas
    pivot_wider(names_from = 'sample', values_from = 'genotype', names_sort = FALSE, ) %>%
    # as_tibble()
    mutate(
      'rs#' =  alleleID, 
      "alleles" = str_replace(SNP.alleles, '>', '/'), 
      chrom = chrom, 
      pos = pos,
      strand = strand, 
      assembly = assembly, 
      center = NA, 
      protLSID = NA, 
      assayLSID = NA, 
      panelLSID = NA, 
      QCcode = NA, 
      .before = 1
    ) %>% 
    # as_tibble()
    dplyr::select(-alleleID, -SNP.alleles, -ref, -alt) %>% 
    relocate(`rs#`, 
             alleles, 
             chrom,
             pos,
             strand,
             assembly,
             center, 
             protLSID, 
             assayLSID,
             panelLSID,
             QCcode, .before = 1
    ) %>% 
    arrange(chrom, pos) %>% 
    mutate(across(.cols = all_of(sample.names), .fns = ~replace_na(.x, 'NN'))) %>% 
    as_tibble()
  
  return(hapmap.table)
  
}

#calcular la distancia modificada de rogers
mrd <- function(mat) {
  require(magrittr)
  if(!is.matrix(mat)){stop('input not a matrix!')}
  mrdmat <- matrix(NA, nrow = ncol(mat), ncol = ncol(mat))
  for (samp in 1:(ncol(mat)-1)) {
    #cat(samp, '\n')
    # mrdmat[,samp] <-  colMeans(    (mat[,samp]  -  mat[,1:ncol(mat)])^2,  na.rm=T) + 
    #                   colMeans(((1- mat[,samp]) - (1- mat[,1:ncol(mat)]))^2, na.rm = TRUE) %>% sqrt()
    
    if(samp != ncol(mat)){
      temp <- colMeans((mat[,samp]  -  mat[,samp:ncol(mat)])^2 + ((1 - mat[,samp]) - (1- mat[,samp:ncol(mat)]))^2   ,  na.rm=T) %>%
        sqrt * 1/sqrt(2) 
      temp <- c(rep(NA, samp -1), temp)
      mrdmat[samp,] <- temp
    }
    
  }
  colnames(mrdmat) <- colnames(mat)
  rownames(mrdmat) <- colnames(mat)
  diag(mrdmat) <- 0
  return(mrdmat)
  #gc()
}


# Convertir media matriz en una matriz completa (i.e. simetrica) 
makeSymm <- function(m) {
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    return(m)
}

# Convertir media matriz a NA (Excluir uno de los triangulos de la matrix (las distancias son, en teoria, simetricas, i.e., AB es igual a BA))
mat.half.na <- function(mat){
  mat[lower.tri(mat)] <- NA
  return(mat)
}


# Pasar una matriz al formato largo
mat_long <- function(x){
  require(dplyr)
  as_tibble(x, rownames = 'row.samples') %>% 
    pivot_longer(-c(row.samples), names_to = "col.samples", values_to = "dist")
}

# Transformar una matriz al formato largo y remover uno de los "triangulos"
prepare.mat <- function(x){
  x <-  x %>% 
    mat.half.na %>% 
    mat_long %>% 
    na.omit
  return(x)
}

# crea un resumen de un objeto genlight y lo retorna como un data frame de una fila
summary.stats.gl.snp <- function(x) {
  if(!class(x) %in% c("genlight", 'dartR') ) stop('Object not genlight or dartR class')
    
  n.samples <- nInd(x)
  n.populations <- nPop(x)
  n.mks <- nLoc(x)
  n.monomorph <- count.monomorphs.snp(x)
  n.all_na <- count.all_na(x)
  
  sample.cr     <- callrate.sample.gl(x) %>% pull(callrate.sample)
  sample.cr.avg <- sample.cr %>% mean(na.rm = TRUE)
  sample.cr.sd  <- sample.cr %>% sd(na.rm = TRUE)
  
  locus.cr      <- callrate.loci.gl(x) %>% pull(callrate.loci)
  locus.cr.avg  <- locus.cr %>% mean(na.rm = TRUE)
  locus.cr.sd   <- locus.cr %>% sd(na.rm = TRUE)
  
  
  df <- tibble(n.samples = n.samples, 
               n.populations = n.populations, 
               n.mks = n.mks, 
               n.monomorphic.loci = n.monomorph, 
               n.polymorphic = n.mks - n.monomorph, 
               n.loci.all.NA = n.all_na, 
               sample.callrate.mean = sample.cr.avg,
               sample.callrate.sd   = sample.cr.sd,
               locus.callrate.mean = locus.cr.avg,
               locus.callrate.sd   = locus.cr.sd)
  return(df)
}

# Contar sitios monomorficos
count.monomorphs.snp <- function(x){
  require(magrittr)
  if(!class(x) %in% c("genlight", 'dartR') ) stop('Object not genlight or dartR class')
  
  #Check columns (i.e markers) if all values are equal.
  x %>% as.matrix() %>% 
    apply(., 2, function(x) all(x == 2, na.rm = TRUE) | all(x == 0, na.rm = TRUE)) %>% 
    sum()
  
}  

# Contar sitios con callrate == 0
count.all_na <- function(x){
  require(magrittr)
  if(!class(x) %in% c("genlight", 'dartR') ) stop('Object not genlight or dartR class')
  
  #Check columns (i.e markers) if all values are missing (NA)
  x %>% as.matrix() %>% 
    apply(., 2, function(x) all(is.na(x))) %>% 
    sum()
  
}

# calcular maf/frecuencia de genotipo menor en SilicoDArTs
maf.silicos <- function(x){
  if(class(x) != "genlight" ) stop('Object not genlight class')
  
  x %>% as.matrix() %>% 
    colMeans(na.rm = TRUE) %>% 
    enframe('Loc.ID', 'maf') %>% 
    mutate(maf = if_else(maf <= 0.5, maf, 1 - maf))
  
}

count.monomorphs.silico <- function(x){
  require(magrittr)
  if(!class(x) %in% c("genlight", 'dartR') ) stop('Object not genlight or dartR class')
  
  #Check columns (i.e markers) if all values are equal.
  x %>% as.matrix() %>% 
    apply(., 2, function(x) all(x == 1, na.rm = TRUE) | all(x == 0, na.rm = TRUE)) %>% 
    sum()
  
}


summary.stats.gl.silico <- function(x) {
  if(class(x) != "genlight" ) stop('Object not genlight class')
  
  
  n.samples <- nInd(x)
  n.populations <- nPop(x)
  n.mks <- nLoc(x)
  n.monomorph <- count.monomorphs.silico(x)
  n.all_na <- count.all_na(x)
  
  sample.cr     <- callrate.sample.gl(x) %>% pull(callrate.sample)
  sample.cr.avg <- sample.cr %>% mean(na.rm = TRUE)
  sample.cr.sd  <- sample.cr %>% sd(na.rm = TRUE)
  
  locus.cr      <- callrate.loci.gl(x) %>% pull(callrate.loci)
  locus.cr.avg  <- locus.cr %>% mean(na.rm = TRUE)
  locus.cr.sd   <- locus.cr %>% sd(na.rm = TRUE)
  
  
  df <- tibble(n.samples = n.samples, 
               n.populations = n.populations, 
               n.mks = n.mks, 
               n.monomorphic.loci = n.monomorph, 
               n.polymorphic = n.mks - n.monomorph, 
               n.loci.all.NA = n.all_na, 
               sample.callrate.mean = sample.cr.avg,
               sample.callrate.sd   = sample.cr.sd,
               locus.callrate.mean = locus.cr.avg,
               locus.callrate.sd   = locus.cr.sd)
  return(df)
}

          
# calcular mac/conteo de genotipo menor en SilicoDArTs
calc.mac.gl.silico <- function(x) {
  
  if (!class(x) %in% c('genlight', 'dartR')) {stop('only objects of class genlight or dartR valid as input')}
  
  x <- gl.recalc.metrics(x, verbose = 1)
  
  # A = CallRate x N.individuos
  # B = A * OneRatio ## 
  # C = A * (OneRatio - 1) 
  
  A <- round(x$other$loc.metrics$CallRate * nInd(x), digits = 0)
  B <- round(A * x$other$loc.metrics$OneRatio, digits = 0)
  C <- round(A * (1 - x$other$loc.metrics$OneRatio), digits = 0)
  D <- C + B

  # all(A == D)
    
    # Hallar el conteo alelico menor por marcador
  MinorAlleleCount <- pmin(B, C)
  
  x$other$loc.metrics$MinorAlleleCount <- MinorAlleleCount
  
  return(x)

}



# convert genlight data to the popclust format
genlight2popclust <- function(x, out.name){
  # x <- snp.data.100
  
  
  genotype <- as.matrix(x)
  genotype[is.na(genotype)] <- 3   
  
  # write popClust
  
  require(data.table)
  require(readr)
  require(tibble)
  # Transponer la matriz de datos
  
  t.genotype <- transpose(as.data.table(genotype))
  # t.genotype <- transpose(as.data.table(genotype, keep.rownames = 'sample'), make.names = 'sample' )
  
  sample.names <- colnames(t.genotype) %>% t() %>% as.data.frame()
  
  out.name <- paste0(out.name, '.dat')
  write_delim(sample.names, file = out.name, col_names =  FALSE, delim = ' ', eol = '\n')
   
  # fwrite(t.genotype, file = paste0(out.name, '.geno'), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '')
  write_delim(t.genotype, file = out.name, quote = 'none', col_names = FALSE, delim = '', eol = '\n', append = TRUE)
  
  # all(indNames(x) == colnames(t.genotype))
  
  return(NULL)
  
}


read.popclust.admixture <- function(x) {
  x <- readLines(x)

  colnames.df <- which(str_detect(x, "^Admixture analysis$")) + 3
  start <-which(str_detect(x, "^Admixture analysis$")) + 4
  end   <- which(str_detect(x, "^Pairwise co-assignment probability from admixture analysis")) -4
  
  
  colnames1 <- x[colnames.df] %>% str_remove(' : Inferred ancestry in clusters') %>% str_split(' ') %>% unlist()
  
  colnames2 <- colnames1[-length(colnames1)]
  
  admix.colnames  <- colnames1[length(colnames1)] %>% str_split('~') %>% unlist() %>% as.numeric()
  admix.colnames2 <- admix.colnames[1]:admix.colnames[2] %>% str_c('X', .)
  
  colnames.all <- c(colnames2, admix.colnames2)
  
  admmix.table <- x[start:end]
  
  final.df <- admmix.table %>% 
    str_remove(':') %>% 
    stringr::str_squish() %>% 
    tibble(all =.) %>% 
    separate_wider_delim(cols = all, delim = ' ', names = colnames.all)
    
  final.df <- suppressMessages(type_convert(final.df))

  return(final.df)
}

# Convertir genlight a formato array de LinkImpute para imputar los datos
gl2array <- function(gl, out.name) {
  require(data.table)
  
  x <- as.matrix(gl)
  
  fwrite(x = x, file = out.name, sep = ' ', quote = FALSE, col.names = FALSE, na = '-1')   
}
          

# read multiple files from BLINK results

read.blink.files <- function(blink.dir, l.traits) {
  
  l.t2 <- l.traits %>% str_c(collapse = '|')
  
  reg.pat <- str_c('GAPIT.Association.GWAS_Results.BLINK.', '(', l.t2, ')\\.')
  
  list.blink.files <- list.files(blink.dir, pattern = reg.pat, full.names = TRUE)
  
  ## leer los archivos y ajustar los P.valores
  ## los resultados ya vienen ajustados por fdr, falta agregar bonferroni
  
  variables <- str_extract(list.blink.files, '(?<=BLINK\\.).+(?=\\.csv)') %>% str_remove('.csv')
  
  names(list.blink.files) <- variables
  # list.blink.files %>% enframe() %>% view
  
  pvalues.blink <- vector(mode = 'list', length = 9)
  names(pvalues.blink) <- variables
  
  
  for (i in variables) {
    print(i)
    # leer la tabla
    df.blink <- read_csv(list.blink.files[i], show_col_types = FALSE)
    
    # agregar p.valores corregidos por bonferroni
    df.blink <- df.blink %>% mutate(bonferroni.adjusted.p = p.adjust(P.value, method = 'bonferroni'))
    
    #guardar tabla en lista
    pvalues.blink[[i]] <- df.blink
    
  }
  
  
  
  all.samples.blink <- bind_rows(pvalues.blink, .id = 'variable') %>% 
    # mutate(Chr = as.numeric(Chr)) %>% 
    select(variable, AlleleID = SNP, Chr, Pos, MAF, p.value = P.value, FDR.adjusted.p = `H&B.P.Value`,
           bonferroni.adjusted.p)

  return(all.samples.blink)

}


          
