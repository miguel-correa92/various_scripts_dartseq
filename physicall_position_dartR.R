library(dartR.base)
# Test setup
#load('~/projects/Pathology_25/ANT/target_gl_filt.robj')
#gl <- target_gl_filt
#ref <- '/catalogue/InterspecificCommonBeanDenovoGenomeAssembly/5.Others/references/genomes/Pvulgaris/G19833/v2.1/assembly/Pvulgaris_442_v2.0.fa'
#minimap2_path  <- '/home/scruz/.conda/envs/purge_haplotigs/bin/minimap2'
#samtools_path  <- '/home/scruz/.conda/envs/purge_haplotigs/bin/samtools'


#' gl.get_physical_position
#'
#' An R funtion for UNIX to get the physical position based on a input reference genome (in fasta)
#' using blastn as a engine
#'
#' @param gl genlight object with tag sequence in loc.metrics
#' @param ref_genome path to reference genome fasta file 
#' @param task
#' @param identity_treshold
#' @param evaluf
#' @param overlap_percentage
#' @param bitscore
#' @param threads
#' @param bitscore
#'
#' @return
#' @examples

gl.get_physical_position <- function (gl, ref_genome, minimap_bin, samtools_bin, number_of_threads = 20) 
{
    # VALIDATIONS
   if (!file.exists(ref_genome)){
     cli::cli_abort("`ref_genome` file not found in {ref_genome}")
   }

   if (!file.exists(minimap_bin)){
     cli::cli_abort("`minimap_bin` file not found in {minimap_bin}")
   }

   if (!file.exists(samtools_bin)){
     cli::cli_abort("`samtools_bin` file not found in {samtools_bin}")
   }
   if (!rlang::is_integerish(number_of_threads)){
     cli::cli_abort("`number_of_threads` must be an integer not {number_of_threads}")
   }

  tag_cols <- c('TrimmedSequence', 'AlleleSequence')

  # Check if is a genlight object `gl`
  if (class(gl)[1] == "genlight" | class(gl)[1] == "dartR") {
    # Check if tagcolumn is present in loc.metrics
    if (sum(tag_cols %in% colnames(gl$other$loc.metrics)) == 0) {
      cli::cli_abort("Fatal Error: TrimmedSequence or AlleleSequence column is required!")
    } else {
      tag_seq_column <- ifelse(tag_cols[1] %in% colnames(gl$other$loc.metrics), tag_cols[1], tag_cols[2])
      # Write tags into a fasta format
      unique_tags <- gl@other$loc.metrics %>%
              distinct(clone, !!sym(tag_seq_column))

      seq_id <- paste0(">",c(unique_tags$clone))
      nuc_seq  <- unique_tags[, tag_seq_column]

      fasta.input <- c(rbind(seq_id, as.character(nuc_seq)))
      tmp_path  <- paste0(tempdir(), "/fasta.input")

      write.table(fasta.input, file = tmp_path, 
                quote = FALSE, row.names = FALSE, col.names = FALSE)

      cli::cli_inform("`tags_fasta` saved in {tmp_path}")
    }
  }

  # CHECK if unix
  if (grepl("unix", .Platform$OS.type, ignore.case = TRUE)) {
# MINIMAP2 APPROACH
    # Check if blastn exists
        if(file.exists(minimap_bin)){
                cli::cli_inform("Checking works minimap2 located in: {minimap_bin}")
                system(glue::glue("{minimap_bin} --version", timeout = 2))                      
                tmp_path_minimap_out  <- paste0(tempdir(), "/tags.sorted.bam")
                cli::cli_inform('Bam output stored in {tmp_path_minimap_out}')

                minimap_cmd <- glue::glue("{minimap_bin} -t {number_of_threads} -ax sr",
                                         "{ref_genome} {tmp_path}",
                                         "| {samtools_bin} sort -o {tmp_path_minimap_out}",
                                        .sep = ' ')
                # Execute minimap2 job
                system(minimap_cmd)
                # ADD the MD field
                tmp_path_minimap_md_out  <- paste0(tempdir(), "/tags.sorted.md.bam")
                minimap_cmd <- glue::glue("{samtools_bin} calmd -b {tmp_path_minimap_out} {ref_genome} > {tmp_path_minimap_md_out}")
                system(minimap_cmd)
                # Filter the resulting bam
                tmp_path_filtbam_out  <- paste0(tempdir(), "/tags.sorted.md.filtered.bam")
                samtools_cmd  <- glue::glue("{samtools_bin} view -b -q 20 -F 0x100 -o {tmp_path_filtbam_out} {tmp_path_minimap_md_out}")
                cli::cli_inform('Filtering unique matchs of BAM stored in {tmp_path_filtbam_out}')
                system(samtools_cmd)
                # Get the interesing fields of filtered BAM
                tmp_path_table_out  <- paste0(tempdir(), "/tags.sorted.md.filtered.txt")
                samtools_cmd  <- glue::glue("{samtools_bin} view {tmp_path_filtbam_out} | cut -f1-6,22 > {tmp_path_table_out}")
                cli::cli_inform('Generating input table in {tmp_path_table_out}')
                system(samtools_cmd)
                tab_colnames <- c('clone', 'flags', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'MD')
                mapping_tab <- read.table(tmp_path_table_out, sep = '\t', header = F)
                colnames(mapping_tab) <- tab_colnames
                
                # Merge locmetrics with mapping_tab
                join_tab <- base::merge(gl@other$loc.metrics, mapping_tab, by = 'clone', all.x = T) %>%
                  mutate(allele_test = stringr::str_sub(!!sym(tag_seq_column), SnpPosition+1, SnpPosition+1),
                    query_len = stringr::str_count(!!sym(tag_seq_column))
                  )


                pred_positions <- purrr::pmap(list(join_tab$CIGAR,join_tab$flags, join_tab$query_len, join_tab$POS, join_tab$SnpPosition), query2ref)
                vec_positions <- purrr::map_vec(pred_positions, ~ifelse(is.null(.x), NA, .x))
                join_tab$snp_position  <-  vec_positions

                return(join_tab)
        } else {
                cli::cli_abort("minimap2 executable not found...")
        }
  } else {
        cli::cli_abort("{.Platform$OS.type} not supported only unix!")
  }
}

query2ref <- function(cigar, flags, query_len, ref_start, target_query_pos){
  # Convert zero pos to 1 based pos
  target_query_pos  <- target_query_pos  + 1
  # Parse CIGAR string into lengths and operations
  cigar_tuples <- unlist(regmatches(cigar, gregexpr("\\d+[MIDNSH]", cigar)))
  # Initialize Positions
  q_pos  <- 0
  r_pos  <- ref_start - 1
  if(!is.na(flags)){
    # If tag is mapped reverse
    if(flags == 16){
      target_query_pos  <- query_len - target_query_pos + 1
      cigar_tuples  <- rev(cigar_tuples)
    }
  }
  # Loop over parsed CIGAR elements
  for (element in cigar_tuples){
    len  <- as.numeric(gsub("[MIDNSH]", "", element))
    op  <- gsub("\\d+", "", element)
    switch(op,
      M = {
        if(q_pos + len >= target_query_pos){
          return(r_pos + (target_query_pos - q_pos))
        } else {
          q_pos  <- q_pos + len
          r_pos  <- r_pos + len
        }
      },
      I = {
        if(q_pos + len >= target_query_pos){
          cli::cli_warn("Target Position located over an insertion, reporting left most pos")
          return(r_pos)
        } else {
          q_pos  <- q_pos + len
        }
      },
      D = {
        r_pos  <- r_pos + len
      },
      N = {
        r_pos  <- r_pos + len
      },
      S = {
        if(q_pos + len >= target_query_pos){
          cli::cli_warn("Target Position located over an softclip, reporting left most pos")
          return(r_pos)
        } else {
          q_pos  <- q_pos + len
        }
      },
      H = {
        cli::cli_warn("Found a hardclip not expected returning NA")
        return(NA)
      },

      cli::cli_abort('Position not found in the alignment')
    )
  }
}

