
dart.extract.header <- function(path.report, type = 'DART', make.names.unique = FALSE){
  
  # Usage: Take the path to a dart report and turn the header into a long table of samples
  # Returns a data.frame with all of the information from the header of a dart report
  # Arguments:
  ## path.report (String): Path to a standard dartseq report (Either SNP-one row, SNP-two row or SilicoDArT)
  ## make.names.unique (Logical): Adds a column with the following information joined into a single string:
  ##                              - Order_number
  ##                              - Plate number
  ##                              - Position in plate, row, then column (example: F1) 
  ##                              - Sample name
  ##                              If some names remain equal, then prints a warning with the repeated samples
  
  require(data.table)
  require(dtplyr)
  require(dplyr)
  
  #Read the first lines of the report to extract the header
  head <- fread(file = path.report, check.names = FALSE, nrows = 10, header = FALSE, stringsAsFactors = FALSE)
  
  # The information we require is on the header of the dart report
  nkeep <- sum(head[, 1] == "*") + 1 #how many rows does the header have
  
  
  ## Normally, the header of a DArT report has the following. But always check the metadata file and the data frame itself!
  ## Rows of header vary between providers (DArT PL or SAGA)
  # 1		Order number where sample belongs to - important for multi-orders reports	
  # 2		DArT plate barcode	
  # 3		client plate barcode	
  # 4		well row position	
  # 5		well column position	
  # 6		sample comments	
  # 7		genotype name	
  
  if (type == 'DART') {
    col_names <- c('order_number',
                   'dart_plate_barcode',
                   'client_plate_barcode',
                   'row_pos',
                   'col_pos',
                   'sample_comment',
                   'id')
    
  }
  
  if (type == 'SAGA') {
    
    col_names <- c('plate_pos', 
                   'dart_plate_barcode', 
                   'client_plate_barcode', 
                   'unknown_header_column', 
                   'sample.name', 
                   'id')
    
  }
  
  dat2 <- head[1:nkeep,] %>% #select only the required rows
    remove_rownames() %>%
    lazy_dt() %>% 
    mutate(col_names = col_names) %>%       #create a column with with column names
    #the next to lines reorganize the df
    pivot_longer(cols = -col_names, names_to = 'variable', values_to = 'value') %>% 
    pivot_wider(names_from = col_names, values_from = value) %>% 
    filter(dart_plate_barcode != '*') %>% 
    select(-variable) %>% 
    as.data.table()
  
  if (make.names.unique) {
    
    if (type == 'DART'){
      dat2 <- dat2 %>% 
        mutate(unique.name = str_c(order_number, client_plate_barcode, str_c(row_pos, col_pos, sep = ''), id, sep = '_'))
    }
    
    if (type == 'SAGA'){
      dat2 <- dat2 %>% 
        mutate(unique.name = str_c(client_plate_barcode, plate_pos, id, sep = '_')) 
    }
    
    
    #Verificar que todos los nombres sean unicos
    check.uniques <- length(dat2$unique.name) ==  length(unique(dat2$unique.name))
    
    # Imprimir nombres repetidos
    if (!check.uniques) {
      
      warning('\n\nSome names remain duplicated! \n Printing list of samples\n\n')  
      
      rep.names <- dat2 %>% 
        count(unique.name, sort = TRUE, name = 'times repeated') %>% 
        filter(`times repeated` > 1) %>% 
        as.data.frame()
      
      print(rep.names)
      
    }
  }
  
  return(dat2)  
  
} # fin de dart.extract.header
