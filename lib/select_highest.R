#This function takes a numeric matrix resulting from PMI, GSEA (connectivity), or some other metric and 
#' returns that same matrix with two additional columns: the name of the highest scoring column for each row and
#' the corresponding score

select.highest = function(score.matrix) {

    #Libraries
    suppressMessages(library(tidyverse))
    
    #testingpaths
    #score.matrix = score.matrix
    

    #maxcolumn and value
    score.matrix.highest = data.frame(max_column = colnames(score.matrix)[apply(score.matrix, 1, which.max)], stringsAsFactors = FALSE) %>%
    cbind(score.matrix %>% apply(1, max) %>% as.data.frame() %>% setNames("max_value") %>%
    rownames_to_column(var = "row"), .) %>%
    left_join(score.matrix %>% as.data.frame() %>% rownames_to_column('row'), by = 'row') %>%
    column_to_rownames('row')

return(score.matrix.highest)
}