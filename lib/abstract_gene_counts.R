#This function reads csv with a list of chemical names and a query map and returns a list of dataframes for all combinations of queries in the data set
#this result can be input to two functions: ```PMI_matrix``` which will give you information about the information or strength of the search phrase relative to the total 

abstract_gene_counts = function(query_result_object) {    

    #Libraries
    suppressMessages(library(tidyverse))
    library(parallel)

    #helper functions
    remove.punc = function(x) {
        string.1 = gsub("(?!\\.)[[:punct:]]", "", x, perl=TRUE)
        string.2 = substr(string.1,1,nchar(string.1)-1) #%>% tolower()
        string.3 = str_split(string.2, " ") %>% unlist()
    return(string.3)
    }

    match.genes <- function(sentence, genelist) {
    gene.matches = sentence[sentence %in% genelist]
    return(gene.matches)
    }

    #testingpaths
    #query_result_object = "input/Laura-HTTR/HTTR-full-chemical-query.2021.07.26.RData"
    
    #read query result
    load(paste0(query_result_object))

    #set up genes dictionary - improve this in the future, allow for a list of homolouges in a csv from bioMart - for now use your H-M-R set up
    #read in coding genes
    #also a data dir or lib is needed
    gene.list = suppressMessages(suppressWarnings(read_tsv("~/ipynb.1/Projects/AbstractR/data/proteincodinggenesHGNC.txt"))) %>%
    pull(symbol) %>% unique()
    
    #read in non human gene lists homolouge lists
    nh.gene.list <- suppressMessages(read_csv("~/ipynb.1/Projects/AbstractR/data/genes.list_h.m.r.csv")) %>% filter(human %in% gene.list) %>% distinct(`.keep_all` = FALSE)
    
    #comparative lists
    mouse.rat.vec = c(nh.gene.list$human, nh.gene.list$mouse, nh.gene.list$rat) %>% unique()

    #filter some bad gene names
    bad.genes.q = c('a')
    mouse.rat.vec = mouse.rat.vec[!mouse.rat.vec %in% bad.genes.q]

    #take only pmids and abstract from the query resul
    abstracts = mclapply(seq_along(query.result[[1]]), function(x) {
        query.result[[1]][[x]] %>%
        select(pmid, abstract)
    }, mc.cores = 32) %>% setNames(names(query.result[[1]]))

    #isolate genes from all abstracts
    abstracts.genes.filtered = mclapply(seq_along(abstracts), function(y)
        lapply(1:nrow(abstracts[[y]]), function(x) abstracts[[y]][x,2] %>%
            remove.punc() %>% match.genes(mouse.rat.vec)) %>% unlist() %>% na.omit, mc.cores = 32) %>% setNames(names(abstracts))

                                      
    #remove empty gene matches
    abstracts.genes.filtered = abstracts.genes.filtered[names(lapply(abstracts.genes.filtered, length)[lapply(abstracts.genes.filtered, length) > 0])]

    #count the genes for all queried chemicals with genes mentioned in abstracts
    query.result.gene.counts = mclapply(seq_along(abstracts.genes.filtered), function(x) {
        abstracts.genes.filtered[[x]] %>%
            table() %>%
            data.frame() %>%
            setNames(c('gene', 'count')) %>%
            mutate(gene = as.character(gene)) %>%
                            mutate(rat = nh.gene.list$human[match(gene, nh.gene.list$rat)], mouse = nh.gene.list$human[match(gene, nh.gene.list$mouse)]) %>%
                            mutate(gene = ifelse(is.na(rat), gene, rat), gene = ifelse(is.na(mouse), gene, mouse)) %>%
            select(gene, count) %>%
            group_by(gene) %>%
            summarise(count = sum(count)) %>%
            arrange(desc(count))}, mc.cores = 32) %>% setNames(names(abstracts.genes.filtered))

    #calulate the final counts matrixs for all fo the samples
    abstract.gene.count.matrix = query.result.gene.counts %>% reduce(full_join, by = "gene") %>%
        arrange(gene) %>%
        column_to_rownames('gene') %>%
        t() %>%
        as.data.frame() %>%
        replace(is.na(.), 0) %>%
        rownames_to_column('sample') %>%
        select(-sample) %>%
        mutate(chemical = names(query.result.gene.counts)) %>%
        column_to_rownames('chemical')

return(abstract.gene.count.matrix)
                                        
}