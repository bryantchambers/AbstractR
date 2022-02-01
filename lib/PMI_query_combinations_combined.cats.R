#This function reads a ```query_object``` output by the ```query_result``` function. The result is a list of 
#parwise mututal information scores resulting from 

PMI_query_combinations = function(query_result_object, chemical.list.path, search.min = 0, search.max = 100000, num.cores = 1) {

    #Libraries
    library(tidyverse)
    library(parallel)
    library(magrittr)

    #testingpaths
    #  query_result_object = "input/build/pmi.testing.query.object.RData"
    # chemical.list.path = "input/build/pmi.testing.chem.list.csv"
    # search.min = 5
    # search.max = 100000
    #  num.cores = 15

        #read query result
    load(paste0(query_result_object))
    chems.data = suppressMessages(read_csv(paste0(chemical.list.path))) %>%
        mutate(q.chemical = paste0('"',chemicals,'"'))

#take only pmids and abstract from the query result and parse by/ aggregate by SRP/query
search_matrix.i = mclapply(seq_along(query.result[[1]]), function(x) { query.result[[1]][[x]] %>%
        select(pmid) %>% distinct() %>% setNames(c("pmid")) %>% mutate(query = names(query.result[[1]][x]))}, mc.cores = num.cores) %>%
          do.call(rbind, .) %>%
         left_join(query.result[[2]], 'query') %>%
         mutate(search = query) %>%
        separate(col = search, into = c('q.chemical', 'search'), sep = ' AND  ') %>%
        left_join(chems.data, 'q.chemical') %>%
        mutate(search = gsub('\\"','', search)) %>%
        group_split(SRP)

 # aggregate by SRP and chemcial name to remove non-inique PMIDs and remove negative query phrases

search_matrix = mclapply(seq_along(search_matrix.i), function(y){
     search_matrix.i[[y]] %>% ungroup() %>%
     group_by(SRP, chemicals, pmid) %>%
     summarise(dir.tot = direction %>% unique %>% sort %>% paste(collapse = ", ")) %>% #gets rid of negative querried PMIDs
     filter(!grepl(pattern = "neg", x = dir.tot)) %>%
     ungroup() %>%
     group_by(SRP, chemicals) %>%
     summarise(count = pmid %>% unique() %>% length) %>% #takes the unqiue PMID vector length as n unique searches
     ungroup()
}, mc.cores = num.cores) %>%
     do.call(rbind, .) %>%
     pivot_wider(names_from = SRP, values_from = count, values_fill = 0) %>%
     filter(!is.na(chemicals)) %>%
     column_to_rownames("chemicals") %>%
     mutate(rsum = rowSums(across(where(is.numeric))))


    #filter based on search returns
    search_matrix.wide.f = search_matrix %>%
        filter(rsum >= search.min & rsum <= search.max) %>%
        select(-rsum)


    #build a result matrix
    search_PMI = matrix(nrow = nrow(search_matrix.wide.f), ncol = ncol(search_matrix.wide.f)) %>%
        magrittr::set_rownames(rownames(search_matrix.wide.f)) %>%
        magrittr::set_colnames(colnames(search_matrix.wide.f))

    #Calculate PMI
    for (i in 1:nrow(search_matrix.wide.f)) {
        for (j in 1:ncol(search_matrix.wide.f)){
        search_PMI[i, j] <- (log2(sum(search_matrix.wide.f[i,])/sum(search_matrix.wide.f))
            + log2(sum(search_matrix.wide.f[,j])/sum(search_matrix.wide.f))
            -log2(search_matrix.wide.f[i,j]/sum(search_matrix.wide.f))) *-1
        }
    }
  

    #clean up infs
    search_PMI = search_PMI %>%
    replace(is.infinite(.), 0)

  
    return(search_PMI)
    
}