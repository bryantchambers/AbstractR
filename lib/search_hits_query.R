#This function reads a ```query_object``` output by the ```query_result``` function. The result is a list of 
#parwise mututal information scores resulting from 

search_hits_query = function(query_result_object, chemical.list.path, search.min = 0, search.max = 100000) {

    #Libraries
    library(tidyverse)
    library(parallel)
    library(magrittr)

    #testingpaths
#     query_result_object = "../AbstractR-projects/LINCS/output/query.result_LINCS_2021-07-29.RData"
#     chemical.list.path = "input/chemical-list-LINCS-2021.07.23.csv"
#     search.min = 5
#     search.max = 100000
    

    #read query result
    load(paste0(query_result_object))
    chems.data = suppressMessages(read_csv(paste0(chemical.list.path))) %>%
        mutate(q.chemical = paste0('"',chemicals,'"'))

#take only pmids and abstract from the query result and calc number results
search_matrix = mclapply(seq_along(query.result[[1]]), function(x) { query.result[[1]][[x]] %>%
        pull(pmid) %>% unique() %>% length() }, mc.cores = 10) %>% setNames(names(query.result[[1]])) %>% unlist() %>%
        as.data.frame() %>%
        setNames('search_results') %>%
        rownames_to_column('query') %>%
        left_join(query.result[[2]], 'query') %>%
        mutate(search = query) %>%
        separate(col = search, into = c('q.chemical', 'search'), sep = ' AND  ') %>%
        left_join(chems.data, 'q.chemical') %>%
        mutate(search = gsub('\\"','', search)) %>%
        select(chemicals, SRP, search, direction, search_results)

    #make a wide form matrix
    search_matrix.wide = search_matrix %>%
        select(chemicals, search, search_results) %>%
        pivot_wider(id_cols = chemicals, names_from = search, values_from = search_results) %>%
        replace(is.na(.), 0) %>%
        mutate(rsum = rowSums(across(where(is.numeric))))

    #filter based on search returns
    search_matrix.wide.f = search_matrix.wide %>%
        filter(rsum >= search.min & rsum <= search.max) %>%
        select(-rsum) %>%
        column_to_rownames('chemicals')



  
    return(search_matrix.wide)
    
}