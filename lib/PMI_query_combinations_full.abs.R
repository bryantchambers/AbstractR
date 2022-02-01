#This function reads a ```query_object``` output by the ```query_result``` function. The result is a list of 
#parwise mututal information scores resulting from 

PMI_query_combinations_full.abs = function(query_result_object, chemical.list.path, search.min = 0, search.max = 100000, num.cores = 1) {

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
        pull(pmid) %>% unique() %>% length() }, mc.cores = num.cores) %>% setNames(names(query.result[[1]])) %>% unlist() %>%
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