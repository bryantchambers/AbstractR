#This function reads csv with a list of chemical names and a query map and returns a list of dataframes for all combinations of queries in the data set
#this result can be input to two functions: ```PMI_matrix``` which will give you information about the information or strength of the search phrase relative to the total 

query_result.ids = function(chemical.list.path, query.map.path = NULL, start.n = 1, end.n = NULL, limit_input_to = NULL, num.cores = 15) {

    #Libraries
    suppressMessages(library(tidyverse))
    suppressMessages(library(easyPubMed))
    suppressMessages(library(parallel))

    #testingpaths
    #chemical.list.path = "input/heparg-chems-input.csv"
    #query.map.path = "input/query.map-genes.in.sigs.csv"
    #query.map.path = NULL
    #limit_input_to = NULL

if(!is.null(query.map.path)) { 
    ######################
    #SEARCH WITH QUERY MAP
    
    #read chemical list
    chems = suppressMessages(read_csv(paste0(chemical.list.path)))
    #read query map
    q.map = suppressMessages(read_csv(paste0(query.map.path))) %>%
        mutate(phrase.q = paste0('"',phrase,'"'))
    
    #cutchems
    if(is.null(limit_input_to)) {chems = chems$chemicals[1:length(chems$chemicals)]} else {chems = chems$chemicals[1:limit_input_to]}

    #build partial search vectors
    chems.search = paste0('"', chems,'" AND ')
    searchterms.search = q.map$phrase.q

    #controlable and bindable search grid with directions
    search.combinations = expand.grid(chems.search, searchterms.search, stringsAsFactors = FALSE) %>%
    setNames(c('chemical', 'phrase.q')) %>% 
    left_join(q.map, by = "phrase.q") %>%
    mutate(query = paste(chemical, phrase.q)) %>%
    select(SRP, direction, query)

    #acutal query list
    search.query = search.combinations$query %>% as.list()
    
    } else {
    #######################    
    #SEARCH WITH ONLY CHEMS
        
    #read chemical list      
    chems = suppressMessages(read_csv(paste0(chemical.list.path)))
    
    #cutchems
    if(is.null(limit_input_to)) {chems = chems$chemicals[1:length(chems$chemicals)]} else {chems = chems$chemicals[1:limit_input_to]}

    #build partial search vectors
    chems.search = paste0('"', chems,'"')
        
    #actual query list
    search.query = chems.search
    }

    #build output object
    query.result = vector(length = 3, mode = "list")
     
     ####################################################################
     ####################################################################
     
     #Queries
     
     #ids only
     options(digits.secs = 12)
     #start.n = 1
     #end.n =   length(search.query)
     #num.cores = 25

     t1.s = Sys.time()
     query.result[[1]] = mclapply(seq_along(search.query[start.n:end.n]), function(x){

          Sys.sleep(time = 0.011) #* num.cores)   #if system is killing  try to increase?? network stability at EPAaa

          query.res.i = get_pubmed_ids(search.query[[x]], api_key = 'b3accc43abfb376d63b0c2d2f7d8f984de09')

          ID.list = query.res.i$IdList %>% unlist(use.names = FALSE) %>% unique()

          df.query.i = data.frame(pmid = ID.list)

          if(nrow(df.query.i) < 1) {df.query.i = 'EMPTY'}

          return(df.query.i)


     }, mc.cores = num.cores) %>% setNames(search.query[start.n:end.n]) #change number of cores to keep below 10 so that time out dose not occur

     #timer
     t2.e = Sys.time()
     difftime(time1 = t2.e, time2 = t1.s, units = "secs")


     #clean out empty hits
     query.result[[1]] = query.result[[1]][query.result[[1]] != 'EMPTY']

     #save additional mapping information if a map was used
     if(is.null(query.map.path)) {query.result[[2]] = 'NO MAP USED'} else {query.result[[2]] = search.combinations}
     if(is.null(query.map.path)) {query.result[[3]] = 'NO MAP USED'} else {query.result[[3]] = q.map}
     
     return(query.result)
     
}