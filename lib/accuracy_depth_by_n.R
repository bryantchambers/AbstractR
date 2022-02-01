accuracy_by_depth_n = function(score.data.object, true.colum, entity, top.n = 3){

#build data 
     #scored.data.object = scored.data.object.i
     #top.n = 3
    # true.colum = "chemical.stress.assignment"
    # entity = "chemical.name"

#Prepare accuracy matrix for accuracy function
     chemical.lists = score.data.object %>%
          group_split(.dots = as.symbol(entity))
          #group_split(vars(as.symbol(entity)))
     


#remap true categories
cat.mapping = score.data.object %>%
     select(any_of(c(true.colum, entity)))

#find score dept by 
score.depth.by = mclapply(seq_along(chemical.lists), function(x){
     
     #build an ordered matrix for a given chemcial's ranked scores
     ordered.rank = chemical.lists[[x]] %>% pivot_longer(cols = !(chemical.name | chemical.stress.assignment) , names_to = "SRP", values_to = "score") %>%
          arrange(desc(score))

     #check if the assigned chemicals are in the top N    
     for (i in 1:top.n) {
               var.name <- paste("n.", i)
               ordered.rank = ordered.rank %>%
                    mutate( "by_{{i}}" := paste0(.$SRP[1:i], collapse = " "))#, by2 = paste0(.$SRP[1:2], collapse = " "), by3 = paste0(.$SRP[1:3], collapse = " "))        
          }

     #build the last operable frame
     ordered.rank = ordered.rank %>% 
          head(n=1) %>%
          column_to_rownames("chemical.name")

     #set the category of interest
     cat.oi = ordered.rank$chemical.stress.assignment %>% unique()

     #find the T/F matrix
     ordered.rank %>%
          select(-score, -chemical.stress.assignment, -SRP) %>%
          mutate_all( ~ grepl(pattern = cat.oi, x = .))

     }, mc.cores =10) %>% do.call(rbind, .) * 1

#calculate counts of hits by depth
accuracy.count = score.depth.by %>% as.data.frame() %>%
     rownames_to_column('chemical.name') %>%
     left_join(cat.mapping) %>%
     group_by_at(true.colum) %>%
     summarize_if(is.numeric, sum) %>%
     ungroup() %>%
     column_to_rownames(paste0(true.colum))

#calculate total ents in each row
accuracy.totals = score.depth.by %>% as.data.frame() %>%
     rownames_to_column('chemical.name') %>%
     left_join(cat.mapping) %>%
     group_by_at(true.colum) %>%
     summarize_if(is.numeric, length) %>%
     ungroup() %>%
     column_to_rownames(paste0(true.colum))



#find the weighted accuracy
accuracy = accuracy.count/accuracy.totals * 100

return(accuracy)
}