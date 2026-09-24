if (!require("ontologyIndex")) install.packages("ontologyIndex")
library("ontologyIndex")

if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)

# get relationship names used in the OBO file
obo_relationship_names <- get_relation_names('~/data/external/gene_ontology/go.obo')

go_relation_trees <- list()

parent_child_relations <- data.frame(child_id = character(),
                                     parent_id = character(),
                                     relation_class = character())
ancestors_relations <- data.frame(id = character(),
                                  ancestor = character(),
                                  relation_class = character())

for (relation in obo_relationship_names){
  print(relation)
  go_relation_tree <- get_ontology('~/data/external/gene_ontology/go.obo', propagate_relationships = relation)
  go_relation_trees[[relation]] <- go_relation_tree
  go_relation_df <- as.data.frame(go_relation_tree)
  go_relation_df <- go_relation_df %>%
    filter(obsolete == 'FALSE')
  
  # all IDs at the same level as and upstream of a given ID
  print('Getting ID-ancestor relationships.')
  ancestors_filtered <- go_relation_df[c('id', 'ancestors')]
  rownames(ancestors_filtered) <- 1:nrow(ancestors_filtered)
  
  ancestors_relation_final = data.frame(id = character(),
                                        ancestor = character())
  for (i in 1:nrow(ancestors_filtered)){
    id = ancestors_filtered[i, 'id']
    ancestors = str_split(ancestors_filtered[i, 'ancestors'], '; ')[[1]]
    ancestors_df = data.frame(id = rep(id, length(ancestors)),
                              ancestor = ancestors)
    ancestors_relation_final <- rbind(ancestors_relation_final, ancestors_df)
  }
  ancestors_relation_final$relation_class = relation
  ancestors_relations <- rbind(ancestors_relations, ancestors_relation_final)
}

write_csv(ancestors_relations, '~/data/external/gene_ontology/go_id_ancestor_relations.csv')

# get names and aspects of all GO IDs
go_relation_tree <- get_ontology('http://purl.obolibrary.org/obo/go.obo') 
id_name_df <- data.frame(id = go_relation_tree$id,
                         name = go_relation_tree$name) %>%
  filter(grepl('GO:', id))
# aspects from manually reading the OBO file
go_obo <- data.frame(V1 = readLines('https://purl.obolibrary.org/obo/go.obo')) %>%
  filter(grepl('id: GO', V1) | 
           grepl('namespace: biological_process', V1) |
           grepl('namespace: molecular_function', V1) |
           grepl('namespace: cellular_component', V1)) %>%
  filter(!grepl('alt_id', V1))

go_aspects <- rep('', length = nrow(id_name_df))
for (i in seq_len(nrow(id_name_df))) {
  go_id <- id_name_df$id[i]
  go_obo_idx <- which(go_obo$V1 == paste0('id: ', go_id))
  
  if (length(go_obo_idx) > 1) {
    print(paste(go_id, "appears more than once in go_obo."))
    break
  }
  
  go_obo_idx <- go_obo_idx[1]
  aspect <- gsub('namespace: ', '', go_obo$V1[go_obo_idx + 1])
  go_aspects[i] <- aspect
}
id_name_df$aspect <- go_aspects
write_csv(id_name_df, '~/data/external/gene_ontology/go_id_names_aspects.csv')
