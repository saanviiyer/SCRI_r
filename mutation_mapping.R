list_of_topics <- list()

for(i in 1:min(50, ncol(df))){
  column_list <- as.list(df[, i])
  list_of_topics[[paste0("topic", i)]] <- column_list
}

# Print the names of the lists created
print(names(list_of_topics))

# each topic has list of genes

mutations <- c(
  "CTNNB1", "DDX3X", "SMARCA4", "CREBBP", "KMT2D",
  "MYC", "GABRA5", "GFI1", "GFI1B", "OTX2",
  "BMI1", "SOX2", "NFATC4", "GALNT14", "APC",
  "TP53", "TOP2A", "CDK1", "RRM2", "STMN2",
  "KIF5C", "SYT11", "NME2", "HK2", "PGM5",
  "LRP4", "APCDD1", "JUNB", "EGR1", "DKK2",
  "AXIN2", "WIF1", "ZIC1", "PAX6", "OLIG3",
  "NKD1", "BARHL1", "PDE1C", "PCSK9", "PTF1A",
  "NEUROG1", "ASCL1", "FOXD3", "BRN3A", "CBFA2T2",
  "CBFA2T3", "PRDM6", "UTX"
)

# Iterate through each mutation in the mutations genes
for (mutation in mutations) {
  # Check in which topics the current mutation is present
  is_in_list <- sapply(list_of_topics, function(topic) mutation %in% topic)
  
  # Find the indices of topics containing the mutation
  topics_with_mutation <- which(is_in_list)
  
  # Print the result for the current mutation
  if (length(topics_with_mutation) > 0) {
    topic_names <- paste("Topic", topics_with_mutation, collapse = ", ")
    print(paste("The mutation", mutation, "is present in:", topic_names))
  } else {
    print(paste("The mutation", mutation, "is not present in any topic."))
  }
}
