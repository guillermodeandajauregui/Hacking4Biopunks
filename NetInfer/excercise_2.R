# excercise 2: MI calculation----

library(tidyverse)
library(infotheo)

# read data ----

x <- read_delim(file = "gene_data.txt", delim = "\t")

# make it a matrix 
y <- 
  x |> 
  select(-c(1:2)) |> 
  as.matrix()

rownames(y) <- x$hgnc_symbol

# discretize data ----

## we will convert each row into a discrete variable 

y_list <- split(y, row(y))
names(y_list) <- rownames(y)
  
y_discretized <- 
  lapply(y_list, function(i){
  discretize(X = i, 
             disc = "equalfreq", #extra: try other options
             nbins = round(length(i)^(1/3))) |> #extra: try other values
      unlist() |> as.numeric()
                        }
  )


# calculate MI ----

mi_matrix <- 
  sapply(X = y_discretized, FUN = function(i){
  sapply(X = y_discretized, FUN = function(j){
    mutinformation(X = i, Y = j, method = "emp")
  })
})


# pheatmap::pheatmap(mi_matrix)

# mi_matrix |> 
#   as.data.frame() |> 
#   rownames_to_column(var = "gene") |> 
#   as_tibble() |> 
#   vroom::vroom_write("mi_matrix.txt")
