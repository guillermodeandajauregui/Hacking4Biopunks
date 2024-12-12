# excercise 1: random variables (dice) ----

library(tidyverse)
library(infotheo)

# make two fair dice (red and blue) ----
# Each die has outcomes from 1 to 6

set.seed(725)
{red_dice <- sample(1:6, size = 50000, replace = TRUE)  
blue_dice <- sample(1:6, size = 50000, replace = TRUE) 
}

dice_data <- data.frame(red = red_dice, blue = blue_dice) |> as_tibble()
# plot probability density function (uniform) ---- 

hist(red_dice, col="red")
hist(blue_dice, col = "blue")

# plot heatmap of outcomes

dice_data |> 
  group_by(red, blue) |>  
  tally() |> 
  ungroup() |> 
  mutate(prob = 100*n/sum(n)) |> 
  ggplot() + 
  aes(x = red, y = blue, fill = prob, label = prob)  + 
  geom_tile() + 
  geom_text() #+ 
  #scale_fill_gradient(low = "lightpink", high = "lightpink")

# make mirror dice 

dice_data <- 
  dice_data |> 
  mutate(mirror = blue)

dice_data |> 
  group_by(mirror, blue) |>  
  tally() |> 
  ungroup() |> 
  complete(mirror, blue, fill = list(n=0)) |> 
  mutate(prob = 100*n/sum(n)) |> 
  ggplot() + 
  aes(x = blue, y = mirror, fill = prob, label = prob)  + 
  geom_tile() + 
  geom_text() 

# calculate mutual information ----

mutinformation(X = red_dice, Y = blue_dice)
mutinformation(X = dice_data$mirror, Y = blue_dice)
