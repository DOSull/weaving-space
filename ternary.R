tadd <- function(x, y) {return((x + y) %% 3)}
tmult <- function(x, y) {return((x * y) %% 3)}

m <- matrix(sample(0:2, 9, replace = TRUE), 3, 3)
x <- 0:2
outer(m, x, tadd)

df <- expand.grid(a = 0:2, b = 0:2, c = 0:2)
df <- df %>% 
  mutate(product = tmult(a, tmult(b, c)),
         sum = tadd(a, tadd(b, c))

A <- 0:2 %>%
  outer(0:2, tadd) %>%
  outer(0:2, tadd)

A[1,2,]
