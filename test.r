rwm_chain <- function(N, s,x0) {
  x <- numeric(N)
  x[1]<- x0
  
  for (i in 2:N) {
    x_star <- rnorm(1, mean = x[i - 1], sd = s)
    
    rfunction <- ((1/2) * exp(-abs(x_star))) /
      ((1/2) * exp(-abs(x[i - 1])))
    
    u <- runif(1, 0, 1)
    
    if (u < rfunction) {
      x[i] <- x_star
    } else {
      x[i] <- x[i - 1]
    }
  }
  
  return(x)
}

y <- rwm_chain(10000, 1,0)


# Grid and true density
x_grid <- seq(-8, 8, length.out = 400)
true_density <- 0.5 * exp(-abs(x_grid))

# Histogram (density-scaled)
hist(y, breaks = 50, probability = TRUE,
     col = rgb(0, 0, 1, 0.5),
     border = "white",
     xlab = "x",
     ylab = "Density",
     main = "Random Walk Metropolis Sampling")

# Kernel density estimate
kde <- density(y)
lines(kde, lwd = 2, col = "red")

# True density
lines(x_grid, true_density, lwd = 2, col = "black")

# Legend
legend("topright",
       legend = c("Histogram", "Kernel density", "True density"),
       col = c("blue", "red", "black"),
       lwd = c(NA, 2, 2),
       pch = c(15, NA, NA),
       pt.cex = 2,
       bty = "n")

# Monte Carlo estimates
cat("Sample mean:", mean(y), "\n")
cat("Sample standard deviation:", sd(y), "\n")




set.seed(4567890)

N <- 2000
s <- 0.001
J <- 4
initial_values <- c(-10, -5, 5, 10)



s_values <- seq(0.001, 1, length.out = 30)
Rb_values <- numeric(length(s_values))

for (k in seq_along(s_values)) {
  s <- s_values[k]
  chains <- vector("list", J)
  
  for (j in 1:J) {
    chains[[j]] <- rwm_chain(N, s, initial_values[j])
  }
  
  # Means of each chain (discard first iteration, like chain[1:])
  Mj <- sapply(chains, function(chain) mean(chain[-1]))
  
  # Within-chain variances
  Vj <- sapply(seq_along(chains), function(j) {
    chain <- chains[[j]]
    mean((chain[-1] - Mj[j])^2)
  })
  
  W <- mean(Vj)
  M <- mean(Mj)
  B <- mean((Mj - M)^2)
  
  Rb <- sqrt((B + W) / W)
  Rb_values[k] <- Rb
}



plot(s_values, Rb_values, type = "b", pch = 19,
     xlab = "Step size s",
     ylab = "Rb",
     main = "Rb diagnostic for Random Walk Metropolis")

abline(h = 1.05, col = "red", lty = 2)

legend("topright",
       legend = "Rb = 1.05",
       col = "red",
       lty = 2,
       bty = "n")

print(Rb)
