

random.string <- function(n=1, lenght=12)
  {
    randomString <- c(1:n)                  # initialize vector
    for (i in 1:n)
      {
        randomString[i] <- paste(sample(c(0:20, letters, LETTERS),
                                        lenght, replace=TRUE),
                                 collapse="")
      }
    return(randomString)
  }
