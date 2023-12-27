#!/usr/bin/env Rscript

# charnov's model - individual-based
library("tidyverse")

# slope of male gain curve
n <- 0.2

# initial value of the amount of allocation towards males
r_init <- 0.9

# the population size
popsize <- 5000

# the mutation rate
mu_r <- 0.02

# the standard deviation of the mutational effect size
sdmu_r <- 0.02

# skip output every n generation
skip <- 1000

# max duration of the simulation
max_time <- 10000

wf <- function(r) {
    return(1.0 - r)
}

wm <- function(r) {
    return(r^n)
}


#clamp <- function(vec, min, max) { pmax( min, pmin( vec, max))}

# get a bunch of individuals as a list of r values
individuals <- rep(r_init,length.out=popsize)

# allocate a tibble to track mean and variance in r
data <- as_tibble(cbind(time=1:max_time, mean_r=0.0, val_r=0.0))

clamp <- function(x) { ifelse(x > 1.0, 1.0, ifelse(x < 0.0, 0.0, x))}

mutate <- function(t) {

    # calculate which individuals mutate at rate mu_r 
    # as a boolean vector
    mutate_events = rbinom(n=popsize, prob=mu_r, size=1)

    # obtain number of mutants
    n_mut <- sum(mutate_events)
    
    # convert to boolean 
    mutate_events <- as.logical(mutate_events)

    # calculate mutational effect size for those mutants
    # (continuum of alleles)
    mutational_effects <- rnorm(mean=0,sd= sdmu_r,n=n_mut)

    stopifnot(length(mutational_effects) == length(individuals[mutate_events]))

    # add mutational effects to selected individuals
    new_values <- individuals[mutate_events] + mutational_effects
    
    new_values <- sapply(X=new_values,FUN=clamp)
    
    individuals[mutate_events] <<- new_values
}

reproduce <- function(t) {

    # calculate fitnesses
    wf_vals <- sapply(individuals, wf)
    wm_vals <- sapply(individuals, wm)

    wf_vals <- wf_vals / mean(wf_vals)
    wm_vals <- wm_vals / mean(wm_vals)

    # sample next gen - sample with replacement which reflects competition among gametes
    ova <- sample(x=individuals, size=popsize, replace=T, prob=wf_vals)
    pollen <- sample(x=individuals, size=popsize, replace=T, prob=wm_vals)

    # meiosis
    inherit_mom <- as.logical(rbinom(n=popsize, prob=0.5, size=1))
    inherit_dad <- !inherit_mom

    stopifnot(sum(c(inherit_mom,inherit_dad)) == popsize)

    # obtain the new individuals. We have one big list of gametes 
    # and we need to sample from that. Note that the Fisher condition is met
    # as we have equal numbers of ova and pollen from which we sample
    new_individuals <- sample(c(ova[inherit_mom],pollen[inherit_dad]))

    stopifnot(length(new_individuals) == popsize)

    # build next generation
    individuals <<- new_individuals
}

stats <- function(time_step) {
    mean <- mean(individuals)
    var <- var(individuals)

    data[time_step,] <<- list(time_step, mean, var)
}

# now have them reproduce
for (generation in 1:max_time)
{
    mutate(generation)
    reproduce(generation)
    stats(generation)
    
    if (generation %% skip == 0)
    {
        print(paste0("generation ",generation))
    }
}

ggplot(data=data
        ,mapping=aes(x=time, y=mean_r)) +
    geom_line() +
    ylim(0,1)

ggsave(filename="evolution_r.pdf")


