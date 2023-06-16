# ------------------------------------------------------------------------------
# SIMPLE HELPERS
string_to_integer <- function(x) {
  x_char <- as.character(x)
  x_factor <- as.factor(x_char)
  x_numeric <- as.numeric(x_factor)

  as.integer(x_numeric)
}

# ------------------------------------------------------------------------------
# FORMATTING AND PRINTING
precis_kable <- function(m, ...) {
  precis(m, ...) %>%
    as_tibble(rownames = "coefficient") %>%
    knitr::kable()
}

# ------------------------------------------------------------------------------
# MATH
exponential_1 <- function(age, y0, lambda) {
  y0 * exp(-1 * lambda * age)
}

exponential_delay <- function(age, baseline, lambda, delay) {
  if_else(
    age < delay,
    baseline,
    exponential_1(age - delay, baseline, lambda)
  )
}

# ------------------------------------------------------------------------------
# PLOTTING
plot_data <- function(d, guide = "none") {
  d %>%
    ggplot(aes(x = age, y = volume)) +
    geom_line(aes(color = rosetta, group = rosetta_laterality)) +
    geom_point(aes(color = rosetta)) +
    scale_color_discrete(guide = guide, name = "Individual") +
    scale_x_continuous(breaks = seq(0, 80, by = 10), limits = c(0, NA))
}
# ------------------------------------------------------------------------------
# DRAWING SAMPLES FROM THE MODEL
tidy_samples <- function(m, index_map = NULL, index_names = NULL, ndraws = 10000, pars, fun = extract.samples, ...) {
  # # TASK: If ndraws = NULL, return all draws.
  # if(is.null(ndraws)) {
  #   ndraws <-
  # }

  samples <- fun(m, n = ndraws, pars = pars, ...)

  parameter_names <- names(samples)

  result_list <- list()

  # Get the actual number of draws
  ndraws_actual <- dim(samples[[1]])[[1]]

  index_map_with_draws <-
    index_map %>%
    expand_grid(.draw = 1:ndraws_actual)

  if(!is.null(index_map)) {
    result_list[[1]] <- index_map_with_draws
  }

  for (i in parameter_names) {
    # Pull out the array
    parameter_samples <- samples[[i]]
    # For each variable, check the dimensions of the array
    n_dim <- length(dim(parameter_samples)) - 1

    # First dimension is always the number of draws from the posterior
    col_names <- c(".draw")
    # The second dimension is the d=2 indexing variable
    # I presume that the third dimensions is the d=3 indexing variable
    if (n_dim > 0) {
      # TASK: Check that there are enough index names to match the number of
      #.      indices needed.
      col_names <- c(col_names, index_names[which(i == names(index_names))])
    }

    col_names <- c(col_names, i)

    # Put the array in long format tibble
    # Rename the value with the name of the variable
    result_list[[i]] <-
      reshape2::melt(parameter_samples, value.name = i) %>%
      # For every dimension above one, rename the dimension with the appropriate
      # indexing variable name
      set_names(col_names)
  }

  result_list %>%
    reduce(inner_join) %>%
    as_tibble()
}


