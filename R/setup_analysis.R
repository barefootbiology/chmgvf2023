# Prevent Stan from recompiling programs that are unchanged.
rstan_options(auto_write = TRUE)
# Detect the number of cores available for parallel processing.
options(mc.cores = parallel::detectCores())

set_ulam_cmdstan(TRUE)

theme_set(theme_bw() + theme(panel.grid = element_blank()))

dataset_0 <-
  read_tsv(file.path(here::here(), "data", "chmgvfdata.tsv"), show_col_types = FALSE) %>%
  set_names(tolower) %>%
  unite(pid_laterality, c("pid", "laterality"), sep = "_", remove = FALSE) %>%
  mutate(mutation_label = if_else(arg253stop == 1, "Arg253Stop", "Other")) %>%
  mutate(notes = if_else(is.na(volume) & (logmar >= 2.4), "AssumeZero", as.character(NA))) %>%
  mutate(volume = if_else(is.na(volume) & (notes == "AssumeZero"), 0, volume)) %>%
  rename(acuity = snellen_acuity)

pid_order <-
  dataset_0 %>%
  select(pid) %>%
  distinct() %>%
  mutate(no = str_remove(pid, "P") %>% as.numeric()) %>%
  arrange(no) %>%
  pluck("pid")

dataset <-
  dataset_0 %>%
  mutate(pid = factor(pid, levels = pid_order))

gvf <-
  dataset %>%
  filter(!is.na(volume))

gvf_data <-
  gvf %>%
  # If the variables being indexed are factors, use:
  # rethinking::coerce_index(pid)
  # instead. Use string_to_integer(pid) since these
  # values are characters vectors.
  mutate(S = string_to_integer(pid)) %>%
  mutate(SE = string_to_integer(pid_laterality)) %>%
  mutate(V = volume) %>%
  mutate(T = age) %>%
  mutate(row_id = 1:n())

gvf_list <-
  list(
    S = gvf_data$S,
    SE = gvf_data$SE,
    V = gvf_data$V,
    T = gvf_data$T
  )
