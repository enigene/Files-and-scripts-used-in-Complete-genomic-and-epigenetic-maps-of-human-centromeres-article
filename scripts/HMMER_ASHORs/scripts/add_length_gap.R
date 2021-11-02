#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
input_file <- args[1]
output_dir <- args[2]

suppressMessages(library("tidyverse"))

file_name <- sub("\\.[[:alnum:]]+\\.[[:alnum:]]+$", "", basename(as.character(input_file)))

cn <- c("Contig", "Start", "End", "Name", "Score", "Strand", "Start2", "End2", "Color")
ct <- cols(.default="d", "Contig" = "c", "Name" = "c", "Strand" = "c", "Color" = "c")
t1 <- read_tsv(file(input_file), col_names=cn, col_types=ct)
t1.m <- t1 %>%
  select(-Start2, -End2, -Color) %>%
  mutate(Length = End-Start) %>%
  group_by(Contig) %>%
  mutate(Gap = Start - lag(End)) %>%
  ungroup()

write_tsv(t1.m, path = paste0(output_dir, "/", file_name, ".lengap.tsv"), col_names = FALSE)
