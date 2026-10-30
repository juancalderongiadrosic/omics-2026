pacman::p_load(
  tidyverse,
  data.table,
  ggpubr
)

setwd(this.path::here())

# grab all SQL tables
# help file says use ftp + rsync
rsync_out <- system2(
  command = "rsync",
  args = c(
    "-avzP",
    "rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gtexGeneV8.txt.gz",
    "."),
    stdout = TRUE
)
print(rsync_out)

# we need metadata
rsync_metadata_out <- system2(
  command = "rsync",
  args = c(
    "-avzP",
    "rsync://hgdownload.soe.ucsc.edu/goldenPath/hgFixed/database/gtexTissueV8.txt.gz",
    "."
  ),
  stdout = TRUE
)
print(rsync_metadata_out)

# we want cerebellum stuff
#
# UCSC schema
# cat(readr::read_file("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gtexGeneV8.sql"))
# cat(readr::read_file("https://hgdownload.soe.ucsc.edu/goldenPath/hgFixed/database/gtexTissueV8.sql"))
metadata <- read_tsv(
  "gtexTissueV8.txt.gz",
  col_names = c("id", "name", "description", "organ", "color", "abbrev"),
  show_col_types = TRUE
)
data_raw <- read_tsv(
  "gtexGeneV8.txt.gz",
  col_names = c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "geneId", "geneType", "expCount",
  "expScores"),
  show_col_types = TRUE
)

data <- data_raw %>% 
  select(
    geneId,
    geneType,
    chrom,
    chromStart,
    chromEnd,
    strand,
    expScores
  ) %>% 
  filter(
    geneType == "protein_coding"
  ) %>% 
  separate_wider_delim(
    .,
    cols = expScores,
    delim = ",",
    names = as.character(c(1:54)),
    too_many = "drop"
  ) %>% 
  pivot_longer(
    cols = -c(geneId, geneType, chrom, chromStart, chromEnd, strand),
    names_to = "id"
  ) %>% 
  mutate(
    id = as.numeric(id)
  ) %>% 
  left_join(
    ., metadata, by = c("id")
  )
data

# now filter cerebellum-related
# lets plot size / loci

cb_data <- data %>% 
  filter(grepl("brainCerebellum|brainCortex|brainAmygdala|stomach", name)) %>% 
  mutate(
    locus_kb = (chromEnd - chromStart) / 1000,
    value = as.numeric(value)
  )
cb_data

# lets see if we get a bi-modal dist
p1 <- cb_data %>% 
  ggplot(aes(
    (locus_kb + 1), (value + 1),
    color = name
  )) +
  geom_point() +
  scale_x_continuous(trans = "log") +
  scale_y_continuous(trans = "log") +
  theme_pubr()
p1

p2 <- cb_data %>% 
  ggplot(aes(
    locus_kb, value, color = name,
    fill = name
  )) +
  stat_ellipse(
    type = "norm",
    level = 0.3935,
    geom = "polygon",
    alpha = 0.15,
    linewidth = 0.8
  ) +
  theme_pubr()
p2

ggsave(
  "p1.pdf",
  dpi = 300,
  plot = p1,
  device = cairo_pdf,
  width = 12,
  height = 12,
  units = "in"
)

ggsave(
  "p2.pdf",
  dpi = 300,
  plot = p2,
  device = cairo_pdf,
  width = 12,
  height = 12,
  units = "in"
)


