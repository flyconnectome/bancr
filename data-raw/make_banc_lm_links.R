#!/usr/bin/env Rscript
#
# make_banc_lm_links.R
#
# Build `banc_lm_links` --- a four-column flat record
# (source, gene, sample, ngl_link) of pre-minted Spelunker viewer URLs
# for every light-microscopy layer hosted under
#   gs://lee-lab_brain-and-nerve-cord-fly-connectome/light_level/
# and indexed by the master registry at
#   gs://lee-lab_brain-and-nerve-cord-fly-connectome/light_level/
#                                       registry.json
# Currently covers Kondo et al. 2020 (receptor / peptide / orphan
# GPCR GFP knock-ins); will extend to Deng et al. 2019 once those
# are bridged + uploaded.
#
# Each `ngl_link` is a fresh state POSTed via `bancr::banc_lm_scene()`
# to the BANC state server. Once minted, links are persistent --- the
# state server resolves the same id forever --- so the on-disk CSV is
# a stable shareable index.
#
# Re-run this script to:
#   * pick up new layers added to the registry
#   * re-mint links if the default render settings need updating
#
# Outputs:
#   data-raw/banc_lm_links.csv  --- flat human-readable record
#   data/banc_lm_links.rda      --- packaged tibble (via usethis::use_data)
#
# Usage:
#   cd <bancr root>
#   Rscript data-raw/make_banc_lm_links.R

suppressMessages({
  library(httr)
  library(jsonlite)
  library(tibble)
  library(dplyr)
  library(devtools)
  library(usethis)
})

# Load bancr from the working tree so banc_lm_scene() picks up local
# changes (and so banc_lm_volumes() resolves).
devtools::load_all(".", quiet = TRUE)

# Per-gene render defaults. Ranges are
# `shaderControls.normalized.range` for the LM layer in the rendered
# Spelunker scene; tighter ranges suit dim stains. Falls back to the
# `default` row for genes not listed.
RENDER_DEFAULTS <- tibble::tribble(
  ~gene,         ~range_lo, ~range_hi,
  "default",        5L,       80L,
  "GluRIIA",        2L,       25L,
  "GluRIB",         3L,       30L,
  "Nmdar1",         5L,      150L,
  "Nmdar2",         5L,      150L,
  "mGluR",          5L,      150L,
  "GluClalpha",     5L,      120L
)

# Pull the master registry from GCS.
registry_url <- "https://storage.googleapis.com/lee-lab_brain-and-nerve-cord-fly-connectome/light_level/registry.json"
res <- httr::GET(registry_url)
if (httr::status_code(res) != 200L)
  stop("Could not fetch master registry: ", registry_url,
       " (HTTP ", httr::status_code(res), ")")
reg <- jsonlite::fromJSON(httr::content(res, as = "text", encoding = "UTF-8"),
                          simplifyVector = TRUE)
vols <- reg[["volumes"]]
if (!nrow(vols)) stop("Registry has 0 volumes; nothing to mint.")
cat("registry has", nrow(vols), "volumes\n")

# Map registry dataset id to a human-readable citation source.
DATASET_SOURCE <- list(
  kondo_et_al_2020 = "Kondo et al. 2020"
)
src_for <- function(dataset) {
  s <- DATASET_SOURCE[[dataset]]
  if (is.null(s)) dataset else s
}

range_for <- function(gene) {
  hit <- RENDER_DEFAULTS[RENDER_DEFAULTS$gene == gene, , drop = FALSE]
  if (!nrow(hit))
    hit <- RENDER_DEFAULTS[RENDER_DEFAULTS$gene == "default", , drop = FALSE]
  c(hit$range_lo[[1]], hit$range_hi[[1]])
}

# Mint URLs.
out <- vector("list", nrow(vols))
for (i in seq_len(nrow(vols))) {
  v   <- vols[i, ]
  rng <- range_for(v$gene)
  cat(sprintf("[%2d/%d] %s_%s  range=[%g,%g] ... ",
              i, nrow(vols), v$gene, v$sample, rng[1], rng[2]))
  url <- bancr::banc_lm_scene(
    lm_url     = v$gs_url,
    layer_name = sprintf("%s - %s %s", src_for(v$dataset), v$gene, v$sample),
    range      = rng,
    opacity    = 0.55,
    blend      = "additive",
    volume_rendering = "on",
    shorten    = TRUE
  )
  cat("ok\n")
  out[[i]] <- tibble::tibble(
    source   = src_for(v$dataset),
    gene     = v$gene,
    sample   = v$sample,
    ngl_link = url
  )
}

banc_lm_links <- do.call(rbind, out)
banc_lm_links <- dplyr::arrange(banc_lm_links, source, gene, sample)

# Write flat CSV (data-raw, not packaged) + packaged .rda (data/).
write.csv(banc_lm_links,
          file = "data-raw/banc_lm_links.csv",
          row.names = FALSE)
usethis::use_data(banc_lm_links, overwrite = TRUE)

cat("\nwrote", nrow(banc_lm_links), "rows to data-raw/banc_lm_links.csv",
    "and data/banc_lm_links.rda\n")
