#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
})

# ---------------------------- helpers ----------------------------

sanitize_for_dndscv <- function(dt, label = "") {
  if (nrow(dt) == 0) return(list(dt = dt, dropped_star = 0L, dropped_invalid = 0L))

  dt[, ref := toupper(trimws(as.character(ref)))]
  dt[, mut := toupper(trimws(as.character(mut)))]

  # Drop spanning deletions (ALT='*')
  dropped_star <- sum(dt$mut == "*", na.rm = TRUE)
  if (dropped_star > 0) dt <- dt[mut != "*"]

  valid_nt <- c("A", "C", "G", "T", "N", "-")

  ok_str <- function(s) {
    if (is.na(s) || s == "") return(FALSE)
    chars <- strsplit(s, "", fixed = TRUE)[[1]]
    all(chars %in% valid_nt)
  }

  valid_mask <- vapply(dt$ref, ok_str, logical(1)) & vapply(dt$mut, ok_str, logical(1))
  dropped_invalid <- sum(!valid_mask, na.rm = TRUE)

  if (dropped_star > 0 || dropped_invalid > 0) {
    message(sprintf(
      "[WARN] %s: dropped %d ALT='*' rows and %d invalid rows (non-ACGTN/-).",
      label, dropped_star, dropped_invalid
    ))
  }

  dt <- dt[valid_mask]
  list(dt = dt, dropped_star = as.integer(dropped_star), dropped_invalid = as.integer(dropped_invalid))
}

assert_schema <- function(dt) {
  req <- c("sampleID", "chr", "pos", "ref", "mut")
  missing <- setdiff(req, names(dt))
  if (length(missing) > 0) stop(sprintf("Missing columns: %s", paste(missing, collapse = ", ")))

  pos_ok <- suppressWarnings(!is.na(as.integer(dt$pos)) & as.integer(dt$pos) > 0)
  if (!all(pos_ok)) stop("Invalid POS values detected (must be positive integers).")

  invisible(TRUE)
}

infer_sample_id_from_path <- function(vcf_path) {
  bn <- basename(vcf_path)
  sub("\\.vcf(\\.gz)?$", "", bn)
}

check_bcftools <- function() {
  if (Sys.which("bcftools") == "") {
    stop("bcftools not found in PATH. Activate your conda env (bcftools required).")
  }
}

vcf_to_dt_bcftools <- function(vcf_path, sample_id = NULL) {
  check_bcftools()

  sid <- if (!is.null(sample_id)) sample_id else infer_sample_id_from_path(vcf_path)

  # Query CHROM, POS, REF, ALT from VCF
  cmd <- sprintf("bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\n' %s", shQuote(vcf_path))
  x <- system(cmd, intern = TRUE)

  if (length(x) == 0) {
    return(data.table(sampleID = character(), chr = character(), pos = integer(), ref = character(), mut = character()))
  }

  dt <- fread(text = x, sep = "\t", header = FALSE, col.names = c("chr", "pos", "ref", "mut"))
  dt[, sampleID := sid]
  setcolorder(dt, c("sampleID", "chr", "pos", "ref", "mut"))

  # Split multi-allelic ALT like "A,C" into multiple rows
  dt <- dt[, .(mut = unlist(strsplit(mut, ",", fixed = TRUE))),
           by = .(sampleID, chr, pos, ref)]

  dt
}

read_variant_table <- function(path) {
  if (is.na(path)) {
    return(data.table(sampleID = character(), chr = character(), pos = integer(), ref = character(), mut = character()))
  }
  fread(path)
}

ensure_dir <- function(p) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE)
  invisible(TRUE)
}

# ---------------------------- plotting ----------------------------

write_qc_plots <- function(all_dt, qc_dt, outdir, top_n_chr = 25) {
  plots_dir <- file.path(outdir, "plots")
  ensure_dir(plots_dir)

  qc_long <- melt(
    qc_dt,
    id.vars = "sampleID",
    measure.vars = c("n_snv", "n_indel"),
    variable.name = "type",
    value.name = "count"
  )
  qc_long[, type := fifelse(type == "n_snv", "SNV", "INDEL")]

  # 1) Total variants per sample
  p1 <- ggplot(qc_dt, aes(x = reorder(sampleID, n_variants), y = n_variants)) +
    geom_col() +
    coord_flip() +
    labs(title = "Total variants per sample", x = "Sample", y = "Variant count") +
    theme_minimal(base_size = 12)

  ggsave(filename = file.path(plots_dir, "qc_variants_bar.png"), plot = p1, width = 9, height = 4.5, dpi = 160)

  # 2) SNV vs INDEL stacked
  p2 <- ggplot(qc_long, aes(x = reorder(sampleID, count, FUN = sum), y = count, fill = type)) +
    geom_col() +
    coord_flip() +
    labs(title = "SNV vs INDEL counts per sample", x = "Sample", y = "Count", fill = "Variant type") +
    theme_minimal(base_size = 12)

  ggsave(filename = file.path(plots_dir, "qc_snv_indel_stacked.png"), plot = p2, width = 9, height = 4.5, dpi = 160)

  # 3) SNV fraction
  qc_dt[, snv_fraction := ifelse(n_variants > 0, n_snv / n_variants, NA_real_)]
  p3 <- ggplot(qc_dt, aes(x = reorder(sampleID, snv_fraction), y = snv_fraction)) +
    geom_col() +
    coord_flip() +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = "SNV fraction per sample", x = "Sample", y = "SNV / total variants") +
    theme_minimal(base_size = 12)

  ggsave(filename = file.path(plots_dir, "qc_snv_fraction.png"), plot = p3, width = 9, height = 4.5, dpi = 160)

  # 4) Variants by chromosome (top N)
  chr_counts <- all_dt[, .N, by = chr][order(-N)]
  chr_counts <- chr_counts[1:min(top_n_chr, .N)]
  p4 <- ggplot(chr_counts, aes(x = reorder(chr, N), y = N)) +
    geom_col() +
    coord_flip() +
    labs(title = sprintf("Variants by chromosome (top %d)", nrow(chr_counts)), x = "Chromosome", y = "Variant count") +
    theme_minimal(base_size = 12)

  ggsave(filename = file.path(plots_dir, "qc_variants_by_chr.png"), plot = p4, width = 9, height = 6, dpi = 160)

  message(sprintf("[INFO] Wrote QC plots to %s", plots_dir))
}

# ---------------------------- Shiny scaffold ----------------------------

write_shiny_app <- function(outdir) {
  app_dir <- file.path(outdir, "shiny_app")
  ensure_dir(app_dir)

  app_path <- file.path(app_dir, "app.R")

  # App reads outputs from the same outdir by default.
  app_txt <- '
suppressPackageStartupMessages({
  library(shiny)
  library(data.table)
  library(ggplot2)
})

ui <- fluidPage(
  titlePanel("dNdScv_prep QC Dashboard"),
  sidebarLayout(
    sidebarPanel(
      helpText("Loads qc_summary.csv and dndscv_input.csv produced by dndscv_prep.R"),
      textInput("outdir", "Output directory", value = "."),
      actionButton("reload", "Reload"),
      hr(),
      checkboxInput("show_table", "Show QC table", value = TRUE)
    ),
    mainPanel(
      plotOutput("p_total"),
      plotOutput("p_stack"),
      plotOutput("p_chr"),
      conditionalPanel(
        condition = "input.show_table == true",
        hr(),
        h4("QC summary"),
        tableOutput("qc_table")
      )
    )
  )
)

server <- function(input, output, session) {

  load_data <- reactive({
    input$reload
    isolate({
      qc_path <- file.path(input$outdir, "qc_summary.csv")
      inp_path <- file.path(input$outdir, "dndscv_input.csv")

      if (!file.exists(qc_path)) stop("qc_summary.csv not found in outdir.")
      if (!file.exists(inp_path)) stop("dndscv_input.csv not found in outdir.")

      qc <- fread(qc_path)
      inp <- fread(inp_path)
      list(qc = qc, inp = inp)
    })
  })

  output$p_total <- renderPlot({
    d <- load_data()$qc
    ggplot(d, aes(x = reorder(sampleID, n_variants), y = n_variants)) +
      geom_col() +
      coord_flip() +
      labs(title = "Total variants per sample", x = "Sample", y = "Variant count") +
      theme_minimal(base_size = 12)
  })

  output$p_stack <- renderPlot({
    d <- load_data()$qc
    long <- melt(d, id.vars = "sampleID", measure.vars = c("n_snv", "n_indel"),
                 variable.name = "type", value.name = "count")
    long[, type := ifelse(type == "n_snv", "SNV", "INDEL")]
    ggplot(long, aes(x = reorder(sampleID, count, FUN = sum), y = count, fill = type)) +
      geom_col() +
      coord_flip() +
      labs(title = "SNV vs INDEL counts per sample", x = "Sample", y = "Count", fill = "Variant type") +
      theme_minimal(base_size = 12)
  })

  output$p_chr <- renderPlot({
    inp <- load_data()$inp
    chr_counts <- inp[, .N, by = chr][order(-N)][1:min(25, .N)]
    ggplot(chr_counts, aes(x = reorder(chr, N), y = N)) +
      geom_col() +
      coord_flip() +
      labs(title = "Variants by chromosome (top 25)", x = "Chromosome", y = "Variant count") +
      theme_minimal(base_size = 12)
  })

  output$qc_table <- renderTable({
    load_data()$qc[order(-n_variants)]
  })
}

shinyApp(ui, server)
'

  writeLines(app_txt, con = app_path)
  message(sprintf("[INFO] Wrote Shiny app scaffold to %s", app_path))
  message(sprintf("[INFO] Run it with: R -q -e \"shiny::runApp('%s')\"", app_dir))
}

# ---------------------------- CLI ----------------------------

opt_list <- list(
  make_option(c("--vcf-glob"), type = "character", default = NA,
              help = "Glob for input VCFs (e.g. 'vcf/*.vcf.gz'). If provided, ignores --snv/--indel."),
  make_option(c("--snv"), type = "character", default = NA,
              help = "Path to combined_snv_variants.csv (schema: sampleID,chr,pos,ref,mut)"),
  make_option(c("--indel"), type = "character", default = NA,
              help = "Path to combined_indels_variants.csv (schema: sampleID,chr,pos,ref,mut)"),
  make_option(c("-o", "--outdir"), type = "character", default = "dndscv",
              help = "Output directory (default: %default)"),
  make_option(c("--write-qc"), action = "store_true", default = TRUE,
              help = "Write qc_summary.csv (default: TRUE)"),
  make_option(c("--write-plots"), action = "store_true", default = TRUE,
              help = "Write QC plots to outdir/plots (default: TRUE)"),
  make_option(c("--write-shiny"), action = "store_true", default = FALSE,
              help = "Create a Shiny app scaffold in outdir/shiny_app (default: FALSE)"),
  make_option(c("--top-chr"), type = "integer", default = 25,
              help = "Top N chromosomes to plot (default: %default)")
)

opt <- parse_args(OptionParser(option_list = opt_list))

outdir <- opt$outdir
ensure_dir(outdir)

message("[INFO] Starting dNdScv prep (R, bcftools-backed + visualisation)")

# ---------------------------- load inputs ----------------------------

all_dt <- NULL

if (!is.na(opt$`vcf-glob`)) {
  vcf_files <- Sys.glob(opt$`vcf-glob`)
  if (length(vcf_files) == 0) stop(sprintf("No VCFs found for glob: %s", opt$`vcf-glob`))

  message(sprintf("[INFO] Found %d VCF files", length(vcf_files)))

  dt_list <- lapply(vcf_files, function(v) {
    message(sprintf("[INFO] Reading %s", v))
    vcf_to_dt_bcftools(v)
  })

  all_dt <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)

} else {
  if (is.na(opt$snv) && is.na(opt$indel)) {
    stop("Provide either --vcf-glob OR at least one of --snv/--indel.")
  }

  message("[INFO] Loading variant tables...")

  snv <- read_variant_table(opt$snv)
  indel <- read_variant_table(opt$indel)

  if (!is.na(opt$snv)) assert_schema(snv)
  if (!is.na(opt$indel)) assert_schema(indel)

  snv_s <- sanitize_for_dndscv(copy(snv), label = "SNVs")$dt
  indel_s <- sanitize_for_dndscv(copy(indel), label = "INDELs")$dt

  all_dt <- rbindlist(list(snv_s, indel_s), use.names = TRUE, fill = TRUE)
}

# ---------------------------- sanitize + write ----------------------------

assert_schema(all_dt)

san <- sanitize_for_dndscv(copy(all_dt), label = "ALL")
all_dt <- san$dt

assert_schema(all_dt)

# Write dNdScv input
inp_path <- file.path(outdir, "dndscv_input.csv")
fwrite(all_dt, inp_path)
message(sprintf("[INFO] Wrote dndscv_input.csv with %d rows", nrow(all_dt)))

# QC summary
qc <- all_dt[, .(
  n_variants = .N,
  n_snv = sum(nchar(ref) == 1 & nchar(mut) == 1),
  n_indel = sum(!(nchar(ref) == 1 & nchar(mut) == 1))
), by = sampleID][order(-n_variants)]

if (isTRUE(opt$`write-qc`)) {
  qc_path <- file.path(outdir, "qc_summary.csv")
  fwrite(qc, qc_path)
  message("[INFO] Wrote qc_summary.csv")
}

# Plots
if (isTRUE(opt$`write-plots`)) {
  write_qc_plots(all_dt, qc, outdir, top_n_chr = opt$`top-chr`)
}

# Shiny scaffold
if (isTRUE(opt$`write-shiny`)) {
  write_shiny_app(outdir)
}

message("[INFO] Done")
