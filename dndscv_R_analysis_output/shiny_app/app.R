
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

