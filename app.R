# File: F:/Structure Harvester/app.R
# Version: STRUCTURE DETERMINANT v1.0
# Optimized for Shinylive + Admixture Analysis

# --- Load required packages ---
library(shiny)
library(dplyr)
library(stringr)
library(ggplot2)
library(readr)
library(purrr)
library(DT)
library(tidyr)

# ============================================================
# --- CORE PARSING FUNCTIONS ---
# ============================================================

# Extract K value from file content (search for MAXPOPS)
extract_K <- function(txt) {
  k_pattern <- "MAXPOPS\\s*=\\s*([0-9]+)"
  matched_line <- txt[str_detect(txt, regex(k_pattern, ignore_case = TRUE))]
  if (length(matched_line) > 0) {
    return(as.numeric(str_match(matched_line[1], regex(k_pattern, ignore_case = TRUE))[, 2]))
  }
  return(NA_real_)
}

# Extract estimated Ln probability of data
extract_LnPD <- function(txt) {
  target_lines <- txt[str_detect(txt, regex("Estimated\\s+Ln\\s+Prob\\s+of\\s+Data", ignore_case = TRUE))]
  if (length(target_lines) == 0) return(NA_real_)
  num <- str_extract(target_lines[1], "[-+]?[0-9]*\\.?[0-9]+")
  return(as.numeric(num))
}

# Extract replicate number from file name
extract_replicate <- function(file_name) {
  res <- str_match(basename(file_name), "run_([0-9]+)")[, 2]
  if (is.na(res)) res <- str_extract(file_name, "[0-9]+")
  return(as.numeric(res))
}

# ============================================================
# --- Q-MATRIX EXTRACTION (with genotype names) ---
# ============================================================
extract_q_matrix_named <- function(file_path) {
  txt <- readLines(file_path, warn = FALSE)
  start_line <- which(str_detect(txt, "Inferred ancestry of individuals:")) + 1
  if (length(start_line) == 0 || start_line > length(txt)) return(NULL)
  
  all_empty <- which(txt == "")
  end_line <- all_empty[all_empty > start_line][1] - 1
  if (is.na(end_line) || end_line < start_line) end_line <- length(txt)
  
  raw_data <- txt[start_line:end_line]
  raw_data <- raw_data[str_detect(raw_data, ":")]
  if (length(raw_data) == 0) return(NULL)
  
  genotype_names <- character(length(raw_data))
  q_values_list <- vector("list", length(raw_data))
  
  for (i in seq_along(raw_data)) {
    line <- raw_data[i]
    # STRUCTURE format: "1 G39 (124) : 0.85 0.10 0.05"
    name_match <- str_match(line, "^\\s*\\d+\\s+([^\\s(]+)")
    genotype_names[i] <- if (!is.na(name_match[1, 2])) name_match[1, 2] else paste0("Ind", i)
    
    parts <- str_split(str_trim(line), "\\s+")[[1]]
    colon_pos <- which(parts == ":")
    if (length(colon_pos) > 0) {
      q_values_list[[i]] <- as.numeric(parts[(colon_pos + 1):length(parts)])
    } else {
      q_values_list[[i]] <- rep(NA_real_, 2)
    }
  }
  
  # Remove rows with all-NA values
  valid_idx <- !sapply(q_values_list, function(x) all(is.na(x)))
  q_values_list <- q_values_list[valid_idx]
  genotype_names <- genotype_names[valid_idx]
  
  if (length(q_values_list) == 0) return(NULL)
  
  mat <- do.call(rbind, q_values_list)
  df <- as.data.frame(mat)
  colnames(df) <- paste0("Cluster", seq_len(ncol(df)))
  df$Genotype <- genotype_names
  
  return(df)
}

# ============================================================
# --- ADMIXTURE CLASSIFICATION ---
# ============================================================

# Classify genotypes based on a membership threshold
classify_genotypes <- function(q_df, threshold = 0.60) {
  cluster_cols <- names(q_df)[str_detect(names(q_df), "^Cluster")]
  
  classification <- apply(q_df[, cluster_cols, drop = FALSE], 1, function(q_vals) {
    q_vals <- as.numeric(q_vals)
    max_q <- max(q_vals, na.rm = TRUE)
    if (max_q >= threshold) {
      return(paste0("Cluster ", which.max(q_vals)))
    } else {
      return("Admixed")
    }
  })
  
  q_df$Max_Q <- apply(q_df[, cluster_cols, drop = FALSE], 1, max, na.rm = TRUE)
  q_df$Classification <- classification
  
  return(q_df)
}

# Summarize admixture statistics
summarize_admixture <- function(classified_df, threshold) {
  summary_tbl <- classified_df %>%
    count(Classification) %>%
    rename(Group = Classification, N = n) %>%
    mutate(Percent = round(100 * N / sum(N), 2))
  
  admixed_genotypes <- classified_df %>%
    filter(Classification == "Admixed") %>%
    pull(Genotype)
  
  list(
    summary = summary_tbl,
    admixed_list = admixed_genotypes,
    n_admixed = length(admixed_genotypes),
    threshold = threshold
  )
}

# Generate publication-ready text describing admixture results
generate_admixture_text <- function(classified_df, K, threshold) {
  summ <- summarize_admixture(classified_df, threshold)
  
  n_total <- nrow(classified_df)
  n_admixed <- summ$n_admixed
  admixed_list <- summ$admixed_list
  
  if (n_admixed == 0) {
    admixed_text <- "No genotypes were classified as admixed."
  } else if (n_admixed <= 10) {
    admixed_text <- paste0(
      "Genotypes ", paste(admixed_list, collapse = ", "),
      " were classified as admixed (all membership coefficients < ",
      sprintf("%.2f", threshold), ")."
    )
  } else {
    admixed_text <- paste0(
      "A total of ", n_admixed, " genotypes were classified as admixed ",
      "(all membership coefficients < ", sprintf("%.2f", threshold), "), ",
      "including ", paste(head(admixed_list, 10), collapse = ", "),
      ", and ", n_admixed - 10, " others."
    )
  }
  
  main_text <- paste0(
    "At K = ", K, ", a genotype was assigned to a genetic cluster when its ",
    "membership coefficient (Q) was >= ", sprintf("%.2f", threshold), ". ",
    "Based on this criterion, ", n_total - n_admixed, " of ", n_total,
    " genotypes were assigned to distinct clusters, while ", n_admixed,
    " were classified as admixed. ", admixed_text
  )
  
  return(main_text)
}

# ============================================================
# --- SHARED COLOR PALETTE ---
# ============================================================
get_fill_colors <- function(cluster_names) {
  base_colors <- c(
    "#1F77B4", "#D62728", "#FFD700", "#2CA02C", "#9467BD",
    "#FF7F0E", "#17BECF", "#8C564B", "#E377C2", "#7F7F7F",
    "#BCBD22", "#AEC7E8"
  )
  n <- length(cluster_names)
  if (n > length(base_colors)) {
    fill_colors <- colorRampPalette(base_colors)(n)
  } else {
    fill_colors <- base_colors[seq_len(n)]
  }
  names(fill_colors) <- cluster_names
  return(fill_colors)
}

# ============================================================
# --- SHARED STRUCTURE PLOT ---
# ============================================================
make_structure_plot <- function(q_df, selected_k) {
  cluster_cols <- names(q_df)[str_detect(names(q_df), "^Cluster")]
  
  plot_df <- q_df %>%
    pivot_longer(cols = all_of(cluster_cols),
                 names_to = "Cluster", values_to = "Proportion")
  
  fill_colors <- get_fill_colors(cluster_cols)
  
  ggplot(plot_df, aes(x = factor(Genotype, levels = q_df$Genotype),
                      y = Proportion, fill = Cluster)) +
    geom_bar(stat = "identity", width = 1, color = NA) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_fill_manual(values = fill_colors) +
    labs(
      title = paste("Structure Plot (K =", selected_k, ")"),
      x = "Genotypes",
      y = "Membership Probability"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

# ============================================================
# --- ADMIXTURE PLOT (with admixed marker) ---
# ============================================================
make_admixture_plot <- function(classified_df, K, threshold) {
  cluster_cols <- names(classified_df)[str_detect(names(classified_df), "^Cluster")]
  
  plot_df <- classified_df %>%
    pivot_longer(cols = all_of(cluster_cols),
                 names_to = "Cluster", values_to = "Proportion") %>%
    mutate(Cluster = factor(Cluster, levels = cluster_cols))
  
  fill_colors <- get_fill_colors(cluster_cols)
  
  # Order: first assigned clusters, then admixed
  ordered_genotypes <- classified_df %>%
    arrange(Classification, desc(Max_Q)) %>%
    pull(Genotype)
  
  plot_df$Genotype <- factor(plot_df$Genotype, levels = ordered_genotypes)
  
  admixed_genotypes <- classified_df %>%
    filter(Classification == "Admixed") %>%
    pull(Genotype)
  
  p <- ggplot(plot_df, aes(x = Genotype, y = Proportion, fill = Cluster)) +
    geom_bar(stat = "identity", width = 1, color = NA) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_fill_manual(values = fill_colors) +
    labs(
      title = paste0("Admixture Plot (K = ", K,
                     ", threshold = ", sprintf("%.2f", threshold), ")"),
      x = "Genotypes (ordered by classification)",
      y = "Membership Coefficient (Q)"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  # Dashed separator between assigned and admixed genotypes
  if (length(admixed_genotypes) > 0 && length(admixed_genotypes) < nrow(classified_df)) {
    n_assigned <- nrow(classified_df) - length(admixed_genotypes)
    p <- p + geom_vline(xintercept = n_assigned + 0.5,
                        linetype = "dashed", color = "black", linewidth = 0.7)
  }
  
  return(p)
}

# ============================================================
# --- USER INTERFACE (UI) ---
# ============================================================

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .btn-primary { background-color: #2c3e50; color: white; }
      .well { background-color: #f8f9fa; }
      .admixture-text {
        background-color: #f4f4f4;
        padding: 10px;
        border-radius: 5px;
        font-size: 13px;
        border-left: 4px solid #2c3e50;
      }
    "))
  ),
  titlePanel("STRUCTURE DETERMINANT v1.0 - Optimal K & Admixture Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("zip_file", "1. Upload Results.zip:", accept = ".zip"),
      numericInput("min_K", "Min K:", value = 1, min = 1),
      numericInput("max_K", "Max K:", value = 10, min = 2),
      actionButton("run_analysis", "Run Analysis", class = "btn-primary", style="width: 100%"),
      hr(),
      h4("Plot Settings"),
      uiOutput("select_k_ui"),
      downloadButton("download_structure_plot", "Download Structure Plot"),
      downloadButton("download_summary", "Download Evanno Table"),
      
      hr(),
      h4("Admixture Settings"),
      sliderInput(
        "admixture_threshold",
        "Membership threshold (Q):",
        min = 0.50, max = 1.00, value = 0.60, step = 0.05
      ),
      helpText("Genotypes with max(Q) >= threshold are assigned to that cluster; otherwise classified as Admixed."),
      downloadButton("download_admixture_table", "Download Admixture Table"),
      downloadButton("download_admixture_plot", "Download Admixture Plot")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Evanno Results",
                 br(),
                 h4(textOutput("best_k_text")),
                 DTOutput("evanno_table")),
        
        tabPanel("Structure Plot",
                 br(),
                 plotOutput("structure_plot", height = "400px"),
                 br(),
                 helpText("The plot above shows the genetic membership of individuals for the selected K.")),
        
        tabPanel("Admixture Analysis",
                 br(),
                 h4("Admixture Assessment (Q-matrix based)"),
                 helpText("This module classifies genotypes based on their membership coefficients
                          (Q) from the STRUCTURE Q-matrix. A genotype is assigned to a cluster
                          if its maximum Q >= the user-defined threshold; otherwise, it is
                          classified as Admixed."),
                 
                 fluidRow(
                   column(5,
                          h5("Cluster Composition"),
                          DTOutput("admixture_summary")
                   ),
                   column(7,
                          h5("Publication-ready text:"),
                          div(class = "admixture-text",
                              verbatimTextOutput("admixture_text"))
                   )
                 ),
                 
                 hr(),
                 h4("Admixture Bar Plot"),
                 plotOutput("admixture_plot", height = "450px"),
                 
                 hr(),
                 h4("Full Admixture Table"),
                 DTOutput("admixture_table")),
        
        tabPanel("Plots (Delta K & LnP)",
                 plotOutput("delta_plot"),
                 plotOutput("lnp_plot")),
        
        tabPanel("Raw Data", DTOutput("raw_table"))
      )
    )
  )
)

# ============================================================
# --- SERVER ---
# ============================================================

server <- function(input, output, session) {
  
  # --- Main analysis reactive ---
  results <- eventReactive(input$run_analysis, {
    req(input$zip_file)
    
    td <- tempfile()
    dir.create(td)
    
    tryCatch({
      unzip(input$zip_file$datapath, exdir = td)
    }, error = function(e) {
      showNotification("Error extracting ZIP file. Please upload a valid ZIP archive.", type = "error")
      return(NULL)
    })
    
    files <- list.files(td, full.names = TRUE, recursive = TRUE)
    files <- files[!dir.exists(files)]
    
    raw_list <- map_df(files, function(f) {
      txt <- readLines(f, warn = FALSE)
      k_val <- extract_K(txt)
      lnp <- extract_LnPD(txt)
      if(!is.na(k_val) && !is.na(lnp)) {
        return(data.frame(File = basename(f), Path = f, K = k_val,
                          Replicate = extract_replicate(f), LnP_D = lnp))
      }
      return(NULL)
    })
    
    req(nrow(raw_list) > 0)
    
    raw_list <- raw_list %>% filter(K >= input$min_K, K <= input$max_K)
    
    summary_df <- raw_list %>%
      group_by(K) %>%
      summarise(
        Reps = n(),
        Mean_LnP = mean(LnP_D),
        SD_LnP = sd(LnP_D),
        .groups = "drop"
      ) %>%
      arrange(K)
    
    if(nrow(summary_df) >= 3) {
      summary_df$L_prime <- c(NA, diff(summary_df$Mean_LnP))
      L_double_prime <- abs(diff(summary_df$Mean_LnP, differences = 2))
      summary_df$L_double_prime <- c(NA, L_double_prime, NA)
      summary_df$Delta_K <- summary_df$L_double_prime / summary_df$SD_LnP
    }
    
    list(raw = raw_list, summary = summary_df)
  })
  
  # --- Evanno outputs ---
  output$best_k_text <- renderText({
    req(results())
    res <- results()$summary
    best <- res %>% filter(!is.na(Delta_K)) %>% slice_max(Delta_K, n = 1)
    if(nrow(best) > 0) paste("Best K according to Evanno method is:", best$K)
    else "Not enough data to calculate Delta K."
  })
  
  output$evanno_table <- renderDT({
    req(results())
    datatable(results()$summary, options = list(pageLength = 10)) %>%
      formatRound(columns=c('Mean_LnP', 'SD_LnP', 'L_prime', 'L_double_prime', 'Delta_K'), digits=3)
  })
  
  output$raw_table <- renderDT({
    req(results())
    datatable(results()$raw %>% select(-Path))
  })
  
  output$select_k_ui <- renderUI({
    req(results())
    ks <- sort(unique(results()$raw$K))
    selectInput("selected_k_plot", "Select K for Structure Plot:", choices = ks, selected = ks[1])
  })
  
  # --- Structure plot ---
  output$structure_plot <- renderPlot({
    req(input$selected_k_plot, results())
    
    target_file <- results()$raw %>%
      filter(K == as.numeric(input$selected_k_plot)) %>%
      slice(1) %>%
      pull(Path)
    
    q_df <- extract_q_matrix_named(target_file)
    req(q_df)
    
    make_structure_plot(q_df, input$selected_k_plot)
  })
  
  # --- Delta K & LnP plots ---
  output$delta_plot <- renderPlot({
    req(results())
    ggplot(results()$summary %>% filter(!is.na(Delta_K)), aes(x=K, y=Delta_K)) +
      geom_line() + geom_point(size=3, color="red") + theme_bw() + labs(title="Delta K Plot")
  })
  
  output$lnp_plot <- renderPlot({
    req(results())
    ggplot(results()$summary, aes(x=K, y=Mean_LnP)) +
      geom_line() + geom_point() +
      geom_errorbar(aes(ymin=Mean_LnP-SD_LnP, ymax=Mean_LnP+SD_LnP), width=0.1) +
      theme_bw() + labs(title="Mean LnP(D) Plot")
  })
  
  # --- Classified Q-matrix reactive (shared by admixture outputs) ---
  classified_q <- reactive({
    req(results(), input$selected_k_plot, input$admixture_threshold)
    
    target_file <- results()$raw %>%
      filter(K == as.numeric(input$selected_k_plot)) %>%
      slice(1) %>%
      pull(Path)
    
    q_df <- extract_q_matrix_named(target_file)
    req(q_df)
    
    classify_genotypes(q_df, input$admixture_threshold)
  })
  
  # --- Admixture summary table ---
  output$admixture_summary <- renderDT({
    req(classified_q())
    summ <- summarize_admixture(classified_q(), input$admixture_threshold)
    datatable(summ$summary, options = list(dom = 't'), rownames = FALSE)
  })
  
  # --- Admixture publication-ready text ---
  output$admixture_text <- renderText({
    req(classified_q())
    generate_admixture_text(classified_q(), input$selected_k_plot, input$admixture_threshold)
  })
  
  # --- Admixture bar plot ---
  output$admixture_plot <- renderPlot({
    req(classified_q())
    make_admixture_plot(classified_q(), input$selected_k_plot, input$admixture_threshold)
  })
  
  # --- Full admixture table ---
  output$admixture_table <- renderDT({
    req(classified_q())
    
    df_display <- classified_q()
    cluster_cols <- names(df_display)[str_detect(names(df_display), "^Cluster")]
    df_display[cluster_cols] <- lapply(df_display[cluster_cols], function(x) round(as.numeric(x), 3))
    df_display$Max_Q <- round(as.numeric(df_display$Max_Q), 3)
    
    datatable(
      df_display,
      options = list(pageLength = 15, scrollX = TRUE),
      caption = paste0("Admixture Table (K = ", input$selected_k_plot,
                       ", threshold = ", input$admixture_threshold, ")")
    )
  })
  # --- Download handlers ---
  output$download_summary <- downloadHandler(
    filename = function() { "Evanno_Summary.csv" },
    content = function(file) { write.csv(results()$summary, file, row.names = FALSE) }
  )
  
  output$download_structure_plot <- downloadHandler(
    filename = function() {
      paste0("Structure_Plot_K_", input$selected_k_plot, "_600dpi.png")
    },
    content = function(file) {
      req(input$selected_k_plot, results())
      
      target_file <- results()$raw %>%
        filter(K == as.numeric(input$selected_k_plot)) %>%
        slice(1) %>%
        pull(Path)
      
      q_df <- extract_q_matrix_named(target_file)
      req(q_df)
      
      plt <- make_structure_plot(q_df, input$selected_k_plot)
      
      ggsave(file, plot = plt, width = 12, height = 4,
             units = "in", dpi = 600, bg = "white")
    }
  )
  
  output$download_admixture_table <- downloadHandler(
    filename = function() {
      paste0("Admixture_Table_K", input$selected_k_plot,
             "_threshold", input$admixture_threshold, ".csv")
    },
    content = function(file) {
      req(classified_q())
      write.csv(classified_q(), file, row.names = FALSE)
    }
  )
  
  output$download_admixture_plot <- downloadHandler(
    filename = function() {
      paste0("Admixture_Plot_K", input$selected_k_plot,
             "_threshold", input$admixture_threshold, "_600dpi.png")
    },
    content = function(file) {
      req(classified_q())
      plt <- make_admixture_plot(classified_q(), input$selected_k_plot, input$admixture_threshold)
      ggsave(file, plot = plt, width = 14, height = 5,
             units = "in", dpi = 600, bg = "white")
    }
  )
}

# ============================================================
# --- RUN APPLICATION ---
# ============================================================
shinyApp(ui = ui, server = server)