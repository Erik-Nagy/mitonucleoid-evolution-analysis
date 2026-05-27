library(shiny)
library(r3dmol)
library(stringr)
library(dplyr)
library(readr)

# ---------------------------------------------------------
# 1. Prepare Data
# ---------------------------------------------------------
tsv_file <- "./data/mt_nucleoid_processed_20260414.tsv"

if(file.exists(tsv_file)) {
  protein_data <- read_tsv(tsv_file, show_col_types = FALSE)
} else {
  stop("TSV file not found: ", tsv_file)
}

# ---------------------------------------------------------
# 2. UI
# ---------------------------------------------------------
ui <- fluidPage(
  titlePanel("Mitochondrial P-sites Explorer"),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      selectInput("protein_select", 
                  "Select protein (UniProt ID):", 
                  choices = protein_data$`Uniprot ID`),

      tags$hr(),

      textInput("psites_input", 
                "P-sites (separate by comma):", 
                value = ""),

      actionButton("reset_btn", "Default", class = "btn-warning"),

      tags$hr(),

      radioButtons("color_select", 
                   "Style of protein:",
                   choices = c(
                     "Spectrum (N-terminus -> C-terminus)" = "spectrum",
                     "Secondary structure (Orange=Alpha Helix, Blue=Beta Sheet, Green=Loop)" = "custom_ss",
                     "Gray" = "gray"
                   ),
                   selected = "spectrum"),
    ),
    mainPanel(
      width = 9,
      r3dmolOutput("mol_viewer", height = "800px") 
    )
  )
)

# ---------------------------------------------------------
# 3. LOGIC
# ---------------------------------------------------------
server <- function(input, output, session) {

  default_psites <- reactive({
    req(input$protein_select)
    sites <- protein_data %>% 
      filter(`Uniprot ID` == input$protein_select) %>% 
      pull(`P-site positions`)
    return(as.character(sites))
  })

  observeEvent(input$protein_select, {
    updateTextInput(session, "psites_input", value = default_psites())
  })

  observeEvent(input$reset_btn, {
    updateTextInput(session, "psites_input", value = default_psites())
  })

  parsed_sites <- reactive({
    raw_text <- input$psites_input
    if(is.null(raw_text) || raw_text == "") return(numeric(0))

    clean_text = str_replace_all(raw_text, " ", "")
    split_text = str_split(clean_text, ",")[[1]]
    numeric_sites = suppressWarnings(as.numeric(split_text))
    valid_sites = numeric_sites[!is.na(numeric_sites)]
    
    return(valid_sites)
  })

  output$mol_viewer <- renderR3dmol({
    req(input$protein_select)

    pdb_path <- file.path("./data/proteins", paste0(input$protein_select, ".pdb"))

    if(!file.exists(pdb_path)) {
      showNotification(paste("PDB file not found:", pdb_path), type = "error")
      return(NULL)
    }

    pdb_content <- paste(readLines(pdb_path, warn = FALSE), collapse = "\n")

    viewer <- r3dmol(
      viewer_spec = m_viewer_spec(
        cartoonQuality = 10,
        lowerZoomLimit = 10,
        upperZoomLimit = 350
      )
    ) %>%
      m_add_model(data = pdb_content, format = "pdb") %>%
      m_zoom_to()

    if(input$color_select == "spectrum") {
      viewer <- viewer %>% 
        m_set_style(style = m_style_cartoon(color = "spectrum"))

    } else if(input$color_select == "custom_ss") {

      viewer <- viewer %>% 
        m_set_style(style = m_style_cartoon(color = "#00cc96")) %>%

        m_set_style(
          sel = m_sel(ss = "s"),
          style = m_style_cartoon(color = "#636efa", arrows = TRUE)
        ) %>%

        m_set_style(
          sel = m_sel(ss = "h"), 
          style = m_style_cartoon(color = "#ff7f0e")
        )

    } else if(input$color_select == "gray") {
      viewer <- viewer %>% 
        m_set_style(style = m_style_cartoon(color = "#E0E0E0"))
    }

    sites_to_plot <- parsed_sites()

    if(length(sites_to_plot) > 0) {
      for(site in sites_to_plot) {
        viewer <- viewer %>%
          m_add_style(
            sel = m_sel(resi = site), 
            style = c(
              m_style_stick(),
              m_style_sphere(scale = 0.3)
            )
          )
      }
    }

    viewer
  })
}

shinyApp(ui = ui, server = server)