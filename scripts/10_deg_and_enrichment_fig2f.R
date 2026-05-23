suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(httr)
  library(jsonlite)
})

project_root <- normalizePath(".", winslash = "/")
processed_dir <- file.path(project_root, "data/processed")
figure_dir <- file.path(project_root, "results/figures")
table_dir <- file.path(project_root, "results/tables")
tmp_dir <- file.path(project_root, "tmp")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMPDIR = tmp_dir, TMP = tmp_dir, TEMP = tmp_dir)

set.seed(20260519)

clean_npc <- readRDS(file.path(processed_dir, "GSE165722_clean_npc_subset_final_named.rds"))
DefaultAssay(clean_npc) <- "RNA"

if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
  clean_npc <- SeuratObject::JoinLayers(clean_npc)
}

subtype_levels <- c(
  "ECM/Adh-NPCs",
  "Hom-NPCs",
  "Eff-NPCs",
  "Ht-NPCs",
  "Fibro-NPCs",
  "Fibro-reg NPCs",
  "Hom/Reg-like NPCs"
)

clean_npc$article_group <- factor(clean_npc$article_group, levels = c("mild", "severe"))
clean_npc$npc_subtype_final <- factor(clean_npc$npc_subtype_final, levels = subtype_levels)

message("Running differential expression: severe vs mild inside each NPC subtype...")

deg_list <- list()
for (subtype in subtype_levels) {
  cells <- colnames(clean_npc)[clean_npc$npc_subtype_final == subtype]
  sub_obj <- subset(clean_npc, cells = cells)
  group_counts <- table(sub_obj$article_group)
  message(subtype, ": ", paste(names(group_counts), group_counts, sep = "=", collapse = ", "))

  if (!all(c("mild", "severe") %in% names(group_counts)) ||
      group_counts[["mild"]] < 20 ||
      group_counts[["severe"]] < 20) {
    warning("Skipping ", subtype, " because one group has fewer than 20 cells.")
    next
  }

  Idents(sub_obj) <- "article_group"
  deg <- FindMarkers(
    sub_obj,
    ident.1 = "severe",
    ident.2 = "mild",
    assay = "RNA",
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0.10,
    only.pos = FALSE
  ) %>%
    rownames_to_column("gene") %>%
    mutate(
      npc_subtype_final = subtype,
      comparison = "severe_vs_mild",
      direction = case_when(
        avg_log2FC > 0 ~ "up_in_severe",
        avg_log2FC < 0 ~ "down_in_severe",
        TRUE ~ "no_change"
      )
    ) %>%
    relocate(npc_subtype_final, comparison, gene, direction)

  deg_list[[subtype]] <- deg
}

deg_all <- bind_rows(deg_list)
readr::write_csv(deg_all, file.path(table_dir, "clean_npc_deg_severe_vs_mild_by_subtype_all.csv"))

deg_sig <- deg_all %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) >= 0.25)
readr::write_csv(deg_sig, file.path(table_dir, "clean_npc_deg_severe_vs_mild_by_subtype_sig.csv"))

deg_summary <- deg_all %>%
  mutate(is_sig = p_val_adj < 0.05 & abs(avg_log2FC) >= 0.25) %>%
  group_by(npc_subtype_final, direction) %>%
  summarise(
    n_tested = n(),
    n_sig = sum(is_sig, na.rm = TRUE),
    .groups = "drop"
  )
readr::write_csv(deg_summary, file.path(table_dir, "clean_npc_deg_summary_by_subtype.csv"))

top_deg <- deg_sig %>%
  group_by(npc_subtype_final, direction) %>%
  arrange(p_val_adj, desc(abs(avg_log2FC)), .by_group = TRUE) %>%
  slice_head(n = 30) %>%
  ungroup()
readr::write_csv(top_deg, file.path(table_dir, "clean_npc_deg_top30_sig_by_subtype_direction.csv"))

plot_volcano <- function(deg_df) {
  deg_df %>%
    mutate(
      sig_class = case_when(
        p_val_adj < 0.05 & avg_log2FC >= 0.25 ~ "Up in severe",
        p_val_adj < 0.05 & avg_log2FC <= -0.25 ~ "Down in severe",
        TRUE ~ "Not significant"
      ),
      minus_log10_padj = -log10(pmax(p_val_adj, 1e-300))
    ) %>%
    ggplot(aes(x = avg_log2FC, y = minus_log10_padj, color = sig_class)) +
    geom_point(size = 0.35, alpha = 0.65) +
    facet_wrap(~ npc_subtype_final, scales = "free_y", ncol = 4) +
    scale_color_manual(
      values = c(
        "Up in severe" = "#E15759",
        "Down in severe" = "#4E79A7",
        "Not significant" = "#B8B8B8"
      )
    ) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", linewidth = 0.25, color = "#777777") +
    labs(
      title = "Differential genes in clean NPC subtypes: severe vs mild",
      x = "avg_log2FC (severe / mild)",
      y = "-log10 adjusted p-value",
      color = NULL
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.background = element_rect(fill = "#F2F2F2", color = NA),
      strip.text = element_text(face = "bold"),
      legend.position = "top"
    )
}

p_volcano <- plot_volcano(deg_all)
ggsave(file.path(figure_dir, "28_clean_npc_deg_volcano_by_subtype.png"), p_volcano, width = 13, height = 8, dpi = 300)

query_enrichr_once <- function(gene_vec, query_name) {
  gene_vec <- unique(gene_vec[!is.na(gene_vec) & gene_vec != ""])
  if (length(gene_vec) < 5) {
    warning("Skipping enrichment for ", query_name, ": fewer than 5 genes.")
    return(tibble())
  }

  add_response <- httr::POST(
    "https://maayanlab.cloud/Enrichr/addList",
    body = list(
      list = paste(gene_vec, collapse = "\n"),
      description = query_name
    ),
    encode = "multipart",
    timeout(120)
  )

  if (httr::http_error(add_response)) {
    warning("Enrichr addList failed for ", query_name, ": HTTP ", httr::status_code(add_response))
    return(tibble())
  }

  added <- jsonlite::fromJSON(httr::content(add_response, as = "text", encoding = "UTF-8"))
  user_list_id <- as.character(added[["userListId"]])
  libraries <- c(
    "GO_Biological_Process_2025",
    "KEGG_2021_Human",
    "Reactome_2022"
  )

  map_dfr(libraries, function(library_name) {
    enrich_url <- paste0(
      "https://maayanlab.cloud/Enrichr/enrich?userListId=",
      user_list_id,
      "&backgroundType=",
      URLencode(library_name, reserved = TRUE)
    )
    enrich_response <- httr::GET(enrich_url, timeout(120))
    if (httr::http_error(enrich_response)) {
      warning("Enrichr enrich failed for ", query_name, " / ", library_name, ": HTTP ", httr::status_code(enrich_response))
      return(tibble())
    }

    parsed <- jsonlite::fromJSON(httr::content(enrich_response, as = "text", encoding = "UTF-8"))
    if (is.null(parsed[[library_name]]) || length(parsed[[library_name]]) == 0) {
      return(tibble())
    }

    parsed_table <- parsed[[library_name]]
    if (length(parsed_table) == 0) {
      return(tibble())
    }

    map_dfr(parsed_table, function(row) {
      tibble(
        rank = as.integer(row[[1]]),
        term_name = as.character(row[[2]]),
        p_value = as.numeric(row[[3]]),
        z_score = as.numeric(row[[4]]),
        combined_score = as.numeric(row[[5]]),
        overlap_genes = list(as.character(row[[6]])),
        adjusted_p_value = as.numeric(row[[7]]),
        old_p_value = as.numeric(row[[8]]),
        old_adjusted_p_value = as.numeric(row[[9]])
      )
    }) %>%
      mutate(
        query_name = query_name,
        library = library_name
      )
  })
}

query_enrichr <- function(gene_vec, query_name, max_tries = 3) {
  for (try_id in seq_len(max_tries)) {
    result <- tryCatch(
      query_enrichr_once(gene_vec, query_name),
      error = function(e) {
        warning(
          "Enrichr try ", try_id, " failed for ", query_name, ": ",
          conditionMessage(e)
        )
        tibble()
      }
    )
    if (nrow(result) > 0 || try_id == max_tries) {
      return(result)
    }
    Sys.sleep(3 * try_id)
  }
  tibble()
}

message("Running GO/KEGG/Reactome enrichment with Enrichr...")

up_gene_sets <- deg_sig %>%
  filter(direction == "up_in_severe") %>%
  group_by(npc_subtype_final) %>%
  arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) %>%
  summarise(genes = list(head(unique(gene), 300)), .groups = "drop")

enrich_list <- list()
for (i in seq_len(nrow(up_gene_sets))) {
  subtype <- up_gene_sets$npc_subtype_final[[i]]
  query_name <- paste0(gsub("[^A-Za-z0-9]+", "_", subtype), "_up_in_severe")
  message("Enrichment: ", subtype, " (", length(up_gene_sets$genes[[i]]), " genes)")
  enrich_list[[subtype]] <- query_enrichr(up_gene_sets$genes[[i]], query_name) %>%
    mutate(npc_subtype_final = subtype, direction = "up_in_severe")
  Sys.sleep(0.5)
}

enrich_all <- bind_rows(enrich_list)

if (nrow(enrich_all) > 0) {
  enrich_clean <- enrich_all %>%
    mutate(
      source = case_when(
        library == "GO_Biological_Process_2025" ~ "GO:BP",
        library == "KEGG_2021_Human" ~ "KEGG",
        library == "Reactome_2022" ~ "REAC",
        TRUE ~ library
      ),
      intersection_size = lengths(overlap_genes),
      intersection = map_chr(overlap_genes, ~ paste(.x, collapse = ";"))
    ) %>%
    transmute(
      npc_subtype_final,
      direction,
      source,
      library,
      rank,
      term_name,
      p_value,
      adjusted_p_value,
      z_score,
      combined_score,
      intersection_size,
      intersection
    ) %>%
    arrange(adjusted_p_value, p_value)
} else {
  enrich_clean <- tibble()
}

readr::write_csv(enrich_clean, file.path(table_dir, "clean_npc_enrichr_enrichment_up_in_severe.csv"))

if (nrow(enrich_clean) > 0) {
  top_enrich <- enrich_clean %>%
    filter(adjusted_p_value < 0.10) %>%
    { if (nrow(.) == 0) enrich_clean else . } %>%
    group_by(npc_subtype_final) %>%
    arrange(adjusted_p_value, p_value, .by_group = TRUE) %>%
    slice_head(n = 6) %>%
    ungroup() %>%
    mutate(
      term_name_short = stringr::str_trunc(term_name, width = 58),
      term_label = paste0(term_name_short, " [", source, "]"),
      minus_log10_padj = -log10(pmax(adjusted_p_value, 1e-300))
    )

  readr::write_csv(top_enrich, file.path(table_dir, "clean_npc_enrichr_enrichment_top_terms_for_fig2f.csv"))

  term_levels <- top_enrich %>%
    arrange(source, p_value) %>%
    pull(term_label) %>%
    unique()

  p_enrich <- ggplot(
    top_enrich,
    aes(
      x = npc_subtype_final,
      y = factor(term_label, levels = rev(term_levels)),
      size = intersection_size,
      color = minus_log10_padj
    )
  ) +
    geom_point(alpha = 0.88) +
    scale_color_gradient(low = "#6BAED6", high = "#B2182B") +
    scale_size_continuous(range = c(2.2, 8)) +
    labs(
      title = "Enrichment of severe-upregulated genes in NPC subtypes",
      x = NULL,
      y = NULL,
      color = "-log10(adj. p)",
      size = "Gene count"
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 35, hjust = 1),
      axis.text.y = element_text(size = 8),
      legend.position = "right"
    )

  ggsave(file.path(figure_dir, "29_clean_npc_fig2f_enrichr_enrichment_dotplot.png"), p_enrich, width = 12.5, height = 8.5, dpi = 300)
} else {
  warning("No enrichment result returned. Skipping enrichment dot plot.")
}

unlink(tmp_dir, recursive = TRUE, force = TRUE)
message("DEG and enrichment analysis finished.")
