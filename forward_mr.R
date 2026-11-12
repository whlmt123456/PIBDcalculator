library(VariantAnnotation)
library(gwasglue)
library(TwoSampleMR)
library(data.table)
library(MRPRESSO)
library(dplyr)
library(ggplot2)
library(patchwork)



analysis_root <- Sys.getenv("FATTY_ACID_MR_ROOT", unset = ".")
data_dir <- file.path(analysis_root, "data")
plink_bin <- Sys.getenv("PLINK_BIN", unset = "plink")
bfile <- file.path(analysis_root, "reference", "1000G_EUR")

outdir <- Sys.getenv("FATTY_ACID_MR_OUT", unset=file.path(analysis_root, "results_R1"))
set.seed(20261111)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

plot_dir <- file.path(outdir, "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

combined_plot_dir <- file.path(outdir, "combined_plots")
dir.create(combined_plot_dir, showWarnings = FALSE, recursive = TRUE)

table_dir <- file.path(outdir, "tables")
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

rds_dir <- file.path(outdir, "rds")
dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)



get_outcome_label <- function(outcome_name) {
  case_when(
    outcome_name == "APPENDACUT" ~ "Acute appendicitis",
    outcome_name == "APPENDACUT_COMPLIC" ~ "Acute appendicitis with complications",
    outcome_name == "APPENDOTH" ~ "Other appendicitis",
    TRUE ~ outcome_name
  )
}

get_exposure_label <- function(exposure_name) {
  case_when(
    exposure_name == "saturated" ~ "Saturated fatty acids",
    exposure_name == "Polyunsaturated" ~ "Polyunsaturated fatty acids",
    exposure_name == "Monounsaturated" ~ "Monounsaturated fatty acids",
    exposure_name == "OMEGA6" ~ "Omega-6 fatty acids",
    exposure_name == "OMEGA3" ~ "Omega-3 fatty acids",
    exposure_name == "Monoradio" ~ "Monounsaturated fatty acid ratio",
    exposure_name == "Polyradio" ~ "Polyunsaturated fatty acid ratio",
    exposure_name == "OMEGAradio" ~ "Omega-6/Omega-3 ratio",
    TRUE ~ exposure_name
  )
}



outcome_order <- c(
  "APPENDACUT",
  "APPENDACUT_COMPLIC",
  "APPENDOTH"
)

exposure_order <- c(
  "saturated",
  "Polyunsaturated",
  "Monounsaturated",
  "OMEGA6",
  "OMEGA3",
  "Monoradio",
  "Polyradio",
  "OMEGAradio"
)

outcome_label_map <- c(
  APPENDACUT = "Acute appendicitis",
  APPENDACUT_COMPLIC = "Acute appendicitis with complications",
  APPENDOTH = "Other appendicitis"
)

exposure_label_map <- c(
  saturated = "Saturated fatty acids",
  Polyunsaturated = "Polyunsaturated fatty acids",
  Monounsaturated = "Monounsaturated fatty acids",
  OMEGA6 = "Omega-6 fatty acids",
  OMEGA3 = "Omega-3 fatty acids",
  Monoradio = "Monounsaturated fatty acid ratio",
  Polyradio = "Polyunsaturated fatty acid ratio",
  OMEGAradio = "Omega-6/Omega-3 ratio"
)



exposure_files <- list(
  saturated        = file.path(data_dir, "exposure", "saturated_fatty_acids.gz"),
  Polyunsaturated  = file.path(data_dir, "exposure", "polyunsaturated_fatty_acids.gz"),
  Monounsaturated  = file.path(data_dir, "exposure", "monounsaturated_fatty_acids.gz"),
  OMEGA6           = file.path(data_dir, "exposure", "omega6_fatty_acids.gz"),
  OMEGA3           = file.path(data_dir, "exposure", "omega3_fatty_acids.gz"),
  Monoradio        = file.path(data_dir, "exposure", "monounsaturated_fatty_acid_ratio.gz"),
  Polyradio        = file.path(data_dir, "exposure", "polyunsaturated_fatty_acid_ratio.gz"),
  OMEGAradio       = file.path(data_dir, "exposure", "omega6_to_omega3_ratio.gz")
)

missing_exposure_files <- exposure_files[!file.exists(unlist(exposure_files))]
if (length(missing_exposure_files) > 0) {
  stop("Missing exposure files: ", paste(unlist(missing_exposure_files), collapse = ", "))
}

expected_ids <- c(saturated="ebi-a-GCST90092980",Polyunsaturated="ebi-a-GCST90092939",
  Monounsaturated="ebi-a-GCST90092928",OMEGA6="ebi-a-GCST90092933",OMEGA3="ebi-a-GCST90092931",
  Monoradio="met-d-MUFA_pct",Polyradio="ebi-a-GCST90092941",OMEGAradio="ebi-a-GCST90092934")
observed_ids <- vapply(exposure_files, function(f) VariantAnnotation::samples(VariantAnnotation::scanVcfHeader(f))[1], character(1))
stopifnot(identical(observed_ids[names(expected_ids)],expected_ids))
exposure_vcfs <- lapply(exposure_files, VariantAnnotation::readVcf)



outcome_files <- list(
  APPENDACUT = file.path(data_dir, "outcome", "finngen_R11_K11_APPENDACUT.gz"),
  APPENDACUT_COMPLIC = file.path(data_dir, "outcome", "finngen_R11_APPENDACUT_COMPLIC.gz"),
  APPENDOTH = file.path(data_dir, "outcome", "finngen_R11_K11_APPENDOTH.gz")
)

missing_outcome_files <- outcome_files[!file.exists(unlist(outcome_files))]
if (length(missing_outcome_files) > 0) {
  stop("Missing outcome files: ", paste(unlist(missing_outcome_files), collapse = ", "))
}

outcome_data <- lapply(outcome_files, function(x) {
  as.data.frame(data.table::fread(x, header = TRUE))
})



results <- data.frame()

het_results <- list()
pleio_results <- list()
mr_presso_results <- list()
single_results <- list()
res_single_results <- list()
dat_results <- list()

mr_presso_global_results <- data.frame()

leaveoneout_plot_list <- list()
single_snp_forest_plot_list <- list()
scatter_plot_list <- list()



for (exposure_name in names(exposure_vcfs)) {

  message("======================================")
  message("正在处理暴露：", exposure_name)
  message("======================================")

  vcf_exposure <- exposure_vcfs[[exposure_name]]
  exposure_label <- get_exposure_label(exposure_name)



  exposure_pval_filter <- query_gwas(
    vcf = vcf_exposure,
    pval = 5e-8
  )

  exposure_pval_filter <- gwasvcf_to_TwoSampleMR(
    vcf = exposure_pval_filter
  )

  exposure_pval_filter$exposure_name <- exposure_name
  exposure_pval_filter$exposure_label <- exposure_label

  exposure_pval_filter$exposure <- exposure_label
  exposure_pval_filter$id.exposure <- exposure_name



  exposure_clumped <- clump_data(
    exposure_pval_filter,
    clump_kb = 10000,
    clump_r2 = 0.001,
    pop = "EUR",
    plink_bin = plink_bin,
    bfile = bfile
  )

  exposure_clumped$exposure_name <- exposure_name
  exposure_clumped$exposure_label <- exposure_label
  exposure_clumped$exposure <- exposure_label
  exposure_clumped$id.exposure <- exposure_name

  if (nrow(exposure_clumped) == 0) {
    message("暴露 ", exposure_name, " clump 后无 SNP，跳过。")
    next
  }


  for (outcome_name in names(outcome_data)) {

    message("-------- 正在处理组合：", exposure_name, " × ", outcome_name, " --------")

    combo_name <- paste0(exposure_name, "__", outcome_name)
    outcome_label <- get_outcome_label(outcome_name)

    tryCatch({

      outcome <- outcome_data[[outcome_name]]



      out_dat <- format_data(
        outcome,
        snps = exposure_clumped$SNP,
        type = "outcome",
        snp_col = "rsids",
        beta_col = "beta",
        se_col = "sebeta",
        eaf_col = "af_alt",
        effect_allele_col = "alt",
        other_allele_col = "ref",
        pval_col = "pval",
        ncase_col = "ncase",
        ncontrol_col = "ncontrol",
        samplesize_col = "samplesize"
      )

      out_dat$outcome_name <- outcome_name
      out_dat$outcome_label <- outcome_label
      out_dat$outcome <- outcome_label
      out_dat$id.outcome <- outcome_name

      if (nrow(out_dat) == 0) {
        message(combo_name, " 在结局中未匹配到 SNP，跳过。")
        next
      }



      dat <- harmonise_data(
        exposure_dat = exposure_clumped,
        outcome_dat = out_dat
      )

      dat$exposure_name <- exposure_name
      dat$outcome_name <- outcome_name
      dat$exposure_label <- exposure_label
      dat$outcome_label <- outcome_label
      dat$combo <- combo_name

      dat$exposure <- exposure_label
      dat$outcome <- outcome_label
      dat$id.exposure <- exposure_name
      dat$id.outcome <- outcome_name

      dat <- dat[dat$mr_keep == TRUE, ]

      if (nrow(dat) == 0) {
        message(combo_name, " harmonise 后无可用 SNP，跳过。")
        next
      }

      dat$F_statistic <- (dat$beta.exposure/dat$se.exposure)^2
      dat_results[[combo_name]] <- dat

      write.csv(
        dat,
        file = file.path(table_dir, paste0(combo_name, "_harmonised_dat.csv")),
        row.names = FALSE
      )

      saveRDS(
        dat,
        file = file.path(rds_dir, paste0(combo_name, "_harmonised_dat.rds"))
      )



      het <- mr_heterogeneity(dat)
      het$exposure_name <- exposure_name
      het$outcome_name <- outcome_name
      het$exposure_label <- exposure_label
      het$outcome_label <- outcome_label
      het$combo <- combo_name

      het_results[[combo_name]] <- het

      write.csv(
        het,
        file = file.path(table_dir, paste0(combo_name, "_het.csv")),
        row.names = FALSE
      )



      pleio <- mr_pleiotropy_test(dat)
      pleio$exposure_name <- exposure_name
      pleio$outcome_name <- outcome_name
      pleio$exposure_label <- exposure_label
      pleio$outcome_label <- outcome_label
      pleio$combo <- combo_name

      pleio_results[[combo_name]] <- pleio

      write.csv(
        pleio,
        file = file.path(table_dir, paste0(combo_name, "_pleio.csv")),
        row.names = FALSE
      )



      mr_presso_result <- NULL
      mr_presso_global_pvalue <- NA

      if (nrow(dat) >= 4) {

        mr_presso_result <- tryCatch({

          mr_presso(
            BetaOutcome = "beta.outcome",
            BetaExposure = "beta.exposure",
            SdOutcome = "se.outcome",
            SdExposure = "se.exposure",
            OUTLIERtest = TRUE,
            DISTORTIONtest = TRUE,
            data = dat,
            NbDistribution = 1000,
            SignifThreshold = 0.05
          )

        }, error = function(e) {

          message(combo_name, " MR-PRESSO 运行失败：", e$message)
          NULL

        })

        mr_presso_results[[combo_name]] <- mr_presso_result

        if (!is.null(mr_presso_result)) {

          saveRDS(
            mr_presso_result,
            file = file.path(rds_dir, paste0(combo_name, "_mr_presso_result.rds"))
          )

          mr_presso_global_pvalue <- tryCatch({

            mr_presso_result[["MR-PRESSO results"]][["Global Test"]][["Pvalue"]]

          }, error = function(e) {

            message(combo_name, " MR-PRESSO Global Test Pvalue 提取失败：", e$message)
            NA

          })
        }

      } else {

        message(combo_name, " SNP 数少于 4，跳过 MR-PRESSO。")
        mr_presso_results[[combo_name]] <- NA
        mr_presso_global_pvalue <- NA
      }

      mr_presso_global_temp <- data.frame(
        exposure_name = exposure_name,
        outcome_name = outcome_name,
        exposure_label = exposure_label,
        outcome_label = outcome_label,
        combo = combo_name,
        nsnp = nrow(dat),
        MR_PRESSO_Global_Test_Pvalue = mr_presso_global_pvalue,
        stringsAsFactors = FALSE
      )

      mr_presso_global_results <- rbind(
        mr_presso_global_results,
        mr_presso_global_temp
      )

      write.csv(
        mr_presso_global_results,
        file = file.path(table_dir, "MR_PRESSO_Global_Test_Pvalue_all_combinations.csv"),
        row.names = FALSE
      )



      single <- mr_leaveoneout(dat)
      single$exposure_name <- exposure_name
      single$outcome_name <- outcome_name
      single$exposure_label <- exposure_label
      single$outcome_label <- outcome_label
      single$combo <- combo_name

      if ("exposure" %in% colnames(single)) single$exposure <- exposure_label
      if ("outcome" %in% colnames(single)) single$outcome <- outcome_label

      single_results[[combo_name]] <- single

      write.csv(
        single,
        file = file.path(table_dir, paste0(combo_name, "_leaveoneout.csv")),
        row.names = FALSE
      )

      loo_plot <- mr_leaveoneout_plot(single)[[1]] +
        labs(
          title = NULL,
          subtitle = NULL,
          caption = NULL,
          x = "MR effect estimate",
          y = "Instrument removed"
        ) +
        theme_bw(base_size = 8) +
        theme(
          plot.title = element_blank(),
          plot.subtitle = element_blank(),
          plot.caption = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 7),
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 5),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          strip.text = element_text(size = 8),
          plot.margin = margin(5, 5, 5, 5)
        )

      if (is.null(leaveoneout_plot_list[[outcome_name]])) {
        leaveoneout_plot_list[[outcome_name]] <- list()
      }
      leaveoneout_plot_list[[outcome_name]][[exposure_name]] <- loo_plot



      res_single <- mr_singlesnp(dat)
      res_single$exposure_name <- exposure_name
      res_single$outcome_name <- outcome_name
      res_single$exposure_label <- exposure_label
      res_single$outcome_label <- outcome_label
      res_single$combo <- combo_name

      if ("exposure" %in% colnames(res_single)) res_single$exposure <- exposure_label
      if ("outcome" %in% colnames(res_single)) res_single$outcome <- outcome_label

      res_single_results[[combo_name]] <- res_single

      write.csv(
        res_single,
        file = file.path(table_dir, paste0(combo_name, "_single_snp.csv")),
        row.names = FALSE
      )

      single_forest_plot <- mr_forest_plot(res_single)[[1]] +
        labs(
          title = NULL,
          subtitle = NULL,
          caption = NULL,
          x = "MR effect estimate"
        ) +
        theme_bw(base_size = 8) +
        theme(
          plot.title = element_blank(),
          plot.subtitle = element_blank(),
          plot.caption = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 7),
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 5),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          strip.text = element_text(size = 8),
          plot.margin = margin(5, 5, 5, 5)
        )

      if (is.null(single_snp_forest_plot_list[[outcome_name]])) {
        single_snp_forest_plot_list[[outcome_name]] <- list()
      }
      single_snp_forest_plot_list[[outcome_name]][[exposure_name]] <- single_forest_plot



      res <- mr(dat)

      res$exposure_name <- exposure_name
      res$outcome_name <- outcome_name
      res$exposure_label <- exposure_label
      res$outcome_label <- outcome_label
      res$combo <- combo_name
      res$MR_PRESSO_Global_Test_Pvalue <- mr_presso_global_pvalue

      if ("exposure" %in% colnames(res)) res$exposure <- exposure_label
      if ("outcome" %in% colnames(res)) res$outcome <- outcome_label

      res$OR <- exp(res$b)
      res$OR_lower <- exp(res$b - 1.96 * res$se)
      res$OR_upper <- exp(res$b + 1.96 * res$se)
      res$OR_95CI <- paste0(
        round(res$OR, 3),
        " (",
        round(res$OR_lower, 3),
        "-",
        round(res$OR_upper, 3),
        ")"
      )

      results <- rbind(results, res)

      write.csv(
        res,
        file = file.path(table_dir, paste0(combo_name, "_mr_result.csv")),
        row.names = FALSE
      )



      res_scatter <- res[1:min(3, nrow(res)), ]

      if ("exposure" %in% colnames(res_scatter)) res_scatter$exposure <- exposure_label
      if ("outcome" %in% colnames(res_scatter)) res_scatter$outcome <- outcome_label

      scatter_plot <- mr_scatter_plot(res_scatter, dat)[[1]] +
        labs(
          title = NULL,
          subtitle = NULL,
          caption = NULL,
          x = paste0("SNP effect on ", exposure_label),
          y = paste0("SNP effect on ", outcome_label)
        ) +
        theme_bw(base_size = 9) +
        theme(
          plot.title = element_blank(),
          plot.subtitle = element_blank(),
          plot.caption = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 7),
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          plot.margin = margin(5, 5, 5, 5)
        )

      if (is.null(scatter_plot_list[[outcome_name]])) {
        scatter_plot_list[[outcome_name]] <- list()
      }
      scatter_plot_list[[outcome_name]][[exposure_name]] <- scatter_plot


      message("完成组合：", combo_name)

    }, error = function(e) {

      message("组合 ", combo_name, " 出错：", e$message)

      error_log <- data.frame(
        exposure = exposure_name,
        outcome = outcome_name,
        combo = combo_name,
        error = e$message,
        stringsAsFactors = FALSE
      )

      write.table(
        error_log,
        file = file.path(outdir, "error_log.txt"),
        sep = "\t",
        row.names = FALSE,
        col.names = !file.exists(file.path(outdir, "error_log.txt")),
        append = TRUE,
        quote = FALSE
      )

    })
  }
}



make_combined_plot <- function(
    plot_nested_list,
    outcome_name,
    plot_type_title,
    ncol_value = 2
) {

  plot_list <- plot_nested_list[[outcome_name]]

  if (is.null(plot_list)) {
    message("没有找到 ", outcome_name, " 的 ", plot_type_title, " 图。")
    return(NULL)
  }

  valid_exposures <- intersect(exposure_order, names(plot_list))
  plot_list <- plot_list[valid_exposures]

  if (length(plot_list) == 0) {
    message(outcome_name, " 的 ", plot_type_title, " 图数量为 0。")
    return(NULL)
  }

  p_combined <- wrap_plots(
    plot_list,
    ncol = ncol_value,
    guides = "collect"
  ) +
    plot_annotation(
      title = paste0(outcome_label_map[[outcome_name]], " - ", plot_type_title),
      tag_levels = "A"
    ) &
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 14
      ),
      plot.tag = element_text(
        face = "bold",
        size = 12
      )
    )

  return(p_combined)
}


for (oc in outcome_order) {


  p_loo <- make_combined_plot(
    plot_nested_list = leaveoneout_plot_list,
    outcome_name = oc,
    plot_type_title = "leave-one-out analysis",
    ncol_value = 1
  )

  if (!is.null(p_loo)) {
    pdf(
      file = file.path(combined_plot_dir, paste0(oc, "_leaveoneout_all_exposures.pdf")),
      width = 12,
      height = 45
    )
    print(p_loo)
    dev.off()
  }



  p_single_forest <- make_combined_plot(
    plot_nested_list = single_snp_forest_plot_list,
    outcome_name = oc,
    plot_type_title = "single-SNP forest plot",
    ncol_value = 1
  )

  if (!is.null(p_single_forest)) {
    pdf(
      file = file.path(combined_plot_dir, paste0(oc, "_single_snp_forest_all_exposures.pdf")),
      width = 12,
      height = 45
    )
    print(p_single_forest)
    dev.off()
  }



  p_scatter <- make_combined_plot(
    plot_nested_list = scatter_plot_list,
    outcome_name = oc,
    plot_type_title = "MR scatter plot",
    ncol_value = 2
  )

  if (!is.null(p_scatter)) {
    pdf(
      file = file.path(combined_plot_dir, paste0(oc, "_scatter_all_exposures.pdf")),
      width = 14,
      height = 18
    )
    print(p_scatter)
    dev.off()
  }
}



write.csv(
  results,
  file = file.path(outdir, "all_MR_results.csv"),
  row.names = FALSE
)

saveRDS(
  results,
  file = file.path(outdir, "all_MR_results.rds")
)

saveRDS(
  het_results,
  file = file.path(outdir, "all_het_results_list.rds")
)

saveRDS(
  pleio_results,
  file = file.path(outdir, "all_pleio_results_list.rds")
)

saveRDS(
  mr_presso_results,
  file = file.path(outdir, "all_mr_presso_results_list.rds")
)

saveRDS(
  mr_presso_global_results,
  file = file.path(outdir, "MR_PRESSO_Global_Test_Pvalue_all_combinations.rds")
)

write.csv(
  mr_presso_global_results,
  file = file.path(outdir, "MR_PRESSO_Global_Test_Pvalue_all_combinations.csv"),
  row.names = FALSE
)

saveRDS(
  single_results,
  file = file.path(outdir, "all_leaveoneout_results_list.rds")
)

saveRDS(
  res_single_results,
  file = file.path(outdir, "all_single_snp_results_list.rds")
)

saveRDS(
  dat_results,
  file = file.path(outdir, "all_harmonised_dat_list.rds")
)

writeLines(
  sort(unique(unlist(lapply(dat_results, function(d) {
    as.character(d$SNP[which(d$mr_keep)])
  }), use.names = FALSE))),
  file.path(outdir, "forward_snps.txt")
)



results_selected <- subset(
  results,
  method %in% c(
    "Inverse variance weighted",
    "MR Egger",
    "Weighted median"
  )
)

write.csv(
  results_selected,
  file = file.path(outdir, "all_MR_results_IVW_Egger_WeightedMedian.csv"),
  row.names = FALSE
)

saveRDS(
  results_selected,
  file = file.path(outdir, "all_MR_results_IVW_Egger_WeightedMedian.rds")
)



results_IVW <- subset(
  results,
  method %in% c("Inverse variance weighted")
)

write.csv(
  results_IVW,
  file = file.path(outdir, "all_MR_results_IVW_only.csv"),
  row.names = FALSE
)

head(results_IVW)



cat("所有结果已保存至：\n")
cat(outdir, "\n\n")

cat("表格文件夹：\n")
cat(table_dir, "\n\n")

cat("单图文件夹：\n")
cat(plot_dir, "\n\n")

cat("合并大图文件夹：\n")
cat(combined_plot_dir, "\n\n")

cat("RDS 文件夹：\n")
cat(rds_dir, "\n\n")

cat("合并大图 PDF 数量：\n")
print(length(list.files(combined_plot_dir, pattern = "\\.pdf$")))

cat("合并大图文件列表：\n")
print(list.files(combined_plot_dir, pattern = "\\.pdf$"))

cat("MR-PRESSO Global Test P value 汇总文件：\n")
cat(file.path(outdir, "MR_PRESSO_Global_Test_Pvalue_all_combinations.csv"), "\n")

suppressPackageStartupMessages({
  library(dplyr)
  library(forestplot)
  library(grid)
})

fpDrawTriangleCI <- function(lower_limit,
                             estimate,
                             upper_limit,
                             size,
                             y.offset = 0.5,
                             clr.line,
                             clr.marker,
                             lwd,
                             lty = 1,
                             vertices,
                             vertices.height = 0.1,
                             shapes_gp,
                             shape_coordinates,
                             ...) {
  if (is.na(estimate) || is.na(lower_limit) || is.na(upper_limit)) return(invisible(NULL))

  gp_line <- forestplot:::prGetShapeGp(
    shapes_gp = shapes_gp,
    coords = shape_coordinates,
    object = "lines",
    default = grid::gpar(col = clr.line, lwd = lwd, lty = lty)
  )
  gp_vert <- forestplot:::prGetShapeGp(
    shapes_gp = shapes_gp,
    coords = shape_coordinates,
    object = "vertices",
    default = gp_line
  )
  gp_box <- forestplot:::prGetShapeGp(
    shapes_gp = shapes_gp,
    coords = shape_coordinates,
    object = "box",
    default = grid::gpar(col = clr.marker, fill = clr.marker)
  )

  grid::grid.segments(
    x0 = grid::unit(lower_limit, "native"),
    x1 = grid::unit(upper_limit, "native"),
    y0 = grid::unit(y.offset, "npc"),
    y1 = grid::unit(y.offset, "npc"),
    gp = gp_line
  )

  if (isTRUE(vertices)) {
    h <- grid::unit(vertices.height, "npc")
    grid::grid.segments(
      x0 = grid::unit(lower_limit, "native"),
      x1 = grid::unit(lower_limit, "native"),
      y0 = grid::unit(y.offset, "npc") - h / 2,
      y1 = grid::unit(y.offset, "npc") + h / 2,
      gp = gp_vert
    )
    grid::grid.segments(
      x0 = grid::unit(upper_limit, "native"),
      x1 = grid::unit(upper_limit, "native"),
      y0 = grid::unit(y.offset, "npc") - h / 2,
      y1 = grid::unit(y.offset, "npc") + h / 2,
      gp = gp_vert
    )
  }

  if (!grid::is.unit(size)) size <- grid::unit(size, "snpc")
  grid::grid.points(
    x = grid::unit(estimate, "native"),
    y = grid::unit(y.offset, "npc"),
    pch = 2,
    size = size,
    gp = grid::gpar(col = gp_box$col %||% gp_box$fill %||% "black")
  )

  invisible(NULL)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

if (!exists("outdir")) {
  outdir <- getwd()
}

plot_dir <- file.path(outdir, "forest_plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

results_IVW <- subset(
  results,
  method %in% c("Inverse variance weighted")
)

results_IVW <- results_IVW %>%
  mutate(
    OR = exp(b),
    OR_lower = exp(b - 1.96 * se),
    OR_upper = exp(b + 1.96 * se),
    OR_95CI = paste0(
      sprintf("%.3f", OR),
      " (",
      sprintf("%.3f", OR_lower),
      "-",
      sprintf("%.3f", OR_upper),
      ")"
    ),
    P_value = case_when(
      pval < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", pval)
    ),
    sig05 = ifelse(!is.na(pval) & pval < 0.05, TRUE, FALSE),
    significance = case_when(pval < 0.05/24 ~ "Bonferroni", pval < 0.05 ~ "Nominal", TRUE ~ "Not significant")
  )



results_IVW <- results_IVW %>%
  mutate(
    outcome_label = case_when(
      outcome_name == "APPENDACUT" ~ "Acute appendicitis",
      outcome_name == "APPENDACUT_COMPLIC" ~ "Acute appendicitis with complications",
      outcome_name == "APPENDOTH" ~ "Other appendicitis",
      TRUE ~ outcome_name
    ),
    exposure_label = case_when(
      exposure_name == "saturated" ~ "Saturated fatty acids",
      exposure_name == "Polyunsaturated" ~ "Polyunsaturated fatty acids",
      exposure_name == "Monounsaturated" ~ "Monounsaturated fatty acids",
      exposure_name == "OMEGA6" ~ "Omega-6 fatty acids",
      exposure_name == "OMEGA3" ~ "Omega-3 fatty acids",
      exposure_name == "Monoradio" ~ "Monounsaturated fatty acid ratio",
      exposure_name == "Polyradio" ~ "Polyunsaturated fatty acid ratio",
      exposure_name == "OMEGAradio" ~ "Omega-6/Omega-3 ratio",
      TRUE ~ exposure_name
    )
  )



outcome_order <- c(
  "Acute appendicitis",
  "Acute appendicitis with complications",
  "Other appendicitis"
)

exposure_order <- c(
  "Saturated fatty acids",
  "Monounsaturated fatty acids",
  "Polyunsaturated fatty acids",
  "Omega-3 fatty acids",
  "Omega-6 fatty acids",
  "Monounsaturated fatty acid ratio",
  "Polyunsaturated fatty acid ratio",
  "Omega-6/Omega-3 ratio"
)

results_IVW <- results_IVW %>%
  mutate(
    outcome_label = factor(outcome_label, levels = outcome_order),
    exposure_label = factor(exposure_label, levels = exposure_order)
  ) %>%
  arrange(outcome_label, exposure_label)

plot_list <- list()

for (oc in outcome_order) {

  subdat <- results_IVW %>%
    filter(outcome_label == oc)

  if (nrow(subdat) == 0) next

  group_row <- data.frame(
    Exposure = oc,
    Method = "",
    OR_95CI = "",
    P_value = "",
    OR = NA,
    OR_lower = NA,
    OR_upper = NA,
    sig05 = NA,
    is_summary = TRUE,
    stringsAsFactors = FALSE
  )

  data_rows <- data.frame(
    Exposure = paste0("  ", as.character(subdat$exposure_label)),
    Method = "IVW",
    OR_95CI = subdat$OR_95CI,
    P_value = subdat$P_value,
    OR = subdat$OR,
    OR_lower = subdat$OR_lower,
    OR_upper = subdat$OR_upper,
    sig05 = subdat$sig05,
    is_summary = FALSE,
    stringsAsFactors = FALSE
  )

  plot_list[[oc]] <- rbind(group_row, data_rows)
}

forest_data <- do.call(rbind, plot_list)

tabletext <- cbind(
  c("Outcome / Exposure", forest_data$Exposure),
  c("Method", forest_data$Method),
  c("OR (95% CI)", forest_data$OR_95CI),
  c("P value", forest_data$P_value)
)

mean_value  <- c(NA, forest_data$OR)
lower_value <- c(NA, forest_data$OR_lower)
upper_value <- c(NA, forest_data$OR_upper)

is_summary <- c(TRUE, forest_data$is_summary)

sig_vec <- c(FALSE, ifelse(is.na(forest_data$sig05), FALSE, forest_data$sig05))
sig_vec[is_summary] <- FALSE

col_vec <- ifelse(sig_vec, "blue4", "black")
col_vec[is_summary] <- NA

shapes_gp <- forestplot::fpShapesGp(
  default  = grid::gpar(col = "black", fill = "black"),
  summary  = grid::gpar(col = "black", fill = "black"),
  box      = lapply(seq_along(col_vec), function(i) {
    grid::gpar(
      col  = if (is.na(col_vec[i])) "transparent" else col_vec[i],
      fill = if (is.na(col_vec[i])) "transparent" else col_vec[i]
    )
  }),
  lines    = lapply(seq_along(col_vec), function(i) {
    grid::gpar(col = if (is.na(col_vec[i])) "transparent" else col_vec[i])
  }),
  vertices = lapply(seq_along(col_vec), function(i) {
    grid::gpar(col = if (is.na(col_vec[i])) "transparent" else col_vec[i])
  })
)

fn_list_all <- lapply(seq_along(mean_value), function(i) {
  if (is_summary[i] || is.na(mean_value[i])) {
    forestplot::fpDrawNormalCI
  } else if (sig_vec[i]) {
    fpDrawTriangleCI
  } else {
    forestplot::fpDrawCircleCI
  }
})

fig_all <- forestplot::forestplot(
  labeltext = tabletext,
  mean  = mean_value,
  lower = lower_value,
  upper = upper_value,
  zero = 1,
  boxsize = 0.25,
  graph.pos = 3,
  graphwidth = grid::unit(0.30, "npc"),
  xticks = c(0.7, 0.85, 1, 1.15, 1.3),
  clip = c(0.6, 1.5),
  xlog = FALSE,
  is.summary = is_summary,
  hrzl_lines = list(
    "1" = grid::gpar(lty = 1, lwd = 2),
    "2" = grid::gpar(lty = 1, lwd = 1)
  ),
  txt_gp = forestplot::fpTxtGp(
    label = grid::gpar(cex = 0.72),
    ticks = grid::gpar(cex = 0.80),
    xlab  = grid::gpar(cex = 0.85),
    title = grid::gpar(cex = 1.00)
  ),
  lwd.zero = 1,
  lwd.ci = 1.4,
  lwd.xaxis = 1,
  lty.ci = 1,
  ci.vertices = TRUE,
  ci.vertices.height = 0.12,
  line.margin = grid::unit(5, "mm"),
  colgap = grid::unit(5, "mm"),
  title = "Causal effects of fatty acid-related traits on appendicitis outcomes",
  shapes_gp = shapes_gp,
  fn.ci_norm = fn_list_all,
  col = forestplot::fpColors(zero = "black")
)

fig_all

pdf(
  file = file.path(plot_dir, "IVW_forest_plot_all_outcomes.pdf"),
  width = 10,
  height = 12
)
print(fig_all)
dev.off()

for (oc in outcome_order) {

  subdat <- results_IVW %>%
    filter(outcome_label == oc) %>%
    arrange(OR)

  if (nrow(subdat) == 0) next

  tabletext_sub <- cbind(
    c("Exposure", as.character(subdat$exposure_label)),
    c("Method", rep("IVW", nrow(subdat))),
    c("OR (95% CI)", subdat$OR_95CI),
    c("P value", subdat$P_value)
  )

  mean_sub  <- c(NA, subdat$OR)
  lower_sub <- c(NA, subdat$OR_lower)
  upper_sub <- c(NA, subdat$OR_upper)

  sig_sub <- ifelse(!is.na(subdat$pval) & subdat$pval < 0.05, TRUE, FALSE)
  col_sub <- c(NA, ifelse(sig_sub, "blue4", "black"))

  shapes_gp_sub <- forestplot::fpShapesGp(
    default  = grid::gpar(col = "black", fill = "black"),
    summary  = grid::gpar(col = "black", fill = "black"),
    box      = lapply(seq_along(col_sub), function(i) {
      grid::gpar(
        col  = if (is.na(col_sub[i])) "transparent" else col_sub[i],
        fill = if (is.na(col_sub[i])) "transparent" else col_sub[i]
      )
    }),
    lines    = lapply(seq_along(col_sub), function(i) {
      grid::gpar(col = if (is.na(col_sub[i])) "transparent" else col_sub[i])
    }),
    vertices = lapply(seq_along(col_sub), function(i) {
      grid::gpar(col = if (is.na(col_sub[i])) "transparent" else col_sub[i])
    })
  )

  is_summary_sub <- c(TRUE, rep(FALSE, nrow(subdat)))
  sig_vec_sub <- c(FALSE, sig_sub)

  fn_list_sub <- lapply(seq_along(mean_sub), function(i) {
    if (is_summary_sub[i] || is.na(mean_sub[i])) {
      forestplot::fpDrawNormalCI
    } else if (sig_vec_sub[i]) {
      fpDrawTriangleCI
    } else {
      forestplot::fpDrawCircleCI
    }
  })

  fig_sub <- forestplot::forestplot(
    labeltext = tabletext_sub,
    mean  = mean_sub,
    lower = lower_sub,
    upper = upper_sub,
    zero = 1,
    boxsize = 0.30,
    graph.pos = 3,
    graphwidth = grid::unit(0.32, "npc"),
    xticks = c(0.7, 0.85, 1, 1.15, 1.3),
    clip = c(0.6, 1.5),
    xlog = FALSE,
    is.summary = is_summary_sub,
    hrzl_lines = list(
      "1" = grid::gpar(lty = 1, lwd = 2),
      "2" = grid::gpar(lty = 1, lwd = 1)
    ),
    txt_gp = forestplot::fpTxtGp(
      label = grid::gpar(cex = 0.78),
      ticks = grid::gpar(cex = 0.82),
      xlab  = grid::gpar(cex = 0.85),
      title = grid::gpar(cex = 1.00)
    ),
    lwd.zero = 1,
    lwd.ci = 1.5,
    lwd.xaxis = 1,
    lty.ci = 1,
    ci.vertices = TRUE,
    ci.vertices.height = 0.15,
    line.margin = grid::unit(5, "mm"),
    colgap = grid::unit(5, "mm"),
    title = oc,
    shapes_gp = shapes_gp_sub,
    fn.ci_norm = fn_list_sub,
    col = forestplot::fpColors(zero = "black")
  )

  pdf_name <- dplyr::case_when(
    oc == "Acute appendicitis" ~ "IVW_forest_plot_APPENDACUT.pdf",
    oc == "Acute appendicitis with complications" ~ "IVW_forest_plot_APPENDACUT_COMPLIC.pdf",
    oc == "Other appendicitis" ~ "IVW_forest_plot_APPENDOTH.pdf",
    TRUE ~ paste0("IVW_forest_plot_", gsub(" ", "_", oc), ".pdf")
  )

  pdf(
    file = file.path(plot_dir, pdf_name),
    width = 9,
    height = 7
  )
  print(fig_sub)
  dev.off()
}

write.csv(
  results_IVW,
  file = file.path(plot_dir, "IVW_forest_plot_data.csv"),
  row.names = FALSE
)

list.files(plot_dir)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(outdir, "sessionInfo.txt")
)


extract_p <- function(x) if (length(x)) x[[1]] else NA
strength_rows <- lapply(names(dat_results),function(k) {
 d <- dat_results[[k]]; f <- d$F_statistic
 standardized <- d$exposure_name[1] != "Monoradio"
 data.frame(combo=k,n=nrow(d),minimum_F=min(f),median_F=median(f),maximum_F=max(f),
            approximate_R2=if(standardized) sum(2*d$eaf.exposure*(1-d$eaf.exposure)*d$beta.exposure^2) else NA)
})
write.csv(do.call(rbind,strength_rows),file.path(outdir,"Instrument_strength.csv"),row.names=FALSE)
presso_rows <- lapply(names(dat_results),function(k) {
 x <- mr_presso_results[[k]]; d <- dat_results[[k]]
 if(!is.list(x)) x <- list()
 z <- x[["MR-PRESSO results"]]; m <- x[["Main MR results"]]
 dt <- z[["Distortion Test"]]; ii <- dt[["Outliers Indices"]]
 b <- extract_p(m[["Causal Estimate"]][2]); s <- extract_p(m[["Sd"]][2])
 data.frame(combo=k,global_p=extract_p(z[["Global Test"]][["Pvalue"]]),
  outliers=if(is.numeric(ii))paste(d$SNP[ii],collapse=";")else NA_character_,
  corrected_OR=exp(b),lower=exp(b-1.96*s),upper=exp(b+1.96*s),
  corrected_p=extract_p(m[["P-value"]][2]),distortion_p=extract_p(dt[["Pvalue"]]))
})
write.csv(do.call(rbind,presso_rows),file.path(outdir,"MR_PRESSO_corrected.csv"),row.names=FALSE)

pub <- results_IVW
pub$outcome_label <- factor(as.character(pub$outcome_label),
 levels=c("Acute appendicitis","Acute appendicitis with complications","Other appendicitis"))
pub$significance <- factor(pub$significance,levels=c("Bonferroni","Nominal","Not significant"))
p_pub <- ggplot(pub,aes(x=OR,y=exposure_label,shape=significance))+
 geom_vline(xintercept=1,colour="grey60")+
 geom_segment(aes(x=OR_lower,xend=OR_upper,yend=exposure_label),colour="grey35")+
 geom_point(size=2.5,colour="#143E5A")+
 scale_shape_manual(values=c(Bonferroni=18,Nominal=2,`Not significant`=16),drop=FALSE)+
 scale_x_log10()+facet_grid(outcome_label~.,scales="free_y",space="free_y")+
 labs(x="Odds ratio (95% CI)",y=NULL,shape=NULL)+theme_bw(base_size=10)+
 theme(legend.position="bottom")
ggsave(file.path(outdir,"Figure_2.tif"),p_pub,width=7,height=9,dpi=300,compression="lzw")
