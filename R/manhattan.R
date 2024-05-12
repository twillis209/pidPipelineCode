##' From https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/
##' @importFrom dplyr group_by summarise mutate lag select inner_join filter pull
##' @importFrom data.table data.table
##' @export
munge_gwas_input <- function(gwas_dat) {
  data_cum <- gwas_dat %>% 
    group_by(chr) %>% 
    summarise(max_bp = max(bp)) %>% 
    mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
    select(chr, bp_add)

  gwas_data <- gwas_dat %>% 
    inner_join(data_cum, by = "chr") %>% 
    mutate(bp_cum = as.numeric(bp + bp_add))

  axis_set <- gwas_data %>% 
    group_by(chr) %>% 
    summarize(center = mean(bp_cum))

                                        # Ignore typically large peak from MHC when setting y scale
  ylim <- gwas_data[chr != 6] %>% 
    filter(p == min(p)) %>% 
    mutate(ylim = abs(floor(log10(p))) + 2) %>% 
    pull(ylim)

  list("gwas_data" = gwas_data, "axis_set" = data.table(axis_set), "ylim" = ylim)
}

##' @importFrom ggplot2 ggplot geom_hline geom_point scale_x_continuous scale_y_continuous labs theme element_blank element_text
##' @param x_color_palette per-chromosome palette
##' @importFrom ggtext element_markdown
##' @export
ggmanhattan <- function(munged_gwas_input, x_color_palette, ylim) {
  manhplot <- ggplot(munged_gwas_input$gwas_data, aes(x = bp_cum, y = -log10(p), color = as_factor(chr))) +
    geom_point(size = 0.3) +
    scale_x_continuous(label = munged_gwas_input$axis_set$chr, breaks = munged_gwas_input$axis_set$center) +
    scale_color_manual(values = x_color_palette) +
    scale_y_neglog10(limits = c(1, ylim))+
    labs(x = NULL,
         y = "-log<sub>10</sub>(p)")
}
