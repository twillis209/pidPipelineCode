##' From https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/
##' @importFrom dplyr group_by summarise mutate lag select inner_join filter pull
##' @export
munge_gwas_input <- function(input_data) {
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

  list("gwas_data" = gwas_data, "axis_set" = axis_set, "ylim" = ylim)
}

##' @importFrom ggplot2 ggplot geom_hline geom_point scale_x_continuous scale_y_continuous labs theme element_blank element_text
##' @importFrom ggtext element_markdown
##' @export
ggmanhattan <- function(munged_gwas_input) {
  manhplot <- ggplot(munged_gwas_input$gwas_data, aes(x = bp_cum, y = -log10(p), 
                                    color = as_factor(chr), size = -log10(p))) +
    geom_hline(yintercept = -log10(5e-8), color = "grey40", linetype = "dashed") + 
    geom_hline(yintercept = -log10(1e-5), color = "grey40", linetype = "dashed") + 
    geom_point(size = 0.3) +
    scale_x_continuous(label = munged_gwas_input$axis_set$chr, breaks = munged_gwas_input$axis_set$center) +
    scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = NULL, 
         y = "-log<sub>10</sub>(p)") + 
    theme( 
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.y = element_markdown(),
      axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
    )
}
