#' @importFrom ggplot2 theme_bw theme element_text
#' @export
theme_thesis <- function() {
  theme_bw()+
    theme(
      axis.title = element_text(size=12),
      plot.title=element_text(hjust=0.5, size=12),
      strip.text=element_text(size=10),
      axis.text.x=element_text(size=6, angle=90, color="black"),
      axis.text.y=element_text(size=10, color="black"),
      legend.title=element_text(size=10),
      legend.text=element_text(size=10)
    )
}
