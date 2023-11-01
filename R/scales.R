
#' @importFrom scales trans_new
neglog_trans <- function(base = exp(1)){
  trans_new("neglog",
            transform = function(x) -log(x, base),
            inverse = function(x) base^(-x),
            domain = c(1e-100, Inf)
            )
}

neglog10_trans <- function(){ neglog_trans(base = 10) }

#' @importFrom scales trans_breaks trans_format math_format
#' @importFrom ggplot2 scale_x_continuous
#' @export
scale_x_neglog10 <- function(...){
  scale_x_continuous(..., 
                     trans = neglog10_trans(), 
                     breaks = scales::trans_breaks(function(x) {log10(x)*-1}, function(x){10^(-1*x)}), 
                     labels = scales::trans_format(function(x) {log10(x)*-1}, scales::math_format(.x)))
}

#' @importFrom scales trans_breaks trans_format math_format
#' @importFrom ggplot2 scale_y_continuous
#' @export
scale_y_neglog10 <- function(...){
  scale_y_continuous(..., 
                     trans = neglog10_trans(), 
                     breaks = scales::trans_breaks(function(x) {log10(x)*-1}, function(x){10^(-1*x)}), 
                     labels = scales::trans_format(function(x) {log10(x)*-1}, scales::math_format(.x)))
}
