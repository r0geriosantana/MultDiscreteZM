# Soluciona os falsos positivos do CRAN para variáveis do ggplot2
utils::globalVariables(c(
  "theoretical_quantiles", "sample_quantiles",
  "lower", "upper", "is_outlier", "label"
))
#' Q-Q Plot for Model Residuals with Outlier Detection
#'
#' @description Internal helper function to generate a diagnostic Q-Q plot
#' for randomized quantile residuals, including confidence bands and the Shapiro-Wilk test.
#'
#' @param normal_residuals Numeric vector of residuals to be evaluated.
#'
#' @keywords internal
#' @noRd
qq_plot_zm <- function(normal_residuals) {

  # --- 1. Basic Definitions ---
  # FIX: 'n' must be defined before creating the indices sequence
  n <- length(normal_residuals)
  original_indices <- 1:n

  # --- 2. Shapiro-Wilk Test for Normality ---
  test <- stats::shapiro.test(normal_residuals)
  Wstats <- test$statistic
  pw <- test$p.value

  # --- 3. Build the Data Frame for ggplot ---
  order_map <- order(normal_residuals)
  P <- stats::ppoints(n)

  df_qq <- data.frame(
    original_index = original_indices[order_map],  # Sorted indices
    sample_quantiles = sort(normal_residuals),     # Sorted residuals
    theoretical_quantiles = stats::qnorm(P)        # Theoretical quantiles
  )

  # --- 4. Reference Line Calculations ---
  Q.x <- stats::quantile(df_qq$sample_quantiles, c(0.25, 0.75))
  Q.z <- stats::qnorm(c(0.25, 0.75))
  b <- diff(Q.x) / diff(Q.z)
  coef <- c(Q.x[1] - b * Q.z[1], b)

  # --- 5. Confidence Bands ---
  conf_level <- 0.95
  zz <- stats::qnorm(1 - (1 - conf_level) / 2)

  # Standard Error calculated for each theoretical quantile
  SE <- (coef[2] / stats::dnorm(df_qq$theoretical_quantiles)) * sqrt(P * (1 - P) / n)

  df_qq$fit_value <- coef[1] + coef[2] * df_qq$theoretical_quantiles
  df_qq$upper <- df_qq$fit_value + zz * SE
  df_qq$lower <- df_qq$fit_value - zz * SE

  # --- 6. Outlier Identification ---
  df_qq$is_outlier <- with(df_qq, sample_quantiles > upper | sample_quantiles < lower)
  df_qq$label <- ifelse(df_qq$is_outlier, as.character(df_qq$original_index), "")

  # --- 7. Generate ggplot2 Object ---
  plot_residual <- ggplot2::ggplot(df_qq, ggplot2::aes(x = theoretical_quantiles, y = sample_quantiles)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.1) +
    ggplot2::geom_abline(intercept = coef[1], slope = coef[2], col = 'red4', lwd = 0.7) +
    ggplot2::geom_point(ggplot2::aes(color = is_outlier, shape = is_outlier), size = 2) +
    ggplot2::scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"),
                                labels = c("No Outlier", "Outlier")) +
    ggplot2::scale_shape_manual(values = c("FALSE" = 19, "TRUE" = 8),
                                labels = c("No Outlier", "Outlier")) +
    ggplot2::geom_text(ggplot2::aes(label = label), hjust = -0.5, vjust = 0.5, size = 3, color = "red") +
    ggplot2::labs(
      title = "Q-Q Plot of Model Residuals with Outlier Detection",
      x = "Theoretical Quantiles",
      y = "Sample Residuals",
      caption = paste("Shapiro-Wilk test, statistic =", round(Wstats, 3), ", p-value =", round(pw, 8)),
      color = "Point",
      shape = "Point"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "bottom")

  return(plot_residual)
}
