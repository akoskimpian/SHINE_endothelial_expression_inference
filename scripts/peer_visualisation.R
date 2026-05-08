library(tidyverse)

alpha <- read.csv("path/to/peer_output/Alpha.csv", header = FALSE)
alpha$factor <- 1:nrow(alpha)
colnames(alpha)[1] <- "alpha"

ggplot(alpha, aes(x = factor, y = alpha)) +
  geom_line() +
  geom_point() +
  labs(x = "PEER Factor", y = "Alpha (inverse variance weight)",
       title = "PEER Factor Relevance") +
  theme_bw()

# Plotting 1/Alpha to visualize actual factor variance/relevance
ggplot(alpha, aes(x = factor, y = 1/alpha)) +
  geom_line() +
  geom_point() +
  labs(x = "PEER Factor", y = "Factor Relevance (1/Alpha)",
       title = "PEER Factor Relevance")


# Create the dataframe
alpha_values <- c(277.213, 283.2, 321.972, 342.752, 293.685, 1167.32, 
                  0.00014243, 0.106725, 0.403414, 0.632274, 0.806955, 
                  2.45662, 4.21383, 1.87241, 5.37107, 6.15954, 
                  9.24807, 12.8296, 12.7042, 15.6356, 16.8369)

df <- data.frame(
  Factor = 1:length(alpha_values),
  Alpha = alpha_values
) %>%
  mutate(
    Relevance = 1 / Alpha,
    LogRelevance = log10(Relevance),
    Type = ifelse(Factor <= 6, "Known Covariate", "Hidden Factor")
  )

# Plotting with Log Scale to see everything at once
ggplot(df, aes(x = Factor, y = LogRelevance, fill = Type)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = round(LogRelevance, 1)), vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("Known Covariate" = "#E64B35FF", "Hidden Factor" = "#4DBBD5FF")) +
  labs(
    title = "PEER Factor Relevance (Log10 Scale)",
    x = "Factor Number",
    y = "Log10(1/Alpha)",
    
  ) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    
  ) +
  theme_minimal()
