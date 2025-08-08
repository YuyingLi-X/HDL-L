set.seed(911)

load("ukb_chr22.1_X0_scaled.RData")  # X0 is scaled genotype X 


AOUT04 = AOUT10 = AOUT100=NULL ## alpha  cor
GOUT04 = GOUT10 = GOUT100=NULL ## alpha  cor

ns = ncol(X0)
nrep=100

for (rep in 1:nrep){
u = runif(ns)
 x = rnorm(ns); y= rnorm(ns)
   G1 = c(X0 %*% x); G2 = c(X0 %*% y)
 AOUT100 = c(AOUT100, cor(x,y))
  GOUT100 = c(GOUT100, cor(G1,G2))

 cau.snp = u<0.1  ## equal locs for causal snps
 x10 = ifelse(cau.snp,x,0); y10 = ifelse(cau.snp,y,0)
  G1 = c(X0 %*% x10); G2 = c(X0 %*% y10)
  AOUT10 = c(AOUT10,cor(x10, y10))
  GOUT10 = c(GOUT10, cor(G1,G2))

 cau.snp = u<0.04
 x04 = ifelse(cau.snp,x,0); y04 = ifelse(cau.snp,y,0)
  G1 = c(X0 %*% x04); G2 = c(X0 %*% y04)
  AOUT04 = c(AOUT04,cor(x04, y04))
  GOUT04 = c(GOUT04, cor(G1,G2))

}
summary(AOUT100); summary(AOUT10); summary(AOUT04)
summary(GOUT100); summary(GOUT10); summary(GOUT04)

# Load required libraries
library(ggplot2)
library(reshape2)

# Combine all vectors into a data frame
df <- data.frame(
  AOUT100 = AOUT100,
  GOUT100 = GOUT100,
  AOUT10  = AOUT10,
  GOUT10  = GOUT10,
  AOUT04  = AOUT04,
  GOUT04  = GOUT04
)

# Melt to long format
df_long <- melt(df, variable.name = "Measure", value.name = "Correlation")

# Define fixed x-axis order
df_long$Measure <- factor(df_long$Measure,
  levels = c("AOUT100", "GOUT100", "AOUT10", "GOUT10", "AOUT04", "GOUT04"))

# Label causal proportion
df_long$Causal <- factor(ifelse(grepl("100", df_long$Measure), "100%",
                          ifelse(grepl("10", df_long$Measure), "10%", "4%")),
                          levels = c("100%", "10%", "4%"))

# Human-readable x-axis labels with math formatting
x_labels <- c(
  "AOUT100" = expression(cor(alpha[1], alpha[2])),
  "GOUT100" = expression(cor(G[1], G[2])),
  "AOUT10"  = expression(cor(alpha[1], alpha[2])),
  "GOUT10"  = expression(cor(G[1], G[2])),
  "AOUT04"  = expression(cor(alpha[1], alpha[2])),
  "GOUT04"  = expression(cor(G[1], G[2]))
)

# Color mapping
color_map <- c("100%" = "#0072B2", "10%" = "#E69F00", "4%" = "#D55E00")

# Plot
p <- ggplot(df_long, aes(x = Measure, y = Correlation, fill = Causal)) +
  geom_boxplot(outlier.size = 0.8, width = 0.6) +
  scale_fill_manual(values = color_map) +
  scale_x_discrete(labels = x_labels) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    title = expression("True " * cor(alpha[1], alpha[2]) * " = 0"),
    x = "", y = "Realized Correlation", fill = "Proportion of causal variants"
  ) + 
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(size = 12, color = "black", angle = 15, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    legend.position = "top"
  )

# Display it
print(p)




