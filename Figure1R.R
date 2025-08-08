set.seed(911)

load("ukb_chr22.1_X0_scaled.RData")  # X0 is scaled genotype X 

# Parameters
N1 <- 10000  # cohort 1
N2 <- 10000  # cohort 2
M <- ncol(X0)    # SNPs
r_g <- 0  # target genetic correlation
sigma_1 <- 1
sigma_2 <- 2
n_sim <- 100

# Storage
cor_G1G2_list <- numeric(n_sim)
cor_alpha_list <- numeric(n_sim)

cor_G1G2_sp_list <- numeric(n_sim)
cor_alpha_sp_list <- numeric(n_sim)

cor_G1G2_sp2_list <- numeric(n_sim)
cor_alpha_sp2_list <- numeric(n_sim)

# 1. Full genotype matrix and split
X_full <- X0[1:(N1 + N2),]
X1 <- X_full[1:N1, ]
X2 <- X_full[(N1 + 1):(N1 + N2), ]

require(svMisc)
for (sim in 1:n_sim) {

  # 2. Simulate alpha_1 and a starting alpha_2
  alpha_1 <- rnorm(M, sd = sigma_1)
  alpha_2 <- rnorm(M, sd = sigma_2)

  M1 <- round(M*.1)
  idx <- sample(1:M, M1)
  alpha_1_sp <- rnorm(M1, sd = sigma_1)
  alpha_2_sp <- rnorm(M1, sd = sigma_2)

  M2 <- round(M*.04)
  idx2 <- sample(1:M, M2)
  alpha_1_sp2 <- rnorm(M2, sd = sigma_1)
  alpha_2_sp2 <- rnorm(M2, sd = sigma_2)

  # 3. Create G1 and G2
  G1_full <- X_full %*% alpha_1
  G2_full <- X_full %*% alpha_2

  G1_full_sp <- X_full[,idx] %*% alpha_1_sp
  G2_full_sp <- X_full[,idx] %*% alpha_2_sp

  G1_full_sp2 <- X_full[,idx2] %*% alpha_1_sp2
  G2_full_sp2 <- X_full[,idx2] %*% alpha_2_sp2

  # 4. Regress G2 on G1, get residual R
  beta_hat <- as.numeric(crossprod(G1_full, G2_full)) / as.numeric(crossprod(G1_full))
  G2_full_hat <- beta_hat * G1_full
  R_full <- as.vector(G2_full - G2_full_hat)  # R ⟂ G1

  beta_hat_sp <- as.numeric(crossprod(G1_full_sp, G2_full_sp)) / as.numeric(crossprod(G1_full_sp))
  G2_full_sp_hat <- beta_hat_sp * G1_full_sp
  R_full_sp <- as.vector(G2_full_sp - G2_full_sp_hat)  # R ⟂ G1

  beta_hat_sp2 <- as.numeric(crossprod(G1_full_sp2, G2_full_sp2)) / as.numeric(crossprod(G1_full_sp2))
  G2_full_sp2_hat <- beta_hat_sp2 * G1_full_sp2
  R_full_sp2 <- as.vector(G2_full_sp2 - G2_full_sp2_hat)  # R ⟂ G1

  # 5. Construct G2 in cohort 2 with desired correlation to G1
  G1_std <- as.vector(scale(G1_full))
  R_std  <- as.vector(scale(R_full))
  G2_full <- r_g * G1_std + sqrt(1 - r_g^2) * R_std
  G2 <- G2_full[(N1 + 1):(N1 + N2)]

  G1_std_sp <- as.vector(scale(G1_full_sp))
  R_std_sp  <- as.vector(scale(R_full_sp))
  G2_full_sp <- r_g * G1_std_sp + sqrt(1 - r_g^2) * R_std_sp
  G2_sp <- G2_full_sp[(N1 + 1):(N1 + N2)]

  G1_std_sp2 <- as.vector(scale(G1_full_sp2))
  R_std_sp2  <- as.vector(scale(R_full_sp2))
  G2_full_sp2 <- r_g * G1_std_sp2 + sqrt(1 - r_g^2) * R_std_sp2
  G2_sp2 <- G2_full_sp2[(N1 + 1):(N1 + N2)]

  # 6. Regress G2 on X2 to get alpha_2
  alpha_2 <- coef(lm(G2 ~ X2 - 1))

  alpha_2_sp <- coef(lm(G2_sp ~ X2[,idx] - 1))

  alpha_2_sp2 <- coef(lm(G2_sp2 ~ X2[,idx2] - 1))

  # 7. Compute correlations
  cor_G1G2 <- cor(G1_full, G2_full)  # match indices arbitrarily
  cor_alpha <- cor(alpha_1, alpha_2, use = 'pairwise.complete.obs')

  cor_G1G2_sp <- cor(G1_full_sp, G2_full_sp)  # match indices arbitrarily
  cor_alpha_sp <- cor(alpha_1_sp, alpha_2_sp)

  cor_G1G2_sp2 <- cor(G1_full_sp2, G2_full_sp2)  # match indices arbitrarily
  cor_alpha_sp2 <- cor(alpha_1_sp2, alpha_2_sp2)

  cor_G1G2_list[sim] <- cor_G1G2
  cor_alpha_list[sim] <- cor_alpha

  cor_G1G2_sp_list[sim] <- cor_G1G2_sp
  cor_alpha_sp_list[sim] <- cor_alpha_sp

  cor_G1G2_sp2_list[sim] <- cor_G1G2_sp2
  cor_alpha_sp2_list[sim] <- cor_alpha_sp2

  progress(sim/n_sim*100)
}

# Plot results
df <- data.frame(
  #cor_G1_R = cor_G1R_list,
  cor_G1_G2 = cor_G1G2_list,
  cor_alpha = cor_alpha_list,

  cor_G1_G2_sp = cor_G1G2_sp_list,
  cor_alpha_sp = cor_alpha_sp_list,

  cor_G1_G2_sp2 = cor_G1G2_sp2_list,
  cor_alpha_sp2 = cor_alpha_sp2_list
)

# Load ggplot2 for visualizations
library(ggplot2)
library(reshape2)

df_melted <- melt(df)

# Add a group label based on the variable name
df_melted$group <- rep(c(
  "100%", 
  "10%",
  "4%"), each = n_sim*2
)

df_melted <- df_melted[df_melted$variable %in% c("cor_G1_G2", "cor_alpha", "cor_alpha_sp", "cor_alpha_sp2"),]

# Set factor levels to control order
df_melted$variable <- factor(df_melted$variable, levels = c(
  "cor_G1_G2", "cor_alpha", "cor_G1_G2_sp", "cor_alpha_sp", "cor_G1_G2_sp2", "cor_alpha_sp2"
))

# Human-readable x-axis labels with math formatting
x_labels <- c(
  "cor_G1_G2" = expression(cor(G[1], G[2])),
  "cor_alpha"  = expression(cor(alpha[1], alpha[2])),
  "cor_alpha_sp"  = expression(cor(alpha[1], alpha[2])),
  "cor_alpha_sp2"  = expression(cor(alpha[1], alpha[2]))
)

# Set desired legend order
df_melted$group <- factor(df_melted$group, levels = c("100%", "10%", "4%"))

# Plot
p2 <- ggplot(df_melted, aes(x = variable, y = value, fill = group)) +
  geom_boxplot() +
  geom_hline(yintercept = r_g, linetype = "dashed", color = "gray") +
  scale_fill_manual(values = c("100%" = "#0072B2", "10%" = "#E69F00", "4%" = "#D55E00")) +
  labs(
    title = expression("Forced " * cor(G[1], G[2]) * " = 0"),
    x = "", y = "Realized Correlation", fill = "Proportion of causal variants"
  )  +
  scale_x_discrete(labels = x_labels) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(size = 11),
    legend.position = "top"
  )

# Print plots
print(p2)


