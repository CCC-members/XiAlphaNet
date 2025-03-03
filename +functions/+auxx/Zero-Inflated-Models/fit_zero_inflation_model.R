library(glmmTMB)
data <- read.csv('temp_data.csv')
model <- glmmTMB(y ~ x + I(x^2) + (1 | group), data = data, ziformula = ~x + I(x^2), family = nbinom2)
coeffs <- coef(model)
conditional_coeffs <- coeffs$cond
zero_infl_coeffs <- coeffs$zi
write.table(conditional_coeffs, file = 'conditional_coeffs.csv', row.names = FALSE, col.names = TRUE)
write.table(zero_infl_coeffs, file = 'zero_infl_coeffs.csv', row.names = FALSE, col.names = TRUE)
