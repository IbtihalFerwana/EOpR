#####################################################################################################
# This script was written by Ibtihal Ferwana, for any inquiries please contact iferwna2@illinois.edu
#####################################################################################################



source("algorithm.R")

# first row is the treated row
application = 'prop99_abadie_reshaped'
df = read.csv(paste("data/",application,".csv",sep=""), stringsAsFactors = F, check.names=FALSE)
pre_treat_end_col_name = 19 # can be table index or column name, e.g. 17 or "1971" for basque
if(is.character(pre_treat_end_col_name)==TRUE)
{
  pre_treat_end = which(colnames(df)==pre_treat_end_col_name)
} else
{
  pre_treat_end = pre_treat_end_col_name
}
u = df[1:1,] # treated unit
u = as.numeric(u)
A = df[2:dim(df)[1],] # the matrix without treated unit
xx = seq(1,pre_treat_end-1) # the indicies of pre-treatment

# optimize for lambda
regs = seq(0.01,0.99,0.01)
vals = c()
for(reg in regs){
  eopr_output = get_ubar(A, u, xx, reg)
  ubar = eopr_output$ubar
  val = sqrt(mean((ubar[xx]-u[xx])^2))
  vals = append(vals, val)
}
# EOpR effects
min_reg = which.min(vals)
eopr_output = get_ubar(A, u, xx, regs[min_reg])
write.csv(eopr_output$ubar, file = paste("results/",application,"_estimates.csv"))

# EOpR errors
worst_bounds = error_bars(eopr_output$PHI, xx, eopr_output$R, eopr_output$Q, A, eopr_output$ubar)
max_error= worst_bounds$max_error

# SC effects
sc_mode = synth_control_est(A,u, xx)
sc_eff =sc_mode$effects

# plot trajectories
pdf(paste("figures/",application,"causal_effects_vis.pdf"), width = 10, height =8)
par(mar=c(10, 4, 4, 4))
plot(colnames(df), u, type="l", xlab = "Year", ylab = "Outcome", col = "black", ylim=c(min(eopr_output$ubar)-1,max(eopr_output$ubar)+1))
lines(colnames(df), eopr_output$ubar, type="l", col = "blue", lty=2, lwd=1.5)
lines(colnames(df), sc_eff, type="l", col = "orange", lty=2, lwd=1.5)
abline(v= pre_treat_end_col_name, col="gray", lwd=1.5, lty=4)
legend("topleft", legend = c("Treated unit", "EOpR","Synthetic Control",  "Intervention"), ncol =1, col=c("black", "blue","orange", "gray"), lty=c(1,2,3,4))
dev.off()

# plot EOpR with worst case lines
pdf(paste("figures/",application,"_error_bars_vis.pdf"), width = 10, height =8)
par(mar=c(10, 4, 4, 4))
plot(colnames(df), u, type="l", xlab = "Year", ylab = "Outcome", lty=1, col = "black", ylim=c(min(ubar-max_error)-1,max(ubar+max_error)+1))
lines(colnames(df), eopr_output$ubar, type="l", col = "blue", lty=2,lwd=1.5)
lines(colnames(df), eopr_output$ubar+max_error, type="l", col = "black", lty=3, lwd=1.3)
lines(colnames(df), eopr_output$ubar-max_error, type="l", col = "black", lty=3, lwd=1.3)
abline(v= pre_treat_end_col_name, col="gray", lwd=1.5, lty=4)
legend("topleft", legend = c("Treated unit", "EOpR", "worst case","Intervention"), ncol =1, col=c("black", "blue", "black", "gray"), lty=c(1,2,3,4))

dev.off()