setwd("~/Documents/Robbins Rotation")
library(runjags)
library(MASS)
library(kinship2)
library(ggplot2)
nrep = 3
seed = 13
set.seed(seed)

replicate = 1:nrep

sample.ped = rbind(c(1,1,0,0,1), c(1,2,0,0,2), c(1,3,0,0,2), c(1,4,1,2,1), 
                   c(1,5,1,3,2), c(1,6,4,5,1), c(1,7,4,5,2), c(1,8,6,7,1))
colnames(sample.ped) = c("ped","id", "father", "mother", "sex")
sample.ped = as.data.frame(sample.ped)
sample.ped
pedAll <- pedigree(id=sample.ped$id, dadid=sample.ped$father, 
                   momid=sample.ped$mother, sex = sample.ped$sex)
plot.pedigree(pedAll, cex = 0.5)
GRM = kinship(pedAll)#[3:7,3:7]
ngeno = dim(GRM)[1]
genotype = 1:ngeno
# variance defines genotype differences in parameter values
# note here that c is implicitly 0
R = diag(c(100,10,10))
df = 4
set.seed(seed)
sigma = solve(rWishart(n = 1, df = df, Sigma = R)[,,1])
sigmaGRM = sigma %x% GRM
mean = c(-0.5,0.5,20)
mvmean = c(rep(mean[1],ngeno), rep(mean[2],ngeno), rep(mean[3],ngeno))
set.seed(seed)
params = mvrnorm(1, mvmean, sigmaGRM)

p = matrix(nrow = ngeno, ncol = 3)
for (i in 1:3) {
  p[,i] = params[(((i-1) * ngeno) + 1) : (i * ngeno)]
}

paramnames = c("b", "d", "e")
colnames(p) = paramnames
data = as.data.frame(cbind(genotype, p))

alldata = c()
time = seq(0, 40,5)
var_y = 0.001
sigma_y = diag(var_y, length(time)) 
for (i in 1:ngeno) {
  subdata = data[data$genotype == i,] # pick the genotype-specific b,d,e
  mu_arr = (subdata$d)/(1 + exp(subdata$b*(time-subdata$e))) # growth curve as mean
  set.seed(seed)
  y = mvrnorm(nrep, mu_arr, sigma_y) # n replicates
  for (j in 1:nrep) {
    alldata = rbind(alldata, cbind(rep(i, length(time)), rep(j, length(time)),time, y[j,]))
  }
}

alldata = as.data.frame(alldata)
colnames(alldata) = c("genotype", "rep", "time", "y")
alldata$rep = as.character(alldata$rep)

## ----echo = FALSE-------------------------------------------------------------
print("covariance matrix (sigma) for growth curve parameters")
sigma
for (i in 1:3){
  print(paste("realized variance for",paramnames[i], ":", var(p[,i])))
}


# library(ggplot2)
plot = ggplot() + geom_point(data = alldata, aes(x = time, y = y, color = rep))
plot = plot + facet_wrap(~genotype, nrow = 2)
plot

nrT = length(time)

for (geno in genotype) {
  ally = c()
  for (rep in replicate) {
    ally = c(ally, list(alldata[alldata$genotype == geno & alldata$rep == rep,]$y))
  }
  ally = do.call(rbind, ally)
  if(geno == 1) {
    res = array(ally, dim = c(nrep, length(time), 1))
  } else {
    res = array(c(res, ally), dim = c(nrep, length(time), dim(res)[3] + 1))
  }
}


## ----JAGS---------------------------------------------------------------------------
R = diag(c(100,100,100)) # uninformative scaled matrix
df = 4

jagsData <- list("Y"=res,"N"=ngeno,"nrT"=nrT,"time"=time, "nRep" = nrep, "G" = GRM, "R" = R, "df" = df, "d" = 3)

model_string <- "model {
  # dimensions of Y matrix = nRep x nrT x N
  for (i in 1:N) { # loop over genotypes
    for (t in 1:nrT) { # loop over time points
      for (j in 1:nRep) {
         Y[j, t, i] ~ dnorm(params[2*i] / (1 + exp(params[1*i] * (time[t] - params[3*i]))), 1/pow(sd,2))    
      }
    }
  }
  sd ~ dunif(0, 100)

  for (i in 1:d) {
    for (j in 1:d) {
      for (k in 1:N) {
        for (l in 1:N) {
          TG[(i-1)*N+k,(j-1)*N+l] <- TAU[i,j]*G[k,l] # Kronecker product
        }
      }
    }
  }
  C <- cholt(TG) # using new function within JAGS
  
  params <- C[,] %*% z[] + mu[]
  TAU ~ dwish(R, df)
  for (i in 1:(N*d)) { z[i] ~ dnorm(0,1) }
  for (j in 1:(N*d)) {
    mu[j] ~ dnorm(0,0.01)
  }
  
  VCOV <- inverse(TAU)
  vars[1] <- VCOV[1,1]; vars[2] <- VCOV[2,2]; vars[3] <- VCOV[3,3]
  cov[1] <- VCOV[1,2]; cov[2] <- VCOV[1,3]; cov[3] <- VCOV[2,3]
}"


## ----run----------------------------------------------------------------------------
nadapt = 100000
nupdate = 100000
nimplement = 50000
parameters = c("params", "sd", "vars", "cov", "mu", "C")
samples <- run.jags(model_string, data = jagsData, monitor=parameters, 
                    method="rjags", modules="runjags", n.chains = 1,
                    adapt = nadapt, burnin = nupdate, sample = nimplement)
s = as.data.frame(as.matrix(as.mcmc.list(samples)))


# model <- jags.model(textConnection(model_string),
#                     data=jagsData, n.chains=1, n.adapt = nadapt, inits = list(.RNG.name = "base::Wichmann-Hill",.RNG.seed = seed))
# 
# update(model, n.iter=nupdate)
# samples <- coda.samples(model, variable.names=parameters, n.iter=nimplement)
# s = as.data.frame(as.matrix(samples))
# head(s[,1:5])


## ----plot param estimates, echo = FALSE---------------------------------------------
total = dim(s)[1]; thin = rep(total/1000, 3)
print(paste("thinning to use every", total/1000))
p = c("b", "d", "e")
print(ngeno)
for (i in 1:length(p)) {
  param_index = grep("param", colnames(s))
  param = colnames(s)[param_index]
  index = param_index[grep(paste0(i,"\\]"),param)]
  subset = as.data.frame(s[,index])
  
  keep = seq(1, total,thin[i])
  subset = as.data.frame(subset[keep,])
  # print(dim(subset))
  subset = as.data.frame(rowMeans(subset))
  
  print(dim(subset))
  colnames(subset) = c("mean")
  subset$iteration = 1:(total/thin[i])
  subset$param = p[i]
  if(i == 1) {
    graphdata = subset
  }
  else {
    graphdata = rbind(graphdata, subset)
  }
}
graphdata$iteration = as.numeric(graphdata$iteration)

plot = ggplot(graphdata) +
  geom_histogram(aes(mean), bins = 20,fill="lightblue", colour="black") +
  facet_wrap(~param, scales = "free") +
  ggtitle(paste("Adapt =", nadapt, "Update =", nupdate, "Run =", nimplement))
ggsave("param_hist_covar_param.png", plot, dpi = 300, height = 5, width = 7, unit = 'in')

plot = ggplot(graphdata) +
  geom_point(aes(iteration, mean)) +
  facet_wrap(~param, scales = "free", nrow = 3) +
  ggtitle(paste("Adapt =", nadapt, "Update =", nupdate, "Run =", nimplement))
ggsave("param_trace_covar_param.png", plot, dpi = 300, height = 5, width = 5, unit = 'in')

## ----plot var estimates, echo = FALSE-----------------------------------------------
total = dim(s)[1]; thin = rep(total/1000, 3)
print(paste("thinning to use every", total/1000))
p = c("var(b)", "var(d)", "var(e)")
index = grep("var", colnames(s))
for (i in 1:length(p)) {
  subset = as.data.frame(s[,index[i]])
  # print(colnames(s)[index[i]])
  
  keep = seq(1, total,thin[i])
  subset = as.data.frame(subset[keep, ])
  colnames(subset) = c("var")
  
  subset$iteration = row.names(subset)
  subset$param = p[i]
  if(i == 1) {
    graphdata = subset
  }
  else {
    graphdata = rbind(graphdata, subset)
  }
}
graphdata$iteration = as.numeric(graphdata$iteration)
plot = ggplot(graphdata) +
  geom_histogram(aes(var), bins = 20,fill="lightblue", colour="black") +
  facet_wrap(~param, scales = "free") +
  ggtitle(paste("Adapt =", nadapt, "Update =", nupdate, "Run =", nimplement))
ggsave("var_hist_covar_param.png", plot, dpi = 300, height = 5, width = 7, unit = 'in')

plot = ggplot(graphdata) +
  geom_point(aes(iteration, var)) +
  facet_wrap(~param, scales = "free", nrow = 3) +
  ggtitle(paste("Adapt =", nadapt, "Update =", nupdate, "Run =", nimplement))
ggsave("var_trace_covar_param.png", plot, dpi = 300, height = 5, width = 5, unit = 'in')


## ----plot sd estimates, echo = FALSE------------------------------------------------
total = dim(s)[1]; thin = total/1000
print(paste("thinning to use every", total/1000))
sd = as.data.frame(s$sd)
keep = seq(1, total,thin)
sd = as.data.frame(sd[keep, 1])
colnames(sd) = "sd"
sd$iteration = row.names(sd)
sd$iteration = as.numeric(sd$iteration)
plot = ggplot(sd) +
  geom_point(aes(iteration, sd)) +
  ggtitle(paste("Adapt =", nadapt, "Update =", nupdate, "Run =", nimplement))
ggsave("sd_trace_covar_param.png", plot, dpi = 300, height = 5, width = 5, unit = 'in')

plot = ggplot(sd) +
  geom_histogram(aes(sd), bins = 20,fill="lightblue", colour="black") +
  ggtitle(paste("Adapt =", nadapt, "Update =", nupdate, "Run =", nimplement))
ggsave("sd_hist_covar_param.png", plot, dpi = 300, height = 5, width = 5, unit = 'in')



