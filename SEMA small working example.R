#################################
##small working example of SEMA, fitting a random-intercept model
#################################

source("SEMA functions small working example.R")

dataset<-gen.data(N=10000, j=1000, mu=10, tau=10, sigma=5) #note: tau and sigma are sd not var

result<-list()

for(i in 1:nrow(dataset))
{
	result<-fitSema(id=dataset$id[i], obs=dataset$obs[i], res=result,r=i)
}


###########################################
##result now contains:
result$global$n
result$global$J
result$global$T1
result$global$T2
result$global$T3

#dataframe: id, mean_y_j, hatmu_j, hatnu_j, mean_ysq, n_j
result$individual 

#vector with mu, tau^2, and sigma^2
result$parameters
##########################################