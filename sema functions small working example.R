############################
##Streaming EM Approximation functions
############################

require("plyr")


############################
##generate data
############################
gen.data<-function(N=10000, j=1000, mu=10, tau=10, sigma=5)
{
	true<-rnorm(j, mu, tau)
	dataset<-data.frame(id=numeric(0), obs=numeric(0), true=numeric(0))
	for(d in 1:j)
	{
	dataset<-rbind(dataset, data.frame(id=d, obs=rnorm((N/j), true[d], sigma),true[d]))
	}
	new.order<-sample(1:N)
	dataset<-dataset[new.order,]
	return(dataset)
}
#########################
##update average
#########################
update_average<-function(oldmean=1, obs=5, n=2)
{	
	return(oldmean+(obs-oldmean)/n)
}
########################
##reliablity 
########################
estimate_rho<-function(tausq=5, sigmasq=10, nj=2)
{
	return(tausq/(tausq+(sigmasq/nj)))
}
#########################
##update hatmu_j
########################
estimate_mu_j<-function(rho=1, average=10, mu=10)
{
	return(rho*average+((1-rho)*mu))
}
#######################
##update hatnu_j
#######################
estimate_nu_j<-function(rho=1,tausq=5)
{
	return(tausq*(1-rho))
}
######################
##update CDSS gamma
######################
update_T1<-function(T1=10, mu_j=5)
{
	T1	<- T1+mu_j
	return(sum(T1))
}
######################
##update CDSS tau^2
######################
update_T2<-function(T2=10, mu_j=5, nu_j=1)
{
	T2	<- T2+(mu_j**2+nu_j)
	return(sum(T2))
}
######################
##update CDSS sigma^2
######################
update_T3<-function(T3=10, mu_j=5, nu_j=1, ysq=5, y=2, nj=2)
{
	T3	<- T3+(ysq - 2*y*mu_j + mu_j**2 + nu_j)*nj
	return(sum(T3))
}
############################################################################
##SEMA for random intercept models, without adaptations
############################################################################

fitSema<-function(id=1, obs=10, res=list(),r=1)
{
	if(r==1){res				<-list(global = list(n=0, J=0,T1=0,T2=0,T3=0),
								individual=data.frame(id=numeric(0), mean_y_j=numeric(0), 
								hatmu_j = numeric(0),n_j=numeric(0), hatnu_j=numeric(0),mean_ysq=numeric(0)),
								parameters=c(mu=0,tausq=1,sigmasq=1))}
	T1						<-res$global$T1			
	T2						<-res$global$T2		
	T3						<-res$global$T3		
	
	tausq						<-as.numeric(res$parameters[2])
	sigmasq					<-as.numeric(res$parameters[3])
	if (res$global$n==0) 			{mu_hat<-obs}			#set mu	
	else 						{mu_hat <- as.numeric(res$parameters[1])}

	res$global$n				<- res$global$n+1
	if(sum(res[[2]]$id==id)==0)		#add new person to storage
		{
			res$global$J		<- res$global$J+1
			res[[2]]			<- rbind(res[[2]], data.frame(id=id, mean_y_j=obs, 
								hatmu_j = 0, n_j=0, hatnu_j=0,mean_ysq=obs**2))						
		}		

	id_row					<-res[[2]][res[[2]]$id==id,]	#select working row
	
	#############################
	##subtract previous contributions from the CDSS
	#############################	
	T1						<-T1-id_row$hatmu_j		
	T2						<-T2-(id_row$hatmu_j**2+id_row$hatnu_j)
	T3						<-T3-(id_row$mean_ysq-2*id_row$mean_y_j*id_row$hatmu_j+id_row$hatmu_j**2+id_row$hatnu_j)*id_row$n_j
	#######################
	##update individual parameters
	#######################	
	id_row$n_j 					<- id_row$n_j+1		
	id_row$mean_y_j				<- update_average(oldmean=id_row$mean_y_j, obs=obs, n=id_row$n_j)
	id_row$mean_ysq				<- update_average(oldmean=id_row$mean_ysq, obs=obs**2, n=id_row$n_j)

	#######################
	##individual parameters
	#######################
	rho						<- estimate_rho(tausq=tausq, sigmasq=sigmasq, nj=id_row$n_j)
	id_row$hatmu_j				<- estimate_mu_j(rho=rho, average=id_row$mean_y_j, mu=mu_hat)
	id_row$hatnu_j				<- estimate_nu_j(rho=rho,tausq=tausq)
		
	#######################
	##model parameters
	#######################
	res$global$T1				<- update_T1(T1=T1, mu_j=id_row$hatmu_j)
	res$global$T2				<- update_T2(T2=T2, mu_j=id_row$hatmu_j, nu_j=id_row$hatnu_j)
	res$global$T3				<- update_T3(T3=T3, mu_j=id_row$hatmu_j, nu_j=id_row$hatnu_j, 
								ysq=id_row$mean_ysq, y=id_row$mean_y_j,  nj=id_row$n_j)


	res[[2]][res[[2]]$id==id,]		<- id_row	#store updated information

	res$parameters[1]				<- res$global$T1/res$global$J
	res$parameters[2]				<- res$global$T2/res$global$J-res$parameters[1]**2
	res$parameters[3]				<- res$global$T3/res$global$n
	
	return(res)

}
