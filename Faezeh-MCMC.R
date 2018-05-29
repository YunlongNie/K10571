# Ben's Function:

library(dplyr)
library(ggplot2)

gFuns= function(t,k, sd_noise=0){
  if (t<0|t>10) {
    #cat('Warning: give a t between 0 and 10\n') 
    value = NA
  } else {
    
    if (k==1|k==2|k==3){
      klevel = 3
      mlevel = 10
      # set the upper limit for t 
      xlimt = 10
      # generate a data.frame with slope and intercept for each k and m level 
      data_km = expand.grid(k0=1:klevel,m0=1:mlevel)
      
      #set.seed(1);data_km$a = rnorm(nrow(data_km));data_km$b= rnorm(nrow(data_km))
      data_km = structure(list(k0 = c(1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 1L, 
      2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L,3L, 1L, 2L, 3L),
                               
      m0 = c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L, 4L, 4L, 4L, 5L, 5L, 5L, 6L, 6L, 6L,
      7L, 7L, 7L, 8L, 8L, 8L, 9L, 9L,9L, 10L, 10L, 10L), 
      a = c(-0.626453810742332, 0.183643324222082,-0.835628612410047, 1.59528080213779, 0.329507771815361, -0.820468384118015, 
      0.487429052428485, 0.738324705129217, 0.575781351653492, -0.305388387156356, 
      1.51178116845085, 0.389843236411431, -0.621240580541804, -2.2146998871775, 
      1.12493091814311, -0.0449336090152309, -0.0161902630989461, 0.943836210685299, 
      0.821221195098089, 0.593901321217509, 0.918977371608218, 0.782136300731067, 
      0.0745649833651906, -1.98935169586337, 0.61982574789471, -0.0561287395290008, 
      -0.155795506705329, -1.47075238389927, -0.47815005510862, 0.417941560199702),
      b = c(1.35867955152904, -0.102787727342996, 0.387671611559369, 
     -0.0538050405829051, -1.37705955682861, -0.41499456329968, -0.394289953710349, 
      -0.0593133967111857, 1.10002537198388, 0.763175748457544, -0.164523596253587, 
      -0.253361680136508, 0.696963375404737, 0.556663198673657, -0.68875569454952, 
     -0.70749515696212, 0.36458196213683, 0.768532924515416, -0.112346212150228, 
      0.881107726454215, 0.398105880367068, -0.612026393250771, 0.341119691424425, 
      -1.12936309608079, 1.43302370170104, 1.98039989850586, -0.367221476466509, 
      -1.04413462631653, 0.569719627442413, -0.135054603880824)), 
     .Names = c("k0", "m0", "a", "b"), out.attrs = structure(list(dim = structure(c(3L, 
      10L), .Names = c("k0", "m0")), dimnames = structure(list(k0 = c("k0=1", 
     "k0=2", "k0=3"), m0 = c("m0= 1", "m0= 2", "m0= 3", "m0= 4", "m0= 5", 
     "m0= 6", "m0= 7", "m0= 8", "m0= 9", "m0=10")), .Names = c("k0","m0"))), 
     .Names = c("dim", "dimnames")), row.names = c(NA, -30L), class = "data.frame")
      
      temp = data_km%>%dplyr::filter(k0==k)
      value = temp$a*t+temp$b+rnorm(length(mlevel), sd=sd_noise)
    } else {
      #cat('Warning: give a k from 1, 2, 3\n')
      value = NA
    }
  }
  return(value)
}
# End of Ben's code
#############################################################################

t = 0.2
tmax = 10
k = 3
sigma_noise = 0.3

#yobs = gFuns(t=0.2,k=3,sd_noise=sigma_noise)

yobs=c(0.9410312,0.1413971,1.9356670,0.5450923,0.2567158,1.6777855,1.3023867,
       -0.8067481,0.3221048,0.6690190)

##############################################################################

#loglikelihood
log.liklihood.fun=function(yobs,k,t,sigma_noise,n=10){
  
log.likli=-n*log(sigma_noise)-n/2*log(2*pi)-(1/(2*(sigma_noise)^2))*sum((yobs-gFuns(t,k,sd_noise=0))^2)

  return(log.likli)
}

# Prior for t
w1=function(t){
  t_value=dunif(t, min = 0, max = tmax, log = TRUE)
  return(t_value)
} 

# Prior for k
w2=log(1/3)


#joint-log-posterior
log.postFunc=function(yobs,k,t,sigma_noise){
  
  log.post=log.liklihood.fun(yobs,k,t,sigma_noise,n=10)+ w1(t) + w2
  
  return(log.post)
}

#log-posterior of t
log.postFunc_t = function(yobs,k,t,sigma_noise){
  log.post_t=log.liklihood.fun(yobs,k,t,sigma_noise,n=10)+ w1(t)
   
return(log.post_t)
}


############################################################ 
# MCMC
MCMC_Fun = function(yobs,k,sigma_noise,iterations=2500,sigma_stepsize=0.2,t_start=1,burnIn=500){

  proposalfunction <- function(par,sigma_stepsize){
    prop = rnorm(1,mean = par, sd=sigma_stepsize)
    return(prop)
  }

chain = array(dim = c(iterations+1,1))
chain[1,] = t_start
for(i in 1: iterations){
 param=chain[i,]
 
 prop=proposalfunction(param,sigma_stepsize)
 
 alph=(log.postFunc_t(yobs,k,prop,sigma_noise)-log.postFunc_t(yobs,k,param,sigma_noise))
  if(is.na(alph)){alph=-Inf}
   if (log(runif(1)) < alph){

      param=prop
      }
  chain[i+1,] = param 
}

return(chain[burnIn:iterations])

}

tsample_1 = MCMC_Fun(yobs,k=1,sigma_noise)
tsample_2 = MCMC_Fun(yobs,k=2,sigma_noise)
tsample_3 = MCMC_Fun(yobs,k=3,sigma_noise)

acceptance_1 =1-mean(duplicated(tsample_1))
acceptance_1
acceptance_2 =1-mean(duplicated(tsample_2))
acceptance_2
acceptance_3 =1-mean(duplicated(tsample_3))
acceptance_3

hist(tsample_1,nclass=30, main="t-samples for k=1", xlab="True value of t = blue line" )
abline(v = mean(tsample_1),col="red")
abline(v = t, col="blue" )

plot(tsample_1, type = "l", xlab="True value = blue line" , main = "t-samples for k=1", )
abline(h=mean(tsample_1) , col="red")
abline(h = t, col="blue" )



hist(tsample_2,nclass=30 , main="t-samples for k=2", xlab="True value of t = blue line" )
abline(v = mean(tsample_2),col="red")
abline(v = t, col="blue" )

plot(tsample_2, type = "l", xlab="True value = blue line" , main = "t-samples for k=2", )
abline(h=mean(tsample_2) , col="red")
abline(h = t, col="blue" )



hist(tsample_3,nclass=30 , main="t-samples for k=3", xlab="True value of t = blue line" )
abline(v = mean(tsample_3),col="red")
abline(v = t, col="blue" )

plot(tsample_3, type = "l", xlab="True value = blue line" , main = "t-samples for k=3", )
abline(h=mean(tsample_3) , col="red")
abline(h = t, col="blue" )

################################################
#Luyao'codes

# Normalizing constant
lalpha= function(t){
  return(50) 
}

lpostk1 = function(t){
return(lalpha(t)+log.postFunc(yobs,k=1,t,sigma_noise))
}

lpostk2 = function(t){
  return(lalpha(t)+log.postFunc(yobs,k=2,t,sigma_noise))
}

lpostk3 = function(t){
  return(lalpha(t)+log.postFunc(yobs,k=3,t,sigma_noise))
}


E1Q2 = mean(exp(as.numeric(lapply(tsample_1,lpostk2))))
E2Q1 = mean(exp(as.numeric(lapply(tsample_2,lpostk1))))

E1Q3 = mean(exp(as.numeric(lapply(tsample_1,lpostk3))))
E3Q1 = mean(exp(as.numeric(lapply(tsample_3,lpostk1))))

E3Q2 = mean(exp(as.numeric(lapply(tsample_3,lpostk2))))
E2Q3 = mean(exp(as.numeric(lapply(tsample_2,lpostk3))))


unnorm_prob = c(1,E1Q2/E2Q1,E1Q3/E3Q1)
norm_prob = unnorm_prob/sum(unnorm_prob)
barplot(norm_prob,ylab="posterior prob.",names.arg=c("k=1", "k=2", "k=3"),
        main=substitute(paste(sigma,"=",noise),list(noise = 0.3)))

############################################ THE END ####################################





