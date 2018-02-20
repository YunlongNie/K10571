#  simulation one 
# all the g's are lines 

# set the levels for k and m
klevel = 3
mlevel = 10
# set the upper limit for t 
xlimt = 10
# generate a data.frame with slope and intercept for each k and m level 
data_km = expand.grid(k=1:klevel,m=1:mlevel)
data_km$a = rnorm(nrow(data_km))
data_km$b= rnorm(nrow(data_km))


# I create a data frame in order to plot all the g functions for each m 

flim = function(data){
t = seq(0,xlimt,len=1000) # eval each line at 1000 equal spaced time points 
data.frame(value = data$a*t+data$b,t=t)
}
library(ggplot2)
limd = data_km%>%dplyr::group_by(m,k)%>%do(flim(.))%>%data.frame
ggplot(limd)+geom_line(aes(x=t,y=value,color=factor(k)))+facet_wrap(~m, scales="free")

## set the true value for k and t 
ttrue = 5
ktrue = 1
sigma_noise = 1

# fytrue is a function with input as the true k and t and output the true value of y
fytrue= function(ttrue,ktrue){
temp = data_km%>%dplyr::filter(k==ktrue)
temp$a*ttrue+temp$b
}

# obtain the observed value of y with some noise 
yobs = fytrue(ttrue,ktrue)+rnorm(length(mlevel), sd=sigma_noise)

# loss function with input as any t and k and output as the square differences between the observed and fitted value (a length of m vector)

loss = function(y,t,k){
	
	(y-fytrue(t,k))^2

}

##  losst give the mean square differences or errors
losst = function(t,k) {
	mean((loss(yobs,t,k))^2)
}

## we can evaluate the loss function surface

settings = expand.grid(t=seq(0,10,len=20), k=1:3)
res = mapply(t=settings$t,k=settings$k, losst)
settings$res = res

head(settings)
ggplot(settings,aes(x=t,group=factor(k),color=factor(k),y=res))+geom_line()

#  given a different value of t 
ttrue = 0.1
ktrue = 1
ytrue = fytrue(ttrue,ktrue)+rnorm(length(mlevel))
res = mapply(t=settings$t,k=settings$k, losst)
settings$res = res
ggplot(settings,aes(x=t,group=factor(k),color=factor(k),y=res))+geom_line()
ggplot(settings%>%dplyr::filter(t<2),aes(x=t,group=factor(k),color=factor(k),y=res))+geom_line()


