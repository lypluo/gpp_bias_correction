f_hardening <- function(temp, par){
  #changing the function according to the logistic function in weikipeidia:
  #https://en.wikipedia.org/wiki/Logistic_function
  xx_ori <- temp # * ppfd
  xx <- (-1)*par["b"] * (xx_ori - par["a"])
  yy <- 1 / (1 + exp(xx))
  plot(xx_ori,yy,xlab="Ta",ylab="par_stress")
  return(yy)
}

##
temp<-seq(-40,70,0.5)

#
par1<-c("a"=5,"b"=0.2)
f_hardening(temp,par1)

par1<-c("a"=5,"b"=0.1)
f_hardening(temp,par1)

par1<-c("a"=-5,"b"=0.1)
f_hardening(temp,par1)
