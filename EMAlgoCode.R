require(mvtnorm)
set.seed(1)

######################
####IMPORT DATASET####
######################
swissbank <- read.csv("~/swiss-bank.dat", sep="")


########################
####RANDOM PARTITION####
########################
###training set/testing set
rsamp=array(dim=170)
n=1
while(n<=170){
  x=floor(runif(1,min=0.005,max=1.005)*200)
  if(! x%in%rsamp){
    rsamp[n]=x
    n=n+1
  }
}

train.data<-swissbank[rsamp,-1]
test.data<-swissbank[-rsamp,]

###population 1/population 2
pop1=array(dim=85)
n=1
while(n<=85){
  x=floor(runif(1,min=0.005,max=1.005)*170)
  if(! x%in%pop1){
    pop1[n]=x
    n=n+1
  }
}

population1<-train.data[pop1,]
population2<-train.data[-pop1,]




##########################################
####INITIALISATION OF PARAMETER VALUES####
##########################################
p1=.5
mu1=as.numeric(sapply(population1,mean))
mu2=as.numeric(sapply(population2,mean))
sgm1=matrix(cov(as.matrix(population1),y=as.matrix(population1)),6,6)
sgm2=matrix(cov(as.matrix(population2),y=as.matrix(population2)),6,6)

#######################################
####COMPUTING RESPONSIBILITY VECTOR####
#######################################
responsibility.vector<-function(train.data){
  n=length(as.matrix(train.data[,1]))
  r<<-array(dim=n)
  for(i in 1:n){
    f1=dmvnorm(as.numeric(train.data[i,]),mean = mu1,sigma = sgm1)
    f2=dmvnorm(as.numeric(train.data[i,]),mean = mu2,sigma = sgm2)
    r[i]<<-p1*f1/(p1*f1+(1-p1)*f2)
  }
}

####################################
####UPDATION OF PARAMETER VALUES####
####################################
parameters<-function(train.data){
  mu1<<-rep(0,6)
  mu2<<-rep(0,6)
  s=1-r
  n=length(r)
  for(i in 1:n){
    mu1<<-mu1+r[i]*as.numeric(train.data[i,])/sum(r)
    mu2<<-mu2+s[i]*as.numeric(train.data[i,])/sum(s)
  }
  sgm1<<-matrix(rep(0,36),nrow=6,ncol=6)
  sgm2<<-matrix(rep(0,36),nrow=6,ncol=6)
  for(i in 1:n){
    sgm1<<-sgm1+r[i]*(as.numeric(train.data[i,])-mu1)%*%t((as.numeric(train.data[i,])-mu1))/sum(r)
    sgm2<<-sgm2+s[i]*(as.numeric(train.data[i,])-mu2)%*%t((as.numeric(train.data[i,])-mu2))/sum(s)
  }
  p1<<-sum(r)/n
}

######################################
####PSEUDO LOG LIKELIHOOD FUNCTION####
######################################
pseudo.log.likelihood<-function(train.data){
  s=1-r
  q1=1-p1
  t1=0
  t2=0
  for (i in 1:length(as.matrix(train.data[,1]))) {
    lnf1=dmvnorm(as.numeric(train.data[i,]),mean=mu1,sigma=sgm1,log=T)
    lnf2=dmvnorm(as.numeric(train.data[i,]),mean=mu2,sigma=sgm2,log=T)
    t1=t1+r[i]*lnf1
    t2=t2+s[i]*lnf2
  }
  pllf<<-log(p1)*sum(r)+log(q1)*sum(s)+t1+t2
  print(pllf)
}

####################
####EM ALGORITHM####
####################
responsibility.vector(train.data)
pseudo.log.likelihood(train.data)
pllf.array<-array()
pllf.array[1]<-pllf
parameters(train.data)
pseudo.log.likelihood(train.data)
pllf.array[2]<-pllf
k=2
while(TRUE){
  pllf2=pllf
  k=k+1
  responsibility.vector(train.data)
  parameters(train.data)
  pseudo.log.likelihood(train.data)
  pllf.array[k]<-pllf
  if(abs(pllf2-pllf)<0.01)
    break()
}

plot(pllf.array,type = "l",lwd=3,xlab = "Number of iterations", ylab = "Pseudo log-likelihood")

######################
####CLASSIFICATION####
######################
r.t<-array()
for(i in 1:length(test.data[,1])){
  f1=dmvnorm(as.numeric(test.data[i,2:7]),mean = mu1,sigma = sgm1)
  f2=dmvnorm(as.numeric(test.data[i,2:7]),mean = mu2,sigma = sgm2)
  r.t[i]<-p1*f1/(p1*f1+(1-p1)*f2)
}
names(test.data)[1]<-c("status")

s1 = sum(abs(test.data$status-round(r.t)))
s2 = sum(abs(test.data$status-(1-round(r.t))))
if(s1<s2){
  result.vector <- cbind("True class"=test.data$status, "Assigned class"=round(r.t))
  }else{
  result.vector <- cbind("True class"=test.data$status, "Assigned class"=(1-round(r.t)))
}
View(result.vector)
