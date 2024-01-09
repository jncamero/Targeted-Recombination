library(Rglpk)
library(Rcpp)
sourceCpp("tmat.users.cpp")

#Number of allowed recombinations
norec=2

#Number of marker loci to simulate, on a single linkage group
no.markers=10

#Creating a random F1, with "no.markers" rows and 2 columns
pars=array(rbinom(prob=0.5,no.markers*2,1),c(no.markers,2))

#Sample random vector of marker effects from a normal distribution with mean 0 and sd of 0.05
meff=round(rnorm(no.markers,mean=0,sd=.05),2)

#
lgend1=seq(1)
lgend2=seq(10)
###########################
#Building constraint matrix
one=array(0,8)
one[c(1,4)]=1
one[4+c(1,2)]=-1

two=array(0,8)
two[c(2,3)]=1
two[4+c(3,4)]=-1

tmat=conmat(one,two,nrow(pars),lgend1,norec)

const=tmat[[1]]
const.dir=tmat[[2]]
const.rhs=tmat[[3]]

# Input data into integer program solver is phased genotype * marker effects
data=pars*meff

#Building objective function
OF=c()
for(i in 1:nrow(data)){
	addon=c(data[i,1],data[i,2],data[i,2],data[i,1])
			OF=c(OF,addon)
	}

#Running integer program
xa=rep("B",ncol(const))
system.time({result=Rglpk_solve_LP(OF, const, const.dir, const.rhs, types = xa, control = list("verbose" =
FALSE), max = TRUE)})

#x: The decision variable matrix; non-recombinant (p1-p1), recombinant1 (p1-p2), non-recombinant (p2-p2), recombinant2 (p2-p1)
x=c()
res=result$solution
for(i in 1:nrow(data)){
x<-rbind(x,res[1:4])
res=res[-c(1:4)]
}

a<-rowSums(x[,c(1,4)])
b<-rowSums(x[,c(2,3)])

#fin: Matrix indicating optimal haplotypes that maximize genetic value | 2 allowed recombinations
fin=cbind(a,b)

#Organizing optimal recombination points
cnt=0
recpts=c()
for(i in c(2,4)){
cnt=cnt+sum(which(x[,i]==1)%in%lgend1)
#print(which(x[,i]==1)%in%lgend1)
recpts=c(recpts,which(x[,i]==1))
}

#indices of optimal recombination points
targrec=recpts[which(!recpts%in%lgend1)]

