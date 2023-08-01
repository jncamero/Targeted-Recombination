library(Rglpk)
library(Rcpp)

#Number of marker loci to simulate
no.markers=100

#Creating a random F1, with no.markers rows and 2 columns
pars=array(rbinom(prob=0.5,no.markers*2,1),c(no.markers,2))
#Sample random vector of marker effects
meff=rnorm(no.markers,mean=0,sd=.05)

marker_index=1:nrow(pars)

lgend1=seq(1,90,10)
lgend2=seq(10,100,10)
#############################################################
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

#
data=pars*meff

#Objective function
OF=c()
for(i in 1:nrow(data)){
	addon=c(data[i,1],data[i,2],data[i,2],data[i,1])
			OF=c(OF,addon)
	}

#Running integer program
xa=rep("B",ncol(const))
system.time({result=Rglpk_solve_LP(OF, const, const.dir, const.rhs, types = xa, control = list("verbose" =
FALSE), max = TRUE)})

#x: The decision variable matrix; non-recombinant, recombinant1, non-recombinant, recombinant2
x=c()
res=result$solution
for(i in 1:nrow(data)){
x<-rbind(x,res[1:4])
res=res[-c(1:4)]
}

a<-rowSums(x[,c(1,4)])
b<-rowSums(x[,c(2,3)])

#fin: Matrix indicating optimal haplotypes 
fin=cbind(a,b)

#Organizing optimal recombination points
#return the Optimal recombination points
cnt=0
recpts=c()
for(i in c(2,4)){
cnt=cnt+sum(which(x[,i]==1)%in%lgend1)
#print(which(x[,i]==1)%in%lgend1)
recpts=c(recpts,which(x[,i]==1))
}
targrec=recpts[which(!recpts%in%lgend1)]

#recpts: matrix of optimal recombination points, subsetted from map.tog
recpts=map.tog.markers[targrec,]
recpts=recpts[order(recpts$Chrom,recpts$marker),]
