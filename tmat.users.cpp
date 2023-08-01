#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
// n of progeny to create
//rr is the recombination rate between loci
//a is the first parent (or homologous chromosome set) of parent 1
//b is the second parent(or homologous chromosome set) of parent 2

// [[Rcpp::export]]
List conmat(NumericVector one, NumericVector two, int m, NumericVector lgend1,
 int norec) {
   //number of loci
   int w;

   //Output dimensions:
   //A: 2 * Number of loci - 2, recombination or not vector, inflow equals outflow
   //B: Number of loci, only 1 allele inherited at a locus
   //C: 1, Recombination constraint
   NumericMatrix output((m-1)*2 + m + 1,m*4);

   //Sign
   CharacterVector sign(output.nrow());
   //LHS
   NumericVector lhs(output.nrow());

   //Numnber of targetted Recombination constraint
   NumericVector oo=lgend1*4-4;
   NumericVector tt=oo+3;

   //Constraint set A: IN and OUT constraint, (no.loci-1)*2 rows
   int start=0;
   for(int i=0;i<(m-1);i++){
      w=0;
      for(int j=start;j<start+8;j++){
      output(i*2,j)=one(w);
      output(i*2+1,j)=two(w);
      w=w+1;
   }
      lhs(i)=0;
      start=start+4;
}
for(int i=0;i<(m-1)*2;i++){
   sign(i)="==";
}
    // Second set of constraints
   //  No.loci constraints
  //   Starting after the last constraint from constraint set 1
int ii=(m-1)*2;

//Constraint set B: Only 1 allele inherited per locus: no.loci rows //////////////////////////////
start=0;
for(int i=ii;i<ii+m;i++){
      for(int j=start;j<start+4;j++){
      output(i,j)=1;
      }
      sign(i)="==";
      start=start+4;
}


//Constraint C: number of recombinations allowed constraint: 1 row ////////////////////////////////
for(int i=1; i<output.ncol(); i += 2){
   output(output.nrow()-1,i)=1;
}

for(int i=0;i<oo.size();i++){
      for(int j=oo(i);j<=tt(i);j++){
         output(output.nrow()-1,j)=0;
   }
}

//Adding norec sign to sign vector

sign(output.nrow()-1)="<=";

////////////////////////////////////////////////
//LHS
for(int i=0;i<ii;i++){
   lhs(i)=0;
}

for(int i=ii;i<lhs.length();i++){
   lhs(i)=1;
}
   //Number or allowed recombinations
   lhs(lhs.length()-1)=norec;

List L = List::create(output, sign,lhs);

return L;
}