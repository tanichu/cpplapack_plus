#include "cpplapack_plus.h"


#define N 10

int main(){

  scatter();

  dgematrix S(3,3);
  S.identity();
  dcovector m(3),x(3);
  m.zero();
  x(0)=0;x(1)=1;x(2)=2;
  cout <<  Cal_MultiNormLikely_with_precision(x,m,i(S)) << endl;
  cout <<  Cal_MultiNormLikely(x,m,S) << endl;




  getchar();

  dgematrix Sig(2,2); Sig.read("Sig.init.dat");
  dcovector Mu(2); Mu.zero();

  cout << Sig << endl;


  dgematrix X(N,2);

  for(int i=0;i<N;i++){
    vec_set(X,i,t(MultiGaussSampler(Mu,Sig)));
  }

  cout << X << endl;

  dgematrix A = t(X)*X;

  cout << (1.0/double(N)) * A << endl;

  cout << IWishartSampler(N,A) << endl;


  return 0;
}
