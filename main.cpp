#include "cpplapack_plus.h"


int main(){

  scatter();

  cout << gauss_rand() << endl;

  dcovector multi(4);

  multi(0)=1;
  multi(1)=2;
  multi(2)=3;
  multi(3)=4;

  cout << ProbSelect(multi) << endl;


  cout << gennor(1,1) << endl;


  dgematrix A(2,2); A.identity();
  cout << IWishartSampler(3,A) << endl; 
 cout << IWishartSampler(3,10*A) << endl; 
 cout << IWishartSampler(3,100*A) << endl;
 cout << IWishartSampler(3,1000*A) << endl;
 cout << IWishartSampler(3,10000*A) << endl;


  return 0;
}
