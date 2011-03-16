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






  return 0;
}
