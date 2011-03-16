#include "cpplapack_plus.h"
#include <algorithm>

using namespace std;
std::string itos(int n)
{
	std::stringstream str_stream;
	str_stream << n;
	return str_stream.str();
}

std::string ftos(double n)
{
	std::stringstream str_stream;
	str_stream << n;
	return str_stream.str();
}

string ScanTillRet(string str){
	return str.substr(0,str.find("\n"));
}
string ScanTillSpace(string str){
	return str.substr(0,str.find(" "));
}
string ScanTillTab(string str){
	return str.substr(0,str.find("\t"));
}
string ScanTillComma(string str){
	return str.substr(0,str.find(","));
}


string TrimLeft( string arg ) {
	int n = int(arg.find_first_not_of("\n"));
	int t = int(arg.find_first_not_of("\t"));
	int s = int(arg.find_first_not_of(" "));

	for(;(abs(n-s)+abs(n-t)+abs(t-s))!=0;){
	
		arg = arg.substr(min(n,min(t,s))+1,arg.length());
	n = int(arg.find_first_not_of("\n"));
	t = int(arg.find_first_not_of("\t"));
	s = int(arg.find_first_not_of(" "));
	}

	return arg;
}

string TrimRight( string arg ) {

	int n = int(arg.find_last_not_of("\n"));
	int t = int(arg.find_last_not_of("\t"));
	int s = int(arg.find_last_not_of(" "));
	
	while((abs(n-s)+abs(n-t)+abs(t-s))!=0){
	
		arg = arg.substr(0,min(n,min(t,s))+1);
		
		n = int(arg.find_last_not_of("\n"));
		t = int(arg.find_last_not_of("\t"));
		s = int(arg.find_last_not_of(" "));
	}


    return arg;
}



string Trim( string arg){
	return TrimRight(TrimLeft(arg));
}

string DelTillRet(string str){
	return str.substr(str.find_first_of("\n")+1,str.length());
}
string DelTillSpace(string str){
	return str.substr(str.find_first_of(" ")+1,str.length());
}
string DelTillTab(string str){
	return str.substr(str.find_first_of("\t")+1,str.length());
}

string DelTillComma(string str){
	return str.substr(str.find_first_of(",")+1,str.length());
}


vector<string> split(const string& str, const string& delimiter) {
    // delimiter(2 文字以上も可) を空白に置換
    std::string item(str);    
    for(unsigned pos = item.find(delimiter); pos != string::npos; pos = item.find(delimiter, pos)) {
        item.replace(pos, delimiter.size(), " ");
    }
    // 分解
    stringstream buf(item);

    // 読み取り
    std::vector<string> result;
    while(buf >> item) {
        result.push_back(item);
    }

    return result;
}


//csvを読み取る．
vector< vector<string> > csv_reader(const char *filename){
	char ch[10000];
	string str;
	ifstream fin(filename);
	vector< vector<string> > csv;

	while(fin.getline(ch,10000)){//データの読み込み
		str = ch;
		vector<string> strs;
		strs = split(str,",");
		csv.push_back(strs);
	}
	fin.close();

	return csv;
}

//ssvを読み取る．space separated value
vector< vector<string> > ssv_reader(const char *filename){
	char ch[10000];
	string str;
	ifstream fin(filename);
	vector< vector<string> > csv;

	while(fin.getline(ch,10000)){//データの読み込み
		str = ch;
		vector<string> strs;
		strs = split(str," ");
		csv.push_back(strs);
	}
	fin.close();

	return csv;
}	

	
///元suplapack

void scatter(){
      srand( (unsigned)time( NULL ) );
}


double gauss_rand(){
        
        double x1 = double(rand()+1.0)/(RAND_MAX +1.0);
        double x2 = double(rand()+1.0)/(RAND_MAX +1.0);

        double y = sqrt(-2*log(x1))*cos(2.0*M_PI*x2) ;
        return y;
}

double gauss_rand(double ave,double sig){
        return ave + sig*gauss_rand();
}
dcovector gauss_rand(dcovector ave,double sig){

        dcovector out(ave);
        int i;

        for(i=0;i<ave.l;i++){
                out(i) =gauss_rand(ave(i),sig) ;
        }
        return out;
}

dgematrix gauss_rand(dgematrix ave,double sig){

        dgematrix out(ave);
        int i,j;

        for(i=0;i<ave.m;i++){
        for(j=0;j<ave.n;j++){
                out(i,j) =gauss_rand(ave(i,j),sig) ;
		}}
        return out;
}


dcovector gauss_rand(dcovector ave,dcovector sig){

        dcovector out(ave);
        int i;

		for(i=0;i<ave.l;i++){
                out(i) =gauss_rand(ave(i),sig(i)) ;
        }
        return out;
}


double LGF(dcovector x,dcovector m,dcovector s){

        double log_p = 0;
        for(int i=0;i<x.l;i++){
			log_p += -0.5*(x(i)-m(i))*(x(i)-m(i))/(s(i)*s(i)) - log(sqrt(2*M_PI)*s(i)) ;
        }
        return log_p;
}

dcovector ng_calculation(dcovector *m,dcovector *s,int N,dcovector x){
//normalized gaussian
                dcovector lg(N),lg_(N),ng(N);

                double maxlg=0;
                double sumg_=0;

                for(int j=0;j<N;j++){
                       lg(j) =LGF(x,m[j],s[j]);
                       if((maxlg<lg(j))||(j==0)){maxlg=lg(j);}
                }
                for(int j=0;j<N;j++){
                        lg_(j) = lg(j)-maxlg;
                }
                sumg_=0;
                for(int j=0;j<N;j++){
                        sumg_ += exp(lg_(j));
                }
                for(int j=0;j<N;j++){
                        ng(j) = exp(lg_(j))/sumg_;
                }

                return ng;
}

dcovector rbf_calculation(dcovector *m,dcovector *s,int N,dcovector x){
//radial basis function gaussian
                dcovector lg(N),rbf(N);
               //argmaxlg=0;
                for(int j=0;j<N;j++){
                       lg(j) =LGF(x,m[j],s[j]);
                       rbf(j) = exp(lg(j));
                 }
                return rbf;
}

void  each_abs(dcovector &v){
        int i;
        for(i=0;i<v.l;i++){
                v(i) =fabs(v(i));
        }
}

void  each_square(dcovector &v){
        int i;
        for(i=0;i<v.l;i++){
                v(i) =v(i)*v(i);
        }
}

double frobnorm(dgematrix A){
    double norm = 0;     
	int i,j;
    for(i=0;i<A.m;i++){
        for(j=0;j<A.n;j++){
			norm+= A(i,j)*A(i,j);
		}
	}
	return sqrt(norm);
	
}







void restrict(dgematrix& A,double min,double max){
        int i,j;
        for(i=0;i<A.m;i++){
        for(j=0;j<A.n;j++){
                if(A(i,j)<min){A(i,j)=min;}
                if(A(i,j)>max){A(i,j)=max;}
        }}
}

void restrict(dcovector &v, double min,double max){

        int i;
        for(i=0;i<v.l;i++){
                if(v(i)<min){v(i)=min;}
                if(v(i)>max){v(i)=max;}
        }
}

void restrict(drovector &v, double min,double max){

        int i;
        for(i=0;i<v.l;i++){
                if(v(i)<min){v(i)=min;}
                if(v(i)>max){v(i)=max;}
        }
}

double trace(dgematrix A){

        int i;
        double trace=0;
        if(A.n==A.m){
                for(i=0;i<A.m;i++){
                       trace = A(i,i) + trace;
                }
        }
        else{
               cout << "正方行列でないのでトレースをとれません"<<endl;
        }
        return trace;
}

void read_multi_vector(dcovector *x,int num, int dim, char *filename){//一行10000バイトまでです．
	ifstream fin(filename);

	char ch[10000];
	string str;

	for(int i=0;i<num;i++){
		x[i].resize(dim);
		fin.getline(ch,10000);
		str = ch;
		for(int j=0;j<dim;j++){
			string xs = ScanTillSpace(str);
			str = DelTillSpace(str);
			str = TrimLeft(str);
			x[i](j) = atof(xs.c_str());
		}
	}

	fin.close();
}
void read_multi_vector(drovector *x,int num, int dim, char *filename){//一行10000バイトまでです．
	ifstream fin(filename);

	char ch[10000];
	string str;

	for(int i=0;i<num;i++){
		x[i].resize(dim);
		fin.getline(ch,10000);
		str = ch;
		for(int j=0;j<dim;j++){
			string xs = ScanTillSpace(str);
			str = DelTillSpace(str);
			str = TrimLeft(str);
			x[i](j) = atof(xs.c_str());
		}
	}

	fin.close();
}

dgematrix product(dcovector c,drovector r){
	dgematrix A(c.l,r.l),C(c.l,1),R(1,r.l);
	
	for(int i=0;i<c.l;i++){
		C(i,0)=c(i);
	}
	for(int j=0;j<r.l;j++){
		R(0,j)=r(j);
	}
	A = C*R;

	return A;
}



void vec_set(dgematrix &A, int k , dcovector x){
	if(k>=A.n){cout << "k is too learge for matrix A" <<endl;}
	else if(A.m!=x.l){cout << "vec dimension is not match to matrix's dimension."<<endl;}
	else{
		for(int i=0;i<A.m;i++){
			A(i,k)=x(i);
		}
	}
}
void vec_set(dgematrix &A, int k , drovector x){
	if(k>=A.m){cout << "k is too learge for matrix A" <<endl;}
	else if(A.n!=x.l){cout << "vec dimension is not match to matrix's dimension."<<endl;}
	else{
		for(int i=0;i<A.n;i++){
			A(k,i)=x(i);
		}
	}
}

void read_multi_matrix(char *filename,dgematrix *A,int n){
//read matrix files from "./filename/??.dat"
// ?? means indices from 0.
	string str,dir;
	str = filename;
	dir = string("./")+str;
	dir += "/";

	string head;
	head = dir;

	for(int i=0;i<n;i++){
		string file;
		file = head + itos(i)+string(".dat");
		A[i].read(file.c_str());
	}
}

int ProbSelect(dcovector v){
        dcovector p = v;

        double sum=0;
        for(int i=0;i<v.l;i++){
                sum+=v(i);
        }

        p *= (1.0/sum);

        double t = double(rand())/(double(RAND_MAX) + 1.0);
        int k =0;
        for(int i=0;i<v.l;i++){
                if(t>=0){
                        t=t-p(i);
                        if(t<0){k=i;}
                }
        }

        return k;
}

double max_value(dcovector v){
       double max =  v(0);
       for(int i=0;i<v.l;i++){if(v(i)>max){max = v(i);}}
       return max;
}

int argmax(dcovector v){

       int i;
       double max =  v(0);
       int argmax = 0;
       for(i=0;i<v.l;i++){
               if(v(i)>max){max = v(i);argmax=i;}
       }

       return argmax;

}

int argmin(dcovector v){
        int i;
        double min = v(0);
        int argmin = 0;
		for(i=0;i<v.l;i++){
                if(v(i)<min){min = v(i); argmin=i;}
        }

        return argmin;
}

double min_value(dcovector v){
        int i;
        double min = v(0);
        int argmin = 0;
        for(i=0;i<v.l;i++){
                if(v(i)<min){min = v(i); argmin=i;}
        }
        return min;
}

dcovector argmax_vector(dcovector v){
	dcovector w(v);
	w.zero();

	w(argmax(v)) = 1.0;
	return w;

}
dcovector argmin_vector(dcovector v){
	dcovector w(v);
	w.zero();

	w(argmin(v)) = 1.0;
	return w;
}



dcovector e_greedy(dcovector a,double greedy_ratio){
	double other_ratio=(1.0-greedy_ratio)/(double(a.l)-1.0);
	dcovector p(a.l);
	for(int i=0;i<a.l;i++){
		p(i)=other_ratio;
	}

	p(argmax(a)) = greedy_ratio;

	return p;

}

dcovector normalize(dcovector x){
	double sum = 0;
	for(int i=0;i<x.l;i++){
		sum += x(i);
	}

	return (1.0/double(sum))*x;
}

dcovector each_product(dcovector a,dcovector b){
	dcovector c =a;

	for(int i=0;i<a.l;i++){
		c(i) = a(i)*b(i);
	}

	return c;
}

double det(dgematrix q){

	dgematrix V,D;
	dcovector s;

	if(0!=q.dgesvd(s,V,D)){cout << "det error!"<< endl << "q=" <<q << endl;/*getch()*/;}

	double d = 1.0;
	for(int i=0;i<s.l;i++){
		d *= s(i);
	}

	return d;
}

double gaussian(dcovector x, dcovector m,dgematrix q){
	dcovector y  = x-m;
	drovector yt = t(y);
	dgematrix qi = i(q);
	return exp(-0.5*yt*qi*y)/sqrt(pow(2.0*M_PI,(double)x.l)*det(q));
}

double LGF(dcovector x, dcovector m,dgematrix q){
	dcovector y  = x-m;
	drovector yt = t(y);
	dgematrix qi = i(q);

	return -0.5*yt*qi*y-log(sqrt(pow(2.0*M_PI,(double)x.l)*det(q)));
}

double gaussian(dcovector x, dcovector m,dgematrix iq,double detq){
	dcovector y  = x-m;
	drovector yt = t(y);
	return exp(-0.5*yt*iq*y)/sqrt(pow(2.0*M_PI,(double)x.l)*detq);
}
double LGF(dcovector x, dcovector m,dgematrix iq,double detq){
	dcovector y  = x-m;
	drovector yt = t(y);
	return -0.5*yt*iq*y-log(sqrt(pow(2.0*M_PI,(double)x.l)*detq));
}

string csv(dcovector x){
	string out("");
	string comma(",");
	for(int i=0;i<x.l;i++){
		out += comma + ftos(x(i));
	}
	return out;
}


string csv(drovector x){
	string out("");
	string comma(",");
	for(int i=0;i<x.l;i++){
		out += comma + ftos(x(i));
	}
	return out;
}

string csv(dgematrix x){
	string out("");
	string comma(",");
	for(int i=0;i<x.m;i++){
	for(int j=0;j<x.n;j++){
		out += comma + ftos(x(i,j));
	}}
	return out;
}

void self_add(dgematrix &A,dgematrix B){
	for(int i=0;i<A.m;i++){
		for(int j=0; j<A.n;j++){
			A(i,j) += B(i,j);
		}
	}
}
void self_add(dcovector &x,dcovector y){
	for(int i=0;i<x.l;i++){
		x(i) += y(i);
	}
}
void self_set(dgematrix &A,dgematrix B){
	for(int i=0;i<A.m;i++){
		for(int j=0; j<A.n;j++){
			A(i,j) = B(i,j);
		}
	}
}
void self_set(dcovector &x,dcovector y){
	for(int i=0;i<x.l;i++){
		x(i) = y(i);
	}
}

void all_set(dgematrix &A,double value){
	for(int i=0;i<A.m;i++){
		for(int j=0; j<A.n;j++){
			A(i,j) = value;
		}
	}
}

void all_set(dcovector &v,double value){
	for(int i=0;i<v.l;i++){
		v(i) = value;
	}	
}


dgematrix submatrix(dgematrix A, int L, int R, int T, int B){
	dgematrix S(B-T+1,R-L+1);

	for(int i=0;i<B-T+1;i++){
		for(int j=0;j<R-L+1;j++){
			S(i,j)=A(T+i,L+j);		
		}
	}
	return S;
}

dgematrix lean_on (dgematrix A, dgematrix B){
	dgematrix C(A.m,A.n+B.n);

	for(int i=0;i<A.m;i++){
		for(int j=0;j<A.n;j++){
			C(i,j)=A(i,j);
		}
	}
	for(int i=0;i<B.m;i++){
		for(int j=0;j<B.n;j++){
			C(i,A.n+j)=B(i,j);
		}
	}
	return C;

}
dgematrix  put_on (dgematrix A, dgematrix B){
	dgematrix C(A.m+B.m,A.n);

	for(int i=0;i<A.m;i++){
		for(int j=0;j<A.n;j++){
			C(i,j)=A(i,j);
		}
	}
	for(int i=0;i<B.m;i++){
		for(int j=0;j<B.n;j++){
			C(A.m+i,j)=B(i,j);
		}
	}

	return C;
}

dcovector covec_read(dgematrix A, int k){
  dcovector v(A.m);
  for(int i=0; i< A.m;i++){
    v(i)= A(i,k); 
  }
  return v;
}

drovector rovec_read(dgematrix A, int k){
  drovector v(A.n);
  for(int i=0;i<A.n;i++){
    v(i) = A(k,i);
  }
  return v;
}
