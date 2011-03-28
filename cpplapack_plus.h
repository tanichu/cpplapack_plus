/*#pragma once
#pragma comment(lib,"libF77.lib")
#pragma comment(lib,"libI77.lib")
#pragma comment(lib,"blas.lib")
#pragma comment(lib,"clapack.lib")
*/
	
#define CSV_STRING vector< vector< string > >
#define _USE_MATH_DEFINES 
	
	
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

#include <cpplapack.h>//cpplapack
#include <randlib.h> // randlib


#include <time.h>//���Ԋ֐����p�@�����������Ȃ�
//#include <conio.h>//�R���\�[���ňꎞ�~�߂�getch()�𗘗p���邽��
#include <math.h>//���w�֐�
#include <string.h>//�����񏈗��p�̃N���Xstring���p�̂���
#include <fstream>//�t�@�C���X�g���[��
using namespace std;
using namespace CPPL;

//����̕�����ϊ��p
std::string itos(int n); // int -> string
std::string ftos(double n);// float,double -> string
string TrimRight( string arg ); //string �̉E���̃X�y�[�X�Ȃǂ��Ƃ�D
string TrimLeft( string arg );//string �̍����̃X�y�[�X�Ȃǂ��Ƃ�D
string Trim( string arg);//string �̗����̃X�y�[�X�Ȃǂ��Ƃ�D

string ScanTillRet(string str); // Return���̕�������擾
string ScanTillSpace(string str); // space���̕�������擾
string ScanTillTab(string str); // tab���̕�������擾
string ScanTillComma(string str);// comma���̕�������擾
string DelTillRet(string str);// Return���̕����������
string DelTillSpace(string str);// space���̕����������
string DelTillTab(string str);// Return���̕����������
string DelTillComma(string str);// Return���̕����������

//���ɂ��f�t�H���g��
// atof // string -> float
// atoi // string -> int
//�����݂��܂��D�����ƕ�����̑��ݕϊ��Ɏg���Ă��������D

//void DEBUG_TRACE(string x);
vector<string> split(const string& str, const string& delimiter);

vector< vector<string> > csv_reader(const char *filename);
vector< vector<string> > ssv_reader(const char *filename);




double gauss_rand(); //���U�P���ςO�̃K�E�X�����𔭐�
double gauss_rand(double ave,double sig);//����ave,�W���΍�sig�̃K�E�X�����𔭐�
dcovector gauss_rand(dcovector ave,double sig);//�x�N�g���ŕ���ave,�S�����Ɨ��ŕ΍�sig�̃K�E�X����
dcovector gauss_rand(dcovector ave,dcovector sig);//�x�N�g���ŕ���ave,�S�����Ɨ��Ŋe�����΍�sig�̃K�E�X����
dgematrix gauss_rand(dgematrix ave,double sig);//�s��̊e�����ɕ΍�sig�̗������ڂ���

double LGF(dcovector x,dcovector m,dcovector s);//�ΐ��������K�E�X�֐��C���Sm�C�e�����Ɨ��ŕW���΍�s�̃K�E�X�֐���log�ɓ��ꂽ���́D�K�E�X�֐��̂����Z�����ۂȂǂɁC�ΐ���ԓ��Řa�Z�ő�p�����ق������x���o��D
double LGF(dcovector x, dcovector m,dgematrix q); //
double LGF(dcovector x, dcovector m,dgematrix iq,double detq);


dcovector ng_calculation(dcovector *m,dcovector *s,int N,dcovector x);//�����̃K�E�X���z�̒��S�C�W���΍��C�K�E�X���z�̐�N��n�����ƂŐ��K���K�E�X�֐�(NGnet)�̏o�͂��o���D�o�͂͊e��ꂩ��̏o�͒l
dcovector rbf_calculation(dcovector *m,dcovector *s,int N,dcovector x);//��̂q�a�e��

void restrict(dcovector &v, double min,double max);//�x�N�g���̒l��min,max�̊Ԃɐ�������i�O��l�����ȂǂɁj
void restrict(drovector &v, double min,double max);//�x�N�g���̒l��min,max�̊Ԃɐ�������
void restrict(dgematrix& A,double min,double max);//�s��̒l��min,max�̊Ԃɐ�������
void  each_square(dcovector &v);//�e�����ɂ��Ď��悷��D
void  each_abs(dcovector &v);//�e�����ɂ��Đ�Βl��^����D
void scatter();//�����������@���Ԋ֐����Ăяo����rand�֐�������������D�v���O�����̈�Ԏn�߂łƂ肠�����ĂԂׂ��D

double frobnorm(dgematrix A);//�s��̃t���x�j�E�X�m�������v�Z
double trace(dgematrix A);//�s��̃g���[�X���v�Z

void read_multi_vector(dcovector *x,int num, int dim, char *filename);//�t�@�C�����畡���̃x�N�g����ǂݎ��D
void read_multi_vector(drovector *x,int num, int dim, char *filename);//�t�@�C�����畡���̃x�N�g����ǂݎ��D
void read_multi_matrix(char *filename,dgematrix *A,int n);//�t�@�C�����畡���̍s���ǂݎ��D


dgematrix product(dcovector c,drovector r);//��x�N�g���ƍs�x�N�g���̐ςōs��𑢂�D
//dgematrix sup::operator*(dcovector c,drovector r);

void vec_set(dgematrix &A, int k , dcovector x); // �s��̈ꕔ�Ƀx�N�g�����Z�b�g

void vec_set(dgematrix &A, int k , drovector x);
dcovector covec_read(dgematrix A, int k);
drovector rovec_read(dgematrix A, int k);

int ProbSelect(dcovector v);//�e�����̒l���m���ƌ��Ȃ��C���̊m���ɉ����ăT�C�R����U�莟����I������D

double max_value(dcovector v);//�e�����l�̒��̂̍ő�l��Ԃ��D
double min_value(dcovector v);//�e�����l�̒��̍ŏ��l��Ԃ�
int argmax(dcovector v);//�ő�l����������Ԃ�
int argmin(dcovector v);//�ŏ��l����������Ԃ�
dcovector argmax_vector(dcovector v);//
dcovector argmin_vector(dcovector v);
dcovector e_greedy(dcovector a,double greedy_ratio);//ProbSlect��e-greedy�ŁC�ő�l����������1-(N-1)e�̊m���C����e=greedy_ratio�̊m���őI���D


dcovector normalize(dcovector x);//x���e�����̘a���P�ɂȂ�悤�ɐ��K������D
dcovector each_product(dcovector a,dcovector b);//
double det(dgematrix q); //�s�񎮂�Ԃ�
double gaussian(dcovector x, dcovector m,dgematrix q); //
double gaussian(dcovector x, dcovector m,dgematrix iq,double detq);

double norm2(dcovector x);//���m������Ԃ�


//dgematrix inv(dgematrix A);

string csv(dcovector x);//csv�̌`���̕�����ɏo�͂���D fout << csv(X) <<endl;�Ƃ������g����
string csv(drovector x);//����
string csv(dgematrix x);//����

void self_add(dgematrix &A,dgematrix B);//self�V���[�Y�͉��Z�q��p�����C����Ɍv�Z���ʂ�������D
void self_add(dcovector &x,dcovector y);//�����̑��x����͊��҂ł��邩������Ȃ���
void self_set(dgematrix &A,dgematrix B);//����Ȃɗ��p���l�͂Ȃ��D
void self_set(dcovector &x,dcovector y);//


void all_set(dgematrix &x,double value);//�s��̑S������value�ɃZ�b�g
void all_set(dcovector &x,double value);//�s��̑S������value�ɃZ�b�g

dgematrix submatrix(dgematrix A, int L, int R, int T, int B);

dgematrix lean_on (dgematrix A, dgematrix B);
dgematrix put_on (dgematrix A, dgematrix B);

//dgematrix operator = (dcovector v);
