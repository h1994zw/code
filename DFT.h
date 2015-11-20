#pragma once
#include"math.h"
typedef float myType;
struct plur
{
	myType real, imag;
};
void re_DFT(float an[], float bn[], float ak[], float bk[], float n_length, float N, float choose);
void conv(float yn[], float hn[], float xn[], int L1, int L2);
void DFT_FFT(float AR[], float AI[], int m, int n, int choose);
void circle_conv(float ai[], float ar[], float bi[], float br[], int length_a, int length_b, float nr[], float ni[]);
void line_back(int xn[], int n, int len);
void back_cal(int V[], int n);
plur mul_plu(plur plu_a, plur plu_b);
void DFT_MY(plur plu[], int n,int N,int choose);
void DFT_FFT_1(plur plu[], int m, int n, int choose);
bool circleConv(plur pluA[], plur pluB[], int lenA, int lenB, int L);