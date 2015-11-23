#include"stdafx.h"
#include "DFT.h"
#include"math.h"
#ifndef PI
#define PI 3.14159265
#endif
/*****************
输入参数：DFT和IDFT变换序列的实部an和虚部bn
n_length 序列长度
N n点变换
choose 1.....DFT
       -1....IDFT
输出：ak，bk
******************/
void re_DFT(float an[], float bn[], float ak[], float bk[], float n_length, float N, float choose)
{
	if (choose == 1)
	{
		for (int i = 0; i < N; i++)
		{
			ak[i] = 0;
			bk[i] = 0;
			for (int j = 0; j < n_length; j++)
			{
				ak[i] += an[j] * cos(PI*(float)(2 * i*j) / (float)N) + choose* bn[j] * sin(PI*(float)(2 * i*j) / (float)N);
				bk[i] += bn[j] * cos(PI*(float)(2 * i*j) / (float)N) + choose* an[j] * sin(PI*(float)(2 * i*j) / (float)N);
			}
		}
	}
	if (choose == -1)
	{
		for (int i = 0; i < N; i++)
		{
			ak[i] /= N;
			bk[i] /= N;
		}
	}
}
/*
卷积：输入hn*xn 
输出yn
*/
void conv(float yn[],float hn[],float xn[],int L1,int L2)
{
	float num = 0;
	for (int i = 0; i < L1 + L2 - 1; i++)
	{
		yn[i] = 0;
	}
	for (int i = 0; i < L1 + L2 - 1; i++)
	{
		if (i < L1)
		{
			for (int j = 0; j <= i; j++)
			{
				yn[i] += hn[j] * xn[i - j];
				//yn[0]=hn[0]*xn[0]
				//yn[1]=hn[0]*xn[1]+hn[1]*xn[0];
				//yn[2]=hn[0]*xn[2]+hn[1]*xn[1]+hn[2]*xn[0];
				//yn[3]=hn[0]*xn[3]+hn[1]*xn[2]+hn[2]*xn[1]+hn[3]*xn[0];
			}
		}
		else{
			int k = L1 - 1;
			for (int j = 1; j <= L1 + L2 - 1 - i; j++, k--)
			{
				yn[i] += xn[k] * hn[i - k];
			}
		}
	}
}
/*
**功能：求解DFT
**入口参数：AR  AI  m：序列长度   n：n点DFT
**choose=1---FFT	choose=-1------IFFT
**出口参数：AR  AI
**
**/
void DFT_FFT(float AR[], float AI[], int m, int n,int choose)
{//m=5 n=32 choose=1
	int i, i1, j1, j, LH, n1, k, b;
	float tr, ti, p;
	if (choose == -1)
	{
		for (i = 0; i < n; i++)
			AI[i] = 0 - AI[i];
	}
	LH = (int)n / 2;
	j = LH;
	n1 = n - 2;//30
	for (i = 1; i <= n1; i++)//1-30
	{
		if (i<j)
		{
			i1 = i;
			j1 = j;
			tr = AR[i];
			ti = AI[i];
			AR[i1] = AR[j1];
			AI[i1] = AI[j1];
			AR[j1] = tr;
			AI[j1] = ti;
		}
		k = LH;
		while (j >= k)
		{
			j = j - k;
			k = (int)k / 2;
		}
		j = j + k;
	}//
	for (i = 1; i <= m; i++)//12345
	{
		b = (int)pow(2.0, i - 1);//012-124
		for (j = 0; j <= b - 1; j++)//0-0  0-1  0-3
		{//j=0-01-0123
			p = (float)(pow(2.0, m - i)*j*2.0*PI / (float)n);//蝶形因子
			//p=0-02-0123  
			for (k = j; k <= n - 1;)
			{//j=0 p=0 b=1
			//k=0
				tr = AR[k + b] * cos(p) + AI[k + b] * sin(p);
				ti = AI[k + b] * cos(p) - AR[k + b] * sin(p);
				AR[k + b] = AR[k] - tr;
				AI[k + b] = AI[k] - ti;
				AR[k] = AR[k] + tr;
				AI[k] = AI[k] + ti;
				k += b * 2;
			}
		}
	}
	if (choose == -1)
	{
		for (i = 0; i < n; i++)
		{
			AR[i] = AR[i] /n;
			AI[i] = 0-AI[i]/n ;
		}
	}
}
void circle_conv(float ai[],float ar[],float bi[],float br[],int length_a,int length_b,float nr[],float ni[])
{
	int length = length_a + length_b-1;
	float *xi = new float[64];
	float *xr = new float[64];
	float *yi = new float[64];
	float *yr = new float[64];
	//for (int i = length_a; i < length; i++)
	//{
	//	xi[i] = 0;
	//	xr[i] = 0;
	//}
	for (int i = 0; i < length; i++)
	{
		if (i < length_a)
		{
			xi[i] = ai[i];
			xr[i] = ar[i];
		}
		else
		{
			xr[i] = 0;
			xi[i] = 0;
		}
		if (i < length_b)
		{
			yi[i] = bi[i];
			yr[i] = br[i];
		}
		else
		{
			yi[i] = 0;
			yr[i] = 0;
		}
	}
	DFT_FFT(xr, xi, length, 64, 1);
	DFT_FFT(yr, yi, length, 64, 1);
//	float *ni = new float[64];
//	float *nr = new float[64];
	for (int i = 0; i < length; i++)
	{
		//(a+bi)*(x+yi)=ax-by  +(bx+ay)i
		nr[i] = xr[i] * yr[i] - xi[i] * yi[i];
		ni[i] = xr[i] * yi[i] + xi[i] * yr[i];
	}
	DFT_FFT(nr, ni, length, 64, -1);
	//int length;
	//if (length_a>length_b)
	//{
	//	float *xi = new float[length_a];
	//	float *xr = new float[length_a];
	//	for (int i = 0; i < length_a; i++)
	//	{
	//		if (i >= length_b)
	//		{
	//			xi[i] = 0;
	//			xr[i] = 0;
	//		}
	//		else
	//		{
	//			xi[i] = ai[i];
	//			xr[i] = ar[i];
	//		}
	//	}
	//	length = length_a;
	//}
	//else if (length_a < length_b)
	//{
	//	float *xi = new float[length_b];
	//	float *xr = new float[length_b];
	//	for (int i = 0; i < length_b; i++)
	//	{
	//		if (i >= length_a)
	//		{
	//			xi[i] = 0;
	//			xr[i] = 0;
	//		}
	//		else
	//		{
	//			xi[i] = bi[i];
	//			xr[i] = br[i];
	//		}
	//	}
	//	length = length_b;
	//}
}
/*
*功能：移形换位
*输入：V[]——需要倒序的序列 n：2的n次方为序列长度
*输出：V[]——倒序后的序列
*/
void back_cal(plur V[], int n)
{
	int len = pow(2.0, n), j = 0, temp;
	j = len / 2;
	int next = len / 2, now = 0;
	//cout << V[0] << "\n";
	for (int i = 1; i<len; i++)
	{
		next = len / 2;
		while (now >= next)//8  4
		{
			now -= next;//0
			next /= 2;//4
		}
		now = now + next;//4
		if (now>i)
		{
			temp = V[i].real;
			V[i].real = V[now].real;
			V[now].real = temp;
			temp = V[i].imag;
			V[i].imag = V[now].imag;
			V[now].imag = temp;
		}
		//V[i]=now;
		//cout<<now<<","<<V[i]/4<<endl;
	}
}
/*
*函数功能：复数乘法
*输入参数：plu_a=ar+jai  plu_b=br+jbi
*返回值：ar*br-ai*bi  ar*bi+ai*br
*/
plur mul_plu(plur plu_a, plur plu_b)
{
	plur plu;
	plu.real = plu_a.real*plu_b.real - plu_a.imag*plu_b.imag;
	plu.imag = plu_a.real*plu_b.imag + plu_a.imag*plu_b.real;
	return plu;
}
/*
*蝶形因子
*/
plur buffer(float k,float N)
{
	plur plu;
	plu.real = cos(2.0 * PI*k / N);
	plu.imag = 0.0 - sin(2.0*PI*k/N);
	return plu;
}

/*
*复数加减法
*/
plur plus(plur p1, plur p2,float cho)
{
	plur p;
	p.real = p1.real + cho*p2.real;
	p.imag = p1.imag + cho*p2.imag;
	return p;
}
plur equal(plur p1)
{
	plur p;
	p.real = p1.real;
	p.imag = p1.imag;
	return p;
}
/*
*函数功能：FFT
*输入参数：序列，长度n，N点FFT，choose=-1取IFFT
*
*/
void DFT_MY(plur plu[],int n,int N,int choose)//(n=5,N=32)
{
	//倒序
	int len = pow(2.0, n);
	int tmp = 0;
	float tmp1 = 0.0, real, imag;
	plur save,save1;
	if (choose == -1)
	{
		for (int i = 0; i < N; i++)
			plu[i].imag =  0.00- plu[i].imag;
	}
	back_cal(plu, n);
	for (int i = 1; i <= n; i++)//12345
	{
		tmp = (int)pow(2.0, i - 1);
		for (int j = 0; j <= tmp - 1; j++)
		{
			tmp1 = (float)(pow(2.0, n - i)*(float)j*2.0*PI / (float)N);
			for (int k = j; k <= N - 1;)
			{
				real = plu[k + tmp].real*cos(tmp1) + plu[k + tmp].imag*sin(tmp1);
				imag = plu[k + tmp].imag*cos(tmp1) - plu[k + tmp].real*sin(tmp1);
				plu[k + tmp].real = plu[k].real - real;
				plu[k + tmp].imag = plu[k].imag - imag;
				plu[k].real = plu[k].real + real;
				plu[k].imag = plu[k].imag + imag;
				//save = mul_plu(plu[k + tmp], save1);
				//plu[k + tmp] = plus(plu[k], save, -1.0);
				//plu[k] = plus(plu[k], save, 1);
				k += tmp * 2;
			}
		}
		//tmp = pow(2.0, i);//1.2.4
		//tmp1 = len / 2 / tmp;//
		//for (int j = 0; j < tmp1; j++)//8  //4  2
		//{
		//	for (int k = 0; k < tmp; k++)//1  202   40123
		//	{
		//		//save = plus(plu[j * 2], mul_plu(buffer(k, len), plu[j * 2 + tmp]), 1);
		//		//i=0 tmp=1 tmp1=8   i=1 tmp=2 tmp1=4  i=2,tmp=4,tmp1=2
		//		//j=0 j=1   0123 048 12  
		//		//0---   01

		//		//save = mul_plu(buffer(pow(2.0, n - i - 1)*k, len), plu[j * 2 * tmp + tmp]);
		//		//save1 = plu[j * 2 * tmp + tmp];
		//		//plu[j * 2 * tmp + tmp] = plus(save, plu[j * 2], 1);
		//		//plu[j * 2 * tmp + tmp] = plus(save1, save, -1);
		//	}
		//}
	}
	if (choose == -1)
	{
		for (int i = 0; i < N; i++)
		{
			plu[i].real /= (float)N;
			plu[i].imag =0- plu[i].imag/(float)N;
		}
	}
}
void DFT_FFT_1(plur plu[], int m, int n, int choose)//float AR[], float AI[]
{//m=5 n=32 choose=1
	int i, i1, j1, j, LH, n1, k, b;
	float tr, ti, p;
	plur pl,p2;
	if (choose == -1)
	{
		for (i = 0; i < n; i++)
			plu[i].imag = 0.00 - plu[i].imag;
			//AI[i] = 0 - AI[i];
	}
	back_cal(plu, m);
	for (i = 1; i <= m; i++)//12345
	{
		b = (int)pow(2.0, i - 1);//012-124
		for (j = 0; j <= b - 1; j++)//0-0  0-1  0-3
		{//j=0-01-0123
			p = (float)(pow(2.0, m - i)*j*2.0*PI / (float)n);//蝶形因子
			//p=0-02-0123  
			for (k = j; k <= n - 1;)
			{//j=0 p=0 b=1
				//k=0
				//p2 = mul_plu(pl, plu[k + b]);
				p2.real = plu[k + b].real*cos(p) + plu[k + b].imag*sin(p);
				p2.imag = plu[k + b].imag*cos(p) - plu[k + b].real*sin(p);
				//tr = AR[k + b] * cos(p) + AI[k + b] * sin(p);
				//ti = AI[k + b] * cos(p) - AR[k + b] * sin(p);
				//plu[k + b] = plus(plu[k], p2,-1.0);
				//plu[k] = plus(plu[k], p2, 1.0);
				plu[k + b].real = plu[k].real - p2.real;
				plu[k + b].imag = plu[k].imag - p2.imag;
				//AR[k + b] = AR[k] - tr;
				//AI[k + b] = AI[k] - ti;
				//AR[k] = AR[k] + tr;
				//AI[k] = AI[k] + ti;
				plu[k].real = plu[k].real + p2.real;
				plu[k].imag = plu[k].imag + p2.imag;
				k += b * 2;
			}
		}
	}
	if (choose == -1)
	{
		for (i = 0; i < n; i++)
		{
			plu[i].real /= (float)n;
			plu[i].imag = 0.0 - plu[i].imag / (float)n;
			//AR[i] = AR[i] / n;
			//AI[i] = 0 - AI[i] / n;
		}
	}
}
int max_AB(int A, int B)
{
	if (A > B)return A;
	else return B;
}
/*
*函数功能：循环卷积
*输入:pluA  长度lenA  pluB  长度lenB  循环卷积长度L
*输出：pluA卷积后的序列
*/
bool circleConv(plur pluA[], plur pluB[], int lenA, int lenB, int L)
{
	if (L < max_AB(lenA, lenB))
		return FALSE;
	int n_fir=0;
	while (L>pow(2.0, n_fir))n_fir++;
	int n = pow(2.0, n_fir);
	float real, imag;
	float *AI = new float[n];
	float *AR = new float[n];
	float *BI = new float[n];
	float *BR = new float[n];
	for (int i = 0; i < n; i++)
	{
		if (i < lenA)
		{
			AI[i] = pluA[i].imag;
			AR[i] = pluA[i].real;
		}
		else
		{
			AI[i] = 0;
			AR[i] = 0;
		}
		if (i < lenB)
		{
			BI[i] = pluB[i].imag;
			BR[i] = pluB[i].real;
		}
		else
		{
			BI[i] = 0;
			BR[i] = 0;
		}
	}
	DFT_FFT(AR, AI, n_fir, n, 1); 
	DFT_FFT(BR, BI, n_fir, n, 1);
	for (int i = 0; i < n; i++)
	{
		real = AR[i] * BR[i] - AI[i] * BI[i];
		imag = AR[i] * BI[i] + AI[i] * BR[i];
		AR[i] = real;
		AI[i] = imag;
	}
	DFT_FFT(AR, AI, n_fir, n, -1);
	for (int i = 0; i < L; i++)
	{
		pluA[i].real = AR[i];
		pluA[i].imag = AI[i];
	}
	return TRUE;
}