// dit-fft.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include<conio.h>
#include"math.h"
#include <graphics.h>
#define PI 3.1415926
#define MAX 512
float AR[MAX];
float AI[MAX];
float A[MAX];
float AP[MAX];

void fun(unsigned char cho, int N);
void DFT_FFT(float AR[], float AI[], int m, int n)
{
	int i, i1, j1, j, LH, n1, k, b;
	float tr, ti, p;
	LH = (int)n / 2;
	j = LH;
	n1 = n - 2;
	for (i = 1; i <= n1; i++)
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

	}
	for (i = 1; i <= m; i++)
	{
		b = (int)pow(2.0, i - 1);
		for (j = 0; j <= b - 1; j++)
		{
			p = (float)(pow(2.0, m - i)*j*2.0*PI / (float)n);
			for (k = j; k <= n - 1;)
			{
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
}

int _tmain(int argc, _TCHAR* argv[])
{
	int gdriver, gmode;
	int i, j, k, M=5, N;
	//printf("input M:");
	//scanf_s("%d", &M);
	N = (int)pow(2.0, M);
	//printf("N=%d\n", N);
	//printf("\n input k:");
	//scanf_s("%d", &k);
	fun(1, N);
	//printf("the fun is %d", k);
	//printf("the values:\n");
	for (j = 0; j<N; j++)
	{
		if (AI[j] == 0.0)A[j] = AR[j];
		else A[j] = sqrt(pow(AR[j], 2) + pow(AI[j], 2));
		AP[j] = A[j];
		//printf("AP[%d]=%f", j, A[j]);
		if (j % 2 == 0)printf("\n");
	}

	DFT_FFT(AR, AI, M, N);
	printf("the value of re:\n");
	for (j = 0; j<N; j++)
	{
		A[j] = sqrt(AR[j] * AR[j] + AI[j] * AI[j]);
		printf("A[%d]=%f+j%f\n", j, AR[j], AI[j]);
	}
	
	initgraph(640, 480);
	setbkcolor(0);
	setlinecolor(WHITE);
		int m, x, y, L;
		x = 50; y = 200; i = 10;
		for (L = 0; L<N; L++)
		line(x + L*i, y, x + L*i, y - A[L] * 5);
		x = 50; y = 400; i = 10;
		line(30, 400, 609, 400);
		for (m = 0; m<N; m++)
		{
			line(x + m*i, y, x + m*i, y - AP[m] * 20);
		}
	
	getchar();
	closegraph();
	return 0;
}
void fun(unsigned char cho, int N)
{
	unsigned int i = 0;
	switch (cho)
	{
	case 0:
		for (i = 0; i<N; i++)
		{
			AI[i] = 0.0;
			if (i <= 3)AR[i] = 1.0;
			else AR[i] = 0.0;
		}
		break;
	case 1:
		for (i = 0; i<N; i++)
		{
			if (i <= 3)AR[i] = i + 1;
			else if (i <= 7)AR[i] = 8 - i;
			else AR[i] = 0.0;
			AI[i] = 0;
		}
		break;
	case 2:
		for (i = 0; i<N; i++)
		{
			if (i <= 3)AR[i] = 4 - i;
			else if (i <= 7)AR[i] = i - 3;
			else AR[i] = 0;
			AI[i] = 0;
		}
		break;
	case 3:
		for (i = 0; i<N; i++)
		{
			AR[i] = cos(PI / 4.00*(float)i);
			AI[i] = 0;
		}
		break;
	case 4:
		for (i = 0; i<N; i++)
		{
			AR[i] = sin(PI / 8.000*(float)i);
			AI[i] = 0;
		}
		break;
	case 5:
		for (i = 0; i<N; i++)
		{
			AR[i] = cos(8.0*PI*(float)i / 64.00) + cos(16.0*PI*(float)i / 64.00) + cos(24.0*PI*(float)i / 64.0);
			AI[i] = 0;
		}
		break;
	case 6:
		for (i = 0; i<N; i++)
		{

		}
		break;
	}
}
