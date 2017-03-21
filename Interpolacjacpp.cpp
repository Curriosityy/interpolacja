// Interpolacjacpp.cpp: Okreœla punkt wejœcia dla aplikacji konsoli.
//

#include <iostream>
#include <conio.h>
using namespace std;
int silnia(int i)
{
	int sil=1;
	for (int j = 1; j <= i; j++)
	{
			sil *= j;
	}
	return sil;
}
void rozniceProgrwsywne(int tabx[], int taby[], double *rozProgTab[], int ilosc_w)
{
	for (int i = 0; i < ilosc_w; i++)
	{
		for (int j = 0; j < ilosc_w-i; j++)
		{
			if (i == 0)
			{
				rozProgTab[0][j] = taby[j];
			}
			else
			{
				rozProgTab[i][j] = rozProgTab[i-1][j + 1] - rozProgTab[i-1][j];
			}
		}
	}
}
double interpolacjaNewtona2(int tabx[], int taby[], double *rozProgTab[], int ilosc_w,double x_szukany)
{
	double wx = rozProgTab[0][0];
	double h = tabx[1] - tabx[0];
	for (int i = 1; i < ilosc_w; i++)
	{
		double wj = 1;
		for (int j = 0; j <= i-1; j++)
		{
			wj *= (x_szukany - tabx[j]);
		}
		double pod = (silnia(i)*pow(h, i));
		wx += (wj*(rozProgTab[i][0]/pod));
	}
	return wx;
}
void ilorazyRoznicowe(int tabx[], int taby[], double *tabiloczyny[], int ilosc_w)
{
	for (int i = 0; i < ilosc_w; i++)
	{
		for (int j = 0; j < ilosc_w-i; j++)
		{
			if (i == 0)
			{
				tabiloczyny[0][j]= taby[j];
			}
			else
			{
				tabiloczyny[i][j] = (tabiloczyny[i - 1][j + 1] - tabiloczyny[i - 1][j]) / (tabx[j+i] - tabx[j]);
			}
		}
	}
}


double interpolacjaNewtona1(int tabx[], int taby[], double *tabiloczyny[], int ilosc_w,double x_szukany)
{
	double wx = 0;
	double wj = 1;
	wx = taby[0];
	for (int i = 1; i < ilosc_w;i++)
	{
		wj = 1;
		for (int j = 0; j <= i-1; j++)
		{
			wj *= (x_szukany - tabx[j]);
		}
		wx += (tabiloczyny[i][0] * wj);
	}
	return wx;
}


double interpolacjaLagrangea(int ilosc_w, int tabx[], int taby[],double x_szukany)
{
	double wx = 0;
	for (int i = 0; i < ilosc_w; i++)
	{
		double li = 1;
		for (int j = 0; j < ilosc_w; j++)
		{
			if (i != j)
				li *= (x_szukany - tabx[j]) / (tabx[i] - tabx[j]);
		}
		wx += taby[i]*li;
	}

	return wx;
}

int main()
{
	int ilosc_w;
	double x_szukany;
	cout << "podaj ilosc wezlow" << endl;
	cin >> ilosc_w;
	int *tabx = new int[ilosc_w];
	int *taby = new int[ilosc_w];
	double **tabiloczyny = new double*[ilosc_w];
	for (int i = 0; i < ilosc_w; i++)
	{
		tabiloczyny[i] = new double[ilosc_w];
	}
	double **rozProgTab = new double*[ilosc_w];
	for (int i = 0; i < ilosc_w; i++)
	{
		rozProgTab[i] = new double[ilosc_w];
	}
	for (int i = 0; i < ilosc_w; i++)
	{
		cin >> tabx[i];
		cin >> taby[i];
	}
	cout << "podaj dla jakiego x znalezc" << endl;
	cin >> x_szukany;
	cout << "Lagrange'a dla x " << x_szukany << " y to " << interpolacjaLagrangea(ilosc_w, tabx, taby, x_szukany) << endl;
	ilorazyRoznicowe(tabx,taby,tabiloczyny,ilosc_w);
	cout << "Newtona ilorazy dla x " << x_szukany << " y to " << interpolacjaNewtona1(tabx, taby, tabiloczyny, ilosc_w, x_szukany) << endl;
	rozniceProgrwsywne(tabx, taby, rozProgTab, ilosc_w);
	cout << "Newtona progrestwna dla x " << x_szukany << " y to " << interpolacjaNewtona2(tabx, taby, rozProgTab, ilosc_w, x_szukany) << endl;
	getch();
	return EXIT_SUCCESS;
}

