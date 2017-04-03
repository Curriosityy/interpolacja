
// Interpolacjacpp.cpp: Okreœla punkt wejœcia dla aplikacji konsoli.
//


// Eliminacja Gaussa
// Data: 15.02.2010
// (C)2012 mgr Jerzy Wa³aszek
//-----------------------------

#include <iostream>
#include <conio.h>
#include <cmath>
using namespace std;


const double eps = 1e-12; // sta³a przybli¿enia zera

						  // Funkcja realizuje algorytm eliminacji Gaussa
						  //---------------------------------------------
bool gauss(int n, double ** AB, double * X)
{

	int i, j, k;
	double m, s;

	// eliminacja wspó³czynników

	for (i = 0; i < n - 1; i++)
	{
		for (j = i + 1; j < n; j++)
		{
			if (fabs(AB[i][i]) < eps) return false;
			m = -AB[j][i] / AB[i][i];
			for (k = i + 1; k <= n; k++)
				AB[j][k] += m * AB[i][k];
		}
	}

	// wyliczanie niewiadomych

	for (i = n - 1; i >= 0; i--)
	{
		s = AB[i][n];
		for (j = n - 1; j >= i + 1; j--)
			s -= AB[i][j] * X[j];
		if (fabs(AB[i][i]) < eps) return false;
		X[i] = s / AB[i][i];
	}
	return true;
}
//------------------------------------------------------------//
int silnia(int i)
{
	int sil = 1;
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
		for (int j = 0; j < ilosc_w - i; j++)
		{
			if (i == 0)
			{
				rozProgTab[0][j] = taby[j];
			}
			else
			{
				rozProgTab[i][j] = rozProgTab[i - 1][j + 1] - rozProgTab[i - 1][j];
			}
		}
	}
}
double interpolacjaNewtona2(int tabx[], int taby[], double *rozProgTab[], int ilosc_w, double x_szukany)
{
	double wx = rozProgTab[0][0];
	double h = tabx[1] - tabx[0];
	for (int i = 1; i < ilosc_w; i++)
	{
		double wj = 1;
		for (int j = 0; j <= i - 1; j++)
		{
			wj *= (x_szukany - tabx[j]);
		}
		double pod = (silnia(i)*pow(h, i));
		wx += (wj*(rozProgTab[i][0] / pod));
	}
	return wx;
}
void ilorazyRoznicowe(int tabx[], int taby[], double *tabiloczyny[], int ilosc_w)
{
	for (int i = 0; i < ilosc_w; i++)
	{
		for (int j = 0; j < ilosc_w - i; j++)
		{
			if (i == 0)
			{
				tabiloczyny[0][j] = taby[j];
			}
			else
			{
				tabiloczyny[i][j] = (tabiloczyny[i - 1][j + 1] - tabiloczyny[i - 1][j]) / (tabx[j + i] - tabx[j]);
			}
		}
	}
}


double interpolacjaNewtona1(int tabx[], int taby[], double *tabiloczyny[], int ilosc_w, double x_szukany)
{
	double wx = 0;
	double wj = 1;
	wx = taby[0];
	for (int i = 1; i < ilosc_w; i++)
	{
		wj = 1;
		for (int j = 0; j <= i - 1; j++)
		{
			wj *= (x_szukany - tabx[j]);
		}
		wx += (tabiloczyny[i][0] * wj);
	}
	return wx;
}


double interpolacjaLagrangea(int ilosc_w, int tabx[], int taby[], double x_szukany)
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
		cout << taby[i] << "*" << li<<"="<< taby[i] * li <<"\n";
		wx += taby[i] * li;
	}

	return wx;
}
void tworzenieTabeliDoGaussa(int tabx[], int taby[], int *tabPoch, double *AB[], int wielkoscTab)
{

	int poch = 0;
	for (int i = 0; i < wielkoscTab; i++)
	{
		if (i < wielkoscTab - 2)
		{
			for (int j = 0; j <= wielkoscTab; j++)
			{
				if (j == wielkoscTab)
				{
					AB[i][j] = taby[i];
				}
				else if (j < 4)
				{
					AB[i][j] = pow(tabx[i], j);
				}else
				{
					AB[i][j] = 0;
				}
			}
		}
		else
		{
			for (int j = 0; j <= wielkoscTab; j++)
			{
				if (j == wielkoscTab)
				{
					AB[i][j] = tabPoch[poch];
					poch++;
				}
				else
				if (j < 4)
				{
					AB[i][j] = j*pow(tabx[poch*(wielkoscTab - 3)], j - 1);
				}
				else
				{
					AB[i][j] = 0;
				}
			}
		}
	}

	for (int i = 1; i < wielkoscTab; i++)
	{
		if (i < wielkoscTab -2) {
			for (int j = wielkoscTab - (wielkoscTab - 4); j < wielkoscTab; j++)
			{
				if(j-2<=i)
				AB[i][j] = pow(tabx[i]-tabx[j-3], 3);
			}
		}else
		{
			for (int j = wielkoscTab - (wielkoscTab - 4); j < wielkoscTab; j++)
			{
				if(i=wielkoscTab-1)
				AB[i][j]=3*pow(tabx[(wielkoscTab - 3)]-tabx[j-3],2);
			}
		}
	}
}
double interpolacjaSklejanie(double x_szukany, double X[],int ilosc_w,int tabx[])
{
	double wynik = 0;
	for (int i = 0; i < 4;i++)
	{
		wynik += (pow(x_szukany, i) * X[i]);
	}
	for (int i = 1; i < ilosc_w; i++)
	{
		if (x_szukany > tabx[i]) {
			wynik += X[3 + i]*pow((x_szukany-tabx[i]),3);
		}
		else break;
	}
	return wynik;
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
	ilorazyRoznicowe(tabx, taby, tabiloczyny, ilosc_w);
	cout << "Newtona ilorazy dla x " << x_szukany << " y to " << interpolacjaNewtona1(tabx, taby, tabiloczyny, ilosc_w, x_szukany) << endl;
	rozniceProgrwsywne(tabx, taby, rozProgTab, ilosc_w);
	cout << "Newtona progrestwna dla x " << x_szukany << " y to " << interpolacjaNewtona2(tabx, taby, rozProgTab, ilosc_w, x_szukany) << endl;

	//-------------funkcje sklejane------------------//


	// tworzymy macierze AB i X
	int n = ilosc_w + 2;
	double **AB = new double *[n];
	double *X = new double[n];
	for (int i = 0; i < n; i++)
		AB[i] = new double[n + 1];
	int tabelePoch[2];
	for (int i = 0; i < 2; i++)
	{
		cout << "pierwsza pochodna\n";
		cin >> tabelePoch[i];
	}
	tworzenieTabeliDoGaussa(tabx, taby, tabelePoch, AB, n);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j <= n; j++)
		{
			cout<<AB[i][j] << " ";
		}
		cout << "\n";
	}
	gauss(n, AB, X);
	for (int i = 0; i < n; i++)
	{
		cout << X[i] << " ";
	}
	cout << "\nMetoda sklejania dla x " << x_szukany << " y to " << interpolacjaSklejanie(x_szukany,X,ilosc_w,tabx)<< endl;
	getch();
	return EXIT_SUCCESS;
}

