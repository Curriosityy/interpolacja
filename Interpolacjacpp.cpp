// Interpolacjacpp.cpp: Okreœla punkt wejœcia dla aplikacji konsoli.
//

#include <iostream>
#include <conio.h>
using namespace std;


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
	for (int i = 0; i < ilosc_w; i++)
	{
		cin >> tabx[i];
		cin >> taby[i];
	}
	cout << "podaj dla jakiego x znalezc" << endl;
	cin >> x_szukany;
	cout<< "dla x "<<x_szukany<<"y to "<<interpolacjaLagrangea(ilosc_w, tabx, taby,x_szukany);
	getch();
	return 0;
}

