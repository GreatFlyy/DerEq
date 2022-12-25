#include<iostream>
#include"Matrix.h"
#include"Reader.h"
#include<functional>
#include<math.h>

using namespace std;

double dot(vector<double>& a, vector<double>& b)
{
	int n = min(a.size(), b.size());
	double result = 0;
	for (int i = 0; i < n; i++)
	{
		result += a[i] * b[i];
	}
	return result;
}

double step(double y, double t, double h, vector<vector<double>>& coeffs, const function<double(double, double)>& f)
{
	int n = coeffs[0].size();
	vector<double> k(n, 0);
	double sum = 0;

	for (int i = 0; i < n; i++)
	{
		k[i] = f(t + coeffs[1][i] * h, y + dot(k, coeffs[i + 2]) * h);
		sum += coeffs[0][i] * k[i];
	}

	return (y + sum * h);
}

// coeffs[0] = b[]; coeffs[1] = c[]; coeffs[2-s+3] - a[][]

double sqrt1(double a)
{
	if (a >= 0)
		return sqrt(a);
	else
		return -sqrt(-a);
}



void pol2dec(vector<vector<double>>& in)
{
	for (auto& el : in)
	{
		double phi = el[0];
		double r = el[1];
		el[0] = r * cos(phi);
		el[1] = r * sin(phi);
	}
}

int main()
{
	string pathout = "C:/Users/GreatFly/Desktop/DerEq/out.txt";
	string pathoutbmp = "C:/Users/GreatFly/Desktop/DerEq/out.bmp";
	
	const double PI = 3.14159265359;

	double r0 = 10000;
	double t0 = 0;
	double h = 0.0001;
	int N = 100000;

	double M = 100;
	double b = 519.5;

	vector<vector<double>> table;
	table.push_back(vector<double>{1./6., 1./3., 1./3., 1./6.});//b
	table.push_back(vector<double>{0., 1./2., 1./2., 1.});//c
	table.push_back(vector<double>{});
	table.push_back(vector<double>{1./2.});
	table.push_back(vector<double>{0., 1./2.});
	table.push_back(vector<double>{0., 0., 1.});//a

	vector<vector<double>> output(N, vector<double>(2, 0));

	output[0][0] = t0;
	output[0][1] = r0;

	function<double(double, double)> equation = [b, M](double t, double y)
	{
		return -y * y * sqrt(1 / (b * b) - ((1 - 2 * M / y) / (y * y)));
	};	

	for (int i = 1; i < N; i++)
	{
		output[i][0] = output[i - 1][0] + h;
		output[i][1] = step(output[i - 1][1], output[i - 1][0], h, table, equation);
		if (isnan(output[i][1]))
		{
			equation = [b, M](double t, double y)
			{
				return y * y * sqrt(1 / (b * b) - ((1 - 2 * M / y) / (y * y)));
			};
			i--;
		}
		if (output[i][1] >= r0)
		{
			i = N;
		}
	}

	pol2dec(output);

	bmp image(1000, 1000, 24);
	int scale = 5000;
	int xb = 0;
	int yb = 0;


	double k = acos(b / r0);
	for (double i = 0; i < PI/2; i+=0.01)
	{
		image.pixarr[max(min((int)(b * (1/cos(i - k)) * sin(-i) / scale * image.heightpx + image.heightpx / 2 + yb), image.heightpx - 1), 0)][max(min((int)(b * (1 / cos(i - k)) * cos(i) / scale * image.widthpx + image.widthpx / 2 + xb), image.widthpx - 1), 0)] = vector<BYTE>{(BYTE)0, (BYTE)0, (BYTE)255};
	}

	for (int i = 0; i < 50; i++)
	{
		image.pixarr[max(min((int)(2 * M * sin(i * 2 * PI / 50) / scale * image.heightpx + image.heightpx / 2 + yb), image.heightpx - 1), 0)][max(min((int)(2 * M * cos(i * 2 * PI / 50) / scale * image.widthpx + image.widthpx / 2 + xb), image.widthpx - 1), 0)] = vector<BYTE>{ (BYTE)255, (BYTE)0, (BYTE)0 };
	}
	for (int i = 0; i < 50; i++)
	{
		image.pixarr[max(min((int)(3 * M * sin(i * 2 * PI / 50) / scale * image.heightpx + image.heightpx / 2 + yb), image.heightpx - 1), 0)][max(min((int)(3 * M * cos(i * 2 * PI / 50) / scale * image.widthpx + image.widthpx / 2 + xb), image.widthpx - 1), 0)] = vector<BYTE>{ (BYTE)0, (BYTE)255, (BYTE)0 };
	}

	for (auto& el : output)
	{
		image.pixarr[max(min((int)(-el[1] / scale * image.heightpx + image.heightpx / 2 + yb), image.heightpx-1), 0)][max(min((int)(el[0] / scale * image.widthpx + image.widthpx / 2 + xb), image.widthpx - 1), 0)] = vector<BYTE>{ (BYTE)255, (BYTE)255, (BYTE)255 };
	}
	
	//read_txt rt(pathout);

	//rt.print(output, 4);

	read_bmp rb(pathoutbmp);
	rb.print(image);

	

	return 0;
}