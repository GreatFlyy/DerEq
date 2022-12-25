#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "Matrix.h"

using namespace std;

typedef unsigned char BYTE;

class read_txt
{
public:
	fstream fs;

	read_txt(const string& path);
	read_txt(const string& path, int);
	void open(const string& path);
	void open(const string& path, int);

	~read_txt();
	void close();
	
	vector<vector<double>> read_txt_rows(int nrows, int ncols);
	vector<vector<double>> read_txt_rows();

	vector<vector<double>> read_txt_cols(int ncols, int nrows);
	vector<vector<double>> read_txt_cols(int ncols);

	void print(const vector<vector<double>> &text, int n);

private:

	void print_double(double a, int n);
};

class read_dat
{
public:
	fstream fs;

	read_dat(const string& path);

	~read_dat();
	void close();

	vector<unsigned char> read_char();
	vector<int> read_int();

};

class bmp
{
public:
	int widthpx;
	int heightpx; //>0 bottom->top; <0 top->bottom

	int bpp;

	vector<vector<vector<BYTE>>> pixarr;

	bmp(int width, int height, int bpp);
};

int bytes2int(vector<char>& x);
int bytes2int(vector<char>& x, int beg, int end);

void int2bytes(vector<char>& file, int x);
void int2bytes(vector<char>& file, int beg, int end, int x);

class read_bmp
{
public:
	fstream fs;

	read_bmp(const string& path);

	~read_bmp();
	void close();

	bmp read();

	void print(bmp& image);
};

