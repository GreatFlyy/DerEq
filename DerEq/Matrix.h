#pragma once
#include<math.h>

class Quat;
class Comp;

class Vec2;
class Vec3;
class Vec4;

class Matrix
{
public:
	int m;
	int n;
	float** mtrx;

	Matrix();

	Matrix(int m, float E);

	Matrix(int m, int n, int E);

	Matrix(Vec2* arr, int n);
	Matrix(Vec3* arr, int n);
	Matrix(Vec4* arr, int n);

	Matrix(const Matrix& A);

	Matrix(Vec2 i, Vec2 j);
	Matrix(Vec3 i, Vec3 j, Vec3 k);
	Matrix(Vec4 i, Vec4 j, Vec4 k, Vec4 w);

	Matrix(const Vec3& axis, float angle);
	Matrix(const Vec3& axisang);

	Matrix(const Quat& q);
	Matrix(const Comp& z);

	~Matrix();

	Matrix operator=(const Matrix& A);

	Matrix operator+(Matrix& B)const;

	Matrix operator*(float a)const;
	Matrix operator*(const Matrix& B)const;

	Matrix operator/(float a)const;

	float* operator[](int i) const;

	Matrix operator-()const;
	Matrix operator-(const Matrix& B)const;

	Matrix Transp()const;

	float minor(int I, int J)const;

	float algdop(int I, int J)const;

	float det()const;

	Matrix matalgdop()const;

	Matrix inv()const;

	Vec2 operator*(const Vec2& v)const;
	
	Vec3 operator*(const Vec3& v)const;	

	Vec4 operator*(const Vec4& v)const;

	float* arr()const;
};

class Vec2
{
public:
	float x;
	float y;	

	Vec2();

	Vec2(float x, float y);

	Vec2(const Vec2& A);

	Vec2(const Comp& z);

	Vec2 operator=(const Vec2& a);

	Vec2 operator+(const Vec2& a)const;

	Vec2 operator*(float c)const;

	Vec2 operator/(float c)const;

	Vec2 operator-()const;
	Vec2 operator-(const Vec2& a)const;

	float lenght()const;

	Vec2 Normalise()const;

	Matrix mtrx()const;
};

class Vec3
{
public:
	float x;
	float y;
	float z;	

	Vec3();

	Vec3(float x, float y, float z);

	Vec3(const Vec3& a);

	Vec3(const Quat& q);

	Vec3 operator=(const Vec3& a);

	Vec3 operator+(const Vec3& a)const;

	Vec3 operator*(float c)const;

	Vec3 operator/(float c)const;

	Vec3 operator-()const;
	Vec3 operator-(const Vec3& a)const;

	float lenght()const;

	Vec3 Normalise()const;

	Matrix mtrx()const;
};

class Vec4
{
public:
	float x;
	float y;
	float z;
	float w;

	Vec4();

	Vec4(float x, float y, float z, float w);

	Vec4(Vec3 vec, float w);

	Vec4(const Vec4& a);

	Vec4(const Quat& q);

	Vec4 operator=(const Vec4& a);

	Vec4 operator+(const Vec4& a)const;

	Vec4 operator*(float c)const;

	Vec4 operator/(float c)const;

	Vec4 operator-()const;
	Vec4 operator-(const Vec4& a)const;

	float lenght()const;

	Vec4 Normalise()const;

	Matrix mtrx()const;
};

float dtprdct(const Vec2& a, const Vec2& b);

float dtprdct(const Vec3& a, const Vec3& b);

float dtprdct(const Vec4& a, const Vec4& b);

Vec3 crssprdct(const Vec3& a, const Vec3& b);

class Comp
{
public:
	double r;
	double i;

	Comp();

	Comp(double r);

	Comp(double r, double i);

	Comp(Vec2& v);

	Comp operator=(const Comp& z);
	Comp operator=(const double a);

	Comp operator+(const Comp& z)const;

	Comp operator*(const Comp& z)const;
	Comp operator*(const double a)const;

	Comp operator-()const;
	Comp operator-(const Comp& z)const;

	Comp operator/(const Comp& z)const;
	Comp operator/(const double a)const;

	Comp conj()const;

	double normsq()const;
	double lenght()const;

	Comp inv()const;

	Comp Normalize()const;

	double phase()const;
	

};

Comp exp(const Comp& z);

class Quat
{
public:
	float r;
	float i;
	float j;
	float k;

	Quat();

	Quat(float r);

	Quat(const Vec3& vec);

	Quat(const Vec2& vec);

	Quat(float r, const Vec3& vec);

	Quat(float r, float i, float j, float k);

	Quat(const Vec4& q);

	Quat(const Vec3& axis, float angle);

	Quat operator=(const Quat& q);

	Quat operator+(const Quat& q)const;

	Vec3 vec()const;

	Quat operator*(const Quat& q)const;
	Quat operator*(float a)const;

	Quat operator/(const Quat& q)const;
	Quat operator/(float a)const;

	Quat operator-()const;
	Quat operator-(const Quat& q)const;	

	Quat conj()const;

	float normsq()const;

	float lenght()const;

	Quat Normalise()const;

	Quat inv()const;
};

template <typename T>
T round(T a, int n)
{
	a = a * pow(10, n);
	a = round(a);
	a = a / pow(10, n);
	return a;
}
