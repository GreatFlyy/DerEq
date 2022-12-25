#include "Matrix.h"
#include <math.h>

Matrix::Matrix()
{
	this->mtrx = new float* [1];
	this->mtrx[0] = new float[1];
	this->m = 1;
	this->n = 1;
	this->mtrx[0][0] = 0;
}

Matrix::Matrix(int m, float E)
{
	this->mtrx = new float* [m];
	this->m = m;
	this->n = m;
	for (int i = 0; i < m; i++)
	{
		this->mtrx[i] = new float[m];
		for (int j = 0; j < m; j++)
		{
			if (i == j)
				this->mtrx[i][j] = E;
			else
				this->mtrx[i][j] = 0;
		}
	}
}

Matrix::Matrix(int m, int n, int E)
{
	this->mtrx = new float* [m];
	this->m = m;
	this->n = n;
	for (int i = 0; i < m; i++)
	{
		this->mtrx[i] = new float[n];
		for (int j = 0; j < n; j++)
		{
			this->mtrx[i][j] = E;
		}
	}
}

Matrix::Matrix(Vec2* arr, int n)
{
	this->mtrx = new float* [2];
	this->m = 2;
	this->n = n;

	this->mtrx[0] = new float[n];
	this->mtrx[1] = new float[n];
	for (int j = 0; j < n; j++)
	{
		this->mtrx[0][j] = arr[j].x;
		this->mtrx[1][j] = arr[j].y;
	}
}

Matrix::Matrix(Vec3* arr, int n)
{
	this->mtrx = new float* [3];
	this->m = 3;
	this->n = n;
	this->mtrx[0] = new float[n];
	this->mtrx[1] = new float[n];
	this->mtrx[2] = new float[n];
	for (int j = 0; j < n; j++)
	{
		this->mtrx[0][j] = arr[j].x;
		this->mtrx[1][j] = arr[j].y;
		this->mtrx[2][j] = arr[j].z;
	}
}

Matrix::Matrix(Vec4* arr, int n)
{
	this->mtrx = new float* [4];
	this->m = 4;
	this->n = n;
	this->mtrx[0] = new float[n];
	this->mtrx[1] = new float[n];
	this->mtrx[2] = new float[n];
	this->mtrx[3] = new float[n];
	for (int j = 0; j < n; j++)
	{
		this->mtrx[0][j] = arr[j].x;
		this->mtrx[1][j] = arr[j].y;
		this->mtrx[2][j] = arr[j].z;
		this->mtrx[3][j] = arr[j].w;
	}
}

Matrix::Matrix(const Matrix& A)
{
	this->m = A.m;
	this->n = A.n;
	this->mtrx = new float* [m];
	for (int i = 0; i < m; i++)
	{
		this->mtrx[i] = new float[n];
		for (int j = 0; j < n; j++)
		{
			this->mtrx[i][j] = A.mtrx[i][j];
		}
	}
}

Matrix::Matrix(Vec2 i, Vec2 j)
{
	this->mtrx = new float*[2];
	this->m = 2;
	this->n = 2;
	for (int a = 0; a < 2; a++)
	{
		this->mtrx[a] = new float[2];
	}

	this->mtrx[0][0] = i.x; this->mtrx[0][1] = j.x;
	this->mtrx[1][0] = i.y; this->mtrx[1][1] = j.y;
}

Matrix::Matrix(Vec3 i, Vec3 j, Vec3 k)
{
	this->mtrx = new float* [3];
	this->m = 3;
	this->n = 3;
	for (int a = 0; a < 3; a++)
	{
		this->mtrx[a] = new float[3];
	}

	this->mtrx[0][0] = i.x; this->mtrx[0][1] = j.x; this->mtrx[0][2] = k.x;
	this->mtrx[1][0] = i.y; this->mtrx[1][1] = j.y; this->mtrx[1][2] = k.y;
	this->mtrx[2][0] = i.z; this->mtrx[2][1] = j.z; this->mtrx[2][2] = k.z;
}

Matrix::Matrix(Vec4 i, Vec4 j, Vec4 k, Vec4 w)
{
	this->mtrx = new float* [4];
	this->m = 4;
	this->n = 4;
	for (int a = 0; a < 4; a++)
	{
		this->mtrx[a] = new float[4];
	}

	this->mtrx[0][0] = i.x; this->mtrx[0][1] = j.x; this->mtrx[0][2] = k.x; this->mtrx[0][3] = w.x;
	this->mtrx[1][0] = i.y; this->mtrx[1][1] = j.y; this->mtrx[1][2] = k.y; this->mtrx[1][3] = w.y;
	this->mtrx[2][0] = i.z; this->mtrx[2][1] = j.z; this->mtrx[2][2] = k.z; this->mtrx[2][3] = w.z;
	this->mtrx[3][0] = i.w; this->mtrx[3][1] = j.w; this->mtrx[3][2] = k.w; this->mtrx[3][3] = w.w;
}

Matrix::Matrix(const Quat& q)
{
	this->mtrx = new float* [3];
	this->m = 3;
	this->n = 3;
	for (int a = 0; a < 3; a++)
	{
		this->mtrx[a] = new float[3];
	}

	this->mtrx[0][0] = 1 - 2 * q.j * q.j - 2 * q.k * q.k; this->mtrx[0][1] = 2 * q.i * q.j - 2 * q.k * q.r; this->mtrx[0][2] = 2 * q.k * q.i + 2 * q.j * q.r;
	this->mtrx[1][0] = 2 * q.i * q.j + 2 * q.k * q.r; this->mtrx[1][1] = 1 - 2 * q.i * q.i - 2 * q.k * q.k; this->mtrx[1][2] = 2 * q.j * q.k - 2 * q.i * q.r;
	this->mtrx[2][0] = 2 * q.i * q.k - 2 * q.j * q.r; this->mtrx[2][1] = 2 * q.j * q.k + 2 * q.i * q.r; this->mtrx[2][2] = 1 - 2 * q.i * q.i - 2 * q.j * q.j;
}

Matrix::Matrix(const Comp& z)
{
	this->mtrx = new float* [2];
	this->m = 2;
	this->n = 2;
	for (int a = 0; a < 3; a++)
	{
		this->mtrx[a] = new float[2];
	}

	this->mtrx[0][0] = z.r; this->mtrx[0][1] = -z.i;
	this->mtrx[1][0] = z.i; this->mtrx[1][1] = z.r;
}


Matrix::~Matrix()
{
	for (int i = 0; i < m; i++)
	{
		delete[] this->mtrx[i];
	}
	delete[] this->mtrx;
}

Matrix Matrix::operator=(const Matrix& A)
{
	this->m = A.m;
	this->n = A.n;
	this->mtrx = new float* [m];
	for (int i = 0; i < m; i++)
	{
		this->mtrx[i] = new float[n];
		for (int j = 0; j < n; j++)
		{
			this->mtrx[i][j] = A.mtrx[i][j];
		}
	}

	return *this;
}

Matrix Matrix::operator+(Matrix& B)const
{
	if (this->m == B.m && this->n == B.n)
	{
		Matrix result(B.m, B.n, 0);
		for (int i = 0; i < B.m; i++)
		{
			for (int j = 0; j < B.n; j++)
			{
				result.mtrx[i][j] = this->mtrx[i][j] + B.mtrx[i][j];
			}
		}
		return result;
	}
	else
		throw "no size";
}

Matrix Matrix::operator*(float a)const
{
	Matrix result(*this);
	for (int i = 0; i < this->m; i++)
	{
		for (int j = 0; j < this->n; j++)
		{
			result.mtrx[i][j] *= a;
		}
	}
	return result;
}

Matrix Matrix::operator*(const Matrix& B)const
{
	if (this->n == B.m)
	{
		Matrix result(this->m, B.n, 0);
		for (int i = 0; i < result.m; i++)
		{
			for (int j = 0; j < result.n; j++)
			{
				for (int k = 0; k < B.m; k++)
				{
					result.mtrx[i][j] += this->mtrx[i][k] * B.mtrx[k][j];
				}
			}
		}
		return result;
	}
	else
		throw "no size";
}

Matrix Matrix::operator/(float a)const
{
	return *this * (1 / a);
}

float* Matrix::operator[](int i) const
{
	return (this->mtrx[i - 1] - 1);
}

Matrix Matrix::operator-()const
{
	return (*this * (-1.0));
}

Matrix Matrix::operator-(const Matrix& B)const
{
	if (this->m == B.m && this->n == B.n)
	{
		Matrix result(B.m, B.n, 0);
		for (int i = 0; i < B.m; i++)
		{
			for (int j = 0; j < B.n; j++)
			{
				result.mtrx[i][j] = this->mtrx[i][j] - B.mtrx[i][j];
			}
		}
		return result;
	}
	else
		throw "no size";
}

Matrix Matrix::Transp()const
{
	Matrix result(this->n, this->m, 0);
	for (int i = 0; i < this->n; i++)
	{
		for (int j = 0; j < this->m; j++)
		{
			result.mtrx[i][j] = this->mtrx[j][i];
		}
	}
	return result;
}

float Matrix::minor(int I, int J)const
{
	if (this->m == this->n)
	{
		Matrix result(this->m - 1, this->n - 1, 0);

		I -= 1;
		J -= 1;

		for (int i = 0; i < m; i++)
		{
			if (i < I)
			{
				for (int j = 0; j < n; j++)
				{
					if (j < J)
					{
						result.mtrx[i][j] = this->mtrx[i][j];
					}
					else if (j > J)
					{
						result.mtrx[i][j - 1] = this->mtrx[i][j];
					}
				}
			}
			else if (i > I)
			{
				for (int j = 0; j < n; j++)
				{
					if (j < J)
					{
						result.mtrx[i - 1][j] = this->mtrx[i][j];
					}
					else if (j > J)
					{
						result.mtrx[i - 1][j - 1] = this->mtrx[i][j];
					}
				}
			}
		}
		return result.det();
	}
	else
		throw "no size";
}

float Matrix::algdop(int I, int J)const
{
	if ((I + J) % 2)
	{
		return -(this->minor(I, J));
	}
	else
	{
		return this->minor(I, J);
	}
}

float Matrix::det()const
{
	if (this->m == this->n)
	{
		if (this->m == 1)
		{
			return this->mtrx[0][0];
		}
		else if (this->m == 2)
		{
			return ((this->mtrx[0][0] * this->mtrx[1][1]) - (this->mtrx[0][1] * this->mtrx[1][0]));
		}
		else
		{
			float result = 0;
			for (int j = 0; j < this->n; j++)
			{
				result += this->mtrx[0][j] * this->algdop(1, j + 1);
			}
			return result;
		}
	}
	else
		throw "no size";
}

Matrix Matrix::matalgdop()const
{
	Matrix result(this->m, 0);
	for (int i = 0; i < this->m; i++)
	{
		for (int j = 0; j < this->m; j++)
		{
			result.mtrx[i][j] = this->algdop(i + 1, j + 1);
		}
	}
	return result;
}

Matrix Matrix::inv()const
{
	if (this->det() != 0)
	{
		return ((this->matalgdop()).Transp() * (1 / this->det()));
	}
	else
	{
		throw "virozhden";
	}
}

float* Matrix::arr()const
{
	float* arr = new float[this->m * this->n];
	for (int i = 0; i < this->m; i++)
	{
		for (int j = 0; j < this->n; j++)
		{
			arr[i * this->n + j] = this->mtrx[i][j];
		}
	}
	return arr;
}


Vec2::Vec2()
{
	this->x = 0;
	this->y = 0;
}

Vec2::Vec2(float x, float y)
{
	this->x = x;
	this->y = y;
}

Vec2::Vec2(const Vec2& A)
{
	this->x = A.x;
	this->y = A.y;
}

Vec2::Vec2(const Comp& z)
{
	this->x = z.r;
	this->y = z.i;
}

Vec2 Vec2::operator=(const Vec2& a)
{
	this->x = a.x;
	this->y = a.y;
	return *this;
}

Vec2 Vec2::operator+(const Vec2& a)const
{
	return Vec2(this->x + a.x, this->y + a.y);
}

Vec2 Vec2::operator*(float c) const
{
	return Vec2(this->x * c, this->y * c);
}

Vec2 Vec2::operator/(float c)const
{
	return Vec2(this->x / c, this->y / c);
}

Vec2 Vec2::operator-()const
{
	return *this * (-1);
}

Vec2 Vec2::operator-(const Vec2& a)const
{
	return Vec2(this->x - a.x, this->y - a.y);
}

float Vec2::lenght() const
{
	return (float)sqrt(this->x * this->x + this->y * this->y);
}

Vec2 Vec2::Normalise()const
{
	if (this->lenght())
	{
		return (*this * (1 / this->lenght()));
	}
	else
	{
		return *this;
	}
}

Matrix Vec2::mtrx()const
{
	Matrix result(2, 1, 0);
	result[1][1] = this->x;
	result[2][1] = this->y;
	return result;
}



Vec3::Vec3()
{
	this->x = 0;
	this->y = 0;
	this->z = 0;
}

Vec3::Vec3(float x, float y, float z)
{
	this->x = x;
	this->y = y;
	this->z = z;
}

Vec3::Vec3(const Vec3& a)
{
	this->x = a.x;
	this->y = a.y;
	this->z = a.z;
}

Vec3::Vec3(const Quat& q)
{
	this->x = q.i;
	this->y = q.j;
	this->z = q.k;
}

Vec3 Vec3::operator=(const Vec3& a)
{
	this->x = a.x;
	this->y = a.y;
	this->z = a.z;
	return *this;
}

Vec3 Vec3::operator+(const Vec3& a)const
{
	return Vec3(this->x + a.x, this->y + a.y, this->z + a.z);
}

Vec3 Vec3::operator*(float c)const
{
	return Vec3(this->x * c, this->y * c, this->z * c);
}

Vec3 Vec3::operator/(float c)const
{
	return Vec3(this->x / c, this->y / c, this->z / c);
}

Vec3 Vec3::operator-()const
{
	return *this * (-1);
}

Vec3 Vec3::operator-(const Vec3& a)const
{
	return Vec3(this->x - a.x, this->y - a.y, this->z - a.z);
}

float Vec3::lenght() const
{
	return (float)sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
}

Vec3 Vec3::Normalise() const
{
	if (this->lenght())
	{
		return (*this * (1 / this->lenght()));
	}
	else
	{
		return *this;
	}
}

Matrix Vec3::mtrx()const
{
	Matrix result(3, 1, 0);
	result[1][1] = this->x;
	result[2][1] = this->y;
	result[3][1] = this->z;
	return result;
}



Vec4::Vec4()
{
	this->x = 0;
	this->y = 0;
	this->z = 0;
	this->w = 0;
}

Vec4::Vec4(float x, float y, float z, float w)
{
	this->x = x;
	this->y = y;
	this->z = z;
	this->w = w;
}

Vec4::Vec4(Vec3 vec, float w)
{
	this->x = vec.x;
	this->y = vec.y;
	this->z = vec.z;
	this->w = w;
}

Vec4::Vec4(const Vec4& a)
{
	this->x = a.x;
	this->y = a.y;
	this->z = a.z;
	this->w = a.w;
}

Vec4::Vec4(const Quat& q)
{
	this->w = q.r;
	this->x = q.i;
	this->y = q.j;
	this->z = q.k;
}

Vec4 Vec4::operator=(const Vec4& a)
{
	this->x = a.x;
	this->y = a.y;
	this->z = a.z;
	this->w = a.w;
	return *this;
}

Vec4 Vec4::operator+(const Vec4& a)const
{
	return Vec4(this->x + a.x, this->y + a.y, this->z + a.z, this->w + a.w);
}

Vec4 Vec4::operator*(float c) const
{
	return Vec4(this->x * c, this->y * c, this->z * c, this->w * c);
}

Vec4 Vec4::operator/(float c)const
{
	return Vec4(this->x / c, this->y / c, this->z / c, this->w / c);
}

Vec4 Vec4::operator-()const
{
	return *this * (-1);
}

Vec4 Vec4::operator-(const Vec4& a)const
{
	return Vec4(this->x - a.x, this->y - a.y, this->z - a.z, this->w + w);
}

float Vec4::lenght() const
{
	return (float)sqrt(this->x * this->x + this->y * this->y + this->z * this->z + this->w * this->w);
}

Vec4 Vec4::Normalise() const
{
	if (this->lenght())
	{
		return (*this * (1 / this->lenght()));
	}
	else
	{
		return *this;
	}
}

Matrix Vec4::mtrx()const
{
	Matrix result(4, 1, 0);
	result[1][1] = this->x;
	result[2][1] = this->y;
	result[3][1] = this->z;
	result[4][1] = this->w;
	return result;
}



Vec2 Matrix::operator*(const Vec2& v)const
{
	return Vec2((*this * v.mtrx())[1][1], (*this * v.mtrx())[2][1]);
}

Vec3 Matrix::operator*(const Vec3& v)const
{
	return Vec3((*this * v.mtrx())[1][1], (*this * v.mtrx())[2][1], (*this * v.mtrx())[3][1]);
}

Vec4 Matrix::operator*(const Vec4& v)const
{
	return Vec4((*this * v.mtrx())[1][1], (*this * v.mtrx())[2][1], (*this * v.mtrx())[3][1], (*this * v.mtrx())[4][1]);
}

Matrix::Matrix(const Vec3& axis, float angle)
{
	this->mtrx = new float* [3];
	this->m = 3;
	this->n = 3;
	for (int i = 0; i < m; i++)
	{
		this->mtrx[i] = new float[3];		
	}

	this->mtrx[0][0] = cos(angle) + (1 - cos(angle)) * axis.x * axis.x;
	this->mtrx[0][1] = (1 - cos(angle)) * axis.x * axis.y - sin(angle) * axis.z;
	this->mtrx[0][2] = (1 - cos(angle)) * axis.x * axis.z + sin(angle) * axis.y;

	this->mtrx[1][0] = (1 - cos(angle)) * axis.y * axis.x + sin(angle) * axis.z;
	this->mtrx[1][1] = cos(angle) + (1 - cos(angle)) * axis.y * axis.y;
	this->mtrx[1][2] = (1 - cos(angle)) * axis.y * axis.z - sin(angle) * axis.x;

	this->mtrx[2][0] = (1 - cos(angle)) * axis.z * axis.x - sin(angle) * axis.y;
	this->mtrx[2][1] = (1 - cos(angle)) * axis.z * axis.y + sin(angle) * axis.x;
	this->mtrx[2][2] = cos(angle) + (1 - cos(angle)) * axis.z * axis.z;
}

Matrix::Matrix(const Vec3& axisang)
{
	Vec3 axis = axisang.Normalise();
	float angle = axisang.lenght();

	this->mtrx = new float* [3];
	this->m = 3;
	this->n = 3;
	for (int i = 0; i < m; i++)
	{
		this->mtrx[i] = new float[3];
	}

	this->mtrx[0][0] = cos(angle) + (1 - cos(angle)) * axis.x * axis.x;
	this->mtrx[0][1] = (1 - cos(angle)) * axis.x * axis.y - sin(angle) * axis.z;
	this->mtrx[0][2] = (1 - cos(angle)) * axis.x * axis.z + sin(angle) * axis.y;

	this->mtrx[1][0] = (1 - cos(angle)) * axis.y * axis.z + sin(angle) * axis.z;
	this->mtrx[1][1] = cos(angle) + (1 - cos(angle)) * axis.y * axis.y;
	this->mtrx[1][2] = (1 - cos(angle)) * axis.y * axis.z - sin(angle) * axis.x;

	this->mtrx[2][0] = (1 - cos(angle)) * axis.z * axis.x - sin(angle) * axis.y;
	this->mtrx[2][1] = (1 - cos(angle)) * axis.z * axis.y + sin(angle) * axis.x;
	this->mtrx[2][2] = cos(angle) + (1 - cos(angle)) * axis.z * axis.z;
}

float dtprdct(const Vec2& a, const Vec2& b)
{
	return(a.x * b.x + a.y * b.y);
}

float dtprdct(const Vec3& a, const Vec3& b)
{
	return(a.x * b.x + a.y * b.y + a.z * b.z);
}

float dtprdct(const Vec4& a, const Vec4& b)
{
	return (a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w);
}

Vec3 crssprdct(const Vec3& a, const Vec3& b)
{
	return Vec3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}



Comp::Comp()
{
	this->r = 0;
	this->i = 0;
}

Comp::Comp(double r)
{
	this->r = r;
	this->i = 0;
}

Comp::Comp(double r, double i)
{
	this->r = r;
	this->i = i;
}

Comp::Comp(Vec2& v)
{
	this->r = v.x;
	this->i = v.y;
}

Comp Comp::operator=(const Comp& z)
{
	this->r = z.r;
	this->i = z.i;
	return *this;
}
Comp Comp::operator=(const double a)
{
	this->r = a;
	this->i = 0;
	return *this;
}

Comp Comp::operator+(const Comp& z)const
{
	return Comp(this->r + z.r, this->i + z.i);
}

Comp Comp::operator*(const Comp& z)const
{
	return Comp(this->r * z.r - this->i * z.i, this->r * z.i + this->i * z.r);
}
Comp Comp::operator*(const double a)const
{
	return Comp(this->r * a, this->i * a);
}

Comp Comp::operator-()const
{
	return Comp(-this->r, -this->i);
}
Comp Comp::operator-(const Comp& z)const
{
	return Comp(this->r - z.r, this->i - z.i);
}

Comp Comp::operator/(const Comp& z)const
{
	return(*this * this->inv());
}
Comp Comp::operator/(const double a)const
{
	return Comp(this->r / a, this->i / a);
}

Comp Comp::conj()const
{
	return Comp(this->r, -this->i);
}

double Comp::normsq()const
{
	return (this->r * this->r + this->i * this->i);
}
double Comp::lenght()const
{
	return sqrt(this->normsq());
}

Comp Comp::inv()const
{
	return this->conj() / this->normsq();
}

Comp Comp::Normalize()const
{
	if (!this->lenght())
		return *this / this->lenght();
	else
		return *this;
}

double Comp::phase()const
{
	return atan2(this->i, this->r);
}



Comp exp(const Comp& z)
{
	return Comp(cos(z.i), sin(z.i)) * exp(z.r);
}

Comp pow(const Comp& z, double power)
{
	double r = z.lenght();
	double ang = atan2(z.i, z.r);

	r = pow(r, power);
	ang = ang * power;

	return Comp(r * cos(ang), r * sin(ang));
}


Quat::Quat()
{
	this->r = 0;
	this->i = 0;
	this->j = 0;
	this->k = 0;
}

Quat::Quat(float r)
{
	this->r = r;
	this->i = 0;
	this->j = 0;
	this->k = 0;
}

Quat::Quat(const Vec3& vec)
{
	this->r = 0;
	this->i = vec.x;
	this->j = vec.y;
	this->k = vec.z;
}

Quat::Quat(const Vec2& vec)
{
	this->r = 0;
	this->i = vec.x;
	this->j = vec.y;
	this->k = 0;
}

Quat::Quat(float r, const Vec3& vec)
{
	this->r = r;
	this->i = vec.x;
	this->j = vec.y;
	this->k = vec.z;
}

Quat::Quat(float r, float i, float j, float k)
{
	this->r = r;
	this->i = i;
	this->j = j;
	this->k = k;
}

Quat::Quat(const Vec4& q)
{
	this->r = q.w;
	this->i = q.x;
	this->j = q.y;
	this->k = q.z;
}

Quat::Quat(const Vec3& axis, float angle)
{
	Vec3 ax = axis.Normalise() * sin(angle / 2);
	this->r = cos(angle / 2);
	this->i = ax.x;
	this->j = ax.y;
	this->k = ax.z;
}

Quat Quat::operator=(const Quat& q)
{
	this->r = q.r;
	this->i = q.i;
	this->j = q.j;
	this->k = q.k;
	return *this;
}

Quat Quat::operator+(const Quat& q)const
{
	return Quat(this->r + q.r, this->i + q.i, this->j + q.j, this->k + q.k);
}

Vec3 Quat::vec()const
{
	return Vec3(this->i, this->j, this->k);
}

Quat Quat::operator*(const Quat& q)const
{
	float rr;
	Vec3 rvec;
	rr = this->r * q.r - dtprdct(this->vec(), q.vec());
	rvec = this->vec() * q.r + q.vec() * this->r + crssprdct(this->vec(), q.vec());
	return Quat(rr, rvec);
}

Quat Quat::operator*(float a)const
{
	return Quat(this->r * a, this->i * a, this->j * a, this->k * a);
}

Quat Quat::operator/(float a)const
{
	return Quat(this->r / a, this->i / a, this->j / a, this->k / a);
}

Quat Quat::operator/(const Quat& q)const
{
	return *this * this->inv();
}

Quat Quat::operator-()const
{
	return Quat(-this->r, -this->i, -this->j, -this->k);
}

Quat Quat::operator-(const Quat& q)const
{
	return Quat(this->r - q.r, this->i - q.i, this->j - q.j, this->k - q.k);
}

Quat Quat::conj()const
{
	return Quat(this->r, -this->i, -this->j, -this->k);
}

float Quat::normsq()const
{
	return(this->r * this->r + this->i * this->i + this->j * this->j + this->k * this->k);
}

float Quat::lenght()const
{
	return sqrt(this->normsq());
}

Quat Quat::Normalise()const
{
	if (this->normsq())
	{
		return *this * (1 / this->lenght());
	}
	else
	{
		return *this;
	}
}

Quat Quat::inv()const
{
	return this->conj() / this->normsq();
}

