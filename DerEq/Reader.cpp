#include "Reader.h"


read_txt::read_txt(const string& path)
{
	fs.open(path, std::fstream::in | std::fstream::out);
	if (!fs.is_open())
	{
		std::cout << "read_txt_open_error" << std::endl;
	}
}

read_txt::read_txt(const string& path, int)
{
	fs.open(path, std::fstream::in | std::fstream::out | std::fstream::app);
	if (!fs.is_open())
	{
		std::cout << "read_txt_open_error" << std::endl;
	}
}


read_txt::~read_txt()
{
	fs.close();
}


void read_txt::open(const string& path)
{
	fs.open(path, std::fstream::in | std::fstream::out);
	if (!fs.is_open())
	{
		std::cout << "read_txt_open_error" << std::endl;
	}
}

void read_txt::open(const string& path, int)
{
	fs.open(path, std::fstream::in | std::fstream::out | std::fstream::app);
	if (!fs.is_open())
	{
		std::cout << "read_txt_open_error" << std::endl;
	}
}


void read_txt::close()
{
	fs.close();
}

vector<vector<double>> read_txt::read_txt_rows(int nrows, int ncols)
{
	vector<vector<double>> result;
	result.reserve(nrows);
	result.push_back(vector<double>(ncols, 0));

	int iColVec = 0;
	int iRowVec = 0;	
	int isign = 1;
	int iPart = 1;
	bool isComma = 0;

	char box;
	//fs.get(box);
	//fs.get(box);
	//fs.get(box);
	while (fs.get(box))
	{
		if (box != '\0' && box != '\r')
		{
			if (box != '\t' && box != '\n' && box != '-')
			{
				if (box != ',')
				{
					if (!isComma)
					{
						result[iColVec][iRowVec] = result[iColVec][iRowVec] * 10 + (float)(box - '0');
					}
					else
					{
						result[iColVec][iRowVec] = result[iColVec][iRowVec] + (float)(box - '0') / pow(10, iPart);
						iPart++;
					}
				}
				else
				{
					isComma = 1;
				}
			}
			else if (box == '-')
			{
				isign = -1;
			}
			else if (box == '\n')
			{
				result[iColVec][iRowVec] *= isign;
				isign = 1;
				++iColVec;
				result.push_back(vector<double>(1, 0));
				result[iColVec].reserve(ncols);
				iRowVec = 0;
				iPart = 1;
				isComma = 0;
			}
			else if (box == '\t')
			{
				result[iColVec][iRowVec] *= isign;
				isign = 1;
				++iRowVec;
				result[iColVec].push_back(0);
				iPart = 1;
				isComma = 0;
			}
		}
	}

	result.pop_back();

	return result;
}

vector<vector<double>> read_txt::read_txt_rows()
{
	int nrows = 1;
	int ncols = 1;
	vector<vector<double>> result;
	result.reserve(nrows);
	result.push_back(vector<double>(ncols, 0));

	int iColVec = 0;
	int iRowVec = 0;
	int iPart = 1;
	int isign = 1;
	bool isComma = 0;

	char box;
	//fs.get(box);
	//fs.get(box);
	while (fs.get(box))
	{
		if (box != '\0' && box != '\r')
		{
			if (box != '\t' && box != '\n' && box != '-')
			{
				if (box != ',')
				{
					if (!isComma)
					{
						result[iColVec][iRowVec] = result[iColVec][iRowVec] * 10 + (float)(box - '0');
					}
					else
					{
						result[iColVec][iRowVec] = result[iColVec][iRowVec] + (float)(box - '0') / pow(10, iPart);
						iPart++;
					}
				}
				else
				{
					isComma = 1;
				}
			}
			else if (box == '-')
			{
				isign = -1;
			}
			else if (box == '\n')
			{
				result[iColVec][iRowVec] *= isign;
				isign = 1;
				++iColVec;
				result.push_back(vector<double>(1, 0));
				result[iColVec].reserve(ncols);
				iRowVec = 0;
				iPart = 1;
				isComma = 0;
			}
			else if (box == '\t')
			{
				result[iColVec][iRowVec] *= isign;
				isign = 1;
				++iRowVec;
				result[iColVec].push_back(0);
				iPart = 1;
				isComma = 0;
			}
		}
	}

	result.pop_back();

	return result;
}


vector<vector<double>> read_txt::read_txt_cols(int ncols, int nrows)
{
	vector<vector<double>> result(ncols, vector<double>(nrows, 0));
	

	int iColVec = 0;
	int iRowVec = 0;
	int isign = 1;
	int iPart = 1;
	bool isComma = 0;

	char box;
	//fs.get(box);
	//fs.get(box);
	//fs.get(box);
	while (fs.get(box))
	{
		if (box != '\0' && box != '\r')
		{
			if (box != '\t' && box != '\n' && box != '-')
			{
				if (box != ',')
				{
					if (!isComma)
					{
						result[iRowVec][iColVec] = result[iRowVec][iColVec] * 10 + (float)(box - '0');
					}
					else
					{
						result[iRowVec][iColVec] = result[iRowVec][iColVec] + (float)(box - '0') / pow(10, iPart);
						iPart++;
					}
				}
				else
				{
					isComma = 1;
				}
			}
			else if (box == '-')
			{
				isign = -1;
			}
			else if (box == '\n')
			{
				result[iRowVec][iColVec] *= isign;
				isign = 1;
				++iColVec;
				iRowVec = 0;
				iPart = 1;
				isComma = 0;

				
			}
			else if (box == '\t')
			{
				result[iRowVec][iColVec] *= isign;
				isign = 1;
				++iRowVec;
				iPart = 1;
				isComma = 0;
			}
		}
	}


	return result;
}

vector<vector<double>> read_txt::read_txt_cols(int ncols)
{
	vector<vector<double>> result(ncols, vector<double>(1, 0));


	int iColVec = 0;
	int iRowVec = 0;
	int isign = 1;
	int iPart = 1;
	bool isComma = 0;

	char box;
	//fs.get(box);
	//fs.get(box);
	//fs.get(box);
	while (fs.get(box))
	{
		if (box != '\0' && box != '\r')
		{
			if (box != '\t' && box != '\n' && box != '-')
			{
				if (box != ',')
				{
					if (!isComma)
					{
						result[iRowVec][iColVec] = result[iRowVec][iColVec] * 10 + (float)(box - '0');
					}
					else
					{
						result[iRowVec][iColVec] = result[iRowVec][iColVec] + (float)(box - '0') / pow(10, iPart);
						iPart++;
					}
				}
				else
				{
					isComma = 1;
				}
			}
			else if (box == '-')
			{
				isign = -1;
			}
			else if (box == '\n')
			{
				result[iRowVec][iColVec] *= isign;
				isign = 1;
				++iColVec;
				iRowVec = 0;
				iPart = 1;
				isComma = 0;

				for (auto &el : result)
				{
					el.push_back(0);
				}
			}
			else if (box == '\t')
			{
				result[iRowVec][iColVec] *= isign;
				isign = 1;
				++iRowVec;
				iPart = 1;
				isComma = 0;
			}
		}
	}

	for (auto& el : result)
	{
		el.pop_back();
	}

	return result;
}


//void read_vecs::print_double(double a, int n)
//{
//	a = a * pow(10, n);
//	a = round(a);
//	vector<int> num;
//	int i = 0;
//	while (a != 0)
//	{
//		num.push_back((int)a % 10);
//		a = (a - num[i]) / 10;
//		++i;
//	}
//	for (int k = num.size() - 1; k >= 0; --k)
//	{
//		if (k != n - 1)
//		{
//			fs << num[k];
//		}
//		else
//		{
//			fs << ",";
//			fs << num[k];
//		}
//	}
//}

void read_txt::print_double(double a, int n)
{
	string num;
	num = to_string(a);
	for (auto& el : num)
	{
		if (el == '.')
			el = ',';
	}
	fs << num;
}

void read_txt::print(const vector<vector<double>> &text, int n)
{
	for (auto i = text.begin(); i != text.end(); ++i)
	{
		print_double((*(*i).begin()), n);
		for (auto j = (*i).begin()+1; j != (*i).end(); ++j)
		{
			fs << "\t";
			print_double(*j, n);			
		}
		fs <<"\n";
	}
}






read_dat::read_dat(const string& path)
{
	fs.open(path, std::fstream::in | std::ios::binary);
	if (!fs.is_open())
	{
		std::cout << "read_dat_open_error" << std::endl;
	}
}

read_dat::~read_dat()
{
	fs.close();
}
void read_dat::close()
{
	fs.close();
}

vector<unsigned char> read_dat::read_char()
{
	vector<char> result;

	std::streampos fsSize;

	fs.seekg(0, std::ios::end);
	fsSize = fs.tellg();
	fs.seekg(0, std::ios::beg);

	vector<unsigned char> res(fsSize);
	fs.read((char*)&res[0], fsSize);

	/*unsigned char box;

	while (fs.get(box))
	{
		result.push_back(box);
	}*/

	return res;
}




bmp::bmp(int widthpx, int heightpx, int bpp)
{
	this->widthpx = widthpx;
	this->heightpx = heightpx;
	this->bpp = bpp;

	this->pixarr = vector<vector<vector<BYTE>>>(heightpx, vector<vector<BYTE>>(widthpx, vector<BYTE>(bpp / 8, 0)));
}



read_bmp::read_bmp(const string& path)
{
	fs.open(path, std::fstream::in | std::fstream::out | std::ios::binary);
	if (!fs.is_open())
	{
		std::cout << "read_bmp_open_error" << std::endl;
	}
}

read_bmp::~read_bmp()
{
	fs.close();
}

void read_bmp::close()
{
	fs.close();
}


int bytes2int(vector<char>& x)
{
	int res = 0;
	for (int i = 0; i < x.size(); i++)
	{
		res += ((int)x[i]) << (i * 8);
	}
	return res;
}
int bytes2int(vector<char>& x, int beg, int end)
{
	int res = 0;
	for (int i = 0; i <= end - beg; i++)
	{
		res += ((int)(BYTE)x[beg+i]) << (i * 8);
	}
	return res;
}


bmp read_bmp::read()
{
	std::streampos fsize;

	fs.seekg(0, std::ios::end);
	fsize = fs.tellg();
	fs.seekg(0, std::ios::beg);

	vector<char> file(fsize);

	fs.read(&file[0], fsize);

	int dibheadersize = bytes2int(file, 14, 17);

	int width = bytes2int(file, 18, 21);
	int height = bytes2int(file, 22, 25);
	int bpp = bytes2int(file, 28, 29);
	int imsize = bytes2int(file, 34, 37);

	bmp image(width, height, bpp);
	//int RowSize = ndatabytes + (ndatabytes % 4 ? 4 - ndatabytes % 4 : 0);
	int RowSize = ((int)((width * bpp + 31) / 32)) * 4;

	if (height > 0)
	{
		int nrow = 0;
		for (int j = image.heightpx - 1; j >= 0; j--)
		{
			for (int i = 0; i < image.widthpx; i++)
			{
				image.pixarr[j][i][0] = (BYTE)file[14 + dibheadersize + nrow * RowSize + i * bpp / 8 + 2];
				image.pixarr[j][i][1] = (BYTE)file[14 + dibheadersize + nrow * RowSize + i * bpp / 8 + 1];
				image.pixarr[j][i][2] = (BYTE)file[14 + dibheadersize + nrow * RowSize + i * bpp / 8 + 0];
			}
			nrow++;
		}
	}
	else
	{
		for (int j = 0; j < -image.heightpx; j++)
		{
			for (int i = 0; i < image.widthpx; i++)
			{
				image.pixarr[j][i][0] = (BYTE)file[14 + dibheadersize + j * RowSize + i * bpp / 8 + 2];
				image.pixarr[j][i][1] = (BYTE)file[14 + dibheadersize + j * RowSize + i * bpp / 8 + 1];
				image.pixarr[j][i][2] = (BYTE)file[14 + dibheadersize + j * RowSize + i * bpp / 8 + 0];
			}
		}
	}	

	return image;
}

void int2bytes(vector<char>& file, int x)
{
	for (int i = 0; i < file.size(); i++)
	{
		file[i] = x % 256;
		x = x >> 8;
	}
}
void int2bytes(vector<char>& file, int beg, int end, int x)
{
	for (int i = beg; i <= end; i++)
	{
		file[i] = x % 256;
		x = x >> 8;
	}
}

void read_bmp::print(bmp& image)
{
	vector<char> file(54 + (((int)((image.widthpx * image.bpp + 31)/ 32)) * 4) * image.heightpx, 0);

	//header
	file[0] = 'B'; file[1] = 'M';
	int2bytes(file, 2, 5, file.size());
	int2bytes(file, 10, 13, 54);

	//dib header
	int2bytes(file, 14, 17, 40);
	int2bytes(file, 18, 21, image.widthpx);
	int2bytes(file, 22, 25, image.heightpx);
	int2bytes(file, 26, 27, 1);
	int2bytes(file, 28, 29, image.bpp);
	int2bytes(file, 34, 37, file.size() - 54);
	//int2bytes(file, 38, 41, 2835);
	//int2bytes(file, 42, 45, 2835); //resolution in px/meter

	//data
	int RowSize = ((int)((image.widthpx * image.bpp + 31)/ 32)) * 4;

	if (image.heightpx > 0)
	{
		int nrow = 0;
		for (int j = image.heightpx - 1; j >= 0; j--)
		{
			for (int i = 0; i < image.widthpx; i++)
			{
				file[54 + nrow * RowSize + i * image.bpp / 8 + 2] = image.pixarr[j][i][0];
				file[54 + nrow * RowSize + i * image.bpp / 8 + 1] = image.pixarr[j][i][1];
				file[54 + nrow * RowSize + i * image.bpp / 8 + 0] = image.pixarr[j][i][2];
			}
			nrow++;
		}
	}
	else
	{
		for (int j = 0; j < -image.heightpx; j++)
		{
			for (int i = 0; i < image.widthpx; i++)
			{
				file[54 + j * RowSize + i * image.bpp / 8 + 2] = image.pixarr[j][i][0];
				file[54 + j * RowSize + i * image.bpp / 8 + 1] = image.pixarr[j][i][1];
				file[54 + j * RowSize + i * image.bpp / 8 + 0] = image.pixarr[j][i][2];
			}
		}
	}

	fs.write(&file[0], file.size());
}