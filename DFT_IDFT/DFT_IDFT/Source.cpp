#include "vector.h"
#include "DFT.h"
#include <chrono>
#include <cmath>
#include <math.h>
#define Eps 0.001

int Sign(double Val)
{
	if (abs(Val) < Eps)
		return 0;
	if (Val > 0.0)
		return 1;
	else
		return -1;
}

int main()
{
	VectorZ vector(512);
	DiscreteTransformation DFT;
	vector.type = cos2pwj;
	vector.w = 100;
	vector.A = -1;
	vector.A0 = -1;
	vector.fi = -PI;
	vector.FillVector();
	


	//auto start = std::chrono::high_resolution_clock::now();

	//DFT.RealizeDFT(vector.Z, vector.DFT_Data);
	//auto end = std::chrono::high_resolution_clock::now();
	//std::chrono::duration<double> diff = end - start;
	//std::cout << "Time DFT " << diff.count() << std::endl;

	//auto start1 = std::chrono::high_resolution_clock::now();

	DFT.RealizeFFT(vector.Z, vector.FFT_Data);
	

	//auto end1 = std::chrono::high_resolution_clock::now();
	//std::chrono::duration<double> diff1 = end1 - start1;
	//std::cout << "Time FFT" << diff1.count() << std::endl;

	for (int i = 0; i < vector.FFT_Data.size(); i++)
	{
		double A = sqrt(vector.FFT_Data[i].real() * vector.FFT_Data[i].real() + vector.FFT_Data[i].imag() * vector.FFT_Data[i].imag());
		vector.ampl.push_back(A);
	}

	for (int i = 0; i < vector.FFT_Data.size(); i++)
	{
		double fi = 0;
		if (vector.FFT_Data[i].real() > 0)
		{
			fi = atan(vector.FFT_Data[i].imag() / vector.FFT_Data[i].real());
			vector.f.push_back(fi);
		}

		if (vector.FFT_Data[i].real() == 0)
		{
			fi = Sign(vector.FFT_Data[i].imag()) * PI / 2;
			vector.f.push_back(fi);
		}

		if (vector.FFT_Data[i].real() < 0 && vector.FFT_Data[i].imag()>=0)
		{
			fi = atan(vector.FFT_Data[i].imag() / vector.FFT_Data[i].real()) + PI;
			vector.f.push_back(fi);
		}

		if (vector.FFT_Data[i].real() < 0 && vector.FFT_Data[i].imag() < 0)
		{
			fi = atan(vector.FFT_Data[i].imag() / vector.FFT_Data[i].real()) - PI;
			vector.f.push_back(fi);
		}
	}

	for (int i = 0; i < vector.FFT_Data.size(); i++)
	{
		if (vector.ampl[i] > 1)
		{
			std::cout << vector.Z[i]._Val[0] << "    " << vector.FFT_Data[i]._Val[0] << "    " << vector.FFT_Data[i]._Val[1] << "    " << vector.ampl[i] << "    " << vector.f[i] << "\n";
		}
	}

	//for (int i = 0; i < vector.DFT_Data.size(); i++)
	//{
	//	std::cout << vector.Z[i]._Val[0] << "    " << vector.DFT_Data[i]._Val[0] << "+ " << vector.DFT_Data[i]._Val[1] << "i\n";
	//}
}