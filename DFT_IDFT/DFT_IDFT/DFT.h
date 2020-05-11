#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>
#include "vector.h"

class DiscreteTransformation
{
public:
	void RealizeDFT(const std::vector<std::complex<double>>& data, std::vector<std::complex<double>>& result)
	{
		int M = data.size();
		for (int i = 0; i < M; i++)
		{
			for (int k = 0; k < M; k++)
			{
				
			}
		} 
	}

	void RealizeFFT(const std::vector<std::complex<double>>& data, std::vector<std::complex<double>>& result)
	{
		int N = data.size(), M = N / 2;
		result.clear(); result.resize(N);
		std::complex<double> Exp, U, V;

		for (int m = 0; m < M; m++)
		{
			U._Val[0] = 0.0; U._Val[1] = 0.0;
			V._Val[0] = 0.0; V._Val[1] = 0.0;
			for (int n = 0; n < M; n++)
			{
				Exp._Val[0] = cos(-2.0 * PI * m * n / M);
				Exp._Val[1] = sin(-2.0 * PI * m * n / M);
				U += data[2 * n] * Exp;
				V += data[2 * n + 1] * Exp;
			}

			Exp._Val[0] = cos(-2.0 * PI * m / N);
			Exp._Val[1] = sin(-2.0 * PI * m / N);
			result[m] = U + Exp * V;
			result[m + M] = U - Exp * V;
		}
	}
	void RealizeIDFT()
	{

	}
	void RealizeIFFT(const std::vector<std::complex<double>>& data, std::vector<std::complex<double>>& result)
	{

		int N = data.size();
		RealizeFFT(data, result);
		std::complex<double> Val;
		for (int i = 1; i <= N / 2; i++)
		{
			Val = result[i];
			result[i] = result[N - i] / double(N);
			result[N - i] = Val / double(N);
		}

		result[0] /= double(N);
	}
};