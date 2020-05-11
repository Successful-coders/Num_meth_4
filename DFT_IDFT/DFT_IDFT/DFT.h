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
		int N = data.size();
		int K = N;
		result.resize(N);

		std::complex<double> sum;

		for (int k = 0; k < K; k++)
		{
			for (int n = 0; n < N; n++)
			{
				double realPart = cos((2 * PI / N) * k * n);
				double imagPart = sin((2 * PI / N) * k * n);
				std::complex<double> w(realPart, imagPart);
				sum += data[n] * w;
			}
			result.push_back(sum);
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
	void RealizeIDFT(const std::vector<std::complex<double>>& data, std::vector<std::complex<double>>& result)
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