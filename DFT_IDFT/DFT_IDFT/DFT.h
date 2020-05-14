#pragma once
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
		std::complex<double> sum;
		for (int k = 0; k < K; k++)
		{
			sum = std::complex<double>(0, 0);
			for (int n = 0; n < N; n++)
			{
				double realPart = cos((2.0 * PI * k * n / N) );
				double imagPart = sin((2.0 * PI * k * n / N) );
				std::complex<double> w(realPart, -imagPart);
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
		std::vector<std::vector<std::complex<double>>> W;
		int N = data.size();
		result.clear(); result.resize(N);
		W.resize(N);
		for (int i = 0; i < N; i++)
		{
			W[i].resize(N);
			W[i][0] = 1.0;
			W[0][i] = 1.0;
		}
		for (int i = 1; i < N; i++)
			for (int j = i; j < N; j++)
			{
				W[i][j]._Val[0] = W[j][i]._Val[0] = cos(-2.0 * PI * i * j / N);
				W[i][j]._Val[1] = W[j][i]._Val[1] = -sin(-2.0 * PI * i * j / N);
			}

		for (int i = 0; i < N; i++)
		{
			result[i]._Val[0] = 0.0;
			result[i]._Val[1] = 0.0;
			for (int j = 0; j < N; j++)
				result[i] += data[j] * W[i][j];
			result[i] /= double(N);
		}
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