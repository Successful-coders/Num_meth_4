#include "Operators.h"

namespace Com_Methods
{
	void Operators::Cyclic_Shift(int k, const std::vector<std::complex<double>>& Data, std::vector<std::complex<double>>& Result)
	{
		int N = Data.size();
		Result.clear(); Result.resize(N);

		for (int i = 0; i < N; i++)
		{
			if (i - k < 0)
				Result[i] = Data[i - k + N];
			else
				Result[i] = Data[i - k];
		}
	}

	void Operators::Downsampling_Operator(int l, const std::vector<std::complex<double>>& Data, std::vector<std::complex<double>>& Result)
	{
		int Val = int(pow(2, l));
		int N = Data.size() / Val;
		Result.clear();
		Result.resize(N);
		for (int i = 0; i < N; i++)
		{
			Result[i] = Data[i * Val];
		}
	}


	void Operators::Upsampling_Operator(int l, const std::vector<std::complex<double>>& Data, std::vector<std::complex<double>>& Result)
	{
		int Val = int(pow(2, l));

		int N = Data.size() * Val;
		Result.clear();
		Result.resize(N);

		for (int i = 0; i < N; i++)
		{
			if (i % Val == 0) Result[i] = Data[i / Val];
			else Result[i] = 0;
		}
	}

	std::complex<double> Operators::Dot_Product(const std::vector<std::complex<double>>& Vec1, const std::vector<std::complex<double>>& Vec2)
	{
		int N = Vec1.size();
		std::complex<double> Result(0, 0), Val;
		for (int i = 0; i < N; i++)
		{
			Val._Val[0] = Vec2[i].real();
			Val._Val[1] = -Vec2[i].imag();
			Result += Vec1[i] * Val;
		}
		return Result;
	}
}