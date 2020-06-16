#pragma once
#ifndef Discrete_Fourier_Transformation_h
#define Discrete_Fourier_Transformation_h

#include <vector>
#include <complex>

namespace Com_Methods
{
	class Discrete_Fourier_Transformation
	{
	public:

		void FFT(const std::vector<std::complex<double>>& Data, std::vector<std::complex<double>>& Result);

		void IFFT(const std::vector<std::complex<double>>& Data, std::vector<std::complex<double>>& Result);

		void Convolution(const std::vector<std::complex<double>>& Vec1,
						 const std::vector<std::complex<double>>& Vec2,
						 std::vector<std::complex<double>>& Result);
	};
}
#endif // !Discrete_Fourier_Transformation_h

