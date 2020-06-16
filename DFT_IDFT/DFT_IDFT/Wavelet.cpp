#include "Wavelet.h"

class Wavelet
{
private:
	std::vector<std::complex<double>> U, V;
	std::vector<std::vector<std::complex<double>>> f, g;

public:
	enum class Basis_Type
	{
		Complex_Shannon = 1,
		Dobeshi
	};


	Wavelet(int Count_Data, Basis_Type Type);

private:
	void Wavelet_Filters_System(int Stages);

	void Wavelet_Basis(int Stage,
		std::vector<std::vector<std::complex<double>>>& Psi,
		std::vector<std::vector<std::complex<double>>>& Fi);

public:
	void Analysis_Phase(int Stage,
		const std::vector<std::complex<double>>& Data,
		std::vector<std::complex<double>>& Koef_Psi,
		std::vector<std::complex<double>>& Koef_Fi);


	void Synthesis_Phase(int Stage,
		const std::vector<std::complex<double>>& Koef_Psi,
		const std::vector<std::complex<double>>& Koef_Fi,
		std::vector<std::complex<double>>& P,
		std::vector<std::complex<double>>& Q,
		std::vector<std::complex<double>>& Recovery);
};
