#include "vector.h"
#include "DFT.h"
#include "Wavelet_Analysis.h"

int main()
{
	//число отчсётов
	int N = 512;

	//число этапов
	int Stage = 4;

	//сигнал и его коеффициенты разложения по вейвлет-базиксу
	std::vector<std::complex<double>> Z(N), Koef_Psi, Koef_Fi;

	//частичное восстановление, уточняющие данные и восстановленный сигнал
	std::vector<std::complex<double>> P, Q, Z_Rec;

	//формирование исходного сигнала
	for (int n = 0; n < N; n++)
	{
		if (n >= 128 && n <= 255) Z[n]._Val[0] = sin(fabs(pow(n - 128, 1.7)) / 128.0);
		if (n >= 384 && n <= 447) Z[n]._Val[0] = sin(fabs(pow(n - 128, 2.0)) / 128.0);
	}

	//выбор базисной системы
	Com_Methods::Wavelet_Analysis Wavelet_Test(N, Com_Methods::Wavelet_Analysis::Basis_Type::Complex_Shannon);

	//Фаза анализа данных
	Wavelet_Test.Analysis_Phase(Stage, Z, Koef_Psi, Koef_Fi);

	//Фаза восстановления
	Wavelet_Test.Synthesis_Phase(Stage, Koef_Psi, Koef_Fi, P, Q, Z_Rec);
	
	//Печать: исходный и восстановленный сигнал
	int SETW = 22;
	std::cout << std::left << std::setw(SETW) << "Number" << std::setw(SETW) << "Z" << std::setw(SETW) << "Z_Recovery" << std::endl;
	for (int i = 0; i < Z.size(); i++)
	{
		std::cout << std::left << std::setw(SETW) << i
			<< std::setw(SETW) << Z[i].real()
			<< std::setw(SETW) << P[i].real() << std::endl;
	}

	//Печать: вейвлет-коэфы
	std::cout << std::left << std::setw(SETW) << "Number" << std::setw(SETW) << "Psi" << std::setw(SETW) << "Fi" << std::endl;
	for (int i = 0; i < Koef_Fi.size(); i++)
	{
		std::cout << std::left << std::setw(SETW) << pow(2, Stage) * i
			<< std::setw(SETW) << sqrt(pow(Koef_Psi[i].real(), 2) + pow(Koef_Psi[i].imag(), 2))
			<< std::setw(SETW) << sqrt(pow(Koef_Fi[i].real(), 2) + pow(Koef_Fi[i].imag(), 2)) << std::endl;
	}
}