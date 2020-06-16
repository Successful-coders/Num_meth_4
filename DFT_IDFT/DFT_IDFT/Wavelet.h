#pragma once
#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>
#include "vector.h"
#include "Operators.h"


Wavelet::Wavelet(int Count_Data, Basis_Type Type)
{
	int N = Count_Data;
	U.clear(); U.resize(N);
	V.clear(); V.resize(N);

	switch (Type)
	{
	case Basis_Type::Complex_Shannon:
	{
		if (N % 4 != 0) throw std::exception("\nError in Complex Shannon Basis: N % 4 != 0 ...\n");

		std::complex<double> Val;

		U[0] = V[0] = 1.0 / sqrt(2.0);

		for (int i = 1; i < N; i++)
		{
			Val._Val[0] = sqrt(2.0) / N * cos(PI * i / N) * sin(PI * i / 2.0) / sin(PI * i / N);
			Val._Val[1] = -sqrt(2.0) / N * sin(PI * i / N) * sin(PI * i / 2.0) / sin(PI * i / N);
			U[i] = Val;
			V[i] = pow(-1, i) * Val;
		}
		break;
	}
	case Basis_Type::Dobeshi:
	{
		float a = 1 - sqrt(10);
		float b = 1 + sqrt(10);
		float c = sqrt(5 + 2 * sqrt(10));
		float d = sqrt(2) / 32;

		std::complex<double> Val;

		Val._Val[0] = d * (b + c);
		Val._Val[1] = 0;
		U[0] = Val;

		Val._Val[0] = d * (2 * a + 3 * b + 3 * c);
		U[1] = Val;

		Val._Val[0] = d * (6 * a + 4 * b + 2 * c);
		U[2] = Val;

		Val._Val[0] = d * (6 * a + 4 * b - 2 * c);
		U[3] = Val;

		Val._Val[0] = d * (2 * a + 3 * b - 3 * c);
		U[4] = Val;

		Val._Val[0] = d * (b - c);
		U[5] = Val;

		Val._Val[0] = 0;

		for (int i = 6; i < N; i++)
		{
			U[i] = Val;
		}

		V[0] = -U[1];
		V[1] = U[0];
		V[N - 1] = U[2];
		V[N - 2] = -U[3];
		V[N - 3] = U[4];
		V[N - 4] = -U[5];

		for (int i = 2; i < N - 4; i++)
		{
			V[i] = Val;
		}
	}
	}
}

void Wavelet::Wavelet_Filters_System(int Stages)
{

	Operators Operator;

	int N = U.size();

	std::vector<std::vector<std::complex<double>>> U_Filters(Stages), V_Filters(Stages);

	U_Filters[0] = U;
	V_Filters[0] = V;

	for (int i = 1; i < Stages; i++)
	{
		int Count_Elements = N / int(pow(2, i));

		U_Filters[i].resize(Count_Elements);
		V_Filters[i].resize(Count_Elements);

		for (int n = 0; n < Count_Elements; n++)
		{
			int Max_Index = int(pow(2, i));
			for (int k = 0; k < Max_Index; k++)
			{
				U_Filters[i][n] += U_Filters[0][n + k * N / Max_Index];
				V_Filters[i][n] += V_Filters[0][n + k * N / Max_Index];
			}
		}
	}

	Com_Methods::Discrete_Fourier_Transformation DFT;

	std::vector<std::complex<double>> U_u, U_v;

	f.resize(Stages); g.resize(Stages);
	f[0] = V_Filters[0];
	g[0] = U_Filters[0];
	for (int i = 1; i < Stages; i++)
	{
		Operator.Upsampling_Operator(i, U_Filters[i], U_u);
		Operator.Upsampling_Operator(i, V_Filters[i], U_v);
		DFT.Convolution(g[i - 1], U_v, f[i]);
		DFT.Convolution(g[i - 1], U_u, g[i]);
	}
}

void Wavelet::Wavelet_Basis(int Stage, std::vector<std::vector<std::complex<double>>>& Psi,	std::vector<std::vector<std::complex<double>>>& Fi)
{
	Operators Operator;

	int Count_Data = U.size();

	int Count_Bas_Elements = Count_Data / int(pow(2, Stage));

	if (g.size() < Stage) Wavelet_Filters_System(Stage + 1);

	Psi.resize(Count_Bas_Elements); Fi.resize(Count_Bas_Elements);

	for (int i = 0; i < Count_Bas_Elements; i++)
	{
		int Index = int(pow(2, Stage)) * i;
		Operator.Cyclic_Shift(Index, f[Stage - 1], Psi[i]);
		Operator.Cyclic_Shift(Index, g[Stage - 1], Fi[i]);
	}
}

void  Wavelet::Analysis_Phase(int Stage, const std::vector<std::complex<double>>& Data, std::vector<std::complex<double>>& Koef_Psi, std::vector<std::complex<double>>& Koef_Fi)
{
	Operators Operator;

	std::vector<std::vector<std::complex<double>>> Psi, Fi;

	Wavelet_Basis(Stage, Psi, Fi);

	int Count_Bas_Elements = Psi.size();

	for (int Bas_Element_Index = 0; Bas_Element_Index < Count_Bas_Elements; Bas_Element_Index++)
	{
		Koef_Psi.push_back(Operator.Dot_Product(Data, Psi[Bas_Element_Index]));
		Koef_Fi.push_back(Operator.Dot_Product(Data, Fi[Bas_Element_Index]));
	}
}


void  Wavelet::Synthesis_Phase(int Stage, const std::vector<std::complex<double>>& Koef_Psi, const std::vector<std::complex<double>>& Koef_Fi, std::vector<std::complex<double>>& P,
	std::vector<std::complex<double>>& Q,std::vector<std::complex<double>>& Recovery)
{
	std::vector<std::vector<std::complex<double>>> Psi, Fi;

	Wavelet_Basis(Stage, Psi, Fi);

	int Count_Bas_Elements = Psi.size();

	int Count_Data = U.size();

	for (int Data_Index = 0; Data_Index < Count_Data; Data_Index++)
	{
		std::complex<double> P_Stage(0, 0), Q_Stage(0, 0);

		for (int Bas_Element_Index = 0; Bas_Element_Index < Count_Bas_Elements; Bas_Element_Index++)
		{
			P_Stage += Koef_Fi[Bas_Element_Index] * Fi[Bas_Element_Index][Data_Index];
			Q_Stage += Koef_Psi[Bas_Element_Index] * Psi[Bas_Element_Index][Data_Index];

		}
		P.push_back(P_Stage);
		Q.push_back(Q_Stage);
		Recovery.push_back(P[Data_Index] + Q[Data_Index]);
	}
}
