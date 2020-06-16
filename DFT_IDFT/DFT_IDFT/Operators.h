#pragma once
#ifndef Complex_Operators_h
#define Complex_Operators_h

#include <complex>
#include <vector>
#include "CONST.h"

namespace Com_Methods
{
	class Operators
	{
	public:
		//�������� ������
		//k - �������� ������, Data - ������� ������, Result - ������ ��� ���������
		void Cyclic_Shift(int k, const std::vector<std::complex<double>>& Data, std::vector<std::complex<double>>& Result);

		//�������� ��������� �������
		//l - ���������� �������, Data - ��������, Result - �����
		void Downsampling_Operator(int l, const std::vector<std::complex<double>>& Data, std::vector<std::complex<double>>& Result);

		//�������� ����������� �������
		//l - ���������� �������, Data - ��������, Result - �����
		void Upsampling_Operator(int l, const std::vector<std::complex<double>>& Data, std::vector<std::complex<double>>& Result);

		//��������� ������������ ���� ��������
	    std::complex<double> Dot_Product(const std::vector<std::complex<double>>& Vec1, const std::vector<std::complex<double>>& Vec2);
	};
}

#endif // ! Complex_Operators_h

