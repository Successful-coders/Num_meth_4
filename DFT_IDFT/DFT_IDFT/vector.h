#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>
#define PI 3.1415926535897932

enum FUNCTYPE {
	cos2pwj,//������� 2
	Acos2pj,//������� 3
	cos2pjplusfi,//������� 4
	A0pluscos2pj,//������� 5
	cos2pplus001//������� 6
};

class VectorZ
{
public:
	int N;
	FUNCTYPE type;
	int w; //�������
	int A; // ��������� ��� 3 �������
	int fi;// ��������� ��� 4 �������
	int A0; // ��������� ��� 5 �������

	//Z - ������, DFT_Data - ������ DFT, IDFT_Data - �������� DFT
	std::vector<std::complex<double>> Z, DFT_Data, FFT_Data, IDFT_Data, IFFT_Data;

	VectorZ(int N);
	void FillVector();
};