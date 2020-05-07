#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>
#define PI 3.1415926535897932

enum FUNCTYPE {
	cos2pwj,//задание 2
	Acos2pj,//задание 3
	cos2pjplusfi,//задание 4
	A0pluscos2pj,//задание 5
	cos2pplus001//задание 6
};

class VectorZ
{
public:
	int N;
	FUNCTYPE type;
	int w; //частота
	int A; // константа для 3 задания
	int fi;// константа для 4 задания
	int A0; // константа для 5 задания

	//Z - данные, DFT_Data - прямое DFT, IDFT_Data - обратное DFT
	std::vector<std::complex<double>> Z, DFT_Data, FFT_Data, IDFT_Data, IFFT_Data;

	VectorZ(int N);
	void FillVector();
};