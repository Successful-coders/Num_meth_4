#include "vector.h"
#include "DFT.h"

int main()
{
	VectorZ vector(256);
	DiscreteTransformation DFT;
	vector.type = cos2pwj;
	vector.w = 100;
	vector.A = -1;
	vector.A0 = -1;
	vector.fi = -PI;
	vector.FillVector();
	
	DFT.RealizeDFT(vector.Z, vector.DFT_Data);
	DFT.RealizeFFT(vector.Z, vector.FFT_Data);

	for (int i = 0; i < vector.FFT_Data.size(); i++)
	{
		std::cout << vector.FFT_Data[i]._Val[0] << "+ " << vector.FFT_Data[i]._Val[1] << "i\n";
	}
	std::cout <<"\n";
	std::cout <<"\n";
	std::cout <<"\n";
	std::cout <<"\n";

	for (int i = 0; i < vector.DFT_Data.size(); i++)
	{
		std::cout << vector.DFT_Data[i]._Val[0] << "+ " << vector.DFT_Data[i]._Val[1] << "i\n";
	}
}