#include "vector.h"
#include "DFT.h"
#include <chrono>
int main()
{
	VectorZ vector(512);
	DiscreteTransformation DFT;
	vector.type = cos2pwj;
	vector.w = 100;
	vector.A = -1;
	vector.A0 = -1;
	vector.fi = -PI;
	vector.FillVector();
	
	//auto start = std::chrono::high_resolution_clock::now();

	DFT.RealizeDFT(vector.Z, vector.DFT_Data);
	//auto end = std::chrono::high_resolution_clock::now();
	//std::chrono::duration<double> diff = end - start;
	//std::cout << "Time DFT " << diff.count() << std::endl;

	//auto start1 = std::chrono::high_resolution_clock::now();

	DFT.RealizeFFT(vector.Z, vector.FFT_Data);
	
	//auto end1 = std::chrono::high_resolution_clock::now();
	//std::chrono::duration<double> diff1 = end1 - start1;
	//std::cout << "Time FFT" << diff1.count() << std::endl;


	for (int i = 0; i < vector.FFT_Data.size(); i++)
	{
		std::cout << vector.Z[i]._Val[0] << "    " <<vector.FFT_Data[i]._Val[0] << "    " << vector.FFT_Data[i]._Val[1] << "i\n";
	}

	//for (int i = 0; i < vector.DFT_Data.size(); i++)
	//{
	//	std::cout << vector.Z[i]._Val[0] << "    " << vector.DFT_Data[i]._Val[0] << "+ " << vector.DFT_Data[i]._Val[1] << "i\n";
	//}
}