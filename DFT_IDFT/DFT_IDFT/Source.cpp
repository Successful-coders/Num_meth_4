#include "vector.h"
#include "DFT.h"

int main()
{
	VectorZ vector(8);
	DiscreteTransformation DFT;
	vector.type = cos2pwj;
	vector.w = 100;
	vector.A = -1;
	vector.A0 = -1;
	vector.fi = -PI;
	vector.FillVector();

	for (int i = 0; i < vector.DFT_Data.size(); i++)
	{
		//std::cout << vector.DFT_Data[i]._Val[0] << " + " << vector.DFT_Data[i]._Val[1] << "i\n";
		std::cout << vector.Z[i] << "i\n";
	}
	std::cout << "\n\n";

	DFT.RealizeDFT(vector.Z, vector.DFT_Data);
	for (int i = 0; i < vector.DFT_Data.size(); i++)
	{
		std::cout << vector.DFT_Data[i]._Val[0] << " + " << vector.DFT_Data[i]._Val[1] << "i\n";
	}
	std::cout << "\n\n";

	DFT.RealizeIDFT(vector.Z, vector.DFT_Data);
	for (int i = 0; i < vector.DFT_Data.size(); i++)
	{
		std::cout << vector.DFT_Data[i]._Val[0] << " + " << vector.DFT_Data[i]._Val[1] << "i\n";
	}
	std::cout << "\n\n";

	for (int i = 0; i < vector.DFT_Data.size(); i++)
	{
		//std::cout << vector.DFT_Data[i]._Val[0] << " + " << vector.DFT_Data[i]._Val[1] << "i\n";
		std::cout << vector.Z[i] << "i\n";
	}
	std::cout << "\n\n";
}