#include "vector.h"

VectorZ::VectorZ(int N1)
{
	N = N1;
	Z.resize(N1);
	DFT_Data.resize(0);
	FFT_Data.resize(0);
	IDFT_Data.resize(0);
	IFFT_Data.resize(0);

}
void VectorZ::FillVector()
{
	if (type == cos2pwj)
	{
		
		for (int i = 0; i < N; i++)
		{
			Z[i]._Val[0] = cos(2.0 * PI * i * w / N);
			Z[i]._Val[1] = 0;
		}
	}

	if (type == Acos2pj)
	{
		for (int i = 0; i < N; i++)
		{
			Z[i]._Val[0] = A*cos(2.0 * PI * i / N);
			Z[i]._Val[1] = 0;
		}
	}

	if (type == cos2pjplusfi)
	{
		for (int i = 0; i < N; i++)
		{
			Z[i]._Val[0] = cos(2.0 * PI * i / N + fi);
			Z[i]._Val[1] = 0;
		}
	}

	if (type == A0pluscos2pj)
	{
		for (int i = 0; i < N; i++)
		{
			Z[i]._Val[0] = A0 + cos(2.0 * PI * i / N);
			Z[i]._Val[1] = 0;
		}
	}

	if (type == cos2pplus001)
	{
		for (int i = 0; i < N; i++)
		{
			Z[i]._Val[0] = cos(2.0 * PI * i / N) + 0.01 * cos(2.0 * PI * i * w/ N);
			Z[i]._Val[1] = 0;
		}
	}
}