#pragma once

#include <vector>

#include "IterationMatrix.h"

namespace classical {

	namespace utils {

		int CombinationsNR(int n, int r);
		int PermutationsNR(int n, int r);
		void CombinationArray(int arr[], int data[], std::vector<std::vector<int>> &indices, int start, int end, int index, int r);
		void CombinationKN(IterationMatrix &matrix, int K, int N);
		void PermutationPairsVector(const std::vector<int> &v, IterationMatrix &indices);

	}

}