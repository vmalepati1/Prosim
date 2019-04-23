#include "IterationTools.h"

#include <algorithm>

#include "String.h"

namespace classical {

	namespace utils {

		int CombinationsNR(int n, int r) {
			if (n < r) return 0;
			if (n == r) return 1;
			if (r == 0) return 1;

			if (r > n / 2) return CombinationsNR(n, n - r);

			long res = 1;

			for (int k = 1; k <= r; ++k)
			{
				res *= n - k + 1;
				res /= k;
			}

			return res;
		}

		int PermutationsNR(int n, int r) {
			if (n < r) return 0;
			if (n == 0 || r == 0) {
				return 1;
			}
			else {
				if (n == r) {
					r = n - 1;
				}
				return (int)lround(((double)n / (double)(n - r)) * exp(lgamma(n) - lgamma(n - r)));
			}
		}

		void CombinationArray(int arr[], int data[], std::vector<std::vector<int>> &indices, int start, int end, int index, int r) {
			if (index == r)
			{
				std::vector<int> indexRow;

				for (int j = 0; j < r; j++)
					indexRow.push_back(data[j]);

				indices.push_back(indexRow);
				return;
			}

			for (int i = start; i <= end && end - i + 1 >= r - index; i++)
			{
				data[index] = arr[i];
				CombinationArray(arr, data, indices, i + 1, end, index + 1, r);
			}
		}

		void CombinationKN(IterationMatrix &matrix, int K, int N) {
			std::vector<int> indexRow;

			for (int i = 0; i<K; ++i) {
				indexRow.push_back(i);
			}

			int rowCount = 0;

			for (;;) {
				for (int i = 0; i < K; i++) {
					matrix(rowCount , i) = indexRow[i];
				}

				int inc_index = K - 1;
				int index_limit = N - 1;
				while (inc_index >= 0 && indexRow[inc_index] >= index_limit) {
					--inc_index;
					--index_limit;
				}
				if (inc_index < 0) {
					break;
				}

				int val = indexRow[inc_index] + 1;
				for (; inc_index<K; ++inc_index, ++val) {
					indexRow[inc_index] = val;
				}

				rowCount++;
			}
		}

		void PermutationPairsVector(const std::vector<int> &v, IterationMatrix &matrix) {
			int matrixRow = 0;

			for (int i = 0; i < v.size(); i++) {
				for (int j = 0; j < v.size(); j++) {
					if (i == j) continue;
					matrix(matrixRow, 0) = v[i];
					matrix(matrixRow, 1) = v[j];
					matrixRow++;
				}
			}
		}

	}

}