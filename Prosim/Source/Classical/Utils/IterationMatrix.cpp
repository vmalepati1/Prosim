#include "IterationMatrix.h"

namespace classical {

	namespace utils {

		IterationMatrix::IterationMatrix(size_t rows, size_t columns)
			: m_rows(rows),
			m_columns(columns),
			m_data(rows * columns)
		{
		}

		int& IterationMatrix::operator()(size_t i, size_t j)
		{
			return m_data[i * m_columns + j];
		}

		int IterationMatrix::operator()(size_t i, size_t j) const
		{
			return m_data[i * m_columns + j];
		}

	}

}