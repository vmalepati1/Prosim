#pragma once

#include <vector>

namespace classical {

	namespace utils {

		class IterationMatrix
		{
		public:
			IterationMatrix(size_t rows, size_t columns);
			int& operator()(size_t i, size_t j);
			int operator()(size_t i, size_t j) const;

			inline size_t GetRows() const { return m_rows; }
			inline size_t GetColumns() const { return m_columns; }
			inline std::vector<int> &GetDataWritable() { return m_data; }
			inline const std::vector<int> &GetData() const { return m_data; }
		private:
			size_t m_rows;
			size_t m_columns;
			std::vector<int> m_data;
		};

	}

}