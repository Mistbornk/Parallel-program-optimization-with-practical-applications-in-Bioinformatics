#ifndef ROW_MAJOR_MATRIX_H
#define ROW_MAJOR_MATRIX_H
#include <iostream>
#include <vector>
#include <random>
#include <thread>
#include <mutex>
#include <cstddef> 
#include "colMajor.hpp"

template <typename T>
class Column_Major_Matrix;

template <typename T>
class Row_Major_Matrix {
public:	
	size_t rows, cols;
	std::vector<std::vector<T>> all_row;

	// rule of five (six ?)
	Row_Major_Matrix(size_t r, size_t c);
	Row_Major_Matrix(const Row_Major_Matrix& other);
	Row_Major_Matrix<T>& operator=(const Row_Major_Matrix& other);
	Row_Major_Matrix(Row_Major_Matrix&& other) noexcept;
	Row_Major_Matrix<T>& operator=(const Row_Major_Matrix&& other) noexcept;
	~Row_Major_Matrix() { };

	// conversion
	operator Column_Major_Matrix<T>() const;

	// getter function
	const std::vector<T> getRow(int row_idx) const;
	const std::vector<T> getColumn(int col_idx) const;
	size_t rowSize() const {return rows;}
	size_t colSize() const {return cols;}

	// setter function
	void setRow(int row_idx, const std::vector<T>& row);
	void setColumn(int col_idx, const std::vector<T>& col);

	// access element
    T& operator()(std::size_t r, std::size_t c);
    const T& operator()(std::size_t r, std::size_t c) const;

	// matrix multiplication
	Row_Major_Matrix<T> operator*(const Column_Major_Matrix<T>& rhs) const;
	Row_Major_Matrix<T> operator%(const Column_Major_Matrix<T>& rhs) const;

	// equal
	bool operator==(const Row_Major_Matrix<T>& other) const;
	bool operator==(const Column_Major_Matrix<T>& other) const;

	// cout << matrix
	template <typename U>
	friend std::ostream& operator<<(std::ostream& os, const Row_Major_Matrix<U>& matrix);
};


template <typename T>
Row_Major_Matrix<T>::Row_Major_Matrix(size_t r, size_t c) : rows(r), cols(c), all_row(r, std::vector<T>(c)) {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<T> dist(1, 10);
    for (auto& row : all_row) {
        std::generate(row.begin(), row.end(), [&]() { return dist(gen); });
    }
}

template <typename T>
Row_Major_Matrix<T>::Row_Major_Matrix(const Row_Major_Matrix& other) 
	:rows(other.rows), cols(other.cols), all_row(other.all_row) { }

template <typename T>
Row_Major_Matrix<T>& Row_Major_Matrix<T>::operator=(const Row_Major_Matrix& other) {
	if (&this != other) {
		rows = other.rows;
		cols = other.cols;
		all_row = other.all_row;
	}
	return *this;
}

template <typename T>
Row_Major_Matrix<T>::Row_Major_Matrix(Row_Major_Matrix&& other) noexcept
	:rows(other.rows), cols(other.cols), all_row(std::move(other.all_row)) {
	other.rows = 0;
	other.cols = 0;
}

template <typename T>
Row_Major_Matrix<T>& Row_Major_Matrix<T>::operator=(const Row_Major_Matrix&& other) noexcept{
	if (this != &other) {
		all_row.clear();

		all_row = std::move(other.all_row);
		rows = other.rows;
		cols = other.cols;

		other.rows = 0;
		other.cols = 0;
	}
	return *this;
}

template <typename T>
Row_Major_Matrix<T>::operator Column_Major_Matrix<T>() const {
	Column_Major_Matrix<T> converted(rows, cols);

	for (size_t i=0; i<rows; i++) {
		for (size_t j=0; j<cols; j++) {
			converted(i, j) = all_row[i][j];
		}
	}
	return converted;
}

template <typename T>
const std::vector<T> Row_Major_Matrix<T>::getRow(int row_idx) const {
    if (row_idx < 0 || row_idx >= rows) {
        throw std::out_of_range("Row index out of range");
    }
	return all_row[row_idx];
}

template <typename T>
const std::vector<T> Row_Major_Matrix<T>::getColumn(int col_idx) const {
    if (col_idx < 0 || col_idx >= cols) {
        throw std::out_of_range("Column index out of range");
    }
	std::vector<T>col(rows);
	for (size_t row=0; row<rows; row++) {
		col[row] = all_row[row][col_idx];
	}
	return col;
}

template <typename T>
void Row_Major_Matrix<T>::setRow(int row_idx, const std::vector<T>& row) {
    if (row_idx < 0 || row_idx >= rows) {
        throw std::out_of_range("Row index out of range");
    }else if (row.size() != static_cast<size_t>(cols)) {
		throw std::invalid_argument("Row size does not match the matrix's column count");
	}
	all_row[row_idx] = row;
}

template <typename T>
void Row_Major_Matrix<T>::setColumn(int col_idx, const std::vector<T>& col) {
    if (col_idx < 0 || col_idx >= cols) {
        throw std::out_of_range("Column index out of range");
    }else if (col.size() != static_cast<size_t>(rows)) {
		throw std::invalid_argument("Column size does not match the matrix's row count");
	}
	for (size_t row=0; row<rows; row++) {
		all_row[row][col_idx] = col[row];
	}
}

template <typename T>
T& Row_Major_Matrix<T>::operator()(std::size_t r, std::size_t c) {
	if (r<0 || r>=rows || c<0 || c>=cols) {
		throw std::out_of_range("Index out of range");
	}
	return all_row[r][c];
}

template <typename T>
const T& Row_Major_Matrix<T>::operator()(std::size_t r, std::size_t c) const {
	if (r<0 || r>=rows || c<0 || c>=cols) {
		throw std::out_of_range("Index out of range");
	}
	return all_row[r][c];
}

template <typename T>
Row_Major_Matrix<T> Row_Major_Matrix<T>::operator*(const Column_Major_Matrix<T>& rhs) const {
	size_t M = this->rows;
	size_t N = this->cols;
	size_t rhs_N = rhs.rows;
	size_t P = rhs.cols;

	if (N != rhs_N) {
		throw std::runtime_error("Matrix dimension mismatch for multiplication.");
	}

	Row_Major_Matrix<T> result(M, P);
	for (size_t i=0; i<M; i++) {
		for (size_t j=0; j<P; j++) {
			T sum = T();
			for (size_t k=0; k<N; k++) {
				sum += all_row[i][k] * rhs.all_column[k][j];
			}
			result(i, j) = sum;
		}
	}
	return result;
}

template <typename T>
Row_Major_Matrix<T> Row_Major_Matrix<T>::operator%(const Column_Major_Matrix<T>& rhs) const {
	size_t M = this->rows;
	size_t N = this->cols;
	size_t rhs_N = rhs.rows;
	size_t P = rhs.cols;
	// check dimension
	if (N != rhs_N) {
		throw std::runtime_error("Matrix dimension mismatch for multiplication.");
	}
	// initial parameter
	Row_Major_Matrix<T> result(M, P);
	std::vector<std::thread> threads;
	const int num_threads = 10;
	std::mutex mtx;
	// define a work lambda for thread
	auto work = [&](const int thread_id) {
		for (size_t i=thread_id; i<M; i+=num_threads) {
			for (size_t j=0; j<P; j++) {
				T sum = T();
				for (size_t k=0; k<N; k++) {
					sum += all_row[i][k] * rhs.all_column[k][j];
				}
				result(i, j) = sum;
			}
		}		
	};
	// add threads
	for (int t=0; t<num_threads; t++) {
		threads.emplace_back(work, t);
	}
	// wait threads
	for (auto& t : threads) {
		t.join();
	}
	return result;
}

template <typename T>
bool Row_Major_Matrix<T>::operator==(const Row_Major_Matrix<T>& other) const {
	if (rows != other.rows || cols != other.cols) return false;
	for (size_t i=0; i<rows; i++) {
		if (all_row[i] != other.all_row[i]) return false;
	}
	return true;
}

template <typename T>
bool Row_Major_Matrix<T>::operator==(const Column_Major_Matrix<T>& other) const {
	if (rows != other.rows || cols != other.cols) return false;
	for (size_t i=0; i<rows; i++) {
		if (all_row[i] != other.getRow(i)) return false;
	}
	return true;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const Row_Major_Matrix<T>& matrix) {
    size_t rows = matrix.rows;
    size_t cols = matrix.cols;

    os << "Row_Major_Matrix (" << rows << " x " << cols << "):\n";
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            os << matrix(i, j) << " ";
        }
        os << "\n";
    }
    return os;
}

#endif