#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <random>
#include "sypy.h"

namespace sypy {

void print_array(const std::vector<double>& vec) {
    for (const auto& i : vec) std::cout << i << " ";
    std::cout << std::endl;
}

std::vector<double> arange(double start, double stop, double step) {
    std::vector<double> values;
    for (double value = start; value < stop; value += step) values.push_back(value);
    return values;
}

std::vector<double> linspace(double start, double stop, int num) {
    std::vector<double> result;
    if (num == 0) return result;
    if (num == 1) { result.push_back(start); return result; }
    double step = (stop - start) / (num - 1);
    for (int i = 0; i < num; ++i) result.push_back(start + i * step);
    return result;
}

double dot(const std::vector<double>& v1, const std::vector<double>& v2) {
    if (v1.size() != v2.size()) throw std::invalid_argument("Vectors must have the same size for dot product.");
    return std::inner_product(v1.begin(), v1.end(), v2.begin(), 0.0);
}

void print_matrix(const std::vector<std::vector<double>>& vec) {
    if (vec.empty()) return;
    for (const auto& row : vec) {
        for (const auto& elem : row) std::cout << elem << " ";
        std::cout << std::endl;
    }
}

std::vector<std::vector<double>> input_matrix(int num_row, int num_col) {
    if (num_row <= 0 || num_col <= 0) throw std::invalid_argument("Matrix dimensions must be positive.");
    std::vector<std::vector<double>> vec(num_row, std::vector<double>(num_col));
    for (int i = 0; i < num_row; ++i) for (int j = 0; j < num_col; ++j) std::cin >> vec[i][j];
    return vec;
}

std::vector<std::vector<double>> gen_random_matrix(int num_row, int num_col) {
    if (num_row <= 0 || num_col <= 0) throw std::invalid_argument("Matrix dimensions must be positive.");
    std::vector<std::vector<double>> vec(num_row, std::vector<double>(num_col));
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    for (int i = 0; i < num_row; ++i) for (int j = 0; j < num_col; ++j) vec[i][j] = dis(gen);
    return vec;
}

void print_header(const std::vector<std::vector<double>>& vec, int a) {
    if (vec.empty()) return;
    for (int i = 0; i < a && i < vec.size(); ++i) {
        for (int j = 0; j < vec[0].size(); ++j) std::cout << vec[i][j] << " ";
        std::cout << std::endl;
    }
}

std::vector<std::vector<double>> add_matrix(const std::vector<std::vector<double>>& vec1, const std::vector<std::vector<double>>& vec2) {
    if (vec1.size() != vec2.size() || vec1[0].size() != vec2[0].size()) throw std::invalid_argument("Matrices must have the same dimensions for addition.");
    std::vector<std::vector<double>> result(vec1.size(), std::vector<double>(vec1[0].size()));
    for (size_t i = 0; i < vec1.size(); ++i) for (size_t j = 0; j < vec1[0].size(); ++j) result[i][j] = vec1[i][j] + vec2[i][j];
    return result;
}

std::vector<std::vector<double>> sub_matrix(const std::vector<std::vector<double>>& vec1, const std::vector<std::vector<double>>& vec2) {
    if (vec1.size() != vec2.size() || vec1[0].size() != vec2[0].size()) throw std::invalid_argument("Matrices must have the same dimensions for subtraction.");
    std::vector<std::vector<double>> result(vec1.size(), std::vector<double>(vec1[0].size()));
    for (size_t i = 0; i < vec1.size(); ++i) for (size_t j = 0; j < vec1[0].size(); ++j) result[i][j] = vec1[i][j] - vec2[i][j];
    return result;
}

std::vector<std::vector<double>> transpose_matrix(const std::vector<std::vector<double>>& vec) {
    if (vec.empty()) return {};
    std::vector<std::vector<double>> result(vec[0].size(), std::vector<double>(vec.size()));
    for (size_t i = 0; i < vec.size(); ++i) for (size_t j = 0; j < vec[0].size(); ++j) result[j][i] = vec[i][j];
    return result;
}

std::vector<std::vector<double>> scalar_mul(const std::vector<std::vector<double>>& vec, double scalar) {
    if (vec.empty()) return {};
    std::vector<std::vector<double>> result(vec.size(), std::vector<double>(vec[0].size()));
    for (size_t i = 0; i < vec.size(); ++i) for (size_t j = 0; j < vec[0].size(); ++j) result[i][j] = vec[i][j] * scalar;
    return result;
}

long long size_of_matrix(const std::vector<std::vector<double>>& vec) {
    if (vec.empty()) return 0;
    return vec.size() * vec[0].size();
}

std::pair<int, int> shape_of_matrix(const std::vector<std::vector<double>>& vec) {
    if (vec.empty()) { return {0, 0}; }
    return {static_cast<int>(vec.size()), static_cast<int>(vec[0].size())};
}

std::vector<std::vector<double>> zeros(int row, int col) {
    if (row <= 0 || col <= 0) throw std::invalid_argument("Matrix dimensions must be positive.");
    return std::vector<std::vector<double>>(row, std::vector<double>(col, 0));
}

std::vector<std::vector<double>> ones(int row, int col) {
    if (row <= 0 || col <= 0) throw std::invalid_argument("Matrix dimensions must be positive.");
    return std::vector<std::vector<double>>(row, std::vector<double>(col, 1));
}

bool is_square(const std::vector<std::vector<double>>& vec) {
    if (vec.empty()) { return true; } // 0x0 is square
    return vec.size() == vec[0].size();
}

std::vector<double> get_row(const std::vector<std::vector<double>>& vec, int num_row) {
    if (num_row < 0 || num_row >= vec.size()) throw std::out_of_range("Row index out of range.");
    return vec[num_row];
}

std::vector<double> get_col(const std::vector<std::vector<double>>& vec, int num_col) {
    if (vec.empty() || num_col < 0 || num_col >= vec[0].size()) throw std::out_of_range("Column index out of range.");
    std::vector<double> col_vec;
    for (size_t i = 0; i < vec.size(); ++i) col_vec.push_back(vec[i][num_col]);
    return col_vec;
}

double determinant(const std::vector<std::vector<double>>& vec) {
    if (vec.empty()) return 1;
    if (vec.size() != vec[0].size()) throw std::invalid_argument("Matrix must be square to calculate determinant.");
    int n = vec.size();
    if (n == 0) return 1;
    if (n == 1) return vec[0][0];

    std::vector<std::vector<double>> temp = vec;

    double det = 1.0;

    for (int i = 0; i < n; ++i) {
        int pivot = i;
        for (int j = i + 1; j < n; ++j) {
            if (std::abs(temp[j][i]) > std::abs(temp[pivot][i])) {
                pivot = j;
            }
        }
        if (pivot != i) {
            std::swap(temp[i], temp[pivot]);
            det *= -1.0;
        }

        if (std::abs(temp[i][i]) < 1e-9) {
            return 0.0;
        }

        det *= temp[i][i];
        for (int j = i + 1; j < n; ++j) {
            temp[i][j] /= temp[i][i];
        }
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                double mult = temp[j][i];
                for (int k = i; k < n; ++k) {
                    temp[j][k] -= mult * temp[i][k];
                }
            }
        }
    }

    return det;
}

std::vector<std::vector<double>> hadamard_product(const std::vector<std::vector<double>>& vec1, const std::vector<std::vector<double>>& vec2) {
    if (vec1.size() != vec2.size() || vec1[0].size() != vec2[0].size()) throw std::invalid_argument("Matrices must have the same dimensions for hadamard product.");
    std::vector<std::vector<double>> result(vec1.size(), std::vector<double>(vec1[0].size()));
    for (size_t i = 0; i < vec1.size(); ++i) for (size_t j = 0; j < vec1[0].size(); ++j) result[i][j] = vec1[i][j] * vec2[i][j];
    return result;
}

double trace(const std::vector<std::vector<double>>& vec) {
    if (vec.size() != vec[0].size()) throw std::invalid_argument("Matrix must be square to calculate trace.");
    double result = 0;
    for (size_t i = 0; i < vec.size(); ++i) result += vec[i][i];
    return result;
}

void swap_rows(std::vector<std::vector<double>>& vec, int row1, int row2) {
    if (row1 < 0 || row1 >= vec.size() || row2 < 0 || row2 >= vec.size()) throw std::out_of_range("Row index out of range.");
    std::swap(vec[row1], vec[row2]);
}

void swap_cols(std::vector<std::vector<double>>& vec, int col1, int col2) {
    if (vec.empty() || col1 < 0 || col1 >= vec[0].size() || col2 < 0 || col2 >= vec[0].size()) throw std::out_of_range("Column index out of range.");
    for (size_t i = 0; i < vec.size(); ++i) std::swap(vec[i][col1], vec[i][col2]);
}

std::vector<std::vector<double>> make_unit_matrix(int num_row, int num_col) {
    if (num_row <= 0 || num_col <= 0) throw std::invalid_argument("Matrix dimensions must be positive.");
    std::vector<std::vector<double>> vec(num_row, std::vector<double>(num_col, 0));
    for (int i = 0; i < std::min(num_row, num_col); ++i) vec[i][i] = 1;
    return vec;
}

std::vector<double> get_diagonal_elements(const std::vector<std::vector<double>>& vec) {
    if (vec.empty()) return {};
    std::vector<double> diagonals;
    for (size_t i = 0; i < std::min(vec.size(), vec[0].size()); ++i) diagonals.push_back(vec[i][i]);
    return diagonals;
}

std::vector<std::vector<double>> multiply_matrices(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B) {
    if (A.empty() || B.empty() || A[0].size() != B.size()) throw std::invalid_argument("Matrix dimensions are not compatible for multiplication.");
    std::vector<std::vector<double>> result(A.size(), std::vector<double>(B[0].size(), 0));
    for (size_t i = 0; i < A.size(); ++i) for (size_t j = 0; j < B[0].size(); ++j) for (size_t k = 0; k < A[0].size(); ++k) result[i][j] += A[i][k] * B[k][j];
    return result;
}

std::vector<double> row_wise_sum(const std::vector<std::vector<double>>& matrix) {
    if (matrix.empty()) return {};
    std::vector<double> sums;
    for (const auto& row : matrix) sums.push_back(std::accumulate(row.begin(), row.end(), 0.0));
    return sums;
}

std::vector<std::vector<double>> reshape(const std::vector<double>& vec, int rows, int cols) {
    if (vec.size() != rows * cols) throw std::invalid_argument("Total number of elements must remain the same for reshape.");
    std::vector<std::vector<double>> result(rows, std::vector<double>(cols));
    for (int i = 0; i < rows; ++i) for (int j = 0; j < cols; ++j) result[i][j] = vec[i * cols + j];
    return result;
}

} 
