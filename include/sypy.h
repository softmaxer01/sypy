#ifndef SYPY_H
#define SYPY_H
#include <vector>
#include <utility> // For std::pair

namespace sypy {

void print_array(const std::vector<double>& vec);
std::vector<double> arange(double start, double stop, double step);
std::vector<double> linspace(double start, double stop, int num);
double dot(const std::vector<double>& v1, const std::vector<double>& v2);
void print_matrix(const std::vector<std::vector<double>>& vec);
std::vector<std::vector<double>> input_matrix(int num_row, int num_col);
std::vector<std::vector<double>> gen_random_matrix(int num_row, int num_col);
void print_header(const std::vector<std::vector<double>>& vec, int a);
std::vector<std::vector<double>> add_matrix(const std::vector<std::vector<double>>& vec1, const std::vector<std::vector<double>>& vec2);
std::vector<std::vector<double>> sub_matrix(const std::vector<std::vector<double>>& vec1, const std::vector<std::vector<double>>& vec2);
std::vector<std::vector<double>> transpose_matrix(const std::vector<std::vector<double>>& vec);
std::vector<std::vector<double>> scalar_mul(const std::vector<std::vector<double>>& vec, double scalar);
long long size_of_matrix(const std::vector<std::vector<double>>& vec);
std::pair<int, int> shape_of_matrix(const std::vector<std::vector<double>>& vec);
std::vector<std::vector<double>> zeros(int row, int col);
std::vector<std::vector<double>> ones(int row, int col);
bool is_square(const std::vector<std::vector<double>>& vec);
std::vector<double> get_row(const std::vector<std::vector<double>>& vec, int num_row);
std::vector<double> get_col(const std::vector<std::vector<double>>& vec, int num_col);
double determinant(const std::vector<std::vector<double>>& vec);
std::vector<std::vector<double>> hadamard_product(const std::vector<std::vector<double>>& vec1, const std::vector<std::vector<double>>& vec2);
double trace(const std::vector<std::vector<double>>& vec);
void swap_rows(std::vector<std::vector<double>>& vec, int row1, int row2);
void swap_cols(std::vector<std::vector<double>>& vec, int col1, int col2);
std::vector<std::vector<double>> make_unit_matrix(int num_row, int num_col);
std::vector<double> get_diagonal_elements(const std::vector<std::vector<double>>& vec);
std::vector<std::vector<double>> multiply_matrices(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B);
std::vector<double> row_wise_sum(const std::vector<std::vector<double>>& matrix);
std::vector<std::vector<double>> reshape(const std::vector<double>& vec, int rows, int cols);

}

#endif
