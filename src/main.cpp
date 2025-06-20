#include <iostream>
#include "sypy.h"

int main() {
    std::vector<double> arr1 = {1.0, 2.0, 3.0, 4.0, 5.0};
    sypy::print_array(arr1);
    std::vector<double> arr2 = {6.0, 7.0, 8.0, 9.0, 10.0};
    std::cout << "Dot product: " << sypy::dot(arr1, arr2) << std::endl;
    auto mat1 = sypy::zeros(2, 3);
    sypy::print_matrix(mat1);
    auto mat2 = sypy::ones(2, 3);
    sypy::print_matrix(mat2);
    auto mat3 = sypy::add_matrix(mat1, mat2);
    sypy::print_matrix(mat3);
    auto mat4 = sypy::gen_random_matrix(3, 2);
    sypy::print_matrix(mat4);
    auto mat5 = sypy::multiply_matrices(mat3, mat4);
    sypy::print_matrix(mat5);
    auto shape = sypy::shape_of_matrix(mat5);
    std::cout << "Shape: " << shape.first << "x" << shape.second << std::endl;
    std::cout << "Trace: " << sypy::trace(mat5) << std::endl;
    auto mat6 = sypy::gen_random_matrix(2, 2);
    sypy::print_matrix(mat6);
    std::cout << "Determinant: " << sypy::determinant(mat6) << std::endl;
    std::vector<double> v = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    auto mat7 = sypy::reshape(v, 2, 3);
    sypy::print_matrix(mat7);
    auto lin = sypy::linspace(0.0, 1.0, 5);
    sypy::print_array(lin);
    auto ara = sypy::arange(0, 10, 2);
    sypy::print_array(ara);
    return 0;
} 
