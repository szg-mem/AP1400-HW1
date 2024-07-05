#include "hw1.h"
#include <random>
#include <iomanip>
#include <iostream>

namespace algebra {
    Matrix zeros(size_t n, size_t m) {
        return std::vector<std::vector<double>>(n, std::vector<double>(m, 0.0));
    }

    Matrix ones(size_t n, size_t m) {
        return std::vector<std::vector<double>>(n, std::vector<double>(m, 1.0));
    }

    Matrix random(size_t n, size_t m, double min, double max) {
        if (min > max) {
            throw std::logic_error("min num grater than max num");
        }
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(min, max);

        Matrix matrix = std::vector<std::vector<double>>(n, std::vector<double>(m, 0.0));
        for (size_t row = 0; row < n; ++row) {
            for (size_t col = 0; col < m; ++col) {
                matrix[row][col] = dis(gen);
            }
        }

        return matrix;
    }

    void show(const Matrix &matrix) {
        const size_t Row = matrix.size();
        const size_t Col = matrix[0].size();

        for (size_t row = 0; row < Row; ++row) {
            for (size_t col = 0; col < Col; ++col) {
                std::cout << std::fixed << std::setprecision(3) << matrix[row][col];
                std::cout << " ";
            }
            std::cout << std::endl;
        }
    }

    Matrix multiply(const Matrix &matrix, double c) {
        const size_t Row = matrix.size();
        const size_t Col = matrix[0].size();
        Matrix result = zeros(Row, Col);

        for (size_t row = 0; row < Row; ++row) {
            for (size_t col = 0; col < Col; ++col) {
                result[row][col] = c * matrix[row][col];
            }
        }

        return result;
    }

    Matrix multiply(const Matrix &matrix1, const Matrix &matrix2) {
        if (matrix1.empty()) {
            if (matrix2.empty()) {
                return Matrix{};
            }
            throw std::logic_error("dimensions are not match");
        }
        if (matrix1[0].size() != matrix2.size()) {
            throw std::logic_error("dimensions are not match");
        }

        const size_t a = matrix1.size();
        const size_t b = matrix1[0].size();
        const size_t c = matrix2[0].size();

        Matrix result = zeros(a, c);
        for (size_t r = 0; r < a; ++r) {
            for (size_t l = 0; l < c; ++l) {
                double product = 0;
                for (size_t count = 0; count < b; ++count) {
                    product += matrix1[r][count] * matrix2[count][l];
                }
                result[r][l] = product;
            }
        }

        return result;
    }

    Matrix sum(const Matrix &matrix, double c) {
        if (matrix.empty()) {
            return Matrix{};
        }
        const size_t Row = matrix.size();
        const size_t Col = matrix[0].size();

        Matrix result = zeros(Row, Col);

        for (size_t row = 0; row < Row; ++row) {
            for (size_t col = 0; col < Col; ++col) {
                result[row][col] = c + matrix[row][col];
            }
        }

        return result;
    }

    Matrix sum(const Matrix &matrix1, const Matrix &matrix2) {
        if (matrix1.empty()) {
            if (matrix2.empty()) {
                return Matrix{};
            }
            throw std::logic_error("dimensions are not match");
        }
        if (matrix2.empty() or matrix1.size() != matrix2.size() or matrix1[0].size() != matrix2[0].size()) {
            throw std::logic_error("dimensions are not match");
        }

        const size_t Row = matrix1.size();
        const size_t Col = matrix1[0].size();

        Matrix result = std::vector<std::vector<double>>(Row, std::vector<double>(Col, 0.0));

        for (size_t row = 0; row < Row; ++row) {
            for (size_t col = 0; col < Col; ++col) {
                result[row][col] = matrix1[row][col] + matrix2[row][col];
            }
        }

        return result;
    }

    Matrix transpose(const Matrix &matrix) {
        if (matrix.empty()) {
            return Matrix{};
        }
        const size_t Row = matrix.size();
        const size_t Col = matrix[0].size();

        Matrix result = zeros(Col, Row);
        for (size_t row = 0; row < Row; ++row) {
            for (size_t col = 0; col < Col; ++col) {
                result[col][row] = matrix[row][col];
            }
        }

        return result;
    }

    Matrix minor(const Matrix &matrix, size_t n, size_t m) {
        if (n >= matrix.size() or m >= matrix[0].size()) {
            throw std::logic_error("row/col is not match");
        }

        const size_t Row = matrix.size();
        const size_t Col = matrix[0].size();
        Matrix result = zeros(Row - 1, Col - 1);

        for (size_t row = 0; row < n; ++row) {
            for (size_t col = 0; col < m; ++col) {
                result[row][col] = matrix[row][col];
            }
        }

        for (size_t row = 0; row < n; ++row) {
            for (size_t col = m + 1; col < Col; ++col) {
                result[row][col - 1] = matrix[row][col];
            }
        }

        for (size_t row = n + 1; row < Row; ++row) {
            for (size_t col = 0; col < m; ++col) {
                result[row - 1][col] = matrix[row][col];
            }
        }

        for (size_t row = n + 1; row < Row; ++row) {
            for (size_t col = m + 1; col < Col; ++col) {
                result[row - 1][col - 1] = matrix[row][col];
            }
        }

        return result;
    }

    double determinant(const Matrix &matrix) {
        if (matrix.empty()) {
            return 1;
        }
        const size_t Row = matrix.size();
        const size_t Col = matrix[0].size();
        if (Row != Col) {
            throw std::logic_error("non-square matrices have no determinant");
        }
        const size_t N = Row;
        if (N == 1) {
            return matrix[0][0];
        }
        if (N == 2) {
            return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        }
        double result = 0;
        int t = 1;
        for (size_t count = 0; count < N; ++count) {
            result += t * matrix[0][count] * determinant(minor(matrix, 0, count));
            t *= -1;
        }
        return result;
    }

    Matrix inverse(const Matrix &matrix) {
        size_t n = matrix.size();
        if (n == 0) {
            return Matrix{};
        }
        if (matrix[0].size() != n or determinant(matrix) == 0) {
            throw std::logic_error("Matrix must be square.");
        }

        Matrix AugmentedMatrix = matrix;
        for (size_t i = 0; i < n; ++i) {
            AugmentedMatrix[i].resize(2 * n, 0);
            AugmentedMatrix[i][n + i] = 1;
        }

        for (size_t i = 0; i < n; ++i) {
            if (AugmentedMatrix[i][i] == 0.0) {
                for (size_t k = i + 1; k < n; ++k) {
                    if (AugmentedMatrix[k][i] != 0.0) {
                        std::swap(AugmentedMatrix[i], AugmentedMatrix[k]);
                        break;
                    }
                }
            }

            double pivot = AugmentedMatrix[i][i];
            for (size_t j = i; j < 2 * n; ++j) {
                AugmentedMatrix[i][j] /= pivot;
            }

            for (size_t k = 0; k < n; ++k) {
                if (k != i) {
                    double factor = AugmentedMatrix[k][i];
                    for (size_t j = i; j < 2 * n; ++j) {
                        AugmentedMatrix[k][j] -= factor * AugmentedMatrix[i][j];
                    }
                }
            }
        }

        Matrix result(n, std::vector<double>(n));
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                result[i][j] = AugmentedMatrix[i][n + j];
            }
        }

        return result;
    }

    Matrix concatenate(const Matrix &matrix1, const Matrix &matrix2, int axis) {
        if (matrix1.empty() or matrix2.empty()) {
            return Matrix{};
        }
        if (axis == 0) {
            if (matrix1[0].size() != matrix2[0].size()) {
                throw std::logic_error("Matrix dimensions cannot match.");
            }
        } else {
            if (matrix1.size() != matrix2.size()) {
                throw std::logic_error("Matrix dimensions cannot match.");
            }
        }
        size_t Row = matrix1.size();
        size_t Col = matrix1[0].size();
        Row += (axis xor 1) * matrix2.size();
        Col += (axis xor 0) * matrix2[0].size();

        Matrix result = zeros(Row, Col);
        for (size_t r = 0; r < matrix1.size(); ++r) {
            for (size_t c = 0; c < matrix1[0].size(); ++c) {
                result[r][c] = matrix1[r][c];
            }
        }

        size_t row_start = axis == 0 ? matrix1.size() : 0;
        size_t col_start = axis == 0 ? 0 : matrix1[0].size();
        for (size_t r = row_start; r < Row; ++r) {
            for (size_t c = col_start; c < Col; ++c) {
                result[r][c] = matrix2[r - row_start][c - col_start];
            }
        }

        return result;
    }

    Matrix ero_swap(const Matrix &matrix, size_t r1, size_t r2) {
        if (r1 >= matrix.size() or r2 >= matrix.size()) {
            throw std::logic_error("out of range.");
        }
        Matrix result = matrix;
        std::swap(result[r1], result[r2]);
        return result;
    }

    Matrix ero_multiply(const Matrix &matrix, size_t r, double c) {
        if (r >= matrix.size()) {
            throw std::logic_error("out of range.");
        }
        Matrix result = matrix;
        for (auto &elem: result[r]) {
            elem *= c;
        }
        return result;
    }

    Matrix ero_sum(const Matrix &matrix, size_t r1, double c, size_t r2) {
        if (r1 >= matrix.size() or r2 >= matrix.size()) {
            throw std::logic_error("out of range.");
        }
        Matrix result = matrix;
        const size_t Col = matrix[0].size();
        for (size_t col = 0; col < Col; ++col) {
            result[r2][col] += matrix[r1][col] * c;
        }
        return result;
    }

    Matrix upper_triangular(const Matrix &matrix) {
        if (matrix.empty()) {
            return matrix;
        }
        const size_t Row = matrix.size();
        const size_t Col = matrix[0].size();
        if (Row != Col) {
            throw std::logic_error("non-square matrices have no upper triangular form");
        }

        Matrix result = matrix;
        for (size_t i = 0; i < Row; ++i) {
            for (size_t j = 0; j < i; ++j) {
                double n = result[i][j];
                double m = result[j][j];
                if (n != 0) {
                    if (m == 0) {
                        std::swap(result[i], result[j]);
                        continue;
                    }
                    double factor = n / m;
                    for (size_t index = j; index < Col; ++index) {
                        result[i][index] -= factor * result[j][index];
                    }
                }
            }
        }
        return result;
    }
}