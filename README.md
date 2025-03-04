# Conjugate Gradient Method for Solving Linear Systems

This repository contains a C++ implementation of the **Conjugate Gradient (CG) algorithm** for solving systems of linear equations of the form **Ax = b**. The program reads the coefficient matrix `A` and the constant vector `b` from a file named `Ax=b.txt` and computes the solution vector `x` using the CG method.

---

## Table of Contents
1. [Introduction](#introduction)
2. [Features](#features)
3. [Usage](#usage)
4. [Input Format](#input-format)
5. [Dependencies](#dependencies)
6. [License](#license)

---

## Introduction

The Conjugate Gradient method is an iterative algorithm used to solve systems of linear equations, particularly for large, sparse, and symmetric positive-definite matrices. This implementation uses a **10th-order Runge-Kutta integration** for high precision and includes functions for matrix inversion, determinant calculation, and permutation generation.

---

## Features

- **Matrix Inversion**: Computes the inverse of the coefficient matrix `A`.
- **Determinant Calculation**: Calculates the determinant of `A` using permutations and the Levi-Civita symbol.
- **Conjugate Gradient Algorithm**: Solves the system `Ax = b` iteratively until convergence.

---

## Usage

1. **Compile and run the code**:
```bash
g++ -o conjugate_gradient conjugate_gradient.cpp -std=c++11 -O3
```
```sh
./conjugate_gradient
```
---

## Input Format

The input file Ax=b.txt should contain the coefficient matrix A and the constant vector b in the following format:
```sh
a11 a12 a13 ... a1n b1
a21 a22 a23 ... a2n b2
...
an1 an2 an3 ... ann bn
```

---

## Dependencies

C++ Compiler: The code requires a C++ compiler (e.g., g++).

Standard Libraries: The program uses the following C++ standard libraries:

- <iostream> for input/output.
- <fstream> for file handling.
- <cmath> for mathematical functions.

---

## License
This code is released under the MIT License – feel free to use, modify, and distribute with proper attribution.





   
