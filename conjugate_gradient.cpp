//*****************************************************************************************************************************************************************************
// Conjugate Gradient Method for Solving Linear Systems (Ax = b)
// This program implements the Conjugate Gradient algorithm to solve a system of linear equations.
// The matrix A and vector b are read from a file named "Ax=b.txt".
//*****************************************************************************************************************************************************************************

#include<iostream>
#include<fstream>
#include<cmath>

using namespace std;

// Function prototypes
void show(double **, int); // Display a matrix
void initialize(double **&, double **&, double *&, double *&, int &); // Initialize matrices and vectors
void upper_triangular(double **, double *, int); // Convert matrix to upper triangular form
void inverse(double **, double **, int); // Compute the inverse of a matrix
void permute(int *, int *, bool *, int, int); // Generate all permutations of an array
void conjugate_gradient(double **, double **, double *, double *); // Conjugate Gradient algorithm
int levi_civita(int *, int); // Compute the Levi-Civita symbol for a permutation
double determinant(double **, int); // Compute the determinant of a matrix

int main() {
    double **A, **B, *x, *b; // A: Coefficient matrix, B: Inverse of A, x: Solution vector, b: Constant vector
    conjugate_gradient(A, B, b, x); // Solve the system using Conjugate Gradient
    return 0;
}

// Function to initialize matrices and vectors from the input file
void initialize(double **&A, double **&B, double *&b, double *&x, int &n) {
    ifstream read;
    string line;
    int j = 0;

    // Open the file and count the number of lines to determine the size of the matrix
    read.open("Ax=b.txt");
    while (getline(read, line)) j++;
    read.close();

    n = j; // Set the size of the matrix

    // Allocate memory for matrices A and B
    A = new double *[n];
    B = new double *[n];
    for (int i = 0; i < n; i++) {
        A[i] = new double[n];
        B[i] = new double[n];
    }

    // Allocate memory for vectors x and b
    x = new double[n];
    b = new double[n];

    // Initialize the solution vector x to zero
    for (int i = 0; i < n; i++) x[i] = 0;

    // Read the matrix A and vector b from the file
    read.open("Ax=b.txt");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n + 1; j++) {
            if (j == n) read >> b[i]; // Last column is the constant vector b
            else read >> A[i][j]; // First n columns are the coefficient matrix A
        }
    }
    read.close();
}

// Function to convert a matrix to upper triangular form
void upper_triangular(double **A, double *b, int n) {
    double factor;

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            bool is_end_of_loop = false;

            // Handle division by zero by swapping rows
            if (A[i][i] == 0) {
                for (int k = i + 1; k < n; k++) {
                    if (A[k][i] != 0) {
                        b[i] += b[k];
                        for (int z = 0; z < n; z++) A[i][z] += A[k][z];
                        break;
                    }
                    if (k == n - 1) is_end_of_loop = true;
                }
                if (is_end_of_loop == true) continue;
            }

            // Perform row operations to make the matrix upper triangular
            factor = (A[j][i] / A[i][i]);
            b[j] = b[j] - (factor * b[i]);
            for (int h = i; h < n; h++) A[j][h] = A[j][h] - (factor * A[i][h]);
        }
    }
}

// Function to compute the inverse of a matrix
void inverse(double **A, double **B, int n) {
    int m, p;
    double det_A = determinant(A, n), **C;

    // Check if the matrix is invertible
    if (det_A == 0) cout << "Error: 'determinant of matrix A is zero, A is not invertible'" << endl;

    // Allocate memory for the cofactor matrix
    C = new double *[n - 1];
    for (int i = 0; i < n - 1; i++) C[i] = new double[n - 1];

    // Compute the inverse using the adjugate matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            m = 0;
            p = 0;

            // Construct the cofactor matrix
            for (int w = 0; w < n; w++) {
                for (int z = 0; z < n; z++) {
                    if (w != i && z != j) {
                        C[m][p] = A[w][z];
                        p++;
                        if (p == n - 1) {
                            p = 0;
                            m++;
                        }
                    }
                }
            }

            // Compute the inverse using the cofactor matrix and determinant
            B[j][i] = pow(-1, i + j) * determinant(C, n - 1) / det_A;
        }
    }
}

// Function to compute the Levi-Civita symbol for a permutation
int levi_civita(int *e, int n) {
    int temp, permutation = 0, lc;

    // Compute the number of swaps required to sort the permutation
    for (int i = 0; i < n; i++) {
        if (e[i] == i) continue;
        for (int j = 0; j < n; j++) {
            if (e[j] == i) {
                temp = e[i];
                e[i] = i;
                e[j] = temp;
                permutation++;
            }
        }
    }

    // Determine the sign of the permutation
    if (permutation % 2 == 0) lc = 1;
    else lc = -1;

    // Check for repeated elements (Levi-Civita is zero if any two labels are the same)
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (e[i] == e[j]) lc = 0;
        }
    }

    return lc;
}

// Function to compute the determinant of a matrix
double determinant(double **A, int n) {
    fstream read, write;
    double det = 0, a1_an = 1;
    int *e, *e2, j = 0;
    bool *del;

    // Open a file to store permutations
    write.open("permutations.txt", ios::out);

    e = new int[n];
    e2 = new int[n];
    del = new bool[n];

    // Initialize arrays for permutations
    for (int i = 0; i < n; i++) {
        del[i] = false;
        e[i] = i;
    }

    // Generate all permutations and compute the determinant
    permute(e, e2, del, n, 0);

    // Read permutations from the file and compute the determinant
    read.open("permutations.txt", ios::in);
    while (read >> e[j]) {
        if (j == n - 1) {
            for (int i = 0; i < n; i++) a1_an *= A[i][e[i]];
            det += a1_an * levi_civita(e, n);
            j = 0;
            a1_an = 1;
        } else j++;
    }

    read.close();
    return det;
}

// Function to generate all permutations of an array
void permute(int *e, int *e2, bool *del, int n, int l) {
    fstream write("permutations.txt", ios::app);
    string gap = "    ";

    // Recursively generate permutations
    for (int i = 0; i < n; i++) {
        if (del[i] == true) continue;
        e2[l] = e[i];
        if (l == n - 1) {
            for (int k = 0; k < n; k++) write << e2[k] << gap;
            write << endl;
        }
        l++;
        del[i] = true;
        permute(e, e2, del, n, l);
        del[i] = false;
        l--;
    }

    write.close();
}

// Conjugate Gradient algorithm to solve Ax = b
void conjugate_gradient(double **A, double **B, double *b, double *x) {
    double *r, *q, *p, *Ap, a, l, pAp, rq, pp, term = pow(10, -50);
    int n;

    // Initialize matrices and vectors
    initialize(A, B, b, x, n);

    // Allocate memory for residual, search direction, and intermediate vectors
    r = new double[n];
    q = new double[n];
    p = new double[n];
    Ap = new double[n];

    // Compute the inverse of A
    inverse(A, B, n);

    // Initialize residual and search direction
    for (int i = 0; i < n; i++) {
        r[i] = b[i];
        for (int j = 0; j < n; j++) q[i] += B[i][j] * r[j];
        p[i] = q[i];
    }

    // Iterate until convergence
    for (int u = 0;; u++) {
        pAp = 0;
        rq = 0;
        pp = 0;
        l = 0;

        // Compute dot products
        for (int i = 0; i < n; i++) {
            rq += r[i] * q[i];
            q[i] = 0;
            Ap[i] = 0;
        }

        // Compute matrix-vector product Ap
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) Ap[i] += A[i][j] * p[j];
            pAp += p[i] * Ap[i];
        }

        // Compute step size
        a = rq / pAp;

        // Update solution and residual
        for (int i = 0; i < n; i++) {
            x[i] += a * p[i];
            r[i] -= a * Ap[i];
        }

        // Compute new search direction
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) q[i] += B[i][j] * r[j];
            l += r[i] * q[i] / rq;
        }

        for (int i = 0; i < n; i++) p[i] = q[i] + l * p[i];

        // Check for convergence
        for (int i = 0; i < n; i++) pp += p[i] * p[i];
        if (sqrt(fabs(pp)) < term) {
            for (int i = 0; i < n; i++) cout << "x[" << i + 1 << "]= " << x[i] << endl;
            break;
        }
    }
}

// Function to display a matrix
void show(double **A, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << A[i][j] << "        ";
        }
        cout << endl;
    }
}