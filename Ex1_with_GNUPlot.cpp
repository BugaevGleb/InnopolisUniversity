// Gleb Bugaev, g.bugaev@innopolis.university, CS-02

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <cstdio>

#define GNUPLOT_NAME "C:\\gnuplot\\bin\\gnuplot -persist"

using namespace std;

class ColumnVector {
public:

    int n, m;
    vector<vector<double> > matrix;

    ColumnVector() = default;

    ColumnVector(int n, int m, vector<vector<double> > matrix) {
        this->n = n;
        this->m = m;
        this->matrix = matrix;
    }

    ColumnVector(int n, vector<double> b) {
        this->n = n;
        this->m = 1;

        for(int i = 0; i < this->n; ++i) {
            vector<double> vec(this->m, 0);
            this->matrix.push_back(vec);
        }

        for(int i = 0; i < this->n; ++i) {
            for(int j = 0; j < this->m; ++j) {
                this->matrix[i][j] = b[i];
            }
        }

    }

    ~ColumnVector() = default;

    friend istream& operator >> (istream& in, ColumnVector &mat)
    {
        in >> mat.n;
        mat.m = 1;
        for(int i = 0; i < mat.n; ++i) {
            vector<double> vec(mat.m, 0);
            mat.matrix.push_back(vec);
        }
        for(int i = 0; i < mat.n; ++i) {
            for(int j = 0; j < mat.m; ++j) {
                in >> mat.matrix[i][j];
            }
        }
        for(int i = 0; i < mat.n; ++i) {
            for(int j = 0; j < mat.m; ++j) {
            }
        }
        return in;
    }

    friend ostream& operator << (ostream& out, const ColumnVector &mat)
    {
        for(int i = 0; i < mat.n; ++i) {
            out << fixed;
            out << setprecision(4);
            if(abs(mat.matrix[i][0]) < 0.00005) {
                out << 0.00;
            } else {
                out << mat.matrix[i][0];
            }
            for(int j = 1; j < mat.m; ++j) {
                out << fixed;
                out << setprecision(4);
                if(abs(mat.matrix[i][j]) < 0.00005) {
                    out << " " << 0.00;
                } else {
                    out << " " << mat.matrix[i][j];
                }
            }
            out << endl;
        }
        return out;
    }

    ColumnVector& operator= (const ColumnVector& mat) {
        this->n = mat.n;
        this->m = mat.m;
        this->matrix = mat.matrix;
        return *this;
    }

    ColumnVector operator+ (const ColumnVector &mat) {
        vector<vector<double> > sumMatrix;
        for(int i = 0; i < mat.n; ++i) {
            vector<double> vec(mat.m, 0);
            sumMatrix.push_back(vec);
        }
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < m; ++j) {
                sumMatrix[i][j] = this->matrix[i][j] + mat.matrix[i][j];
            }
        }
        return *new ColumnVector(n, m, sumMatrix);
    }

    ColumnVector operator- (const ColumnVector &mat) {
        vector<vector<double> > subMatrix;
        for(int i = 0; i < mat.n; ++i) {
            vector<double> vec(mat.m, 0);
            subMatrix.push_back(vec);
        }
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < m; ++j) {
                subMatrix[i][j] = this->matrix[i][j] - mat.matrix[i][j];
            }
        }
        return *new ColumnVector(n, m, subMatrix);
    }

    ColumnVector operator* (const ColumnVector &mat) {
        vector<vector<double> > mulMatrix;
        for(int i = 0; i < n; ++i) {
            vector<double> vec(mat.m, 0);
            mulMatrix.push_back(vec);
        }
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < mat.m; ++j) {
                for(int k = 0; k < m; ++k) {
                    mulMatrix[i][j] += (this->matrix[i][k] * mat.matrix[k][j]);
                }
            }
        }
        return *new ColumnVector(n, mat.m, mulMatrix);
    }

};

class Matrix {
public:

    int n, m;
    vector<vector<double> > matrix;

    Matrix() = default;

    Matrix(int n, int m, vector<vector<double> > matrix) {
        this->n = n;
        this->m = m;
        this->matrix = matrix;
    }

    Matrix(int n, int m) {
        this->n = n;
        this->m = m;
        matrix = *new vector(n, vector<double>(m, 0));
    }

    Matrix(const Matrix &mat1, const Matrix &mat2) {
        this->n = mat1.n;
        this->m = mat1.m + mat2.m;
        vector<vector<double> > augMatrix;
        for(int i = 0; i < n; ++i) {
            vector<double> vec(m, 0);
            augMatrix.push_back(vec);
        }
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < mat1.m; ++j) {
                augMatrix[i][j] = mat1.matrix[i][j];
            }
            for(int j = 0; j < mat2.m; ++j) {
                augMatrix[i][mat1.m + j] = mat2.matrix[i][j];
            }
        }
        this->matrix = augMatrix;
    }

    ~Matrix() = default;

    friend istream& operator >> (istream& in, Matrix &mat)
    {
        in >> mat.n;
        mat.m = mat.n;
        for(int i = 0; i < mat.n; ++i) {
            vector<double> vec(mat.m, 0);
            mat.matrix.push_back(vec);
        }
        for(int i = 0; i < mat.n; ++i) {
            for(int j = 0; j < mat.m; ++j) {
                in >> mat.matrix[i][j];
            }
        }
        for(int i = 0; i < mat.n; ++i) {
            for(int j = 0; j < mat.m; ++j) {
            }
        }
        return in;
    }

    friend ostream& operator << (ostream& out, const Matrix &mat)
    {
        for(int i = 0; i < mat.n; ++i) {
            out << fixed;
            out << setprecision(4);
            if(abs(mat.matrix[i][0]) < 0.00005) {
                out << 0.00;
            } else {
                out << mat.matrix[i][0];
            }
            for(int j = 1; j < mat.m; ++j) {
                out << fixed;
                out << setprecision(4);
                if(abs(mat.matrix[i][j]) < 0.00005) {
                    out << " " << 0.00;
                } else {
                    out << " " << mat.matrix[i][j];
                }
            }
            out << endl;
        }
        return out;
    }

    Matrix& operator= (const Matrix& mat) {
        this->n = mat.n;
        this->m = mat.m;
        this->matrix = mat.matrix;
        return *this;
    }

    Matrix operator+ (const Matrix &mat) {
        vector<vector<double> > sumMatrix;
        for(int i = 0; i < mat.n; ++i) {
            vector<double> vec(mat.m, 0);
            sumMatrix.push_back(vec);
        }
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < m; ++j) {
                sumMatrix[i][j] = this->matrix[i][j] + mat.matrix[i][j];
            }
        }
        return *new Matrix(n, m, sumMatrix);
    }

    Matrix operator- (const Matrix &mat) {
        vector<vector<double> > subMatrix;
        for(int i = 0; i < mat.n; ++i) {
            vector<double> vec(mat.m, 0);
            subMatrix.push_back(vec);
        }
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < m; ++j) {
                subMatrix[i][j] = this->matrix[i][j] - mat.matrix[i][j];
            }
        }
        return *new Matrix(n, m, subMatrix);
    }

    Matrix operator* (const Matrix &mat) {
        vector<vector<double> > mulMatrix;
        for(int i = 0; i < n; ++i) {
            vector<double> vec(mat.m, 0);
            mulMatrix.push_back(vec);
        }
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < mat.m; ++j) {
                for(int k = 0; k < m; ++k) {
                    mulMatrix[i][j] += (this->matrix[i][k] * mat.matrix[k][j]);
                }
            }
        }
        return *new Matrix(n, mat.m, mulMatrix);
    }

    ColumnVector operator* (const ColumnVector &mat) {
        vector<vector<double> > mulMatrix;
        for(int i = 0; i < n; ++i) {
            vector<double> vec(mat.m, 0);
            mulMatrix.push_back(vec);
        }
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < mat.m; ++j) {
                for(int k = 0; k < m; ++k) {
                    mulMatrix[i][j] += (this->matrix[i][k] * mat.matrix[k][j]);
                }
            }
        }
        return *new ColumnVector(n, mat.m, mulMatrix);
    }

    Matrix transpose () {
        vector<vector<double> > transposeMatrix;
        for(int i = 0; i < m; ++i) {
            vector<double> vec(n, 0);
            transposeMatrix.push_back(vec);
        }
        for(int i = 0; i < m; ++i) {
            for(int j = 0; j < n; ++j) {
                transposeMatrix[i][j] = this->matrix[j][i];
            }
        }
        return *new Matrix(m, n, transposeMatrix);
    }

    Matrix diagonalNormalize () {
        vector<vector<double> > diagonalMatrix;
        for(int i = 0; i < m; ++i) {
            vector<double> vec(n, 0);
            diagonalMatrix.push_back(vec);
        }
        for(int i = 0; i < n; ++i) {
            diagonalMatrix[i][i] = (1 / this->matrix[i][i]);
        }
        return *new Matrix(n, m, diagonalMatrix);
    }

};

class IdentityMatrix: public Matrix {
public:

    int n;
    vector<vector<double> > matrix;

    IdentityMatrix() = default;

    IdentityMatrix(int n, vector<vector<double> > matrix): Matrix(n, n, matrix) {
        this->n = n;
        this->m = m;
        this->matrix = matrix;
    }

    explicit IdentityMatrix(Matrix& mat) {
        this->n = mat.n;

        vector<vector<double> > idMatrix;
        for(int i = 0; i < mat.n; ++i) {
            vector<double> vec(mat.n, 0);
            idMatrix.push_back(vec);
        }

        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < n; ++j) {
                if(i == j) {
                    idMatrix[i][j] = 1;
                } else {
                    idMatrix[i][j] = 0;
                }
            }
        }

        matrix = idMatrix;

        Matrix::matrix = idMatrix;
        Matrix::n = n;
        Matrix::m = n;

    }

    explicit IdentityMatrix(int nn) {

        n = nn;
        vector<vector<double> > idMatrix;
        for(int i = 0; i < n; ++i) {
            vector<double> vec(n, 0);
            idMatrix.push_back(vec);
        }

        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < n; ++j) {
                if(i == j) {
                    idMatrix[i][j] = 1;
                } else {
                    idMatrix[i][j] = 0;
                }
            }
        }

        matrix = idMatrix;

        Matrix::matrix = idMatrix;
        Matrix::n = n;
        Matrix::m = n;

    }

    ~IdentityMatrix() = default;

    Matrix& operator= (const Matrix& mat) {
        this->n = mat.n;
        this->m = mat.m;
        this->matrix = mat.matrix;
        Matrix::matrix = mat.matrix;
        return *this;
    }

};

class EliminationMatrix: public Matrix {
public:

    int n;
    vector<vector<double> > matrix;

    EliminationMatrix() = default;

    explicit EliminationMatrix(Matrix& mat, int nn, int mm) {
        this->n = mat.n;

        vector<vector<double> > elMatrix;
        for(int i = 0; i < mat.n; ++i) {
            vector<double> vec(mat.n, 0);
            elMatrix.push_back(vec);
        }

        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < n; ++j) {
                if(i == j) {
                    elMatrix[i][j] = 1;
                } else {
                    elMatrix[i][j] = 0;
                }
            }
        }
        nn--;
        mm--;
        elMatrix[nn][mm] = - (mat.matrix[nn][mm] / mat.matrix[mm][mm]);

        matrix = elMatrix;

        Matrix::matrix = elMatrix;
        Matrix::n = n;
        Matrix::m = n;
    }

    ~EliminationMatrix() = default;

};

class PermutationMatrix: public Matrix {
public:

    int n;
    vector<vector<double> > matrix;

    PermutationMatrix() = default;

    explicit PermutationMatrix(Matrix& mat, int nn, int mm) {
        this->n = mat.n;

        vector<vector<double> > perMatrix;
        for(int i = 0; i < mat.n; ++i) {
            vector<double> vec(mat.n, 0);
            perMatrix.push_back(vec);
        }

        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < n; ++j) {
                if(i == j) {
                    perMatrix[i][j] = 1;
                } else {
                    perMatrix[i][j] = 0;
                }
            }
        }
        nn--;
        mm--;
        swap(perMatrix[nn], perMatrix[mm]);

        matrix = perMatrix;

        Matrix::matrix = perMatrix;
        Matrix::n = n;
        Matrix::m = n;
    }

    ~PermutationMatrix() = default;

};

Matrix findInverse(Matrix& A1) {

    IdentityMatrix A2 = *new IdentityMatrix(A1);
    double maxx = -100000000;
    int indMaxx = -1;

    for(int i = 0; i < A1.n; ++i) {

        maxx = -100000000;
        indMaxx = -1;
        for(int j = i; j < A1.n; ++j) {
            if(abs(A1.matrix[j][i]) > maxx) {
                maxx = abs(A1.matrix[j][i]);
                indMaxx = j;
            }
        }

        if(indMaxx != i) {
            PermutationMatrix P = *new PermutationMatrix(A1, i + 1, indMaxx + 1);
            Matrix B1 = P * A1;
            Matrix B2 = P * A2;
            A1 = B1;
            A2 = B2;
        }

        for(int j = i + 1; j < A1.n; ++j) {
            if(A1.matrix[j][i] != 0) {
                EliminationMatrix E = *new EliminationMatrix(A1, j + 1, i + 1);
                Matrix B1 = E * A1;
                Matrix B2 = E * A2;
                A1 = B1;
                A2 = B2;
            }
        }

    }

    for(int i = A1.n - 1; i >= 0; --i) {

        for(int j = i - 1; j >= 0; --j) {
            if(A1.matrix[j][i] != 0) {
                EliminationMatrix E = *new EliminationMatrix(A1, j + 1, i + 1);
                Matrix B1 = E * A1;
                Matrix B2 = E * A2;
                A1 = B1;
                A2 = B2;
            }
        }

    }

    Matrix D = A1.diagonalNormalize();
    Matrix B1 = D * A1;
    Matrix B2 = D * A2;
    A1 = B1;
    A2 = B2;
    return A2;

}

void fillMatrix(Matrix & A, vector<double> t) {
    for(int i = 0; i < A.m; ++i) {
        for(int j = 0; j < A.n; ++j) {
            A.matrix[j][i] = pow(t[j], i);
        }
    }
}

double findValue(ColumnVector& X, double t) {
    double y = 0;
    for(int i = 0; i < X.n; ++i) {
        y += (X.matrix[i][0] * pow(t, i));
    }
    return y;
}

int main(){

    int m = 0;
    cin >> m;

    vector<double> t, b;
    int tInput, bInput;

    for(int i = 0; i < m; ++i) {
        cin >> tInput >> bInput;
        t.push_back(tInput);
        b.push_back(bInput);
    }

    int n;
    cin >> n;

    Matrix A = *new Matrix(m, n + 1);
    fillMatrix(A, t);
    cout << "A:" << endl;
    cout << A;

    Matrix A_T = A.transpose();

    Matrix First = A_T * A;
    cout << "A_T*A:" << endl;
    cout << First;

    Matrix Second = findInverse(First);
    cout << "(A_T*A)^-1:" << endl;
    cout << Second;

    ColumnVector B = *new ColumnVector(b.size(), b);

    ColumnVector Third = A_T * B;
    cout << "A_T*b:" << endl;
    cout << Third;

    ColumnVector X = Second * Third;
    cout << "x~:" << endl;
    cout << X;

    FILE* pipe = _popen(GNUPLOT_NAME, "w");

    if(pipe != nullptr) {

        fprintf(pipe, "set terminal wxt size 1000,600\n");

        fprintf(pipe, "plot [-40 : 40] [-200 : 200] ");
        fprintf(pipe, "%lf*x**%d", X.matrix[0][0], 0);
        for(int i = 1; i < X.n; ++i) {
            fprintf(pipe, " + %lf*x**%d", X.matrix[i][0], i);
        }
        fprintf(pipe, " , '-' using 1:2 title 'Given Points' with points\n");

        for(int i = 0; i < m; ++i) {
            fprintf(pipe, "%f\t%f\n", t[i], b[i]);
            cout << t[i] << " " << b[i] << endl;
        }

        fprintf(pipe, "e\n");
        fflush(pipe);
        _pclose(pipe);

    }

    return 0;

}
