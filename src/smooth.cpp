#include "smooth.h"
#include <igl/edge_lengths.h>
#include "massmatrix.h"
#include "cotmatrix.h"
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>

void smooth(
    const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &F,
    const Eigen::MatrixXd &G,
    double lambda,
    Eigen::MatrixXd &U)
{

  Eigen::MatrixXd l = Eigen::MatrixXd(F.rows(), 3);
  igl::edge_lengths(V, F, l);
  l = l / 2; // edge_lenghts / 2

  Eigen::DiagonalMatrix<double, Eigen::Dynamic> M = Eigen::DiagonalMatrix<double, Eigen::Dynamic>(V.rows());

  /*
  M =
    matrice diag, la posizione ii è 1/3 dell'area delle facce a cui è associato quel vertice
    per ogni vertice, si associano le aree delle facce in cui esso è incluso
    */

  massmatrix(l, F, M);
  Eigen::MatrixXd dense_M;
  dense_M = Eigen::MatrixXd(M);

  // dbg
  std::cout << "ho fatto massmatrix" << std::endl;

  // cot matrix
  Eigen::SparseMatrix<double> L = Eigen::SparseMatrix<double>(V.rows(), V.rows());
  cotmatrix(l, F, L);
  std::cout << "ho fatto cotmatrix" << std::endl;

  Eigen::MatrixXd dense_L;
  dense_L = Eigen::MatrixXd(L);


  // metodo 2
  Eigen::SparseMatrix<double> sparse_M;
  sparse_M = Eigen::SparseMatrix<double>(M);

  // Solve (M-delta*L) U = M*U
  //  const auto & S = (dense_M - 0.001* L);  original
  const auto &S = (sparse_M - lambda * L);
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver(S);
  assert(solver.info() == Eigen::Success);
  U = solver.solve(sparse_M * U).eval();

  std::cout << "Calcolo svolto..." << std::endl;
}

bool compare(double x, double y, double epsilon = 0.0000001f)
{
  if (fabs(x - y) < epsilon)
    return true; // they are same

  return false; // they are not same
}

/*
* test svolti

- confronto massMatrix con igl_massMatrix
- confronto cotMatrix con igl_cotMatrix


* test mass matrix

    Eigen::SparseMatrix<double> MassMatrix = Eigen::SparseMatrix<double>();
    igl::MassMatrixType type = igl::MASSMATRIX_TYPE_BARYCENTRIC;
    igl::massmatrix(V, F, type, MassMatrix);


    for (int i = 0; i < 10; ++i)
    {
      std::cout << "myFunct: " << M.diagonal()[i] << std::endl;
      std::cout << "LibrFunct: " << MassMatrix.coeff(i, i) << std::endl;
    }



  *test cotmatrix
  Eigen::SparseMatrix<double> L_igl = Eigen::SparseMatrix<double>(V.rows(), V.rows());
  igl::cotmatrix(V, F, L_igl);
  std::cout << "ho fatto igl cotmatrix" << std::endl;

  for (int i = 0; i < L.rows(); i++)
  {
    for (int j = 0; j < L.cols(); j++)
    {
      if (!compare(L.coeff(i, j), L_igl.coeff(i, j), 0.0000001f))
      {
        std::cout << "Elementi diversi: " << L.coeff(i, j) << " // " << L_igl.coeff(i, j);
      }
    }
  }
  std::cout << "Ho finito il confronto" << std::endl;

*/