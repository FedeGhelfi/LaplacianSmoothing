#include "cotmatrix.h"
#include <igl/cotmatrix.h>
#include <cmath>

void cotmatrix(
    const Eigen::MatrixXd &l,
    const Eigen::MatrixXi &F,
    Eigen::SparseMatrix<double> &L)
{
  int first_edge_first_index = 0;
  int first_edge_second_index = 0;
  double first_edge_length = 0;

  int second_edge_first_index = 0;
  int second_edge_second_index = 0;
  double second_edge_length = 0;

  int third_edge_first_index = 0;
  int third_edge_second_index = 0;
  double third_edge_length = 0;

  double cos_alpha;
  double alpha;
  double half_cot_alpha;

  double cos_beta;
  double beta;
  double half_cot_beta;

  double cos_gamma;
  double gamma;
  double half_cot_gamma;

  L.setZero();

  /*
  for each triangle let's find every edges' vertex and
  the edges' length. Now let's find all the angle 
  through the cosine law. 
  */

  for (int i = 0; i < F.rows(); ++i)
  {

    first_edge_first_index = F(i, 0);
    first_edge_second_index = F(i, 1);

    first_edge_length = l(i, 2) * 2.0;


    second_edge_first_index = F(i, 1);
    second_edge_second_index = F(i, 2);

    second_edge_length = l(i, 0) * 2.0;


    third_edge_first_index = F(i, 2);
    third_edge_second_index = F(i, 0);

    third_edge_length = l(i, 1) * 2.0;

    /*
      Find the cosine of the angle opposite the first edge
      using the cosine law: c^2 = a^2 + b^2 - 2abcos(y)
    */

    cos_alpha = (std::pow(first_edge_length, 2) - std::pow(second_edge_length, 2) - std::pow(third_edge_length, 2)) /
                (-2 * second_edge_length * third_edge_length);

    alpha = acos(cos_alpha);
    half_cot_alpha = 0.5 * (1.0 / tan(alpha));

    L.coeffRef(first_edge_first_index, first_edge_second_index) += half_cot_alpha;
    L.coeffRef(first_edge_second_index, first_edge_first_index) += half_cot_alpha;
    L.coeffRef(first_edge_second_index, first_edge_second_index) -= half_cot_alpha;
    L.coeffRef(first_edge_first_index, first_edge_first_index) -= half_cot_alpha;

    /*
      Same thing with the second edge:
       b^2 = a^2 + c^2 - 2accos(y)
    */

    cos_beta = (std::pow(second_edge_length, 2) - std::pow(first_edge_length, 2) - std::pow(third_edge_length, 2)) /
               (-2 * first_edge_length * third_edge_length);

    beta = acos(cos_beta);
    half_cot_beta = 0.5 * (1.0 / tan(beta));

    L.coeffRef(second_edge_first_index, second_edge_second_index) += half_cot_beta;
    L.coeffRef(second_edge_second_index, second_edge_first_index) += half_cot_beta;
    L.coeffRef(second_edge_second_index, second_edge_second_index) -= half_cot_beta;
    L.coeffRef(second_edge_first_index, second_edge_first_index) -= half_cot_beta;

    /*
      same thing with the third edge:
      a^2 = b^2 + c^2 - 2bccos(y)

    */
    cos_gamma = (std::pow(third_edge_length, 2) - std::pow(first_edge_length, 2) - std::pow(second_edge_length, 2)) /
                (-2 * first_edge_length * second_edge_length);

    gamma = acos(cos_gamma);
    half_cot_gamma = 0.5 * (1.0 / tan(gamma));

    L.coeffRef(third_edge_first_index, third_edge_second_index) += half_cot_gamma;
    L.coeffRef(third_edge_second_index, third_edge_first_index) += half_cot_gamma;
    L.coeffRef(third_edge_second_index, third_edge_second_index) -= half_cot_gamma;
    L.coeffRef(third_edge_first_index, third_edge_first_index) -= half_cot_gamma;
  }
}
