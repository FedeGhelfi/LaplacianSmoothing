#include "massmatrix.h"
#include "igl/doublearea.h"
#include <iostream>

void massmatrix(
  const Eigen::MatrixXd &l,     // l =  rows = faces, cols = length/2 edges. #Faces * 3 
  const Eigen::MatrixXi &F,     // F =  rows = faces, cols = vertex of the faces. #Faces * 3 
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> &M)
{

  int num_vertex = F.maxCoeff()+1;   // number of V

  Eigen::MatrixXd matrixArea = Eigen::MatrixXd(F.rows(),1); //  #F * 1: Area of each face


  // fill the area matrix
  for(int i = 0; i < l.rows(); ++i) {
      double semi_p = l(i,0) + l(i,1) + l(i,2);
      double area = sqrt(semi_p*(semi_p - (l(i,0)*2)) * (semi_p - (l(i,1)*2)) * (semi_p - l(i,2)*2));

      matrixArea(i,0) = area;
  }

  // * TEST PER CONFRONTARE: corretto! *
  /*
    verifico se il mio calcolo delle aree di ogni faccia è corretto
    confronto con doubleArea(lenght,A) che mi riempie A con le aree.
  */
/*
  Eigen::MatrixXd matrixArea2 = Eigen::MatrixXd(F.rows(),1);
  igl::doublearea(l*2,matrixArea2);
  matrixArea2 = matrixArea2 / 2.0;

  std::cout << "Scorro i primi 10 elementi Matrix Area di Libigl" << std::endl;
  for (int i = 0; i < 10; ++i) {
    std::cout << matrixArea2(i) << std::endl;
  } 

  std::cout << "Scorro i primi 10 elementi Matrix Area di Fede" << std::endl;
  for (int i = 0; i < 10; ++i) {
    std::cout << matrixArea(i) << std::endl;
  }
*/

  int vertex_a, vertex_b, vertex_c;
  M.setZero(num_vertex);

  for (int i = 0; i < F.rows(); ++i){

    // vertex of the i_th face
     vertex_a = F(i,0);
     vertex_b = F(i,1);
     vertex_c = F(i,2);
    
    // contribute of this face
    double face_area = matrixArea(i);


    M.diagonal()[vertex_a] += (face_area / 3.0);
    M.diagonal()[vertex_b] += (face_area / 3.0);
    M.diagonal()[vertex_c] += (face_area / 3.0);

  }


  /*
  Area of triangle = √[s(s – a)(s – b)(s – c)],
   where s is the semi-perimeter of the triangle, and a, b, c are lengths of the three sides of the triangle.
  */



}

