

#ifndef CUBIC_STYLE_DATA_H
#define CUBIC_STYLE_DATA_H
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <igl/min_quad_with_fixed.h>

/*
Parameters for computation
*/
struct cubic_style_data
{
    double rho_init; // rho is the penalty cond
    double u_init;   //scaled dual variable
    double e_abs;
    double e_rel = 0;
    double mu = 0;
    double t_incr = 0;
    double t_decr = 0;
    double lambda = .0001; // cubeness parameter
    Eigen::MatrixXd bc;
    //cubness factro
    Eigen::SparseMatrix<double> L;                // Cotangeant matrix
    Eigen::SparseMatrix<double> K;                // ARAP Precompute matrix
    Eigen::SparseMatrix<double> D;                // |V| x 3|V| sparse matrix of rim/sopke edge vectors
    Eigen::SparseMatrix<double> W;                // |V| x |V| matrix of such that W(i) is the diagonal of a diagonal matrix of cotangant weights
    Eigen::MatrixXd VertexArea;                   //diagonal of  Mass matrix barycentric area for vertex i
    Eigen::MatrixXd N;                            // |V| x 3 matrix of per-vertex-normals
    Eigen::MatrixXd V;                            // Extra copy in case
    Eigen::MatrixXi F;                            // Extra copy in case
    Eigen::MatrixXd Z;                            // |V| x 3 matrix of minimized z terms using shrinkage from lasoo
    int ADMM_iterations;                          // Number of maximum iterations for ADMM
    igl::min_quad_with_fixed_data<double> q_data; //
    Eigen::VectorXd rho_list;                     // |V| x 1 matrix of the rho values
    Eigen::VectorXd u_list;                       // |V| x 1 of the  of the scaled dual variable
    std::vector<std::vector<int>> vertex_adjacent_faces;
};

#endif