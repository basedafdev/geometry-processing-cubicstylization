
#include "../include/cubic_style_single_iteration.h"
#include "../include/cubic_style_local_step.h"
#include "igl/cotmatrix.h"
#include "igl/massmatrix.h"
#include "igl/min_quad_with_fixed.h"

void cubic_style_single_iteration(
    const cubic_style_data &data,
    const Eigen::MatrixXd &bc,
    Eigen::MatrixXd &U)
{

    Eigen::MatrixXd R;
    cubic_style_local_step(data, U, R);
    Eigen::SparseMatrix<double> L;
    L = data.L;
    Eigen::VectorXd Beq;
    Eigen::MatrixXd B = data.K.transpose() * R;
    Beq.resize(U.rows());
    igl::min_quad_with_fixed_solve(data.q_data, B, bc, Beq, U);
}