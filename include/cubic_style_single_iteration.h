
#ifndef CUBIC_STYLE_SINGLE_ITERATION_H
#define CUBIC_STYLE_SINGLE_ITERATION_H
#include <Eigen/Core>
#include <Eigen/Sparse>

#include "cubic_style_data.h"
// Given precomputed data and current positions of all vertices 'U', conduct a single
// iteration of the local-global solver for minimizing the cubic stylization energy.
// Output the positions of all vertices of the stylized mesh by overwriting 'U'
//
// Inputs:
//   data struct that contains all the precomputed information for cubic stylization
//   U  #V by 3 list of current stylized mesh vertex positions
// Outputs:
//   U  #V by 3 list of new stylized mesh vertex positions
void cubic_style_single_iteration(
    const cubic_style_data &data,
    const Eigen::MatrixXd &bc,
    Eigen::MatrixXd &U);
#endif