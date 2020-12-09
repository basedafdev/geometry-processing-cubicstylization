
#ifndef CUBIC_STYLE_PRECOMPUTE_H
#define CUBIC_STYLE_PRECOMPUTE_H
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>
#include <cubic_style_data.h>
namespace igl
{
    template <typename T>
    struct cubic_style_data;
}
// Precompute data needed to efficiently conduct local-global iterations for
// cubic stylization.
//
// Inputs:
//   V  #V by 3 list of vertex positions
//   F  #F by 3 list of triangle indices into the rows of V
//   data struct that may contain additional input constraints etc
// Outputs:
//   data struct that contains all the precomputed information for cubic stylization
void cubic_style_precompute(
    const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &F,
    const Eigen::VectorXi &b,
    cubic_style_data &data);

#endif