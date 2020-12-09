
#ifndef CUBIC_STYLE_LOCAL_STEP_H
#define CUBIC_STYLE_LOCAL_STEP_H
#include <Eigen/Core>
#include <cubic_style_data.h>
#include <Eigen/Sparse>
namespace igl
{
    template <typename T>
    struct cubic_style_data;
}
// Given the list of vertices and cubic parameters conduct a local step and return
// the rotation matrices pervertex i
//
//
// Inputs:
//   data struct that contains all the precomputed information for cubic stylization
//   U  #V by 3 list of current stylized mesh vertex positions
// Outputs:
//   R  #3V by 3 list of new stylized mesh vertex positions
void cubic_style_local_step(
    const cubic_style_data &data,
    const Eigen::MatrixXd &U,
    Eigen::MatrixXd &R);
#endif