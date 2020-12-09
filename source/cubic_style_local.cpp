#include "../include/cubic_style_local_step.h"
#include <iostream>
#include "igl/find.h"
typedef Eigen::Triplet<double> T;

void cubic_style_local_step(
    const cubic_style_data &data,
    const Eigen::MatrixXd &U,
    Eigen::MatrixXd &R)
{

    double FIND = 0.0000001; //I add this so that igl::find will work with edge vectors where one of the components (x,y,z) is 0
    int iterations = data.ADMM_iterations;
    int n = U.rows();
    R.resize(3 * n, n);
    Eigen::SparseMatrix<double> D = data.D;
    Eigen::SparseMatrix<double> W = data.W;
    Eigen::SparseMatrix<double> Dtilde;
    Dtilde.resize(n, 3 * n);
    std::vector<T> triplets_D;

    //Constructing Dtilde
    for (int i = 0; i < U.rows(); i++)
    {
        std::vector<int> adjacent_faces = data.vertex_adjacent_faces[i];
        int count = 0;
        for (int j = 0; j < adjacent_faces.size(); j++)
        {
            int f_idx = adjacent_faces[j];

            Eigen::VectorXi face = data.F.row(f_idx);
            for (int k = 0; k < 3; k++)
            {
                int v1_idx = face[k % 3];
                int v2_idx = face[(k + 1) % 3];
                Eigen::Vector3d v1 = U.row(v1_idx);
                Eigen::Vector3d v2 = U.row(v2_idx);
                Eigen::Vector3d e = v2 - v1;
                triplets_D.push_back(T(i, count * 3, e[0] + FIND));
                triplets_D.push_back(T(i, count * 3 + 1, e[1] + FIND));
                triplets_D.push_back(T(i, count * 3 + 2, e[2] + FIND));
                count++;
            }
        }
        count = 0;
    }

    Dtilde.setFromTriplets(triplets_D.begin(), triplets_D.end());
    //per vertex we need to find the optimal local rotation
    double rho_init = data.rho_init;
    R.resize(n * 3, 3);
    for (int i = 0; i < n; i++)
    {

        //want to minimize z and R in two D_fferent steps then update

        //our parameters

        Eigen::VectorXd z = Eigen::VectorXd::Zero(3);
        Eigen::VectorXd u = Eigen::VectorXd::Zero(3);
        Eigen::RowVectorXd Dtilde_ = Dtilde.row(i);
        Eigen::RowVectorXd D_ = D.row(i);
        Eigen::RowVectorXd W_ = W.row(i);
        double rho = data.rho_init;

        Eigen::VectorXd D_idx;
        Eigen::VectorXd W_idx;
        Eigen::VectorXd Dtilde_idx;

        igl::find(D_, D_idx);
        igl::find(W_, W_idx);
        igl::find(D_, Dtilde_idx);

        Eigen::MatrixXd Di;                                                     // |N(i)| x 3 edge vectors rest state
        Eigen::MatrixXd Wi = Eigen::MatrixXd::Zero(W_idx.rows(), W_idx.rows()); // |N(i)| x |N(i)| diagonal Matrix
        Eigen::MatrixXd Di_tilde;                                               // |N(i)| x 3 edge vectors deformed state

        Di_tilde.resize(W_idx.rows(), 3);
        Di.resize(W_idx.rows(), 3);
        Wi.resize(W_idx.rows(), W_idx.rows());

        //Getting Wi and Di

        for (int a = 0; a < W_idx.rows(); a++)
        {
            int x_idx = D_idx(a * 3);
            int y_idx = D_idx(a * 3 + 1);
            int z_idx = D_idx(a * 3 + 2);

            Eigen::RowVector3d r_D;
            Eigen::RowVector3d r_Dtilde;
            r_D << D_.coeff(x_idx), D_.coeff(y_idx), D_.coeff(z_idx);
            r_Dtilde << Dtilde_.coeff(x_idx), Dtilde_.coeff(y_idx), Dtilde_.coeff(z_idx);
            double c = W_.coeff(W_idx(a));
            Wi(a, a) = c;
            Di.row(a) = r_D;
            Di_tilde.row(a) = r_Dtilde;
        }

        Eigen::MatrixXd D_bar = Di_tilde;
        //ADMM
        Eigen::MatrixXd Ri;
        for (int k = 0; k < iterations; k++)
        {

            //Compose Mi for SVD
            Eigen::MatrixXd left;
            Eigen::MatrixXd mid = Eigen::MatrixXd::Zero(Di.rows() + 1, Di.rows() + 1);
            Eigen::MatrixXd right;
            left.resize(Di.rows() + 1, 3);
            right.resize(Di.rows() + 1, 3);
            left.block(0, 0, Di.rows(), 3) = Di;
            left.row(Di.rows()) = data.N.row(i);

            right.block(0, 0, Di.rows(), 3) = D_bar;
            right.row(Di.rows()) = z - u;

            mid.block(0, 0, Di.rows(), Di.rows()) = Wi;
            mid(Di.rows(), Di.rows()) = rho;
            /*
            std::cout << left.rows() << "by" << left.cols() << std::endl;
            std::cout << mid.rows() << "by" << mid.cols() << std::endl;
            std::cout << right.rows() << "by" << right.cols() << std::endl;
            */
            Eigen::MatrixXd Mi;
            Mi = left.transpose() * mid * right;

            Eigen::JacobiSVD<Eigen::MatrixXd> svd(Mi, Eigen::ComputeFullU | Eigen::ComputeFullV);

            Ri = svd.matrixV() * svd.matrixU().transpose();
            if (Ri.determinant() <= 0)
            {
                //std::cout << Ri.determinant() << std::endl;
                //std::cout << svd.matrixV() << std::endl;
                //std::cout << svd.matrixU().transpose() << std::endl;

                Eigen::Matrix3d U = svd.matrixU().transpose();
                U.row(1) = U.row(1) * -1;
                Ri = svd.matrixV() * U;
            }
            Eigen::Vector3d shrink_input = Ri * data.N.row(i).transpose() + u;

            double shrink_param = data.lambda * data.VertexArea(i) * (1 / rho);

            //apply shrink operator
            for (int j = 0; j < 3; j++)
            {
                double sign = shrink_input(j);
                Eigen::Vector2d higher;
                higher(0) = abs(sign) - shrink_param;
                higher(1) = 0;
                if (sign < 0)
                {
                    z.coeffRef(j) = -1 * higher.maxCoeff();
                }
                if (sign == 0)
                {
                    z.coeffRef(j) = 0 * higher.maxCoeff();
                }
                if (sign > 0)
                {
                    z.coeffRef(j) = 1 * higher.maxCoeff();
                }
            }

            u = u + Ri * data.N.row(i).transpose() - z;

            double convergence = (Ri * Di.transpose() - Di_tilde.transpose()).norm() + data.lambda * data.VertexArea(i) * z.norm();
            //std::cout << "iteration " << k << " Objective value:" << convergence << std::endl;
        }
        R.block(3 * i, 0, 3, 3) = Ri;
    }
}