

#include "../include/cubic_style_precompute.h"
#include "igl/massmatrix.h"
#include "igl/cotmatrix.h"
#include "igl/arap_rhs.h"
#include "igl/per_vertex_normals.h"
#include "igl/edge_vectors.h"
#include "igl/adjacency_list.h"
#include "igl/vertex_triangle_adjacency.h"
#include "igl/collapse_edge.h"
#include <iostream>
#include "igl/edge_vectors.h"
typedef Eigen::Triplet<double> T;

void cubic_style_precompute(
    const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &F,
    const Eigen::VectorXi &b,
    cubic_style_data &data)
{
    double FIND = 0.0000001; //I add this so that igl::find will work with edge vectors where one of the components (x,y,z) is 0

    Eigen::SparseMatrix<double> L;

    igl::min_quad_with_fixed_data<double> q_data;

    Eigen::SparseMatrix<double> M;
    igl::cotmatrix(V, F, L);
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
    Eigen::MatrixXd normals;
    igl::per_vertex_normals(V, F, normals);
    Eigen::SparseMatrix<double> K;
    igl::arap_rhs(V, F, 3, igl::ARAP_ENERGY_TYPE_SPOKES_AND_RIMS, K);
    Eigen::SparseMatrix<double> aeq;
    igl::min_quad_with_fixed_precompute(L, b, aeq, 0, data.q_data);
    data.lambda = 0.000001;
    data.rho_init = 0.0001;
    data.e_abs = 0.00001;
    data.e_rel = 0.001;
    data.u_init = 10;
    data.t_decr = 2;
    data.t_incr = 2;
    data.ADMM_iterations = 100;
    data.K = K;
    data.L = L;
    data.VertexArea = M.diagonal();
    data.V = V;
    data.F = F;
    data.N = normals;
    std::vector<std::vector<int>> A;
    igl::adjacency_list(F, A, true);
    Eigen::SparseMatrix<double> D;
    D.resize(V.rows(), 3 * V.rows());
    Eigen::SparseMatrix<double> W;
    W.resize(V.rows(), V.rows());
    std::vector<T> triplets_D;
    std::vector<T> triplets_W;
    std::vector<T> adjacencies;
    std::vector<std::vector<int>> vertex_adjacent_faces;
    vertex_adjacent_faces.resize(V.rows());
    //pulled from "vertex_triangle_adjacency.h"
    for (int j = 0; j < F.rows(); j++)
    {
        for (int i = 0; i < 3; i++)
        {
            int v_idx = F(j, i);
            //v_id, i,
            vertex_adjacent_faces[v_idx].push_back(j);
        }
    }
    data.vertex_adjacent_faces = vertex_adjacent_faces;
    int count = 0;
    std::cout << "Step1" << std::endl;

    // Get the edge vectors Spokes and Rims
    for (int i = 0; i < V.rows(); i++)
    {
        std::vector<int> adjacent_faces = vertex_adjacent_faces[i];
        for (int j = 0; j < adjacent_faces.size(); j++)
        {
            int f_idx = adjacent_faces[j];

            Eigen::VectorXi face = F.row(f_idx);
            // Find the edge vectors
            for (int k = 0; k < 3; k++)
            {
                int v1_idx = face[k % 3];
                int v2_idx = face[(k + 1) % 3];

                Eigen::Vector3d v1 = V.row(v1_idx);
                Eigen::Vector3d v2 = V.row(v2_idx);
                Eigen::Vector3d e = v2 - v1;
                triplets_D.push_back(T(i, count * 3, e[0] + FIND));
                triplets_D.push_back(T(i, count * 3 + 1, e[1] + FIND));
                triplets_D.push_back(T(i, count * 3 + 2, e[2] + FIND));
                triplets_W.push_back(T(i, count, L.coeff(v2_idx, v1_idx)));

                count++;
            }
        }

        count = 0;
    }
    /*
    for (int i = 0; i < F.rows(); i++)
    {
        Eigen::Vector3i face = F.row(i);
        //For every vertex on the face
        for (int j = 0; j < 3; j++)
        {
            int v_idx = face(j);

            std::vector<int> adj_vertices = A[v_idx];
            for (int k = 0; k < adj_vertices.size(); k++)
            {
                int adj_v = adj_vertices[k];
                Eigen::RowVector3d e = V.row(adj_v);
                //Constructing D

                triplets_D.push_back(T(v_idx, 3 * adj_v + 0, e[0]));
                triplets_D.push_back(T(v_idx, 3 * adj_v + 1, e[1]));
                triplets_D.push_back(T(v_idx, 3 * adj_v + 2, e[2]));

                //Constructing W

                triplets_W.push_back(T(v_idx, k, L.coeff(v_idx, adj_v)));
            }
        }
    }
    */

    D.setFromTriplets(triplets_D.begin(), triplets_D.end());
    W.setFromTriplets(triplets_W.begin(), triplets_W.end());
    data.D = D;
    data.W = W;
}