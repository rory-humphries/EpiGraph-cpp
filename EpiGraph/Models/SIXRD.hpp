//
// Created by roryh on 19/06/2020.
//

#ifndef EPIGRAPH_SIRX_NETWORK_H
#define EPIGRAPH_SIRX_NETWORK_H

#include <Eigen/Dense>
#include <EpiGraph/EigenUtil/Spectral.hpp>
#include <iostream>
#include <map>

namespace EpiGraph {

enum SixrdId : Eigen::Index { S, I, X, R, D };

enum SixrdParamId : Eigen::Index { beta, c, mu, alpha, kappa };

/**
 * @brief ODE which dercribes the SIXRD model
 *
 * @tparam Mat1 A column or row vector. Inherits from MatrixBase
 * @tparam Mat2 A column or row vector. Inherits from MatrixBase
 * @param x The state vector. A column or row vector of length 5. The  entries
 * index corresponds to 0=S, 1=I, 2=X, 3=R, 4=D.
 * @param dxdt The derivative vector. A column or row vector of length 5. The
 * entries index corresponds to 0=dS/dt, 1=dI/dt, 2=dX/dt, 3=dR/dt, 4=dD/dt.
 * @param params The parameter vector. The  entries index corresponds to
 * 0=beta, 1=c, 2=mu, 3=alpha, 4=kappa.
 */
template <IsVector Mat1, IsVector Mat2>
auto sixrd_ode(const Mat1 &x, Mat1 &dxdt, const Mat2 &params) -> void {
  double beta = params[0];
  double c = params[1];
  double mu = params[2];
  double alpha = params[3];
  double kappa = params[4];

  Eigen::VectorXd n = x.sum();

  dxdt[4] = alpha * x[1] + alpha * x[2];             // D
  dxdt[3] = mu * x[1] + mu * x[2];                   // R
  dxdt[2] = kappa * x[1] - mu * x[2] - alpha * x[2]; // X
  dxdt[1] =
      beta * c * x[0] * x[1] / n - mu * x[1] - alpha * x[1] - kappa * x[1]; // I
  dxdt[0] = -beta * c * x[0] * x[1] / n;                                    // S
}

/**
 * @brief ODE which dercribes the SIXRD model on a network
 *
 * @tparam Mat1
 * @tparam Mat2
 * @tparam Mat3
 * @param xmat
 * @param adj
 * @param params
 * @return Mat1
 */
template <IsMatrix Mat1, IsSparseOrDenseMatrix Mat2, IsMatrix Mat3>
auto sixrd_network_ode(const Mat1 &xmat, const Mat2 &adj, const Mat3 &params)
    -> Mat1 {
  /*
   * xmat is expected to be an n x 5 matrix (columns checked at compile time)
   * where the following column indices represent 0=S, 1=I, 2=X, 3=R, 4=D. The
   * rows represent the nodes of the system.
   *
   * adj is the couplings in the system representing the movements of
   * individual between nodes. If xmat is an n x 5 matrix then adj must be a n
   * x n matrix.
   */
  // SIXRDState_assert(xmat);

  Eigen::VectorXd n = xmat.rowwise().sum();

  Mat2 S = (xmat.col(0).cwiseProduct(n.cwiseInverse())).asDiagonal() * adj;
  Mat2 I = (xmat.col(1).cwiseProduct(n.cwiseInverse())).asDiagonal() * adj;

  Eigen::VectorXd new_s = xmat.col(0) - S * Eigen::VectorXd::Ones(S.rows());

  Eigen::VectorXd new_i = xmat.col(1) +
                          (Eigen::RowVectorXd::Ones(I.rows()) * I).transpose() -
                          (I * Eigen::VectorXd::Ones(I.cols()));

  n += (Eigen::RowVectorXd::Ones(adj.rows()) * adj).transpose() -
       (adj * Eigen::VectorXd::Ones(adj.cols()));

  Eigen::VectorXd idiff_ndiff_inv = new_i.cwiseProduct(n.cwiseInverse());

  Mat1 dxdt(xmat.rows(), xmat.cols());

  dxdt.col(4) = params.col(alpha).cwiseProduct(xmat.col(1) + xmat.col(2)); // D
  dxdt.col(3) = params.col(mu).cwiseProduct(xmat.col(1) + xmat.col(2));    // R
  dxdt.col(2) = params.col(kappa).cwiseProduct(xmat.col(1)) -
                params.col(mu).cwiseProduct(xmat.col(2)) -
                params.col(alpha).cwiseProduct(xmat.col(2)); // X

  dxdt.col(1) =
      params.col(beta).cwiseProduct(
          params.col(c).cwiseProduct(new_s.cwiseProduct(idiff_ndiff_inv))) +
      (S * (params.col(beta).cwiseProduct(params.col(c)).asDiagonal()) *
       idiff_ndiff_inv) -
      params.col(mu).cwiseProduct(xmat.col(1)) -
      params.col(alpha).cwiseProduct(xmat.col(1)) -
      params.col(kappa).cwiseProduct(xmat.col(1)); // I

  dxdt.col(0) =
      -params.col(beta).cwiseProduct(
          params.col(c).cwiseProduct(new_s.cwiseProduct(idiff_ndiff_inv))) -
      (S * (params.col(beta).cwiseProduct(params.col(c)).asDiagonal()) *
       idiff_ndiff_inv); // S

  return dxdt;
}

/**
 * @brief ODE which dercribes the SIXRD model on a network where the parameters
 * are the same for each node
 *
 * @tparam DerivedA inherits from MatrixBase
 * @tparam DerivedB inherits from MatrixBase or SparseMatrixBase
 * @tparam DerivedC inherits from MatrixBase
 * @param x is expected to be an n x 5 matrix (columns checked at compile time)
 * where the following column indices represent 0=S, 1=I, 2=X, 3=R, 4=D. The
 * rows represent the nodes of the system.
 * @param adj is the couplings in the system representing the movements of
 * individuals between nodes. If xmat is an n x 5 matrix then adj must be a n x
 * n matrix.
 * @param sixrd_params
 * @return DerivedA
 */
template <IsMatrix Mat1, IsSparseOrDenseMatrix Mat2, IsVector Mat3>
auto sixrd_network_uniform_params_ode(const Mat1 &x, const Mat2 &adj,
                                      const Mat3 &params) -> Mat1 {

  static_assert(
      std::is_base_of<Eigen::MatrixBase<Mat2>, Mat2>::value ||
          std::is_base_of<Eigen::SparseMatrixBase<Mat2>, Mat2>::value,
      "DerivedB must inherit from either MatrixBase or SparseMatrixBase");

  double beta = params[0];
  double c = params[1];
  double mu = params[2];
  double alpha = params[3];
  double kappa = params[4];

  const Mat2 &adj_d = adj.derived();

  Eigen::VectorXd n = x.rowwise().sum();

  Mat2 S = (x.col(0).cwiseProduct(n.cwiseInverse())).asDiagonal() * adj_d;
  Mat2 I = (x.col(1).cwiseProduct(n.cwiseInverse())).asDiagonal() * adj_d;

  Eigen::VectorXd new_s = x.col(0) - S * Eigen::VectorXd::Ones(S.rows());

  Eigen::VectorXd new_i = x.col(1) +
                          (Eigen::RowVectorXd::Ones(I.rows()) * I).transpose() -
                          (I * Eigen::VectorXd::Ones(I.cols()));

  n += (Eigen::RowVectorXd::Ones(adj_d.rows()) * adj_d).transpose() -
       (adj_d * Eigen::VectorXd::Ones(adj_d.cols()));

  Eigen::VectorXd idiff_ndiff_inv = new_i.cwiseProduct(n.cwiseInverse());

  Mat1 dxdt(x.rows(), x.cols());

  dxdt.col(4) = alpha * x.col(1) + alpha * x.col(2);                 // D
  dxdt.col(3) = mu * x.col(1) + mu * x.col(2);                       // R
  dxdt.col(2) = kappa * x.col(1) - mu * x.col(2) - alpha * x.col(2); // X

  dxdt.col(1) = beta * c * new_s.cwiseProduct(idiff_ndiff_inv) +
                beta * c * S * idiff_ndiff_inv - mu * x.col(1) -
                alpha * x.col(1) - kappa * x.col(1); // I

  dxdt.col(0) = -beta * c * new_s.cwiseProduct(idiff_ndiff_inv) -
                beta * c * S * idiff_ndiff_inv; // S

  return dxdt;
}

/**
 * @brief Get the next generation matrix for the SIXRD model on a network with
 * same parameters for each node.
 *
 * @tparam Mat1
 * @tparam Mat2
 * @tparam Mat3
 * @param xmat
 * @param adj
 * @param params
 * @return Mat2
 */
template <IsMatrix Mat1, IsSparseOrDenseMatrix Mat2, IsVector Mat3>
auto sixrd_network_uniform_params_next_gen_matrix(const Mat1 &xmat,
                                                  const Mat2 &adj,
                                                  const Mat3 &params) -> Mat2 {

  int dim = xmat.rows();

  double alpha = params[SixrdParamId::alpha];
  double beta = params[SixrdParamId::beta];
  double c = params[SixrdParamId::c];
  double mu = params[SixrdParamId::mu];
  double kappa = params[SixrdParamId::kappa];

  const Mat2 &adj_d = adj.derived();

  Eigen::VectorXd n = xmat.rowwise().sum().array();
  Eigen::VectorXd n_inv = n.cwiseInverse();

  Eigen::VectorXd one_vec = Eigen::VectorXd::Ones(dim);
  Eigen::RowVectorXd one_rowvec = Eigen::RowVectorXd::Ones(dim);

  Eigen::Matrix<double, Eigen::Dynamic, 1> S = xmat.col(SixrdId::S);
  Eigen::VectorXd Nprop = (adj_d * one_vec).array() / n.array();

  Eigen::VectorXd n_diff =
      n + (one_rowvec * adj_d).transpose() - (adj_d * one_vec);
  Eigen::VectorXd n_diff_inv = n_diff.cwiseInverse();

  Eigen::VectorXd mat1 = S.array() - S.array() * Nprop.array();

  Mat2 mat2((one_vec.array() - ((adj_d * one_vec).array() / n.array()))
                .matrix()
                .asDiagonal());
  mat2 += Mat2(((n.cwiseInverse().asDiagonal()) * adj_d).transpose());
  mat2 = n_diff_inv.asDiagonal() * mat2;

  Mat2 mat3 = Mat2(adj_d * mat2);
  mat3 = (S.cwiseProduct(n_inv)).asDiagonal() * mat3;
  // DerivedB T = (mat1.asDiagonal() * mat2) + mat3;
  Mat2 T =
      ((beta * c) / (mu + alpha + kappa)) * ((mat1.asDiagonal() * mat2) + mat3);
  // DerivedB T = beta*c*(mat1.asDiagonal() * mat2 + mat3) -
  // DerivedB(one_vec.asDiagonal()*(mu+alpha+kappa));

  return T;
}

/**
 * @brief Computes the reproduction number for the SIXRD model on a network with
 * the same parameters for each node from the spectral radius of the next
 * generation matrix.
 *
 * @tparam Mat1
 * @tparam Mat2
 * @tparam Mat3
 * @param xmat
 * @param adj
 * @param param
 * @return double
 */
template <IsMatrix Mat1, IsSparseOrDenseMatrix Mat2, IsVector Mat3>
auto sixrd_network_uniform_params_r0(const Mat1 &xmat, const Mat2 &adj,
                                     const Mat3 &params) -> double {

  Mat2 T = sixrd_next_gen_matrix(xmat, adj, params);
  return SpectralRadius(T);
}

} // namespace EpiGraph

#endif // EPIGRAPH_SIRX_NETWORK_H
