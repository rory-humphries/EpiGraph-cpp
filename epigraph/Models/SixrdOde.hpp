//
// Created by roryh on 19/06/2020.
//


#ifndef EPIGRAPH_SIRX_NETWORK_H
#define EPIGRAPH_SIRX_NETWORK_H

#include <Eigen/Core>
#include <map>
#include <iostream>

#include <epigraph/Eigen/EigenUtil.hpp>


enum SixrdId : Eigen::Index {
    S, I, X, R, D
};

enum SixrdParamId : Eigen::Index {
    beta, c, mu, alpha, kappa
};


template<typename Derived>
auto SIXRDState_assert(const Eigen::MatrixBase<Derived> &mat) -> void {
    static_assert(Derived::ColsAtCompileTime == 5, "Expected a matrix with 5 cols");
}

template<typename DerivedA, typename DerivedB, typename DerivedC>
auto
net_SIXRD_ode(const Eigen::MatrixBase<DerivedA> &x, const Eigen::EigenBase<DerivedB> &adj_,
              const Eigen::MatrixBase<DerivedC> &sixrd_params) -> DerivedA {
    /*
     * xmat is expected to be an n x 5 matrix (columns checked at compile time) where the following column indices
     * represent 0=S, 1=I, 2=X, 3=R, 4=D. The rows represent the nodes of the system.
     *
     * adj is the couplings in the system representing the movements of individual between nodes. If xmat is an n x 5
     * matrix then adj must be a n x n matrix.
     */
    //col_vector_assert(sixrd_params);
    //SIXRDState_assert(x);

    double beta = sixrd_params[0];//param.beta;
    double c = sixrd_params[1];//param.c;
    double mu = sixrd_params[2];//param.mu;
    double alpha = sixrd_params[3];//param.alpha;
    double kappa = sixrd_params[4];//param.kappa;

    const DerivedB &adj = adj_.derived();

    Eigen::VectorXd n = x.rowwise().sum();

    DerivedB S = (x.col(0).cwiseProduct(n.cwiseInverse())).asDiagonal() * adj;
    DerivedB I = (x.col(1).cwiseProduct(n.cwiseInverse())).asDiagonal() * adj;

    Eigen::VectorXd new_s = x.col(0) - S * Eigen::VectorXd::Ones(S.rows());

    Eigen::VectorXd new_i =
            x.col(1) + (Eigen::RowVectorXd::Ones(I.rows()) * I).transpose() - (I * Eigen::VectorXd::Ones(I.cols()));


    n += (Eigen::RowVectorXd::Ones(adj.rows()) * adj).transpose() - (adj * Eigen::VectorXd::Ones(adj.cols()));

    Eigen::VectorXd idiff_ndiff_inv = new_i.cwiseProduct(n.cwiseInverse());

    DerivedA dxdt(x.rows(), x.cols());

    dxdt.col(4) = alpha * x.col(1) + alpha * x.col(2);// D
    dxdt.col(3) = mu * x.col(1) + mu * x.col(2);// R
    dxdt.col(2) = kappa * x.col(1) - mu * x.col(2) - alpha * x.col(2); // X

    dxdt.col(1) = beta * c * new_s.cwiseProduct(idiff_ndiff_inv) +
                  beta * c * S * idiff_ndiff_inv -
                  mu * x.col(1) -
                  alpha * x.col(1) -
                  kappa * x.col(1);// I

    dxdt.col(0) = -beta * c * new_s.cwiseProduct(idiff_ndiff_inv) -
                  beta * c * S * idiff_ndiff_inv; // S

    return dxdt;
};

template<typename DerivedA, typename DerivedB, typename DerivedC>
auto
net_SIXRD_ode_inhom(const Eigen::MatrixBase<DerivedA> &xmat, const Eigen::EigenBase<DerivedB> &adj_,
                    const Eigen::MatrixBase<DerivedC> &sixrd_params) -> DerivedA {
    /*
     * xmat is expected to be an n x 5 matrix (columns checked at compile time) where the following column indices
     * represent 0=S, 1=I, 2=X, 3=R, 4=D. The rows represent the nodes of the system.
     *
     * adj is the couplings in the system representing the movements of individual between nodes. If xmat is an n x 5
     * matrix then adj must be a n x n matrix.
     */
    //SIXRDState_assert(xmat);

    const DerivedB &adj = adj_.derived();

    Eigen::VectorXd n = xmat.rowwise().sum();

    DerivedB S = (xmat.col(0).cwiseProduct(n.cwiseInverse())).asDiagonal() * adj;
    DerivedB I = (xmat.col(1).cwiseProduct(n.cwiseInverse())).asDiagonal() * adj;

    Eigen::VectorXd new_s = xmat.col(0) - S * Eigen::VectorXd::Ones(S.rows());

    Eigen::VectorXd new_i =
            xmat.col(1) + (Eigen::RowVectorXd::Ones(I.rows()) * I).transpose() - (I * Eigen::VectorXd::Ones(I.cols()));


    n += (Eigen::RowVectorXd::Ones(adj.rows()) * adj).transpose() - (adj * Eigen::VectorXd::Ones(adj.cols()));

    Eigen::VectorXd idiff_ndiff_inv = new_i.cwiseProduct(n.cwiseInverse());

    DerivedA dxdt(xmat.rows(), xmat.cols());

    dxdt.col(4) = sixrd_params.col(alpha).cwiseProduct(xmat.col(1) + xmat.col(2));// D
    dxdt.col(3) = sixrd_params.col(mu).cwiseProduct(xmat.col(1) + xmat.col(2));// R
    dxdt.col(2) = sixrd_params.col(kappa).cwiseProduct(xmat.col(1)) -
                  sixrd_params.col(mu).cwiseProduct(xmat.col(2)) -
                  sixrd_params.col(alpha).cwiseProduct(xmat.col(2)); // X

    dxdt.col(1) = sixrd_params.col(beta).cwiseProduct(
            sixrd_params.col(c).cwiseProduct(new_s.cwiseProduct(idiff_ndiff_inv))) +
                  (S * (sixrd_params.col(beta).cwiseProduct(sixrd_params.col(c)).asDiagonal()) *
                   idiff_ndiff_inv) -
                  sixrd_params.col(mu).cwiseProduct(xmat.col(1)) -
                  sixrd_params.col(alpha).cwiseProduct(xmat.col(1)) -
                  sixrd_params.col(kappa).cwiseProduct(xmat.col(1));// I

    dxdt.col(0) = -sixrd_params.col(beta).cwiseProduct(
            sixrd_params.col(c).cwiseProduct(new_s.cwiseProduct(idiff_ndiff_inv))) -
                  (S * (sixrd_params.col(beta).cwiseProduct(sixrd_params.col(c)).asDiagonal()) *
                   idiff_ndiff_inv); // S

    return dxdt;
};


#endif //EPIGRAPH_SIRX_NETWORK_H
