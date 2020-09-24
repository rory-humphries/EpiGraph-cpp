//
// Created by roryh on 19/06/2020.
//


#ifndef EPIGRAPH_SIRX_NETWORK_H
#define EPIGRAPH_SIRX_NETWORK_H

#include <map>
#include <iostream>
#include "spatial_util.h"
#include "eigen_util.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>
//#include "spectra/include/Spectra/GenEigsSolver.h"
//#include "spectra/include/Spectra/MatOp/SparseGenMatProd.h"


template<typename Scalar>
using SIXRDState = Eigen::Matrix<Scalar, Eigen::Dynamic, 5>;

template<typename Derived>
auto SIXRDState_assert(const Eigen::MatrixBase<Derived> &mat) -> void {
    static_assert(Derived::ColsAtCompileTime == 5, "Expected a matrix with 5 cols");
}

enum SIXRD_idx : Eigen::Index {
    Sidx, Iidx, Xidx, Ridx, Didx
};

enum SIXRD_param_idx : Eigen::Index {
    beta_idx, c_idx, mu_idx, alpha_idx, kappa_idx
};

template<typename Derived>
auto infect_SIXRD_state(Eigen::MatrixBase<Derived> &x, Eigen::Index i, double N = 1) -> void {
    SIXRDState_assert(x);

    double change = x(i, 0) - (0 > x(i, 0) - N ? 0 : x(i, 0) - N);
    x(i, 0) -= change;
    x(i, 1) += change;
}

template<typename Derived>
auto print_SIXRD_totals(Eigen::MatrixBase<Derived> &x) -> void {
    SIXRDState_assert(x);

    Eigen::RowVectorXd op_vec = x.colwise().sum();
    std::cout << "S : " << op_vec[Sidx];
    std::cout << ", I : " << op_vec[Iidx];
    std::cout << ", X : " << op_vec[Xidx];
    std::cout << ", R : " << op_vec[Ridx];
    std::cout << ", D : " << op_vec[Didx] << "\n";
}

template<typename DerivedA, typename DerivedB, typename DerivedC>
auto
net_SIXRD_ode(Eigen::MatrixBase<DerivedA> &xmat, Eigen::EigenBase<DerivedB> &adj_,
              Eigen::MatrixBase<DerivedC> &sixrd_params) -> void {
    /*
     * xmat is expected to be an n x 5 matrix (columns checked at compile time) where the following column indices
     * represent 0=S, 1=I, 2=X, 3=R, 4=D. The rows represent the nodes of the system.
     *
     * adj is the couplings in the system representing the movements of individual between nodes. If xmat is an n x 5
     * matrix then adj must be a n x n matrix.
     */
    col_vector_assert(sixrd_params);
    SIXRDState_assert(xmat);

    double beta = sixrd_params[0];//param.beta;
    double c = sixrd_params[1];//param.c;
    double mu = sixrd_params[2];//param.mu;
    double alpha = sixrd_params[3];//param.alpha;
    double kappa = sixrd_params[4];//param.kappa;

    const DerivedB &adj = adj_.derived();

    Eigen::VectorXd n = xmat.rowwise().sum();

    DerivedB S = (xmat.col(0).cwiseProduct(n.cwiseInverse())).asDiagonal() * adj;
    DerivedB I = (xmat.col(1).cwiseProduct(n.cwiseInverse())).asDiagonal() * adj;

    Eigen::VectorXd new_s = xmat.col(0) - S * Eigen::VectorXd::Ones(S.rows());

    Eigen::VectorXd new_i =
            xmat.col(1) + (Eigen::RowVectorXd::Ones(I.rows()) * I).transpose() - (I * Eigen::VectorXd::Ones(I.cols()));


    n += (Eigen::RowVectorXd::Ones(adj.rows()) * adj).transpose() - (adj * Eigen::VectorXd::Ones(adj.cols()));

    Eigen::VectorXd idiff_ndiff_inv = new_i.cwiseProduct(n.cwiseInverse());

    xmat.col(4) += alpha * xmat.col(1) + alpha * xmat.col(2);// D
    xmat.col(3) += mu * xmat.col(1) + mu * xmat.col(2);// R
    xmat.col(2) += kappa * xmat.col(1) - mu * xmat.col(2) - alpha * xmat.col(2); // X

    xmat.col(1) += beta * c * new_s.cwiseProduct(idiff_ndiff_inv) +
                   beta * c * S * idiff_ndiff_inv -
                   mu * xmat.col(1) -
                   alpha * xmat.col(1) -
                   kappa * xmat.col(1);// I

    xmat.col(0) += -beta * c * new_s.cwiseProduct(idiff_ndiff_inv) -
                   beta * c * S * idiff_ndiff_inv; // S

};


template<typename Derived>
void write_net_SIXRD_state(const Eigen::MatrixBase<Derived> &x, const std::string &id, std::string dir) {
    SIXRDState_assert(x);

    std::ofstream myfile;

    myfile.open(dir + id + ".csv");

    myfile << "S,";
    myfile << "I,";
    myfile << "X,";
    myfile << "R,";
    myfile << "D,";
    myfile << "N\n";

    for (int i = 0; i < x.rows(); i++) {
        myfile << x(i, Sidx) << ",";
        myfile << x(i, Iidx) << ",";
        myfile << x(i, Xidx) << ",";
        myfile << x(i, Ridx) << ",";
        myfile << x(i, Didx) << ",";
        myfile << x.row(i).sum() << "\n";
    }
    myfile.close();
}

template<typename Derived>
void write_net_SIXRD_state_totals(const Eigen::MatrixBase<Derived> &x, const std::string &path_to_file,
                                  bool append_file = false) {

    double sumi = 0;
    double sums = 0;
    double sumr = 0;
    double sumx = 0;
    double sumd = 0;

    for (int i = 0; i < x.rows(); i++) {
        sums += x(i, Sidx);
        sumi += x(i, Iidx);
        sumx += x(i, Xidx);
        sumr += x(i, Ridx);
        sumd += x(i, Didx);
    }

    std::ofstream myfile;

    if (append_file)
        myfile.open(path_to_file, std::ios_base::app);
    else {
        myfile.open(path_to_file);
        myfile << "S,";
        myfile << "I,";
        myfile << "X,";
        myfile << "R,";
        myfile << "D\n";
    }

    myfile << sums << ",";
    myfile << sumi << ",";
    myfile << sumx << ",";
    myfile << sumr << ",";
    myfile << sumd << "\n";

    myfile.close();
}

/*
template<typename DerivedA, typename DerivedB, typename DerivedC>
auto
net_SIXRD_R0(const Eigen::MatrixBase<DerivedA> &xmat, const Eigen::EigenBase<DerivedB> &adj_,
             const Eigen::MatrixBase<DerivedC> &sixrd_params) -> double {

    /*
     * Find the reproduction number for the network SIXRD model
     */
/*
    static_assert(5 == DerivedA::ColsAtCompileTime, "Input matrix must have 5 cols");

    int dim = xmat.rows();

    double beta = sixrd_params[0];
    double c = sixrd_params[1];
    double mu = sixrd_params[2];
    double alpha = sixrd_params[3];
    double kappa = sixrd_params[4];

    const DerivedB &adj = adj_.derived();

    Eigen::VectorXd n = xmat.rowwise().sum();
    Eigen::VectorXd n_inv = n.cwiseInverse();

    DerivedB S = (xmat.col(Sidx).cwiseProduct(n_inv)).asDiagonal() * adj;

    Eigen::VectorXd one_vec = Eigen::VectorXd::Ones(dim);
    Eigen::RowVectorXd one_rowvec = Eigen::RowVectorXd::Ones(dim);

    Eigen::VectorXd n_diff = n + (one_rowvec * adj).transpose() - (adj * one_vec);
    Eigen::VectorXd n_diff_inv = n_diff.cwiseInverse();


    DerivedB mat = (n_diff_inv.asDiagonal() *
                    DerivedB((one_vec - (DerivedB(n_inv.asDiagonal() * adj)) * one_vec).asDiagonal()) +
                    DerivedB(adj.transpose() * n_inv.asDiagonal()));


    DerivedB T = ((beta * c) / (mu + alpha + kappa)) * DerivedB(xmat.col(Sidx).asDiagonal()) *
                 DerivedB(n_inv.asDiagonal());// + (S * mat));

    T.makeCompressed();

    // Construct matrix operation object using the wrapper class
    Spectra::SparseGenMatProd<double> op(T);
    // Construct eigen solver object, requesting the largest
    // (in magnitude, or norm) three eigenvalues
    Spectra::GenEigsSolver<double, Spectra::LARGEST_MAGN, Spectra::SparseGenMatProd<double> > eigs(&op, 1, 6);
    // Initialize and compute
    eigs.init();
    int nconv = eigs.compute();
    // Retrieve results
    Eigen::VectorXcd evalues;
    if (eigs.info() == Spectra::SUCCESSFUL) {
        evalues = eigs.eigenvalues();
        return evalues.cwiseAbs().maxCoeff();
    } else {
        throw std::runtime_error("Could not find eigenvalue");
    }

}

*/

template<typename DerivedA, typename DerivedB>
auto SIXRD_R0(Eigen::MatrixBase<DerivedA> &xmat, Eigen::MatrixBase<DerivedB> &sixrd_params) -> double {
    /*
     * Find the reproduction number for the SIXRD model
     */
    col_vector_assert(xmat);
    col_vector_assert(sixrd_params);

    using Mat = Eigen::Matrix2d;

    double beta = sixrd_params[0];
    double c = sixrd_params[1];
    double mu = sixrd_params[2];
    double alpha = sixrd_params[3];
    double kappa = sixrd_params[4];

    Mat T = Mat::Zero(2, 2);
    T(0, 0) = beta * c * xmat(Sidx) / xmat.sum();

    Mat E = Mat::Zero(2, 2);
    E(0, 0) = -mu - alpha - kappa;
    E(1, 0) = kappa;
    E(1, 1) = -mu - alpha;

    Mat E_inv = E.inverse();

    Mat prod = -T * E_inv;

    Eigen::EigenSolver<Mat> eigensolver;
    eigensolver.compute(prod);
    Eigen::VectorXd eigen_values = eigensolver.eigenvalues().cwiseAbs();

    return eigen_values.maxCoeff();
};

#endif //EPIGRAPH_SIRX_NETWORK_H
