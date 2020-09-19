//
// Created by roryh on 19/06/2020.
//


#ifndef EPIGRAPH_SIRX_NETWORK_H
#define EPIGRAPH_SIRX_NETWORK_H

#include <map>
#include "spatial_util.h"
#include "eigen_util.h"

#include <Eigen/Core>
#include "spectra/include/Spectra/GenEigsSolver.h"
#include "spectra/include/Spectra/MatOp/SparseGenMatProd.h"


template<typename Scalar>
using SIXRD_state = Eigen::Matrix<Scalar, Eigen::Dynamic, 5>;

struct SIXRDParam {
    double beta, c, alpha, mu, kappa;
};

enum SIXRD_idx : Eigen::Index {
    Sidx, Iidx, Xidx, Ridx, Didx
};

template<typename Derived>
auto infect_SIXRD_state(Eigen::MatrixBase<Derived> &x, Eigen::Index i, double N = 1) -> void {
    static_assert(5 == Derived::ColsAtCompileTime, "Input matrix must have 5 cols");

    double change = x(i, 0) - (0 > x(i, 0) - N ? 0 : x(i, 0) - N);
    x(i, 0) -= change;
    x(i, 1) += change;
}

template<typename DerivedA, typename DerivedB>
auto
net_SIXRD_ode(Eigen::MatrixBase<DerivedA> &xmat, Eigen::EigenBase<DerivedB> &adj_, const SIXRDParam param) -> void {
    /*
     * xmat is expected to be an n x 5 matrix (columns checked at compile time) where the following column indices
     * represent 0=S, 1=I, 2=X, 3=R, 4=D. The rows represent the nodes of the system.
     *
     * adj is the couplings in the system representing the movements of individual between nodes. If xmat is an n x 5
     * matrix then adj must be a n x n matrix.
     */

    static_assert(5 == DerivedA::ColsAtCompileTime, "Input matrix must have 5 cols");

    double alpha = param.alpha;
    double beta = param.beta;
    double c = param.c;
    double mu = param.mu;
    double kappa = param.kappa;

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
    static_assert(5 == Derived::ColsAtCompileTime, "Input matrix must have 5 cols");

    std::ofstream myfile;
    std::string newd = dir;

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

template<typename DerivedA, typename DerivedB>
auto
net_SIXRD_R0(const Eigen::MatrixBase<DerivedA> &xmat, const Eigen::EigenBase<DerivedB> &adj_,
             const SIXRDParam param) -> double {

    /*
     * Find the reproduction number for the network SIXRD model
     */

    static_assert(5 == DerivedA::ColsAtCompileTime, "Input matrix must have 5 cols");

    int dim = xmat.rows();

    double alpha = param.alpha;
    double beta = param.beta;
    double c = param.c;
    double mu = param.mu;
    double kappa = param.kappa;

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

template<typename DerivedA>
auto SIXRD_R0(Eigen::MatrixBase<DerivedA> &xmat, const SIXRDParam param) -> double {
    /*
     * Find the reproduction number for the SIXRD model
     */

    static_assert(1 == DerivedA::ColsAtCompileTime, "Input matrix must have 1 cols (must be a vector)");

    using Mat = Eigen::Matrix2d;

    double alpha = param.alpha;
    double beta = param.beta;
    double c = param.c;
    double mu = param.mu;
    double kappa = param.kappa;

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

/*
using state_type = std::vector<std::vector<double>>;
using coupling_type = std::map<std::pair<int, int>, double>;

class SIXRDModel {

public:
    using state_type = std::vector<std::vector<double>>;
    using coupling_type = std::map<std::pair<int, int>, double>;

public:
    state_type x;
    dmat xmat;

public:
    // constructors
    SIXRDModel() {}

    SIXRDModel(int dim) : x(dim, std::vector<double>(5, 0)) {}

    auto get_dim() -> size_t {
        return x.size();
    }

    auto set_state(size_t dim, size_t comp, double val) -> void {
        x[dim][comp] = val;
    }

    auto compertment_total(int comp) -> double {
        double sum = 0;
        for (auto &k: x) {
            sum += k[comp];
        }
        return sum;
    }

    auto infect_vertex(Vertex v, double N = 1) -> void {
        double change = x[v][0] - (0 > x[v][0] - N ? 0 : x[v][0] - N);
        x[v][0] -= change;
        x[v][1] += change;
        x[v][2] = 0;
        x[v][3] = 0;
        x[v][4] = 0;
    }

    // update functions
    template<class TNetwork, class Emap>
    auto update(TNetwork &g, Emap &coup, const SIXRDParam param) -> void {
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();


        double alpha = param.alpha;
        double beta = param.beta;
        double c = param.c;
        double mu = param.mu;
        double kappa = param.kappa;

        std::vector<double> tmp(5, 0);
        state_type dxdt(x.size(), tmp);

        //std::cout << dxdt.size() << x.size() << g.num_vertices();


         * The ode describing the evolution of the network_ptr SIXRD system.
         * The paramater x and dxdt should be of the type std::vector<std::vector<double>> where the inner
         * vector is of size 5 such that
         *
         * x[i][0] and dxdt[i][0] corresponds to the variable S
         * x[i][1] and dxdt[i][1] corresponds to the variable I
         * x[i][2] and dxdt[i][2] corresponds to the variable X
         * x[i][3] and dxdt[i][3] corresponds to the variable R
         * x[i][4] and dxdt[i][4] corresponds to the variable D


        for (int v = 0; v < g.num_vertices(); v++) {

            double new_s = x[v][0];
            double new_i = x[v][1];
            double new_n = x[v][0] + x[v][1] + x[v][2] + x[v][3] + x[v][4];

            double n = x[v][0] + x[v][1] + x[v][2] + x[v][3] + x[v][4];
            auto it_pair = g.out_edges(v);
            for (auto it = it_pair.first; it != it_pair.second; it++) {
                double tmp_n = coup[*it];
                new_n -= tmp_n;
                new_i -= tmp_n * x[v][1] / n;
                new_s -= tmp_n * x[v][0] / n;
            }
            it_pair = g.in_edges(v);
            for (auto it = it_pair.first; it != it_pair.second; it++) {
                double tmp_n = coup[*it];
                auto v_src = it->src;
                double n = x[v_src][0] + x[v_src][1] + x[v_src][2] + x[v_src][3] + x[v_src][4];
                new_n += tmp_n;
                new_i += tmp_n * x[v_src][1] / n;
            }

            it_pair = g.in_edges(v);
            for (auto it = it_pair.first; it != it_pair.second; it++) {
                double tmp_n = coup[*it];
                auto v_src = it->src;
                double n = x[v_src][0] + x[v_src][1] + x[v_src][2] + x[v_src][3] + x[v_src][4];

                double tmp_s = tmp_n * x[v_src][0] / n;

                dxdt[it->src][0] -= beta * c * tmp_s * new_i / new_n;
                dxdt[it->src][1] += beta * c * tmp_s * new_i / new_n;

            }

            dxdt[v][0] -= beta * c * new_s * new_i / new_n; // S
            dxdt[v][1] += beta * c * new_s * new_i / new_n - mu * x[v][1] - alpha * x[v][1] - kappa * x[v][1];// I
            dxdt[v][2] += kappa * x[v][1] - mu * x[v][2] - alpha * x[v][2]; // X
            dxdt[v][3] += mu * x[v][1] + mu * x[v][2];// R
            dxdt[v][4] += alpha * x[v][1] + alpha * x[v][2];// D

        }

        for (int count = 0; count < x.size(); count++) {
            x[count][0] += dxdt[count][0];
            x[count][1] += dxdt[count][1];
            x[count][2] += dxdt[count][2];
            x[count][3] += dxdt[count][3];
            x[count][4] += dxdt[count][4];
        }

        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        std::cout << time_span.count() << " seconds for update\n" << std::endl;


    }

    auto update(const state_type &x, state_type &dxdt, const double t) -> void {

         * The ode describing the evolution of the network_ptr SIXRD system.
         * The paramater x and dxdt should be of the type std::vector<std::vector<double>> where the inner
         * vector is of size 5 such that
         *
         * x[i][0] and dxdt[i][0] corresponds to the variable S
         * x[i][1] and dxdt[i][1] corresponds to the variable I
         * x[i][2] and dxdt[i][2] corresponds to the variable X
         * x[i][3] and dxdt[i][3] corresponds to the variable R
         * x[i][4] and dxdt[i][4] corresponds to the variable D

        for (int v = 0; v < (*g).num_vertices(); v++) {

            double new_s = x[v][0];
            double new_i = x[v][1];
            double new_n = (*g).vprop[v].population;


            auto it_pair = (*g).out_edges(v);
            for (auto it = it_pair.first; it != it_pair.second; it++) {
                double tmp_n = (*g).eprop[*it].population;
                new_n -= tmp_n;
                new_i -= tmp_n * x[v][1] / (*g).vprop[v].population;
                new_s -= tmp_n * x[v][0] / (*g).vprop[v].population;
            }
            it_pair = (*g).in_edges(v);
            for (auto it = it_pair.first; it != it_pair.second; it++) {
                double tmp_n = (*g).eprop[*it].population;
                new_n += tmp_n;
                new_i += tmp_n * x[it->src][1] / (*g).vprop[it->src].population;
            }

            it_pair = (*g).in_edges(v);
            for (auto it = it_pair.first; it != it_pair.second; it++) {
                double tmp_n = (*g).eprop[*it].population;
                double tmp_s = tmp_n * x[it->src][0] / (*g).vprop[it->src].population;

                dxdt[it->src][0] -= (*g).vprop[v].beta * (*g).vprop[v].c * tmp_s * new_i / new_n;
                dxdt[it->src][1] += (*g).vprop[v].beta * (*g).vprop[v].c * tmp_s * new_i / new_n;

            }

            dxdt[v][0] -= (*g).vprop[v].beta * (*g).vprop[v].c * new_s * new_i / new_n; // S
            dxdt[v][1] += (*g).vprop[v].beta * (*g).vprop[v].c * new_s * new_i / new_n -
                          (*g).vprop[v].mu * x[v][1] - (*g).vprop[v].alpha * x[v][1] -
                          (*g).vprop[v].kappa * x[v][1];// I
            dxdt[v][2] += (*g).vprop[v].kappa * x[v][1] - (*g).vprop[v].mu * x[v][2] - (*g).vprop[v].alpha * x[v][2]; // X
            dxdt[v][3] += (*g).vprop[v].mu * x[v][1] + (*g).vprop[v].mu * x[v][2];// R
            dxdt[v][4] += (*g).vprop[v].alpha * x[v][1] + (*g).vprop[v].alpha * x[v][2];// D

        }

    }

    // io
    void write_state(const std::string &id, std::string dir) {
        std::ofstream myfile;
        std::string newd = dir;

        myfile.open(dir + id + ".csv");

        myfile << "S,";
        myfile << "I,";
        myfile << "X,";
        myfile << "R,";
        myfile << "D,";
        myfile << "N\n";

        for (int v = 0; v < x.size(); v++) {
            myfile << x[v][0] << ",";
            myfile << x[v][1] << ",";
            myfile << x[v][2] << ",";
            myfile << x[v][3] << ",";
            myfile << x[v][4] << ",";
            myfile << x[v][0] + x[v][1] + x[v][2] + x[v][3] + x[v][4] << "\n";
        }
        myfile.close();
    }

    void write_compartment_totals(const std::string &path_to_file, bool append_file = false) {

        double sumi = 0;
        double sums = 0;
        double sumr = 0;
        double sumx = 0;
        double sumd = 0;

        for (auto &k: x) {
            sums += k[0];
            sumi += k[1];
            sumx += k[2];
            sumr += k[3];
            sumd += k[4];
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

    void print_compartment_totals() {

        double sumi = 0;
        double sums = 0;
        double sumr = 0;
        double sumx = 0;
        double sumd = 0;

        for (auto &k: x) {
            sums += k[0];
            sumi += k[1];
            sumx += k[2];
            sumr += k[3];
            sumd += k[4];
        }

        std::cout << sums << ", " << sumi << ", " << sumx << ", " << sumr << ", " << sumd << std::endl;
    }

    template<class TNetwork, class EPopMap, class VPopMap>
    void
    next_gen_matrix(TNetwork &g, EPopMap &e_pop_map, VPopMap &v_pop_map, const SIXRDParam param, const std::string &id,
                    std::string dir) {

        std::ofstream myfile;
        std::string newd = dir;

        myfile.open(dir + id + ".csv");

        std::vector<double> Nouts(g.num_vertices(), 0);
        std::vector<double> Nins(g.num_vertices(), 0);

#pragma omp parallel for
        for (int vi = 0; vi < g.num_vertices(); vi++) {
            auto it_pair_vi_out = g.out_edges(vi);
            auto it_beg_vi_out = it_pair_vi_out.first;
            auto it_end_vi_out = it_pair_vi_out.second;
            double Ni_out = 0;

            for (; it_beg_vi_out != it_end_vi_out; it_beg_vi_out++) {
                Ni_out += e_pop_map[*it_beg_vi_out];
            }

            auto it_pair_vi_in = g.in_edges(vi);
            auto it_beg_vi_in = it_pair_vi_in.first;
            auto it_end_vi_in = it_pair_vi_in.second;
            double Ni_in = 0;

            for (; it_beg_vi_in != it_end_vi_in; it_beg_vi_in++) {
                Ni_in += e_pop_map[*it_beg_vi_in];
            }
            Nouts[vi] = Ni_out;
            Nins[vi] = Ni_in;
        }

        double a = (param.beta * param.c) / (param.mu * param.alpha * param.kappa);

        std::vector<double> mat(g.num_vertices() * g.num_vertices());

        for (int vi = 0; vi < g.num_vertices(); vi++) {
            std::cout << vi << std::endl;
            double Ni = v_pop_map[vi];
            double Ni_out = Nouts[vi];
            double Ni_in = Nins[vi];
            double Si = x[vi][0];

            for (int vj = 0; vj < g.num_vertices(); vj++) {

                double Nj = v_pop_map[vj];
                double Nj_out = Nouts[vj];
                double Nj_in = Nins[vj];

                auto ep_ij = g.edge(vi, vj);
                auto ep_ji = g.edge(vj, vi);

                double Nij = 0, Nji = 0;
                if (ep_ij.second)
                    Nij = e_pop_map[ep_ij.first];
                if (ep_ji.second)
                    Nji = e_pop_map[ep_ji.first];

                double b = 0;
                double c = 0;
                // check denominators
                if (Ni != 0 && Nj != 0 && Ni + Ni_in + Ni_out != 0) {
                    if (vi != vj)
                        b = Si * (1 - Nij / Ni) * (Nji / Nj) / (Ni + Ni_in - Ni_out);
                    else
                        b = Si * (1 - Nij / Ni) * (1 - Ni_out / Ni + Nji / Nj) / (Ni + Ni_in - Ni_out);
                }


                auto it_pair_vj_out = g.out_edges(vi);
                auto it_beg_vj_out = it_pair_vj_out.first;
                auto it_end_vj_out = it_pair_vj_out.second;
                for (; it_beg_vj_out != it_end_vj_out; it_beg_vj_out++) {

                    auto vk = it_beg_vj_out->src;

                    double Nk = v_pop_map[vk];
                    double Nk_out = Nouts[vk];
                    double Nk_in = Nins[vk];

                    auto ep_ik = g.edge(vi, vk);
                    auto ep_jk = g.edge(vj, vk);
                    auto ep_kj = g.edge(vk, vj);

                    double Nik = 0, Njk = 0, Nkj = 0;
                    if (ep_ik.second)
                        Nik = e_pop_map[ep_ik.first];
                    if (ep_jk.second)
                        Njk = e_pop_map[ep_jk.first];
                    if (ep_kj.second)
                        Njk = e_pop_map[ep_kj.first];

                    if (Ni != 0 && Nj != 0 && Ni + Ni_in + Ni_out != 0) {
                        if (vj != vk)
                            c += Si * (Nik / Ni) * (Njk / Nj) / (Nk + Nk_in - Nk_out);
                        else
                            c += Si * (Nik / Ni) * (1 - Nj_out / Nj + Njk / Nj) / (Nk + Nk_in - Nk_out);
                    }

                }

                mat[g.num_vertices() * vi + vj] = a * (b + c);

                //double FVinv_ij = a*(b+c);

                //myfile << FVinv_ij;
                //if (vj < g.num_vertices() - 1)
                //   myfile << ",";
                //else
                //    myfile << "\n";
            }

        }
        for (int i = 1; i < mat.size(); i++) {
            myfile << mat[i];

            if (i % g.num_vertices() != 0)
                myfile << ",";
            else
                myfile << "\n";
        }

        myfile.close();
    }



    void write_type_totals(const state_type &x, const std::string& id, std::string dir) {

        std::map<int, std::vector<double>> type_map;
        for (int v = 0; v < (*g).num_vertices(); v++) {
            int type = (*g).vprop[v].type;
            if (type_map.count(type)==0)
                type_map[type].resize(5, 0);
            type_map[type][0] += x[v][0];
            type_map[type][1] += x[v][1];
            type_map[type][2] += x[v][2];
            type_map[type][3] += x[v][3];
            type_map[type][4] += x[v][4];
        }

        std::ofstream myfile;
        myfile.open (dir + id + ".csv");

        myfile << "S,";
        myfile << "I,";
        myfile << "X,";
        myfile << "R,";
        myfile << "D,";
        myfile << "N\n";


        for (auto &k : type_map)
        {
            auto vec = k.second;
            myfile << vec[0] <<",";
            myfile << vec[1] <<",";
            myfile << vec[2] <<",";
            myfile << vec[3] <<",";
            myfile << vec[4] <<",";
            myfile << vec[0]+vec[1]+vec[2]+vec[3]+vec[4] << "\n";
        }
        myfile.close();
    }

    auto print_compartment_totals(const state_type& x) {
        double sumi = 0;
        double sums = 0;
        double sumr = 0;
        double sumx = 0;
        double sumd = 0;
        for (auto &k: x) {
            sums += k[0];
            sumi += k[1];
            sumx += k[2];
            sumr += k[3];
            sumd += k[4];
        }

        std :: cout << "S = " << sums << ", I = " << sumi<< ", R = " << sumr << ", X = " << sumx << ", D = " << sumd << ", N = " << sums+sumi+sumx+sumr+sumd << "\n";
    }



};

*/
#endif //EPIGRAPH_SIRX_NETWORK_H
