//
// Created by roryh on 20/10/2020.
//

#ifndef EPIGRAPH_SIXRD_MODEL_H
#define EPIGRAPH_SIXRD_MODEL_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <epigraph/sixrd_odes.hpp>

class SIXRDNetModelUP {
public:
    using state_type = Eigen::Matrix<double, Eigen::Dynamic, 5>;
    using param_type = Eigen::Matrix<double, 5, 1>;
    using coupling_type = Eigen::SparseMatrix<double>;

    enum StateId : Eigen::Index {
        Sidx, Iidx, Xidx, Ridx, Didx
    };

    enum ParamId : Eigen::Index {
        beta_idx, c_idx, mu_idx, alpha_idx, kappa_idx
    };

    explicit SIXRDNetModelUP(Eigen::Index dim) : m_state(state_type::Zero(dim, 5)),
                                                 m_params(param_type::Zero(dim)),
                                                 m_coupling(dim, dim) {}

    auto state() const -> const state_type & { return m_state; }

    auto params() const -> const param_type & { return m_params; }

    auto coupling() const -> const coupling_type & { return m_coupling; }

    auto set_state(const state_type &x) -> void;

    auto set_compartment(StateId comp, const Eigen::VectorXd &vec) -> void;

    auto set_compartment(StateId comp, Eigen::Index i, double val) -> void;

    auto add_infected(Eigen::Index i, double N = 1) -> void;

    auto set_params(const param_type &params) -> void;

    auto set_params(ParamId i, double val) -> void;

    auto set_coupling(const coupling_type &coup) -> void;

private:
    state_type m_state;
    param_type m_params;
    coupling_type m_coupling;
};

auto derivative(const SIXRDNetModelUP &model) -> SIXRDNetModelUP::state_type;

auto print_totals(const SIXRDNetModelUP &model) -> void;

auto write_state(const SIXRDNetModelUP &model, const std::string &id, std::string dir) -> void;

auto write_state_totals(const SIXRDNetModelUP &model, const std::string &path_to_file,
                        bool append_file) -> void;


#endif //EPIGRAPH_SIXRD_MODEL_H
