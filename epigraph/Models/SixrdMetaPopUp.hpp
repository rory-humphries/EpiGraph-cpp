//
// Created by roryh on 20/10/2020.
//

#ifndef EPIGRAPH_SIXRD_MODEL_H
#define EPIGRAPH_SIXRD_MODEL_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <epigraph/SixrdOde.hpp>

class SIXRDMetaPopUP {
public:
    using state_type = Eigen::Matrix<double, Eigen::Dynamic, 5>;
    using param_type = Eigen::Matrix<double, 5, 1>;
    using coupling_type = Eigen::SparseMatrix<double>;

    enum state_id : Eigen::Index {
        Sidx, Iidx, Xidx, Ridx, Didx
    };

    enum param_id : Eigen::Index {
        beta_idx, c_idx, mu_idx, alpha_idx, kappa_idx
    };

    explicit SIXRDMetaPopUP(Eigen::Index dim) : m_state(state_type::Zero(dim, 5)),
                                                m_params(param_type::Zero(dim)),
                                                m_coupling(dim, dim) {}

    auto state() const -> const state_type & { return m_state; }

    auto params() const -> const param_type & { return m_params; }

    auto coupling() const -> const coupling_type & { return m_coupling; }

    auto set_state(const state_type &x) -> void;

    auto set_compartment(state_id comp, const Eigen::VectorXd &vec) -> void;

    auto set_compartment(state_id comp, Eigen::Index i, double val) -> void;

    auto add_infected(Eigen::Index i, double N = 1) -> void;

    auto set_params(const param_type &params) -> void;

    auto set_params(param_id i, double val) -> void;

    auto set_coupling(const coupling_type &coup) -> void;

private:
    state_type m_state;
    param_type m_params;
    coupling_type m_coupling;
};

auto derivative(const SIXRDMetaPopUP &model) -> SIXRDMetaPopUP::state_type;

auto print_totals(const SIXRDMetaPopUP &model) -> void;

auto write_state(const SIXRDMetaPopUP &model, const std::string &id, std::string dir) -> void;

auto write_state_totals(const SIXRDMetaPopUP &model, const std::string &path_to_file,
                        bool append_file) -> void;


#endif //EPIGRAPH_SIXRD_MODEL_H
