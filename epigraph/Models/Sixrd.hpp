//
// Created by roryh on 22/10/2020.
//

#ifndef EPIGRAPH_CPP_SIXRD_HPP
#define EPIGRAPH_CPP_SIXRD_HPP

#include <epigraph/Core/EpiCompBase.hpp>

class Sixrd;

template<>
struct traits<Sixrd> {
    using state_type = Eigen::Matrix<double, 5, 1>;
    using param_type = Eigen::Matrix<double, 5, 1>;

    enum state_id : Eigen::Index {
        Sidx, Iidx, Xidx, Ridx, Didx
    };

    enum param_id : Eigen::Index {
        beta_idx, c_idx, mu_idx, alpha_idx, kappa_idx
    };
};

class Sixrd : public EpiCompBase<Sixrd> {
public:
    using state_type = typename traits<Sixrd>::state_type;
    using param_type = typename traits<Sixrd>::param_type;
    using state_id = typename traits<Sixrd>::state_id;
    using param_id = typename traits<Sixrd>::param_id;

    explicit Sixrd(Eigen::Index dim) : m_state(state_type::Zero(5, 1)),
                                       m_params(param_type::Zero(5, 1)) {}

    auto state() const -> const state_type &;

    auto params() const -> const param_type &;

    auto set_state(const state_type &x) -> void;

    auto set_compartment(state_id comp, double val) -> void;

    auto set_params(const param_type &params) -> void;

    auto set_params(param_id i, double val) -> void;

private:
    state_type m_state;
    param_type m_params;
};


#endif //EPIGRAPH_CPP_SIXRD_HPP
