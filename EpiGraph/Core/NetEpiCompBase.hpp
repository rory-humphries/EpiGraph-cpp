//
// Created by roryh on 02/11/2020.
//

#ifndef EPIGRAPH_CPP_NETEPICOMPBASE_HPP
#define EPIGRAPH_CPP_NETEPICOMPBASE_HPP

#include <EpiGraph/Core/Util.hpp>

#include <Eigen/Dense>

namespace EpiGraph {
    template<typename Derived>
    class NetEpiCompBase {
    public:
        using state_type = typename traits<Derived>::state_type;
        using param_type = typename traits<Derived>::param_type;
        using coupling_type = typename traits<Derived>::coupling_type;

        auto state() const -> const state_type & { return m_state; }

        auto params() const -> const param_type & { return m_params; }

        auto coupling() const -> const coupling_type & { return m_coupling; }

        auto set_state(const state_type &x) -> void {
            return static_cast<Derived *>(this)->set_state(x);
        }

        auto set_state(Eigen::Index comp, const Eigen::VectorXd &vec) -> void {
            return static_cast<Derived *>(this)->set_state(comp, vec);
        }

        auto set_state(Eigen::Index comp, Eigen::Index i, double val) -> void {
            return static_cast<Derived *>(this)->set_state(i, val);
        }

        auto move_state(Eigen::Index comp1, Eigen::Index comp2, Eigen::Index i, double val) -> void {
            return static_cast<Derived *>(this)->set_state(comp1, comp2, i, val);
        }

        auto set_params(const param_type &params) -> void {
            return static_cast<Derived *>(this)->set_state(params);
        }

        auto set_params(Eigen::Index i, double val) -> void {
            return static_cast<Derived *>(this)->set_state(i, val);
        }

        auto set_params(Eigen::Index i, Eigen::Index j, double val) -> void {
            return static_cast<Derived *>(this)->set_state(i, j, val);
        }

        auto
        set_params(Eigen::Index i, const Eigen::Matrix<typename param_type::Scalar, Eigen::Dynamic, 1> &vec) -> void {
            return static_cast<Derived *>(this)->set_state(i, vec);
        }

        auto set_coupling(const coupling_type &coup) -> void {
            return static_cast<Derived *>(this)->set_state(coup);
        }

    protected:
        state_type m_state;
        param_type m_params;
        coupling_type m_coupling;
    };
}

#endif //EPIGRAPH_CPP_NETEPICOMPBASE_HPP
