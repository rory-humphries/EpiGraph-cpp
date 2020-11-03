//
// Created by roryh on 22/10/2020.
//

#ifndef EPIGRAPH_CPP_EPICOMP_HPP
#define EPIGRAPH_CPP_EPICOMP_HPP

#include <EpiGraph/Core/Util.hpp>

#include <Eigen/Dense>

namespace EpiGraph {
    class EpiComp;

    template<>
    struct traits<EpiComp> {
        using state_scalar_type = double;
        using param_scalar_type = double;
    };

    template<typename Derived, typename StateMat, typename ParamMat>
    class CompParamBase {
    public:
        using state_type = Eigen::Matrix<typename StateMat::Scalar, StateMat::RowsAtCompileTime, StateMat::ColsAtCompileTime>;
        using param_type = Eigen::Matrix<typename StateMat::Scalar, StateMat::RowsAtCompileTime, StateMat::ColsAtCompileTime>;

        auto state() const -> const state_type & {
            return m_state;
        }

        auto params() const -> const state_type & {
            return m_state;
        }

        auto set_state(const state_type &x) -> void {
            static_cast<Derived*>(this)->set_state(x);
        }

        auto set_compartment(Eigen::Index comp, double val) -> void {
            static_cast<Derived*>(this)->set_compartment(comp, val);
        }

        auto set_params(const param_type &params) -> void {
            static_cast<Derived*>(this)->set_params(params);
        }

        auto set_params(Eigen::Index i, double val) -> void {
            static_cast<Derived*>(this)->set_params(i, val);
        }

    private:
        state_type m_state;
        param_type m_params;
    };

    class EpiComp {
    public:
        using state_scalar_type = typename traits<EpiComp>::state_scalar_type;
        using state_type = Eigen::Matrix<state_scalar_type, Eigen::Dynamic, 1>;

        using param_scalar_type = typename traits<EpiComp>::param_scalar_type;
        using param_type = Eigen::Matrix<param_scalar_type, Eigen::Dynamic, 1>;

        explicit EpiComp(Eigen::Index num_comps) : m_state(state_type::Zero(num_comps, 1)),
                                                   m_params(param_type::Zero(num_comps, 1)) {}

        auto state() const -> const state_type & {
            return m_state;
        }

        auto params() const -> const state_type & {
            return m_state;
        }

        auto set_state(const state_type &x) -> void {
            if ((x.array() < 0).any())
                throw std::domain_error("Encountered negative value");
            else
                m_state = x;
        }

        auto set_compartment(Eigen::Index comp, double val) -> void {
            if (val < 0)
                throw std::domain_error("Encountered negative number");
            else
                m_state(comp) = val;
        }

        auto set_params(const EpiComp::param_type &params) -> void {
            if ((params.array() < 0).any() || (params.array() > 1).any())
                throw std::domain_error("Encountered value greater than 1 or less than 0");
            else
                m_params = params;
        }

        auto set_params(Eigen::Index i, double val) -> void {
            if ((val < 0) || (val > 1))
                throw std::domain_error("Encountered value greater than 1 or less than 0");
            else
                m_params(i) = val;
        }

    private:
        state_type m_state;
        param_type m_params;
    };
}

#endif //EPIGRAPH_CPP_EPICOMP_HPP
