//
// Created by roryh on 22/10/2020.
//
#include <epigraph/Models/Sixrd.hpp>

auto Sixrd::state() const -> const state_type & {
    return m_state;
}

auto Sixrd::params() const -> const state_type & {
    return m_state;
}

auto Sixrd::set_state(const Sixrd::state_type &x) -> void {
    if ((x.array() < 0).any())
        throw std::domain_error("Encountered negative value");
    else
        m_state = x;
}

auto Sixrd::set_compartment(Sixrd::state_id comp, double val) -> void {
    if (val < 0)
        throw std::domain_error("Encountered negative number");
    else
        m_state(comp) = val;
}

auto Sixrd::set_params(const Sixrd::param_type &params) -> void {
    if ((params.array() < 0).any() || (params.array() > 1).any())
        throw std::domain_error("Encountered value greater than 1 or less than 0");
    else
        m_params = params;
}

auto Sixrd::set_params(Sixrd::param_id i, double val) -> void {
    if ((val < 0) || (val > 1))
        throw std::domain_error("Encountered value greater than 1 or less than 0");
    else
        m_params(i) = val;
}
