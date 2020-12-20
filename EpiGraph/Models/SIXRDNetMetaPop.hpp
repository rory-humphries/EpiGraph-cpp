//
// Created by roryh on 02/11/2020.
//

#ifndef EPIGRAPH_CPP_SIXRDNETMETAPOP_HPP
#define EPIGRAPH_CPP_SIXRDNETMETAPOP_HPP

#include <EpiGraph/Core/NetMetaPop.hpp>
#include <EpiGraph/EigenUtil/StaticAsserts.hpp>
#include <EpiGraph/Models/SIXRD.hpp>

namespace EpiGraph {

template <IsMatrix StateMat, IsMatrix ParamType, IsSparseOrDenseMatrix CoupType>
class SIXRDNetMetaPop : public NetMetaPop<StateMat, ParamType, CoupType> {
public:
  using Base = NetMetaPop<StateMat, ParamType, CoupType>;
  using state_type = typename Base::state_type;
  using param_type = typename Base::param_type;
  using coupling_type = typename Base::coupling_type;

  explicit SIXRDNetMetaPop(Eigen::Index dim)
      : NetMetaPop<StateMat, ParamType, CoupType>(dim, 5, 5) {}

private:
  using Base::m_coupling;
  using Base::m_params;
  using Base::m_state;
};

template <IsMatrix StateMat, IsVector ParamType, IsSparseOrDenseMatrix CoupType>
auto dXdt(const SIXRDNetMetaPop<StateMat, ParamType, CoupType> &model)
    -> StateMat {
  return sixrd_network_uniform_params_ode(model.state(), model.coupling(),
                                          model.params());
}

template <IsMatrix StateMat, IsMatrix ParamType, IsSparseOrDenseMatrix CoupType>
auto dXdt(const SIXRDNetMetaPop<StateMat, ParamType, CoupType> &model)
    -> StateMat {
  return sixrd_network_ode(model.state(), model.coupling(), model.params());
}

// only defined for uniform parameters
template <IsMatrix StateMat, IsVector ParamType, IsSparseOrDenseMatrix CoupType>
auto next_gen_matrix(
    const SIXRDNetMetaPop<StateMat, ParamType, CoupType> &model) -> CoupType {
  return sixrd_network_uniform_params_next_gen_matrix(
      model.state(), model.coupling(), model.params());
}

// only defined for uniform parameters
template <IsMatrix StateMat, IsVector ParamType, IsSparseOrDenseMatrix CoupType>
auto r0(const SIXRDNetMetaPop<StateMat, ParamType, CoupType> &model) -> double {
  return sixrd_network_uniform_params_r0(model.state(),
                                         Eigen::MatrixXd(model.coupling()),
                                         model.params().transpose());
}

} // namespace EpiGraph

#endif // EPIGRAPH_CPP_SIXRDNETMETAPOP_HPP
