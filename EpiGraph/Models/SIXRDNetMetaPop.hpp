//
// Created by roryh on 02/11/2020.
//

#ifndef EPIGRAPH_CPP_SIXRDNETMETAPOP_HPP
#define EPIGRAPH_CPP_SIXRDNETMETAPOP_HPP

#include <EpiGraph/Core/NetMetaPop.hpp>
#include <EpiGraph/Models/SIXRD.hpp>

namespace EpiGraph {

template <typename StateMat, typename ParamType, typename CoupType>
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

template <typename StateMat, typename ParamType, typename CoupType,
          typename std::enable_if<ParamType::IsVectorAtCompileTime == 1>::type
              * = nullptr>
auto dXdt(const SIXRDNetMetaPop<StateMat, ParamType, CoupType> &model)
    -> StateMat {
  return sixrd_net_ode_uni_params(model.state(), model.coupling(),
                                  model.params());
}

template <typename StateMat, typename ParamType, typename CoupType,
          typename std::enable_if<ParamType::IsVectorAtCompileTime != 1>::type
              * = nullptr>
auto dXdt(const SIXRDNetMetaPop<StateMat, ParamType, CoupType> &model)
    -> StateMat {
  return sixrd_net_ode(model.state(), model.coupling(), model.params());
}

// only defined for uniform parameters
template <typename StateMat, typename ParamType, typename CoupType,
          typename std::enable_if<ParamType::IsVectorAtCompileTime == 1>::type
              * = nullptr>
auto next_gen_matrix(
    const SIXRDNetMetaPop<StateMat, ParamType, CoupType> &model) -> CoupType {
  return sixrd_net_next_gen_matrix(model.state(), model.coupling(),
                                   model.params());
}

// only defined for uniform parameters
template <typename StateMat, typename ParamType, typename CoupType,
          typename std::enable_if<ParamType::IsVectorAtCompileTime == 1>::type
              * = nullptr>
auto r0(const SIXRDNetMetaPop<StateMat, ParamType, CoupType> &model) -> double {
  return sixrd_net_r0(model.state(), Eigen::MatrixXd(model.coupling()),
                      model.params().transpose());
}

} // namespace EpiGraph

#endif // EPIGRAPH_CPP_SIXRDNETMETAPOP_HPP
