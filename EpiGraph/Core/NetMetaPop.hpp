//
// Created by roryh on 20/10/2020.
//

#ifndef EPIGRAPH_SIXRD_MODEL_H
#define EPIGRAPH_SIXRD_MODEL_H

#include <EpiGraph/Core/NetEpiCompBase.hpp>
#include <EpiGraph/Core/Util.hpp>

#include <Eigen/Dense>

namespace EpiGraph {

template <typename StateMat, typename ParamType, typename CoupType>
class NetMetaPop;

template <typename StateMat, typename ParamMat, typename CoupMat>
struct traits<NetMetaPop<StateMat, ParamMat, CoupMat>> {
  using state_type = StateMat;
  using param_type = ParamMat;
  using coupling_type = CoupMat;
};

template <typename StateMat, typename ParamType, typename CoupType>
class NetMetaPop
    : public NetEpiCompBase<NetMetaPop<StateMat, ParamType, CoupType>> {
public:
  using Base = NetEpiCompBase<NetMetaPop<StateMat, ParamType, CoupType>>;

  using state_type = typename Base::state_type;
  using param_type = typename Base::param_type;
  using coupling_type = typename Base::coupling_type;

  NetMetaPop(Eigen::Index dim, Eigen::Index compartments,
             Eigen::Index num_params) {
    m_state = state_type::Zero(dim, compartments);
    m_params = param_type::Zero(
        (param_type::IsVectorAtCompileTime == 1 ? 1 : dim), num_params);
    m_coupling = coupling_type(dim, dim);
  }

  auto set_state(const state_type &x) -> void {
    if ((x.array() < 0).any())
      throw std::domain_error("Encountered negative value");
    else
      m_state = x;
  }

  auto set_state(Eigen::Index comp, const Eigen::VectorXd &vec) -> void {
    if ((vec.array() < 0).any())
      throw std::domain_error("Encountered negative number");
    else
      m_state.col(comp) = vec;
  }

  auto set_state(Eigen::Index comp, Eigen::Index i, double val) -> void {
    if (val < 0)
      throw std::domain_error("Encountered negative number");
    else
      m_state(i, comp) = val;
  }

  auto move_state(Eigen::Index comp1, Eigen::Index comp2, Eigen::Index i,
                  double val) -> void {
    using Scalar = typename state_type::Scalar;
    Scalar change =
        m_state(i, comp1) -
        (comp1 > m_state(i, comp1) - val ? comp1 : m_state(i, comp1) - val);
    m_state(i, comp1) -= change;
    m_state(i, comp2) += change;
  }

  auto set_params(const param_type &params) -> void {
    if ((params.array() < 0).any() || (params.array() > 1).any())
      throw std::domain_error(
          "Encountered value greater than 1 or less than 0");
    else
      m_params = params;
  }

  // specialised for m_params being a vector
  template <typename T = param_type>
  auto set_params(Eigen::Index i, double val) ->
      typename std::enable_if<T::IsVectorAtCompileTime == 1, void>::type {
    if ((val < 0) || (val > 1))
      throw std::domain_error(
          "Encountered value greater than 1 or less than 0");
    else
      m_params[i] = val;
  }

  // specialised for m_params being a matrix
  template <typename T = param_type>
  auto set_params(Eigen::Index i, double val) ->
      typename std::enable_if<T::IsVectorAtCompileTime != 1, void>::type {
    if ((val < 0) || (val > 1))
      throw std::domain_error(
          "Encountered value greater than 1 or less than 0");
    else
      m_params.col(i) =
          Eigen::Matrix<typename param_type::Scalar,
                        param_type::RowsAtCompileTime, 1>::Ones(m_params
                                                                    .rows());
  }

  auto set_params(Eigen::Index i, Eigen::Index j, double val) -> void {
    if ((val < 0) || (val > 1))
      throw std::domain_error(
          "Encountered value greater than 1 or less than 0");
    else
      m_params(j, i) = val;
  }

  auto set_params(
      Eigen::Index i,
      const Eigen::Matrix<typename ParamType::Scalar, Eigen::Dynamic, 1> &vec)
      -> void {
    if ((vec.array() < 0).any() || (vec.array() > 1).any())
      throw std::domain_error(
          "Encountered value greater than 1 or less than 0");
    else
      m_params.col(i) = vec;
  }

  auto set_coupling(const coupling_type &coup) -> void {
    Eigen::VectorXd coup_sum = coup * Eigen::VectorXd::Ones(coup.rows());
    if ((coup_sum.array() > m_state.rowwise().sum().array()).any())
      throw std::domain_error("Coupling row sums greater than population");
    else if (coup.coeffs().minCoeff() < 0)
      throw std::domain_error("Encountered value greater less than 0");
    else
      m_coupling = coup;
  }

protected:
  using Base::m_coupling;
  using Base::m_params;
  using Base::m_state;
};

} // namespace EpiGraph
#endif // EPIGRAPH_SIXRD_MODEL_H
