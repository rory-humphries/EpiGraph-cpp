//
// Created by roryh on 28/10/2020.
//

#ifndef EPIGRAPH_CPP_SIR_HPP
#define EPIGRAPH_CPP_SIR_HPP

#include <EpiGraph/EigenUtil/StaticAsserts.hpp>

namespace EpiGraph {

/**
 * @brief ODE which describes the SIR model. Computes dx/dt in place.
 *
 * @tparam Vec1 A row or column vector which inherits form
 * Eigen::MatrixBase<Vec1>.
 * @tparam Vec2 A row or column vector which inherits form
 * Eigen::MatrixBase<Vec2>.
 * @tparam Vec3 A row or column vector which inherits form
 * Eigen::MatrixBase<Vec3>.
 * @param x A vector of length 3. x[0] = S. x[1] = I. x[2] = R.
 * @param dxdt A vector of length 3. dxdt[0] = dS/dt. dxdt[1] = dI/dt. dxdt[2] =
 * dR/dt.
 * @param p A vector of length 4. p[0] = probability of infection on contact.
 * p[1] = recovery rate. p[2] = birth rate. p[3] = death rate.
 */
template <IsVector Vec1, IsVector Vec2, IsVector Vec3>
auto sir_ode(const Vec1 &x, Vec2 &dxdt, const Vec3 &p) -> void {
  using StateType = Vec1::Scalar;
  using ParamType = Vec3::Scalar;

  const StateType &S = x[0];
  const StateType &I = x[1];
  const StateType &R = x[2];

  const StateType N = x.sum();

  const ParamType beta = p[0];
  const ParamType mu = p[1];
  const ParamType alpha = p[2];
  const ParamType gamma = p[3];

  dxdt[0] = alpha * N - gamma * S - beta * S * I / N;
  dxdt[1] = beta * S * I / N - mu * I - gamma * I;
  dxdt[2] = mu * I - gamma * R;
}

/**
 * @brief ODE which describes the SIR model.
 *
 * @tparam Vec1 A row or column vector which inherits form
 * Eigen::MatrixBase<Vec1>.
 * @tparam Vec2 A row or column vector which inherits form
 * Eigen::MatrixBase<Vec2>.
 * @param x A vector of length 3. x[0] = S. x[1] = I. x[2] = R
 * @param p A vector of length 4. p[0] = probability of infection on contact.
 * p[1] = recovery rate. p[2] = birth rate. p[3] = death rate.
 * @return Vec1 A vector of length 3 representing the derivative dx/dt such
 * that, dxdt[0] = dS/dt. dxdt[1] = dI/dt. dxdt[2] = dR/dt.
 */
template <IsVector Vec1, IsVector Vec2>
auto sir_ode(const Vec1 &x, const Vec2 &p) -> Vec1 {

  Vec1 dxdt(3);

  sir_ode(x, dxdt, p);

  return dxdt;
}
} // namespace EpiGraph
#endif // EPIGRAPH_CPP_SIR_HPP
