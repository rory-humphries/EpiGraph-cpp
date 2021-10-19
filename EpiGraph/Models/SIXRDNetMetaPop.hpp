//
// Created by roryh on 02/11/2020.
//

#ifndef EPIGRAPH_CPP_SIXRDNETMETAPOP_C_HPP
#define EPIGRAPH_CPP_SIXRDNETMETAPOP_C_HPP

#include <EpiGraph/Core/NetMetaPop.hpp>
#include <EpiGraph/EigenUtil/StaticAsserts.hpp>
#include <EpiGraph/Models/SIXRD.hpp>

namespace EpiGraph {

class SIXRDParams {
public:
  double beta;
  double c;
  double mu;
  double alpha;
  double kappa;
};

template <IsMatrix StateMat, IsMatrix ParamType, IsSparseOrDenseMatrix CoupType>
class SIXRDNetMetaPop_c {
public:
  explicit SIXRDNetMetaPop_c(Eigen::Index dim) : m_dim(dim) {
    m_state = StateMat::Zero(dim, 5);
    m_params = ParamType::Zero(dim, 5);
    m_coupling = CoupType(dim, dim);
  }

  auto dim() const -> int { return m_dim; }

  auto state() const -> const StateMat & { return m_state; }

  auto params() const -> const ParamType & { return m_params; }

  auto coupling() const -> const CoupType & { return m_coupling; }

  auto set_state(const StateMat &x) -> void {
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
    using Scalar = typename StateMat::Scalar;
    Scalar change =
        m_state(i, comp1) -
        (comp1 > m_state(i, comp1) - val ? comp1 : m_state(i, comp1) - val);
    m_state(i, comp1) -= change;
    m_state(i, comp2) += change;
  }

  auto set_params(const ParamType &params) -> void {
    if ((params.array() < 0).any() || (params.array() > 1).any())
      throw std::domain_error(
          "Encountered value greater than 1 or less than 0");
    else
      m_params = params;
  }

  auto set_params(Eigen::Index i, double val) -> void {
    if ((val < 0) || (val > 1))
      throw std::domain_error(
          "Encountered value greater than 1 or less than 0");
    else
      m_params.col(i) =
          Eigen::Matrix<typename ParamType::Scalar,
                        ParamType::RowsAtCompileTime, 1>::Ones(m_params.rows());
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

  auto set_coupling(const CoupType &coup) -> void {
    Eigen::VectorXd coup_sum = coup * Eigen::VectorXd::Ones(coup.rows());
    if ((coup_sum.array() > m_state.rowwise().sum().array()).any())
      throw std::domain_error("Coupling row sums greater than population");
    else if (coup.coeffs().minCoeff() < 0)
      throw std::domain_error("Encountered value greater less than 0");
    else
      m_coupling = coup;
  }

  auto dXdt() const -> StateMat {
    /*
     * xmat is expected to be an n x 5 matrix (columns checked at compile time)
     * where the following column indices represent 0=S, 1=I, 2=X, 3=R, 4=D. The
     * rows represent the nodes of the system.
     *
     * adj is the couplings in the system representing the movements of
     * individual between nodes. If xmat is an n x 5 matrix then adj must be a n
     * x n matrix.
     */
    // SIXRDState_assert(xmat);

    Eigen::VectorXd n = m_state.rowwise().sum();

    CoupType S = (m_state.col(0).cwiseProduct(n.cwiseInverse())).asDiagonal() *
                 m_coupling;
    CoupType I = (m_state.col(1).cwiseProduct(n.cwiseInverse())).asDiagonal() *
                 m_coupling;

    Eigen::VectorXd new_s =
        m_state.col(0) - S * Eigen::VectorXd::Ones(S.rows());

    Eigen::VectorXd new_i =
        m_state.col(1) + (Eigen::RowVectorXd::Ones(I.rows()) * I).transpose() -
        (I * Eigen::VectorXd::Ones(I.cols()));

    n +=
        (Eigen::RowVectorXd::Ones(m_coupling.rows()) * m_coupling).transpose() -
        (m_coupling * Eigen::VectorXd::Ones(m_coupling.cols()));

    Eigen::VectorXd idiff_ndiff_inv = new_i.cwiseProduct(n.cwiseInverse());

    StateMat dxdt(m_state.rows(), m_state.cols());

    dxdt.col(4) =
        m_params.col(alpha).cwiseProduct(m_state.col(1) + m_state.col(2)); // D
    dxdt.col(3) =
        m_params.col(mu).cwiseProduct(m_state.col(1) + m_state.col(2)); // R
    dxdt.col(2) = m_params.col(kappa).cwiseProduct(m_state.col(1)) -
                  m_params.col(mu).cwiseProduct(m_state.col(2)) -
                  m_params.col(alpha).cwiseProduct(m_state.col(2)); // X

    dxdt.col(1) =
        m_params.col(beta).cwiseProduct(
            m_params.col(c).cwiseProduct(new_s.cwiseProduct(idiff_ndiff_inv))) +
        (S * (m_params.col(beta).cwiseProduct(m_params.col(c)).asDiagonal()) *
         idiff_ndiff_inv) -
        m_params.col(mu).cwiseProduct(m_state.col(1)) -
        m_params.col(alpha).cwiseProduct(m_state.col(1)) -
        m_params.col(kappa).cwiseProduct(m_state.col(1)); // I

    dxdt.col(0) =
        -m_params.col(beta).cwiseProduct(
            m_params.col(c).cwiseProduct(new_s.cwiseProduct(idiff_ndiff_inv))) -
        (S * (m_params.col(beta).cwiseProduct(m_params.col(c)).asDiagonal()) *
         idiff_ndiff_inv); // S

    return dxdt;
  }

private:
  CoupType m_coupling;
  ParamType m_params;
  StateMat m_state;
  int m_dim;
};

} // namespace EpiGraph

#endif // EPIGRAPH_CPP_SIXRDNETMETAPOP_HPP
