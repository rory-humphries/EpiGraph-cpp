//
// Created by roryh on 02/11/2020.
//

#ifndef EPIGRAPH_CPP_COVID_SIM_HPP
#define EPIGRAPH_CPP_COVID_SIM_HPP

#include <EpiGraph/Core/NetMetaPop.hpp>
#include <EpiGraph/EigenUtil/Accumulate.hpp>
#include <EpiGraph/EigenUtil/LinAlg.hpp>
#include <EpiGraph/EigenUtil/StaticAsserts.hpp>
#include <EpiGraph/Models/SIXRD.hpp>
#include <EpiGraph/Models/SIXRDNetMetaPop.hpp>
#include <EpiGraph/Random/RandomMatrix.hpp>

struct PhaseParams {
  EpiGraph::SIXRDParams sixrd_params;
  double compliance;
  double max_dist;
  int duration;
};

struct CovidSimParams {
  // add constructor to read a toml file
  using Model = EpiGraph::SIXRDNetMetaPop_c<Eigen::MatrixXd, Eigen::MatrixXd,
                                            Eigen::SparseMatrix<double>>;

  int dim;
  int num_phases;

  std::vector<PhaseParams> phase_params_list;

  EpiGraph::RandomAdjMat random_adj_mat;

  Eigen::VectorXd population;
  Eigen::MatrixXd distance_mat;

  std::vector<std::string> phase_order;
  std::vector<std::string> zone;
};

class CovidSim {
 public:
  using Model = EpiGraph::SIXRDNetMetaPop_c<Eigen::MatrixXd, Eigen::MatrixXd,
                                            Eigen::SparseMatrix<double>>;

 private:
  int m_dim;
  int m_num_phases;

  Model m_model;
  EpiGraph::RandomAdjMat random_adj_mat;

  Eigen::VectorXd m_population;
  Eigen::MatrixXd m_distance_mat;

  std::vector<PhaseParams> phase_params_list;

  std::vector<std::string> m_phase_order;
  std::vector<std::string> m_zone;

  std::vector<int> m_cur_phase;
  Eigen::VectorXd m_cur_max_dist;
  Eigen::VectorXd m_cur_compliance;
  Eigen::VectorXi m_cur_duration;

 public:
  CovidSim(CovidSimParams &params) : m_model(params.dim) {
    m_dim = params.dim;
    m_num_phases = params.num_phases;

    // must add check here
    random_adj_mat = params.random_adj_mat;

    assert(params.population.size() == m_dim);
    m_population = params.population;

    Eigen::MatrixXd state = m_model.state();
    state.col(0) = m_population;
    m_model.set_state(state);

    assert(params.distance_mat.rows() == m_dim);
    assert(params.distance_mat.cols() == m_dim);
    m_distance_mat = params.distance_mat;

    // check all the vectors are of the correct length
    phase_params_list = params.phase_params_list;

    assert(params.phase_order.size() == m_num_phases);
    m_phase_order = params.phase_order;

    assert(params.zone.size() == m_dim);
    m_zone = params.zone;

    std::vector<int> m_cur_phase(m_dim);

    m_cur_max_dist =
        Eigen::VectorXd::Ones(m_dim, 1) * phase_params_list[0].max_dist;
    m_cur_compliance =
        Eigen::VectorXd::Ones(m_dim, 1) * phase_params_list[0].compliance;
    m_cur_duration =
        Eigen::VectorXi::Ones(m_dim, 1) * phase_params_list[0].duration;

    m_model.set_params(
        EpiGraph::SixrdParamId::beta,
        Eigen::VectorXd::Ones(m_dim) * phase_params_list[0].sixrd_params.beta);
    m_model.set_params(
        EpiGraph::SixrdParamId::c,
        Eigen::VectorXd::Ones(m_dim) * phase_params_list[0].sixrd_params.c);
    m_model.set_params(
        EpiGraph::SixrdParamId::mu,
        Eigen::VectorXd::Ones(m_dim) * phase_params_list[0].sixrd_params.mu);
    m_model.set_params(
        EpiGraph::SixrdParamId::alpha,
        Eigen::VectorXd::Ones(m_dim) * phase_params_list[0].sixrd_params.alpha);
    m_model.set_params(
        EpiGraph::SixrdParamId::kappa,
        Eigen::VectorXd::Ones(m_dim) * phase_params_list[0].sixrd_params.kappa);
  }

  auto get_model() const -> const Model & { return m_model; }

  auto add_infections(int area_id, int num) -> void {
    m_model.move_state(EpiGraph::SixrdId::S, EpiGraph::SixrdId::I, area_id,
                       num);
  }

  auto update_phase() -> void {
    std::map<std::string, double> zone_I =
        EpiGraph::accumulate_groups(m_model.state().col(1), m_zone);

    // update each nodes current phase
    for (int node = 0; node < m_dim; node++) {
      // decrement duration
      m_cur_duration[node] -= 1;

      if ((zone_I[m_zone[node]] > 1400 && m_cur_phase[node] > 2) ||
          (m_cur_duration[node] == 0 && zone_I[m_zone[node]] > 1400 &&
           m_cur_phase[node] >= 2)) {  // if county gone into lockdown
        m_cur_phase[node] = 2;
      } else if (m_cur_duration[node] == 0 &&
                 (m_cur_phase[node] <
                  m_phase_order.size() - 1)) {  // if node reaches end of phase
        m_cur_phase[node] += 1;
      } else {
        continue;
      }

      // update all parameters if phase has changed
      m_model.set_params(
          EpiGraph::SixrdParamId::beta, node,
          phase_params_list[m_cur_phase[node]].sixrd_params.beta);
      m_model.set_params(EpiGraph::SixrdParamId::c, node,
                         phase_params_list[m_cur_phase[node]].sixrd_params.c);
      m_model.set_params(EpiGraph::SixrdParamId::mu, node,
                         phase_params_list[m_cur_phase[node]].sixrd_params.mu);
      m_model.set_params(
          EpiGraph::SixrdParamId::alpha, node,
          phase_params_list[m_cur_phase[node]].sixrd_params.alpha);
      m_model.set_params(
          EpiGraph::SixrdParamId::kappa, node,
          phase_params_list[m_cur_phase[node]].sixrd_params.kappa);

      m_cur_max_dist[node] = phase_params_list[m_cur_phase[node]].max_dist;
      m_cur_compliance[node] = phase_params_list[m_cur_phase[node]].compliance;
      m_cur_duration[node] = phase_params_list[m_cur_phase[node]].duration;
    }
  }
  auto update_travels() -> void {
    Eigen::MatrixXd max_dist_mat = EpiGraph::general_outer_product(
        m_cur_max_dist, m_cur_max_dist,
        [](double a, double b) { return std::min(a, b); });

    Eigen::VectorXd travel_pop = m_population * 0.6;

    Eigen::MatrixXd tmp = random_adj_mat.gen_sparse_mat(travel_pop);

    Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> bool_mat =
        m_distance_mat.array() < max_dist_mat.array();

    tmp = (bool_mat).select(
        tmp.array(),
        tmp.array() *
            (1 - m_cur_compliance.rowwise().replicate(tmp.cols()).array()));

    m_model.set_coupling(tmp.array().round().matrix().sparseView());
  }
  auto update_state() -> void {
    // Update the state_impl matrix
    Eigen::MatrixXd new_state = m_model.state() + m_model.dXdt();

    Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> eps_mat =
        100000 * new_state.col(EpiGraph::SixrdId::I).array() /
            m_population.array() <
        1.0 / 100000.0;

    Eigen::VectorXd state_diff =
        (eps_mat).select(new_state.col(EpiGraph::SixrdId::I), 0);

    new_state.col(EpiGraph::SixrdId::S) += state_diff;
    new_state.col(EpiGraph::SixrdId::I) -= state_diff;

    m_model.set_state(new_state);
  }

  auto print_compartment_totals() -> void {
    auto comp_vec = m_model.state().colwise().sum();

    std::cout << "\n\nS : " << comp_vec[EpiGraph::SixrdId::S];
    std::cout << "\nI : " << comp_vec[EpiGraph::SixrdId::I];
    std::cout << "\nX : " << comp_vec[EpiGraph::SixrdId::X];
    std::cout << "\nR : " << comp_vec[EpiGraph::SixrdId::R];
    std::cout << "\nD : " << comp_vec[EpiGraph::SixrdId::D];
  }
};

#endif  // EPIGRAPH_CPP_SIXRDNETMETAPOP_HPP
