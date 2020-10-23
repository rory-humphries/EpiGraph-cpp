//
// Created by roryh on 21/10/2020.
//
#include <epigraph/Models/SixrdMetaPop.hpp>
#include <epigraph/SixrdOde.hpp>

#include <iostream>

auto SixrdMetaPop::set_state(const SixrdMetaPop::state_type &x) -> void {
    if ((x.array() < 0).any())
        throw std::domain_error("Encountered negative value");
    else
        m_state = x;
}

auto SixrdMetaPop::set_compartment(SixrdMetaPop::StateId comp,
                                   const Eigen::VectorXd &vec) -> void {
    if ((vec.array() < 0).any())
        throw std::domain_error("Encountered negative number");
    else
        m_state.col(comp) = vec;
}

auto SixrdMetaPop::set_compartment(SixrdMetaPop::StateId comp, Eigen::Index i, double val) -> void {
    if (val < 0)
        throw std::domain_error("Encountered negative number");
    else
        m_state(i, comp) = val;
}

auto SixrdMetaPop::add_infected(Eigen::Index i, double N) -> void {
    double change = m_state(i, 0) - (0 > m_state(i, 0) - N ? 0 : m_state(i, 0) - N);
    m_state(i, 0) -= change;
    m_state(i, 1) += change;
}

auto SixrdMetaPop::set_params(const SixrdMetaPop::param_type &params) -> void {
    if ((params.array() < 0).any() || (params.array() > 1).any())
        throw std::domain_error("Encountered value greater than 1 or less than 0");
    else
        m_params = params;
}

auto SixrdMetaPop::set_params(SixrdMetaPop::ParamId param_id, const Eigen::VectorXd &vec) -> void {
    if ((vec.array() < 0).any() || (vec.array() > 1).any())
        throw std::domain_error("Encountered value greater than 1 or less than 0");
    else
        m_params.col(param_id) = vec;
}


auto SixrdMetaPop::set_params(SixrdMetaPop::ParamId param_id, Eigen::Index i, double val) -> void {
    if ((val < 0) || (val > 1))
        throw std::domain_error("Encountered value greater than 1 or less than 0");
    else
        m_params(param_id, i) = val;
}

auto SixrdMetaPop::set_coupling(const SixrdMetaPop::coupling_type &coup) -> void {
    Eigen::VectorXd coup_sum = coup * Eigen::VectorXd::Ones(coup.rows());
    if ((coup_sum.array() > m_state.rowwise().sum().array()).any())
        throw std::domain_error("Coupling row sums greater than population");
    else if (coup.coeffs().minCoeff() < 0)
        throw std::domain_error("Encountered value greater less than 0");
    else
        m_coupling = coup;
}

auto derivative(const SixrdMetaPop &model) -> SixrdMetaPop::state_type {
    return net_SIXRD_ode_inhom(model.state(), model.coupling(), model.params());
}

auto print_totals(const SixrdMetaPop &model) -> void {
    Eigen::RowVectorXd op_vec = model.state().colwise().sum();
    std::cout << "S : " << op_vec[Sidx];
    std::cout << ", I : " << op_vec[Iidx];
    std::cout << ", X : " << op_vec[Xidx];
    std::cout << ", R : " << op_vec[Ridx];
    std::cout << ", D : " << op_vec[Didx];
}

auto write_state(const SixrdMetaPop &model, const std::string &id, std::string dir) -> void {

    auto &x = model.state();
    std::ofstream myfile;

    myfile.open(dir + id + ".csv");

    myfile << "S,";
    myfile << "I,";
    myfile << "X,";
    myfile << "R,";
    myfile << "D,";
    myfile << "N\n";

    for (int i = 0; i < x.rows(); i++) {
        myfile << x(i, Sidx) << ",";
        myfile << x(i, Iidx) << ",";
        myfile << x(i, Xidx) << ",";
        myfile << x(i, Ridx) << ",";
        myfile << x(i, Didx) << ",";
        myfile << x.row(i).sum() << "\n";
    }
    myfile.close();
}

auto write_state_totals(const SixrdMetaPop &model, const std::string &path_to_file,
                        bool append_file = false) -> void {

    auto &x = model.state();
    double sumi = 0;
    double sums = 0;
    double sumr = 0;
    double sumx = 0;
    double sumd = 0;

    for (int i = 0; i < x.rows(); i++) {
        sums += x(i, Sidx);
        sumi += x(i, Iidx);
        sumx += x(i, Xidx);
        sumr += x(i, Ridx);
        sumd += x(i, Didx);
    }

    std::ofstream myfile;

    if (append_file)
        myfile.open(path_to_file, std::ios_base::app);
    else {
        myfile.open(path_to_file);
        myfile << "S,";
        myfile << "I,";
        myfile << "X,";
        myfile << "R,";
        myfile << "D\n";
    }

    myfile << sums << ",";
    myfile << sumi << ",";
    myfile << sumx << ",";
    myfile << sumr << ",";
    myfile << sumd << "\n";

    myfile.close();
}