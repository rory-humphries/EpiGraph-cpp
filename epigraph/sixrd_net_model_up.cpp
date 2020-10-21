//
// Created by roryh on 17/10/2020.
//

#include <epigraph/sixrd_net_model_up.hpp>
#include <epigraph/sixrd_odes.hpp>

#include <iostream>

auto print_banner() -> void {
    std::cout <<
              " _____       _  ____                 _    \n" <<
              "| ____|_ __ (_)/ ___|_ __ __ _ _ __ | |__ \n" <<
              "|  _| | '_ \\| | |  _| '__/ _` | '_ \\| '_ \\\n" <<
              "| |___| |_) | | |_| | | | (_| | |_) | | | |\n" <<
              "|_____| .__/|_|\\____|_|  \\__,_| .__/|_| |_|\n" <<
              "      |_|                     |_|          \n";
}

auto print_banner2() -> void {
    std::cout << "                                                                               \n"
                 "                       ,,                                        ,,            \n"
                 "`7MM\"\"\"YMM             db   .g8\"\"\"bgd                          `7MM            \n"
                 "  MM    `7                .dP'     `M                            MM            \n"
                 "  MM   d   `7MMpdMAo.`7MM dM'       ` `7Mb,od8 ,6\"Yb. `7MMpdMAo. MMpMMMb.      \n"
                 "  MMmmMM     MM   `Wb  MM MM            MM' \"'8)   MM   MM   `Wb MM    MM      \n"
                 "  MM   Y  ,  MM    M8  MM MM.    `7MMF' MM     ,pm9MM   MM    M8 MM    MM      \n"
                 "  MM     ,M  MM   ,AP  MM `Mb.     MM   MM    8M   MM   MM   ,AP MM    MM      \n"
                 ".JMMmmmmMMM  MMbmmd' .JMML. `\"bmmmdPY .JMML.  `Moo9^Yo. MMbmmd'.JMML  JMML.    \n"
                 "             MM                                         MM                     \n"
                 "           .JMML.                                     .JMML.                   \n"
                 "";
}

auto SIXRDNetModelUP::set_state(const SIXRDNetModelUP::state_type &x) -> void {
    if ((x.array() < 0).any())
        throw std::domain_error("Encountered negative value");
    else
        m_state = x;
}

auto SIXRDNetModelUP::set_compartment(SIXRDNetModelUP::StateId comp,
                                      const Eigen::VectorXd &vec) -> void {
    if ((vec.array() < 0).any())
        throw std::domain_error("Encountered negative number");
    else
        m_state.col(comp) = vec;
}

auto SIXRDNetModelUP::set_compartment(SIXRDNetModelUP::StateId comp, Eigen::Index i, double val) -> void {
    if (val < 0)
        throw std::domain_error("Encountered negative number");
    else
        m_state(i, comp) = val;
}

auto SIXRDNetModelUP::add_infected(Eigen::Index i, double N) -> void {
    double change = m_state(i, 0) - (0 > m_state(i, 0) - N ? 0 : m_state(i, 0) - N);
    m_state(i, 0) -= change;
    m_state(i, 1) += change;
}

auto SIXRDNetModelUP::set_params(const SIXRDNetModelUP::param_type &params) -> void {
    if ((params.array() < 0).any() || (params.array() > 1).any())
        throw std::domain_error("Encountered value greater than 1 or less than 0");
    else
        m_params = params;
}

auto SIXRDNetModelUP::set_params(SIXRDNetModelUP::ParamId i, double val) -> void {
    if ((val < 0) || (val > 1))
        throw std::domain_error("Encountered value greater than 1 or less than 0");
    else
        m_params[i] = val;
}

auto SIXRDNetModelUP::set_coupling(const SIXRDNetModelUP::coupling_type &coup) -> void {
    Eigen::VectorXd coup_sum = coup * Eigen::VectorXd::Ones(coup.rows());
    if ((coup_sum.array() > m_state.rowwise().sum().array()).any())
        throw std::domain_error("Coupling row sums greater than population");
    else if (coup.coeffs().minCoeff() < 0)
        throw std::domain_error("Encountered value greater less than 0");
    else
        m_coupling = coup;
}

auto derivative(const SIXRDNetModelUP &model) -> SIXRDNetModelUP::state_type {
    return net_SIXRD_ode(model.state(), model.coupling(), model.params());
}

auto print_totals(const SIXRDNetModelUP &model) -> void {
    Eigen::RowVectorXd op_vec = model.state().colwise().sum();
    std::cout << "S : " << op_vec[Sidx];
    std::cout << ", I : " << op_vec[Iidx];
    std::cout << ", X : " << op_vec[Xidx];
    std::cout << ", R : " << op_vec[Ridx];
    std::cout << ", D : " << op_vec[Didx];
}

auto write_state(const SIXRDNetModelUP &model, const std::string &id, std::string dir) -> void {

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

auto write_state_totals(const SIXRDNetModelUP &model, const std::string &path_to_file,
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
