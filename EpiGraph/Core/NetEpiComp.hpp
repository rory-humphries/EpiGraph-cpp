//
// Created by roryh on 20/10/2020.
//

#ifndef EPIGRAPH_SIXRD_MODEL_H
#define EPIGRAPH_SIXRD_MODEL_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <EpiGraph/Core/Util.hpp>
#include <EpiGraph/Models/SIXRD.hpp>

namespace EpiGraph {
    template<int ParamType>
    class NetEpiComp;

    template<int ParamType>
    struct traits<NetEpiComp<ParamType>> {

        using state_scalar_type = double;
        using param_scalar_type = double;
        using coupling_scalar_type = double;
        constexpr static int ParaType = ParamType == 1 ? 1 : Eigen::Dynamic;
    };

    template<int ParamType>
    class NetEpiComp {
    public:
        using state_scalar_type = typename traits<NetEpiComp<ParamType>>::state_scalar_type;
        using state_type = Eigen::Matrix<state_scalar_type, Eigen::Dynamic, Eigen::Dynamic>;

        using param_scalar_type = typename traits<NetEpiComp<ParamType>>::param_scalar_type;
        using param_type = Eigen::Matrix<param_scalar_type, traits<NetEpiComp<ParamType>>::ParaType, Eigen::Dynamic>;

        using coupling_scalar_type = typename traits<NetEpiComp<ParamType>>::coupling_scalar_type;
        using coupling_type = Eigen::SparseMatrix<coupling_scalar_type>;

        NetEpiComp(Eigen::Index dim, Eigen::Index compartments) {
            m_state = state_type::Zero(dim, compartments);
            m_params = param_type::Zero((ParamType == 1 ? 1 : dim), compartments);
            m_coupling = coupling_type(dim, dim);
        }

        auto state() const -> const state_type & { return m_state; }

        auto params() const -> const param_type & { return m_params; }

        auto coupling() const -> const coupling_type & { return m_coupling; }

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

        auto add_infected(Eigen::Index i, double N) -> void {
            double change = m_state(i, 0) - (0 > m_state(i, 0) - N ? 0 : m_state(i, 0) - N);
            m_state(i, 0) -= change;
            m_state(i, 1) += change;
        }

        auto set_params(const param_type &params) -> void {
            if ((params.array() < 0).any() || (params.array() > 1).any())
                throw std::domain_error("Encountered value greater than 1 or less than 0");
            else
                m_params = params;
        }

        auto set_params(Eigen::Index i, double val) -> void {
            if ((val < 0) || (val > 1))
                throw std::domain_error("Encountered value greater than 1 or less than 0");
            else
                m_params[i] = val;
        }

        auto set_params(Eigen::Index i, Eigen::Index j, double val) -> void {
            if ((val < 0) || (val > 1))
                throw std::domain_error("Encountered value greater than 1 or less than 0");
            else
                m_params(j, i) = val;
        }

        auto set_params(Eigen::Index i, const Eigen::Matrix<param_scalar_type, Eigen::Dynamic, 1> &vec) -> void {
            if ((vec.array() < 0).any() || (vec.array() > 1).any())
                throw std::domain_error("Encountered value greater than 1 or less than 0");
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

    private:
        state_type m_state;
        param_type m_params;
        coupling_type m_coupling;


    };

    template<int ParamType>
    auto sixrd_meta_pop_ode(const NetEpiComp<ParamType> &model) -> typename NetEpiComp<ParamType>::state_type;

    template<int ParamType>
    auto print_totals(const NetEpiComp<ParamType> &model) -> void {
        Eigen::RowVectorXd op_vec = model.state().colwise().sum();
        std::cout << "S : " << op_vec[SixrdId::S];
        std::cout << ", I : " << op_vec[SixrdId::I];
        std::cout << ", X : " << op_vec[SixrdId::X];
        std::cout << ", R : " << op_vec[SixrdId::R];
        std::cout << ", D : " << op_vec[SixrdId::D];
    }

    template<int ParamType>
    auto write_state(const NetEpiComp<ParamType> &model, const std::string &id, std::string dir) -> void {

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
            myfile << x(i, SixrdId::S) << ",";
            myfile << x(i, SixrdId::I) << ",";
            myfile << x(i, SixrdId::X) << ",";
            myfile << x(i, SixrdId::R) << ",";
            myfile << x(i, SixrdId::D) << ",";
            myfile << x.row(i).sum() << "\n";
        }
        myfile.close();
    }

    template<int ParamType>
    auto write_state_totals(const NetEpiComp<ParamType> &model, const std::string &path_to_file,
                            bool append_file = false) -> void {

        auto &x = model.state();
        double sumi = 0;
        double sums = 0;
        double sumr = 0;
        double sumx = 0;
        double sumd = 0;

        for (int i = 0; i < x.rows(); i++) {
            sums += x(i, SixrdId::S);
            sumi += x(i, SixrdId::I);
            sumx += x(i, SixrdId::X);
            sumr += x(i, SixrdId::R);
            sumd += x(i, SixrdId::D);
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
}
#endif //EPIGRAPH_SIXRD_MODEL_H
