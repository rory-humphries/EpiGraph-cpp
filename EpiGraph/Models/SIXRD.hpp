//
// Created by roryh on 19/06/2020.
//


#ifndef EPIGRAPH_SIRX_NETWORK_H
#define EPIGRAPH_SIRX_NETWORK_H

#include <Eigen/Dense>

#include <map>
#include <iostream>

namespace EpiGraph {

	enum SixrdId : Eigen::Index {
		S, I, X, R, D
	};

	enum SixrdParamId : Eigen::Index {
		beta, c, mu, alpha, kappa
	};

	template<typename DerivedA, typename DerivedB, typename DerivedC>
	auto sixrd_net_ode_uni_params(const Eigen::MatrixBase<DerivedA> &x,
								  const Eigen::EigenBase<DerivedB> &adj_,
								  const Eigen::MatrixBase<DerivedC> &sixrd_params) -> DerivedA {
		/*
		 * xmat is expected to be an n x 5 matrix (columns checked at compile time) where the following column indices
		 * represent 0=S, 1=I, 2=X, 3=R, 4=D. The rows represent the nodes of the system.
		 *
		 * adj is the couplings in the system representing the movements of individual between nodes. If xmat is an n x 5
		 * matrix then adj must be a n x n matrix.
		 */
		//col_vector_assert(sixrd_params);
		//SIXRDState_assert(x);

		double beta = sixrd_params[0];//param.beta;
		double c = sixrd_params[1];//param.c;
		double mu = sixrd_params[2];//param.mu;
		double alpha = sixrd_params[3];//param.alpha;
		double kappa = sixrd_params[4];//param.kappa;

		const DerivedB &adj = adj_.derived();

		Eigen::VectorXd n = x.rowwise().sum();

		DerivedB S = (x.col(0).cwiseProduct(n.cwiseInverse())).asDiagonal() * adj;
		DerivedB I = (x.col(1).cwiseProduct(n.cwiseInverse())).asDiagonal() * adj;

		Eigen::VectorXd new_s = x.col(0) - S * Eigen::VectorXd::Ones(S.rows());

		Eigen::VectorXd new_i = x.col(1) + (Eigen::RowVectorXd::Ones(I.rows()) * I).transpose() -
								(I * Eigen::VectorXd::Ones(I.cols()));


		n += (Eigen::RowVectorXd::Ones(adj.rows()) * adj).transpose() -
			 (adj * Eigen::VectorXd::Ones(adj.cols()));

		Eigen::VectorXd idiff_ndiff_inv = new_i.cwiseProduct(n.cwiseInverse());

		DerivedA dxdt(x.rows(), x.cols());

		dxdt.col(4) = alpha * x.col(1) + alpha * x.col(2);// D
		dxdt.col(3) = mu * x.col(1) + mu * x.col(2);// R
		dxdt.col(2) = kappa * x.col(1) - mu * x.col(2) - alpha * x.col(2); // X

		dxdt.col(1) =
				beta * c * new_s.cwiseProduct(idiff_ndiff_inv) + beta * c * S * idiff_ndiff_inv -
				mu * x.col(1) - alpha * x.col(1) - kappa * x.col(1);// I

		dxdt.col(0) = -beta * c * new_s.cwiseProduct(idiff_ndiff_inv) -
					  beta * c * S * idiff_ndiff_inv; // S

		return dxdt;
	}

	template<typename DerivedA, typename DerivedB, typename DerivedC>
	auto
	sixrd_net_ode(const Eigen::MatrixBase<DerivedA> &xmat, const Eigen::EigenBase<DerivedB> &adj_,
				  const Eigen::MatrixBase<DerivedC> &sixrd_params) -> DerivedA {
		/*
		 * xmat is expected to be an n x 5 matrix (columns checked at compile time) where the following column indices
		 * represent 0=S, 1=I, 2=X, 3=R, 4=D. The rows represent the nodes of the system.
		 *
		 * adj is the couplings in the system representing the movements of individual between nodes. If xmat is an n x 5
		 * matrix then adj must be a n x n matrix.
		 */
		//SIXRDState_assert(xmat);

		const DerivedB &adj = adj_.derived();

		Eigen::VectorXd n = xmat.rowwise().sum();

		DerivedB S = (xmat.col(0).cwiseProduct(n.cwiseInverse())).asDiagonal() * adj;
		DerivedB I = (xmat.col(1).cwiseProduct(n.cwiseInverse())).asDiagonal() * adj;

		Eigen::VectorXd new_s = xmat.col(0) - S * Eigen::VectorXd::Ones(S.rows());

		Eigen::VectorXd new_i = xmat.col(1) + (Eigen::RowVectorXd::Ones(I.rows()) * I).transpose() -
								(I * Eigen::VectorXd::Ones(I.cols()));


		n += (Eigen::RowVectorXd::Ones(adj.rows()) * adj).transpose() -
			 (adj * Eigen::VectorXd::Ones(adj.cols()));

		Eigen::VectorXd idiff_ndiff_inv = new_i.cwiseProduct(n.cwiseInverse());

		DerivedA dxdt(xmat.rows(), xmat.cols());

		dxdt.col(4) = sixrd_params.col(alpha).cwiseProduct(xmat.col(1) + xmat.col(2));// D
		dxdt.col(3) = sixrd_params.col(mu).cwiseProduct(xmat.col(1) + xmat.col(2));// R
		dxdt.col(2) = sixrd_params.col(kappa).cwiseProduct(xmat.col(1)) -
					  sixrd_params.col(mu).cwiseProduct(xmat.col(2)) -
					  sixrd_params.col(alpha).cwiseProduct(xmat.col(2)); // X

		dxdt.col(1) = sixrd_params.col(beta).cwiseProduct(
				sixrd_params.col(c).cwiseProduct(new_s.cwiseProduct(idiff_ndiff_inv))) +
					  (S * (sixrd_params.col(beta).cwiseProduct(sixrd_params.col(c)).asDiagonal()) *
					   idiff_ndiff_inv) - sixrd_params.col(mu).cwiseProduct(xmat.col(1)) -
					  sixrd_params.col(alpha).cwiseProduct(xmat.col(1)) -
					  sixrd_params.col(kappa).cwiseProduct(xmat.col(1));// I

		dxdt.col(0) = -sixrd_params.col(beta).cwiseProduct(
				sixrd_params.col(c).cwiseProduct(new_s.cwiseProduct(idiff_ndiff_inv))) -
					  (S * (sixrd_params.col(beta).cwiseProduct(sixrd_params.col(c)).asDiagonal()) *
					   idiff_ndiff_inv); // S

		return dxdt;
	}

	template<typename DerivedA, typename DerivedB, typename DerivedC>
	auto sixrd_next_gen_matrix(const Eigen::MatrixBase<DerivedA> &xmat,
							   const Eigen::EigenBase<DerivedB> &adj,
							   const Eigen::MatrixBase<DerivedC> &param) -> DerivedB {
		/*
		 * Find the reproduction number for the network SIXRD model
		 */

		vector_assert(param);
		int dim = xmat.rows();

		double alpha = param[SixrdParamId::alpha];
		double beta = param[SixrdParamId::beta];
		double c = param[SixrdParamId::c];
		double mu = param[SixrdParamId::mu];
		double kappa = param[SixrdParamId::kappa];

		const DerivedB &adj_d = adj.derived();

		Eigen::VectorXd n = xmat.rowwise().sum().array();
		for (int i = 0; i < xmat.rows(); i++) {
			if (xmat.row(i).sum() <= 0)
				std::cout << "here" << std::endl;
		}
		Eigen::VectorXd n_inv = n.cwiseInverse();

		DerivedB S = (xmat.col(SixrdId::S).cwiseProduct(n_inv)).asDiagonal() * adj_d;

		Eigen::VectorXd one_vec = Eigen::VectorXd::Ones(dim);
		Eigen::RowVectorXd one_rowvec = Eigen::RowVectorXd::Ones(dim);

		Eigen::VectorXd n_diff = n + (one_rowvec * adj_d).transpose() - (adj_d * one_vec);
		Eigen::VectorXd n_diff_inv = n_diff.cwiseInverse();

		Eigen::VectorXd mat1 = xmat.col(SixrdId::S) - S * one_vec;

		DerivedB mat2(
				(one_vec.array() - ((adj_d * one_vec).array() / n.array())).matrix().asDiagonal());
		mat2 += DerivedB(adj_d.transpose() * (n.cwiseInverse().asDiagonal()));
		mat2 = n_diff_inv.asDiagonal() * mat2;

		DerivedB mat3 = S * mat2;

		DerivedB T = ((beta * c) / (mu + alpha + kappa)) * (mat1.asDiagonal() * mat2 + mat3);

		return T;
	}

	template<typename DerivedA, typename DerivedB, typename DerivedC>
	auto sixrd_net_r0(const Eigen::MatrixBase<DerivedA> &xmat,
					  const Eigen::SparseMatrixBase<DerivedB> &adj,
					  const Eigen::MatrixBase<DerivedC> &param) -> double {

		/*
		 * Find the reproduction number for the network SIXRD model
		 */


		DerivedB T = sixrd_next_gen_matrix(xmat, adj, param);
		return SpectralRadius(T);
	}

	/*
	template<typename DerivedA>
	auto SIXRD_R0(EigenUtil::MatrixBase<DerivedA> &xmat, const SIXRDParam param) -> double {

		 * Find the reproduction number for the SIXRD model


		static_assert(1 == DerivedA::ColsAtCompileTime, "Input matrix must have 1 cols (must be a vector)");

		using Mat = EigenUtil::Matrix2d;

		double alpha = param.alpha;
		double beta = param.beta;
		double c = param.c;
		double mu = param.mu;
		double kappa = param.kappa;

		Mat T = Mat::Zero(2, 2);
		T(0, 0) = beta * c * xmat(Sidx) / xmat.sum();

		Mat E = Mat::Zero(2, 2);
		E(0, 0) = -mu - alpha - kappa;
		E(1, 0) = kappa;
		E(1, 1) = -mu - alpha;

		Mat E_inv = E.inverse();

		Mat prod = -T * E_inv;

		EigenUtil::EigenSolver<Mat> eigensolver;
		eigensolver.compute(prod);
		EigenUtil::VectorXd eigen_values = eigensolver.eigenvalues().cwiseAbs();

		return eigen_values.maxCoeff();
	};

	*/
}

#endif //EPIGRAPH_SIRX_NETWORK_H
