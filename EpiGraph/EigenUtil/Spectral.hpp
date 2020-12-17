//
// Created by roryh on 29/10/2020.
//

#ifndef EPIGRAPH_CPP_SPECTRAL_HPP
#define EPIGRAPH_CPP_SPECTRAL_HPP


#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/MatOp/DenseGenMatProd.h>
#include <iostream>

template<typename Derived>
auto SpectralRadius(const Eigen::MatrixBase<Derived> &mat, int ev = 6) -> double{
	using Scalar = typename Derived::Scalar;
	Spectra::DenseGenMatProd<Scalar> op(mat);
	// Construct eigen solver object, requesting the largest
	// (in magnitude, or norm) three eigenvalues
	Spectra::GenEigsSolver<Scalar, Spectra::LARGEST_MAGN, Spectra::DenseGenMatProd<Scalar> > eigs(
			&op, 1, ev);
	// Initialize and compute
	eigs.init();
	int nconv = eigs.compute();
	// Retrieve results
	Eigen::VectorXcd evalues;
	if (eigs.info() == Spectra::SUCCESSFUL) {
		evalues = eigs.eigenvalues();
		return evalues.cwiseAbs().maxCoeff();
	} else {
		throw std::runtime_error("Could not find eigenvalue");
	}
}

template<typename Derived>
auto SpectralRadius(const Eigen::SparseMatrixBase<Derived> &mat, int ev = 6) -> double {
	using Scalar = typename Derived::Scalar;
	Spectra::SparseGenMatProd<Scalar> op(mat);
	// Construct eigen solver object, requesting the largest
	// (in magnitude, or norm) three eigenvalues
	Spectra::GenEigsSolver<Scalar, Spectra::LARGEST_MAGN, Spectra::SparseGenMatProd<Scalar> > eigs(
			&op, 1, ev);
	// Initialize and compute
	eigs.init();
	int nconv = eigs.compute();
	// Retrieve results
	Eigen::VectorXcd evalues;
	if (eigs.info() == Spectra::SUCCESSFUL) {
		evalues = eigs.eigenvalues();
		return evalues.cwiseAbs().maxCoeff();
	} else {
		throw std::runtime_error("Could not find eigenvalue");
	}
}

#endif //EPIGRAPH_CPP_SPECTRAL_HPP
