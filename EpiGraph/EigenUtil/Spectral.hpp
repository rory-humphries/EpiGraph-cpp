//
// Created by roryh on 29/10/2020.
//

#ifndef EPIGRAPH_CPP_SPECTRAL_HPP
#define EPIGRAPH_CPP_SPECTRAL_HPP

#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/DenseGenMatProd.h>
#include <Spectra/MatOp/SparseGenMatProd.h>

#include <iostream>

#include <EpiGraph/EigenUtil/StaticAsserts.hpp>

namespace EpiGraph {

template <IsMatrix Mat>
auto SpectralRadius(const Mat &mat, int ev = 6) -> double {
  using Scalar = typename Mat::Scalar;
  Spectra::DenseGenMatProd<Scalar> op(mat);
  // Construct eigen solver object, requesting the largest
  // (in magnitude, or norm) three eigenvalues
  Spectra::GenEigsSolver<Scalar, Spectra::LARGEST_MAGN,
                         Spectra::DenseGenMatProd<Scalar>>
      eigs(&op, 1, ev);
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

template <IsSparseMatrix SpMat>
auto SpectralRadius(const SpMat &mat, int ev = 6) -> double {
  using Scalar = typename SpMat::Scalar;
  Spectra::SparseGenMatProd<Scalar> op(mat);
  // Construct eigen solver object, requesting the largest
  // (in magnitude, or norm) three eigenvalues
  Spectra::GenEigsSolver<Scalar, Spectra::LARGEST_MAGN,
                         Spectra::SparseGenMatProd<Scalar>>
      eigs(&op, 1, ev);
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
} // namespace EpiGraph

#endif // EPIGRAPH_CPP_SPECTRAL_HPP
