
#ifndef EPIGRAPH_BFS_HPP
#define EPIGRAPH_BFS_HPP

#include <algorithm>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <EpiGraph/EigenUtil/StaticAsserts.hpp>
#include <iostream>
#include <unordered_set>
#include <set>

namespace EpiGraph
{

  template <IsSparseMatrix Adj>
  auto BFS(const Adj &adj, Eigen::Index v)
      -> Eigen::VectorXi
  {
    Eigen::VectorXi T = Eigen::VectorXi::Zero(adj.cols());
    Eigen::VectorXi L = Eigen::VectorXi::Zero(adj.cols());

    int start = 0;
    int end = 1;
    int z = 1;

    L[start] = v;

    while (start != end)
    {
      for (int j = start; j < end; j++)
      {
        int node = L[j];
        T[node] = 1;
        for (typename Adj::InnerIterator it(adj, node); it; ++it)
        {
          int i = it.row();

          if (T[i] == 0)
          {
            T[i] = 1;
            L[z] = i;
            z++;
          }
        }
      }
      start = end;
      end = z;
    }
    return T;
  }

  template <IsSparseMatrix Adj>
  auto BFS_parallel(const Adj &adj, Eigen::Index v)
      -> Eigen::VectorXi
  {
    Eigen::VectorXi T = Eigen::VectorXi::Zero(adj.cols());

    std::vector<int> L(adj.cols(), 0);

    int d = 1;
    int start = 0;
    int end = 1;

    L[start] = v;
    
    std::vector<std::vector<int>> pvec(8);

    while (start != end)
    {
#pragma omp parallel
      {
        int tid = omp_get_thread_num();
        pvec[tid].clear();
      }
#pragma omp parallel for
      for (int j = start; j < end; j++)
      {

        int node = L[j];
        T[node] = d;

        int tid = omp_get_thread_num();

        for (typename Adj::InnerIterator it(adj, node); it; ++it)
        {
          int i = it.row();
          int c;

          if (T[i] == 0)
          {
#pragma omp atomic capture
            {
              c = T[i];
              T[i] = 1;
            }
            if (c == 0)
            {
              pvec[tid].push_back(i);
            }
          }
        }
      }
     
#pragma omp parallel
      {
        int tid = omp_get_thread_num();
        int offset = 0;
        for (int i = 0; i < tid; i++)
          offset += pvec[i].size();
        std::copy(pvec[tid].begin(), pvec[tid].end(), L.begin() + end + offset);
      }

      start = end;
      for (int i = 0; i < 8; i++)
          end += pvec[i].size();

      d +=1;
    }
    return T;
  }

} // namespace EpiGraph

#endif
