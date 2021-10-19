
#include <EpiGraph/Graph/BFS.hpp>
#include <cstdlib>

void BFS_c(const int *row, const int *outer, int n, int *res, int v) {

  //(*int)std::calloc(n, sizeof(int));
  // Eigen::VectorXi T = Eigen::VectorXi::Zero(n);
  int *L = (int *)std::calloc(n, sizeof(int));

  int start = 0;
  int end = 1;
  int z = 1;

  int *L_curr = L;
  *L_curr = v;

  while (start != end) {

    for (int j = start; j < end; j++) {

      int node = L[j];

      *(res + node) = 1;

      int nnz = *(outer + node + 1) - *(outer + node);
      const int *i = row + *(outer + node);
      // for (int k = 0; k < nnz; k++) {
      const int *i_end = row + *(outer + node) + nnz;
      for (i; i != i_end; i++) {
        // const int *i = row + *(outer + node) + k;
        if (*(res + *i) == 0) {
          *(res + *i) = 1;

          L_curr++;
          *L_curr = *i;

          z++;
        }
      }
    }
    start = end;
    end = z;
  }
  free(L);
}
