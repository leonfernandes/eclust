#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

NumericVector vi, vi0, vj, vj0;

double euclidean_norm(int d)
{
  double norm_z = 0.0;
  for (int i = 0; i < d; ++i)
  {
    norm_z = norm_z + (vi[i] - vj[i]) * (vi[i] - vj[i]);
  }
  return std::pow(norm_z, 0.5);
}

double euclidean_norm0(int d)
{
  double norm_z = 0.0;
  for (int i = 0; i < d; ++i)
  {
    norm_z = norm_z + (vi0[i] - vj0[i]) * (vi0[i] - vj0[i]);
  }
  return std::pow(norm_z, 0.5);
}

double energy_distance_xyh(
    NumericVector x,
    NumericVector y,
    int h,
    int m,
    int n,
    int p,
    double a = 1.0)
{
  if (std::min(m, n) <= h)
  {
    return -1.0;
  }

  int m_new = m - h;
  int n_new = n - h;
  int p_new = (h + 1) * p;

  double sumxy, sumxx, sumyy;
  sumxy = sumxx = sumyy = 0.0;

  for (int i = 0; i < m_new; ++i)
  {
    vi = x[Range(i * p, i * p + p_new - 1)];
    for (int j = 0; j < n_new; ++j)
    {
      vj = y[Range(j * p, j * p + p_new - 1)];
      sumxy += std::pow(euclidean_norm(p_new), a);
    }
  }

  for (int i = 0; i < m_new; ++i)
  {
    vi = x[Range(i * p, i * p + p_new - 1)];
    for (int j = i + 1; j < m_new; ++j)
    {
      vj = x[Range(j * p, j * p + p_new - 1)];
      sumxx += std::pow(euclidean_norm(p_new), a);
    }
  }

  for (int i = 0; i < n_new; ++i)
  {
    vi = y[Range(i * p, i * p + p_new - 1)];
    for (int j = i + 1; j < n_new; ++j)
    {
      vj = y[Range(j * p, j * p + p_new - 1)];
      sumyy += std::pow(euclidean_norm(p_new), a);
    }
  }

  sumxy /= ((double)(m_new * n_new));
  sumxx /= ((double)(m_new * m_new));
  sumyy /= ((double)(n_new * n_new));

  return 2.0 * (sumxy - sumxx - sumyy);
}

double energy_distance_xyh(
    NumericVector x,
    NumericVector y,
    NumericVector h,
    int m,
    int n,
    int p,
    double a = 1.0)
{
  int lag = h.length() - 1;
  if (std::min(m, n) <= lag)
  {
    return -1.0;
  }
  int p_new = 2 * p;

  double sumxy, sumxx, sumyy;
  sumxy = sumxx = sumyy = 0.0;

  int k, l;
  NumericVector::iterator it, start_x, start_y, start_x2, start_y2;

  for (int i = 0; i < m; ++i)
  {
    k = 0;
    start_x = x.begin() + i * p;
    for (it = start_x; it != start_x + p; ++it)
    {
      vi[k] = *it;
      vi0[k] = *it;
      ++k;
    }
    for (int j = 0; j < n; ++j)
    {
      l = 0;
      start_y = y.begin() + j * p;
      for (it = start_y; it != start_y + p; ++it)
      {
        vj[l] = *it;
        vj0[l] = *it;
        ++l;
      }
      for (int lagg = 0; lagg <= lag; ++lagg)
      {
        if (h[lagg] == 0)
        {
          continue;
        }
        if (lagg > std::min(m - i - 1, n - j - 1))
        {
          break;
        }
        if (lagg == 0)
        {
          sumxy += h[0] * std::pow(euclidean_norm0(p_new), a) / ((double)(m * n));
          continue;
        }
        k = p;
        start_x2 = start_x + lagg * p;
        for (it = start_x2; it != start_x2 + p; ++it)
        {
          vi[k] = *it;
          ++k;
        }
        l = p;
        start_y2 = start_y + lagg * p;
        for (it = start_y2; it != start_y2 + p; ++it)
        {
          vj[l] = *it;
          ++l;
        }
        sumxy += h[lagg] * std::pow(euclidean_norm(p_new), a) / ((double)((m - lagg) * (n - lagg)));
      }
    }
  }

  for (int i = 0; i < m; ++i)
  {
    k = 0;
    start_x = x.begin() + i * p;
    for (it = start_x; it != start_x + p; ++it)
    {
      vi[k] = *it;
      vi0[k] = *it;
      ++k;
    }
    for (int j = i + 1; j < m; ++j)
    {
      l = 0;
      start_y = x.begin() + j * p;
      for (it = start_y; it != start_y + p; ++it)
      {
        vj[l] = *it;
        vj0[l] = *it;
        ++l;
      }
      for (int lagg = 0; lagg <= lag; ++lagg)
      {
        if (h[lagg] == 0)
        {
          continue;
        }
        if (lagg > std::min(m - i - 1, m - j - 1))
        {
          break;
        }
        if (lagg == 0)
        {
          sumxx += h[0] * std::pow(euclidean_norm0(p_new), a) / ((double)(m * m));
          continue;
        }
        k = p;
        start_x2 = start_x + lagg * p;
        for (it = start_x2; it != start_x2 + p; ++it)
        {
          vi[k] = *it;
          ++k;
        }
        l = p;
        start_y2 = start_y + lagg * p;
        for (it = start_y2; it != start_y2 + p; ++it)
        {
          vj[l] = *it;
          ++l;
        }
        sumxx += h[lagg] * std::pow(euclidean_norm(p_new), a) / ((double)((m - lagg) * (m - lagg)));
      }
    }
  }

  for (int i = 0; i < n; ++i)
  {
    k = 0;
    start_x = y.begin() + i * p;
    for (it = start_x; it != start_x + p; ++it)
    {
      vi[k] = *it;
      vi0[k] = *it;
      ++k;
    }
    for (int j = i + 1; j < n; ++j)
    {
      l = 0;
      start_y = y.begin() + j * p;
      for (it = start_y; it != start_y + p; ++it)
      {
        vj[l] = *it;
        vj0[l] = *it;
        ++l;
      }
      for (int lagg = 0; lagg <= lag; ++lagg)
      {
        if (h[lagg] == 0)
        {
          continue;
        }
        if (lagg > std::min(n - i - 1, n - j - 1))
        {
          break;
        }
        if (lagg == 0)
        {
          sumyy += h[0] * std::pow(euclidean_norm0(p_new), a) / ((double)(n * n));
          continue;
        }
        k = p;
        start_x2 = start_x + lagg * p;
        for (it = start_x2; it != start_x2 + p; ++it)
        {
          vi[k] = *it;
          ++k;
        }
        l = p;
        start_y2 = start_y + lagg * p;
        for (it = start_y2; it != start_y2 + p; ++it)
        {
          vj[l] = *it;
          ++l;
        }
        sumyy += h[lagg] * std::pow(euclidean_norm(p_new), a) / ((double)((n - lagg) * (n - lagg)));
      }
    }
  }

  return 2.0 * (sumxy - sumxx - sumyy);
}

//  [[Rcpp::export]]
double energy_distance_xy(NumericMatrix x, NumericMatrix y, int h = 0, double a = 1.0)
{
  int px = x.ncol();
  int py = y.ncol();
  if (px != py)
  {
    return -1.0;
  }
  NumericVector v1 = x;
  NumericVector v2 = y;
  return energy_distance_xyh(v1, v2, h, x.nrow(), y.nrow(), px, a);
}

//  [[Rcpp::export]]
NumericMatrix energy_distance_mat(NumericMatrix mat, NumericVector sizes, NumericVector h, double a = 1.0)
{
  int d = sizes.length();
  int p = mat.ncol();
  NumericMatrix ret_mat(d);
  std::vector<int> lindex(d);
  std::vector<int> rindex(d);
  int temp = 0;
  if (h.length() == 1)
  {
    int h0 = (int) h[0];
    NumericVector ti(p * h0 + p), tj(p * h0 + p);
    vi = ti;
    vj = tj;
    for (int i = 0; i < d; ++i)
    {
      lindex[i] = temp;
      temp += sizes[i];
      rindex[i] = temp - 1;
    }
    for (int i = 0; i < d; ++i)
    {
      NumericVector v1 = (NumericMatrix) mat(Range(lindex[i], rindex[i]), _);
      for (int j = i + 1; j < d; ++j)
      {
        NumericVector v2 = (NumericMatrix) mat(Range(lindex[j], rindex[j]), _);
        ret_mat(i, j) = ret_mat(j, i) = energy_distance_xyh(v1, v2, h0, sizes[i], sizes[j], p, a);
      }
    }
  }
  else
  {
    NumericVector ti(2 * p), tj(2 * p), ti0(p), tj0(p);
    vi = ti;
    vj = tj;
    vi0 = ti0;
    vj0 = tj0;
    for (int i = 0; i < d; ++i)
    {
      lindex[i] = temp;
      temp += sizes[i];
      rindex[i] = temp - 1;
    }
    for (int i = 0; i < d; ++i)
    {
      NumericVector v1 = (NumericMatrix) mat(Range(lindex[i], rindex[i]), _);
      for (int j = i + 1; j < d; ++j)
      {
        NumericVector v2 = (NumericMatrix) mat(Range(lindex[j], rindex[j]), _);
        ret_mat(i, j) = ret_mat(j, i) = energy_distance_xyh(v1, v2, h, sizes[i], sizes[j], p, a);
      }
    }
  }

  return ret_mat;
}