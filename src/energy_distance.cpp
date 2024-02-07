#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

double euclidean_norm(NumericVector z)
{
  int d = z.length();
  double norm_z = 0.0;
  for (int i = 0; i < d; ++i)
  {
    norm_z = norm_z + z[i] * z[i];
  }
  return std::pow(norm_z, 0.5);
}

double energy_distance_xyh(
    NumericMatrix x,
    NumericMatrix y,
    int h,
    double a = 1.0)
{
  int m = x.nrow();
  int p = x.ncol();
  int n = y.nrow();
  int p_y = y.ncol();
  if (p != p_y)
  {
    return -1.0;
  }
  if (std::min(m, n) <= h)
  {
    return -1.0;
  }

  NumericVector flat_x = x;
  NumericVector flat_y = y;

  int m_new = m - h;
  int n_new = n - h;
  int p_new = (h + 1) * p;

  double sumxy, sumxx, sumyy;
  sumxy = sumxx = sumyy = 0.0;

  NumericVector vi(p_new);
  NumericVector vj(p_new);

  for (int i = 0; i < m_new; ++i)
  {
    vi = flat_x[Range(i * p, i * p + p_new - 1)];
    for (int j = 0; j < n_new; ++j)
    {
      vj = flat_y[Range(j * p, j * p + p_new - 1)];
      sumxy += std::pow(euclidean_norm(vj - vi), a);
    }
  }

  for (int i = 0; i < m_new; ++i)
  {
    vi = flat_x[Range(i * p, i * p + p_new - 1)];
    for (int j = i + 1; j < m_new; ++j)
    {
      vj = flat_x[Range(j * p, j * p + p_new - 1)];
      sumxx += std::pow(euclidean_norm(vj - vi), a);
    }
  }

  for (int i = 0; i < n_new; ++i)
  {
    vi = flat_y[Range(i * p, i * p + p_new - 1)];
    for (int j = i + 1; j < n_new; ++j)
    {
      vj = flat_y[Range(j * p, j * p + p_new - 1)];
      sumyy += std::pow(euclidean_norm(vj - vi), a);
    }
  }

  sumxy /= ((double)(m_new * n_new));
  sumxx /= ((double)(m_new * m_new));
  sumyy /= ((double)(n_new * n_new));

  return 2.0 * (sumxy - sumxx - sumyy);
}

double energy_distance_xyh(
    NumericMatrix x,
    NumericMatrix y,
    NumericVector h,
    double a = 1.0)
{
  int m = x.nrow();
  int p = x.ncol();
  int n = y.nrow();
  int p_y = y.ncol();
  int lag = h.length() - 1;
  if (p != p_y)
  {
    return -1.0;
  }
  if (std::min(m, n) <= lag)
  {
    return -1.0;
  }

  NumericVector flat_x = x;
  NumericVector flat_y = y;

  int p_new = 2 * p;

  double sumxy, sumxx, sumyy;
  sumxy = sumxx = sumyy = 0.0;

  NumericVector vi(p_new);
  NumericVector vj(p_new);
  NumericVector vi0(p);
  NumericVector vj0(p);

  int k, l;
  NumericVector::iterator it, start_x, start_y, start_x2, start_y2;

  for (int i = 0; i < m; ++i)
  {
    k = 0;
    start_x = flat_x.begin() + i * p;
    for (it = start_x; it != start_x + p; ++it)
    {
      vi[k] = *it;
      vi0[k] = *it;
      ++k;
    }
    for (int j = 0; j < n; ++j)
    {
      l = 0;
      start_y = flat_y.begin() + j * p;
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
          sumxy += h[0] * std::pow(euclidean_norm(vj0 - vi0), a) / ((double)(m * n));
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
        sumxy += h[lagg] * std::pow(euclidean_norm(vj - vi), a) / ((double)((m - lagg) * (n - lagg)));
      }
    }
  }

  for (int i = 0; i < m; ++i)
  {
    k = 0;
    start_x = flat_x.begin() + i * p;
    for (it = start_x; it != start_x + p; ++it)
    {
      vi[k] = *it;
      vi0[k] = *it;
      ++k;
    }
    for (int j = i + 1; j < m; ++j)
    {
      l = 0;
      start_y = flat_x.begin() + j * p;
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
          sumxx += h[0] * std::pow(euclidean_norm(vj0 - vi0), a) / ((double)(m * m));
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
        sumxx += h[lagg] * std::pow(euclidean_norm(vj - vi), a) / ((double)((m - lagg) * (m - lagg)));
      }
    }
  }

  for (int i = 0; i < n; ++i)
  {
    k = 0;
    start_x = flat_y.begin() + i * p;
    for (it = start_x; it != start_x + p; ++it)
    {
      vi[k] = *it;
      vi0[k] = *it;
      ++k;
    }
    for (int j = i + 1; j < n; ++j)
    {
      l = 0;
      start_y = flat_y.begin() + j * p;
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
          sumyy += h[0] * std::pow(euclidean_norm(vj0 - vi0), a) / ((double)(n * n));
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
        sumyy += h[lagg] * std::pow(euclidean_norm(vj - vi), a) / ((double)((n - lagg) * (n - lagg)));
      }
    }
  }

  return 2.0 * (sumxy - sumxx - sumyy);
}

//  [[Rcpp::export]]
double energy_distance_xy(NumericMatrix x, NumericMatrix y, int h = 0, double a = 1.0)
{
  return energy_distance_xyh(x, y, h, a);
}

//  [[Rcpp::export]]
NumericMatrix energy_distance_mat(NumericMatrix mat, NumericVector sizes, NumericVector h, double a = 1.0)
{
  int p = mat.ncol();
  int d = sizes.length();
  NumericMatrix ret_mat(d);
  std::vector<int> lindex(d);
  std::vector<int> rindex(d);
  int temp = 0;
  if (h.length() == 1)
  {
    int h0 = (int) h[0];
    for (int i = 0; i < d; ++i)
    {
      lindex[i] = temp;
      temp += sizes[i];
      rindex[i] = temp - 1;
    }
    for (int i = 0; i < d; ++i)
    {
      NumericMatrix v1(sizes[i], p);
      v1 = mat(Range(lindex[i], rindex[i]), _);
      for (int j = i + 1; j < d; ++j)
      {
        NumericMatrix v2(sizes[j], p);
        v2 = mat(Range(lindex[j], rindex[j]), _);
        ret_mat(i, j) = ret_mat(j, i) = energy_distance_xyh(v1, v2, h0, a);
      }
    }
  }
  else
  {
    for (int i = 0; i < d; ++i)
    {
      lindex[i] = temp;
      temp += sizes[i];
      rindex[i] = temp - 1;
    }
    for (int i = 0; i < d; ++i)
    {
      NumericMatrix v1(sizes[i], p);
      v1 = mat(Range(lindex[i], rindex[i]), _);
      for (int j = i + 1; j < d; ++j)
      {
        NumericMatrix v2(sizes[j], p);
        v2 = mat(Range(lindex[j], rindex[j]), _);
        ret_mat(i, j) = ret_mat(j, i) = energy_distance_xyh(v1, v2, h, a);
      }
    }
  }

  return ret_mat;
}