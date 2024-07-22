#include <vector>
#include <cmath>
#include <complex>

// adapted from: 17-cholesky/matrix.hpp by Thorsten Koch

template <typename T>
T max_norm(std::vector<std::complex<T>> const& vec)
{
    T max_x = 0.0;
    T abs_x;

    for(auto x : vec) {
        abs_x = std::abs(x);
        if (abs_x > max_x)
            max_x = abs_x;
    }

    return max_x;
}

template <typename T>
T two_norm(std::vector<std::complex<T>> const& vec)
{
   T sum = 0.0;
   
   for(auto x : vec)
      sum += std::norm(x);

   return std::sqrt(sum);
}

template <typename T>
std::vector<T> operator-(std::vector<T> const& a, std::vector<T> const& b)
{
   assert(a.size() == b.size());

   std::vector<T> r(a.size(), 0.0);
   
   for(size_t i = 0; i < a.size(); ++i)
      r[i] = a[i] - b[i];

   return r;
}