#include <vector>
#include <algorithm>
#include <numeric>
#include <cassert>

#include "ffts.hpp"

std::vector<int> compute_radices(int size, int option, int threshold) {
    assert( size > 0 );
    assert( 1 <= option and option <= 3 );
    assert( 2 <= threshold );
    
    using std::accumulate;
    using std::vector;
    using std::reverse;

    vector<int> factors{};
    vector<int> radices{};

    int n = size;
    // first compute all factors of n
    for (int i = 2; i * i <= n; ++i){
        while (n % i == 0) {
            n = n / i;
            factors.push_back(i);
        }
    }
    if (n > 1)
        factors.push_back(n);

    assert( accumulate(factors.begin(), factors.end(), 1, std::multiplies<int>()) == size);

    // different choices of radices possible
    // 1. factors in increasing order
    // 2. factors in decreasing order
    // 3. group factors if below threshold (starting from smallest or largest factor?)
    switch (option) {
        case 1:
            return factors;
        case 2:
            reverse(factors.begin(), factors.end());
            return factors;
        case 3:
            int radix_count = 0;
            int radix       = factors[0];
            int factor;

            for (size_t i = 1; i < factors.size(); ++i) {
                factor = factors[i];
                
                if (radix * factor <= threshold) {
                    radix *= factor;
                }
                else {
                    radices.push_back(radix);
                    ++radix_count;
                    radix = factor;
                }
            }
            radices.push_back(radix);
            return radices;
    }

    assert( accumulate(radices.begin(), radices.end(), 1, std::multiplies<int>()) == size);

    return radices;
}

std::vector<std::vector<std::complex<long double>>> precompute_phases(std::vector<int>& radices) 
{

    using std::vector;
    using std::complex;
    using std::polar;
    using std::accumulate;

    constexpr auto PI = 3.14159265358979323846264338327950288419716939937510L;

    vector<vector<complex<long double>>> phase_table(radices.size(), vector<complex<long double>>{});
    
    int     radix_partial_product   = accumulate(radices.begin(), radices.end(), 1, std::multiplies<int>());
    auto    it                      = radices.begin();
    
    for (size_t i = 0; i < radices.size(); ++i)
    {
        ++it;

        phase_table[i]          = vector<complex<long double>>(radix_partial_product, 1.0);
        phase_table[i][1]       = polar(1.0L, -2 * PI / radices[i]);
        
        int j;
        for (j = 2; j < radices[i]; ++j)
            phase_table[i][j] = phase_table[i][j-1] * phase_table[i][1];

        for(; j < radix_partial_product; ++j) 
        {
            vector<int> digits = compute_digits(j, radices.begin(), it);
            
            for (size_t k = 0; k < i; ++k) 
                phase_table[i][j] *= phase_table[k][digits[k]];

        }

        radix_partial_product   /= radices[i];
    }

    return phase_table;
}