#ifndef FFTS_H_
#define FFTS_H_

#include <stack>
#include <vector>
#include <cassert>
#include <complex>
#include <numeric>
#include <algorithm>
#include <iostream>

#include "StridedVector.hpp"

constexpr auto PI = 3.14159265358979323846264338327950288419716939937510L;

std::vector<int> compute_radices(int size, int option, int threshold);

std::vector<std::vector<std::complex<long double>>> precompute_phases(std::vector<int>& radices);

template <typename InputIt>
std::vector<int> compute_digits(int value, InputIt radix_low, InputIt radix_high) 
{

    std::vector<int> digits;

    for (; radix_low != radix_high; ++radix_low)
    {
        digits.push_back(value % (*radix_low));
        value = value / (*radix_low);
    }

    return digits;
}

// utility function to calculate the digit reversed value of a number
template <typename InputIt>
int reverse_digits(int value, InputIt radix_low, InputIt radix_high) {
    
    assert(radix_low != radix_high);

    int j;
    int digit;
    
    j       = value % (*radix_low);
    value   = value / (*radix_low);
    ++radix_low;
    
    for (; radix_low != radix_high; ++radix_low)
    {
    
        digit   = value % (*radix_low);
        value   = value / (*radix_low);
        j       = j * (*radix_low) + digit;
    
    }

    return j;
}

// permute elements of a vector by reversing the digits of the indices
template <typename T>
void permute_by_digit_reversal(std::vector<T>& v, std::vector<int>& radices)
{
    std::vector<int> permutation(v.size(), 0);
    
    for (int i = 0; i < static_cast<int>(v.size()); ++i)
        permutation[reverse_digits(i, radices.rbegin(), radices.rend())] = i;

    // apply permutation
    std::vector<T> w;
    w.reserve(v.size());
    for (size_t i = 0; i < v.size(); i++)
        w.push_back(std::move(v[permutation[i]]));
    v = std::move(w);
}

// Discrete Fourier Transform by matrix multiplication O(n^2) runtime complexity
template <typename complex_t>
std::vector<complex_t> dft_matrix_mult(std::vector<complex_t> &in)
{
    
    int size = static_cast<int>(in.size());
    
    std::vector<complex_t> out;
    out.resize(size);

    complex_t y;
    complex_t phase;
    complex_t phase_step;

    for (int k = 0; k < size; ++k) 
    {
        
        phase_step  = static_cast<complex_t>(std::polar(1.0L, -2 * PI / size * k));
        y           = 0.0;
        
        for (int j = size - 1; j >= 0; --j) 
        {
            
            y = y * phase_step + in[j];
        
        }
        
        out[k] = y;
    
    }

    return out;
}


template <typename complex_t>
void solve_dft_recursive(StridedVector<complex_t>& in, std::vector<int> radices)
{
    assert(radices.size() > 0 and accumulate(radices.begin(), radices.end(), 1, std::multiplies<int>()) == static_cast<int>(in.size()));

    int size = static_cast<int>(in.size());

    complex_t phase_step;

    // base case: column size will not be further reduced
    if (radices.size() == 1)
    {
        // calculate DFT in O(n^2)
        std::vector<complex_t> buffer(size, 0.0);
        complex_t y;
        for (int k = 0; k < size; ++k)
        {
            
            phase_step  = static_cast<complex_t>(std::polar(1.0L, -2 * PI / size * k));
            y           = 0.0;
            
            for (int j = size - 1; j >= 0; --j) {
            
                y = y * phase_step + in[j];
            
            }
            
            buffer[k] = y;
        }
        for (int k = 0; k < size; ++k)
            in[k] = buffer[k];
    }
    else
    {
        // Interpret in as matrix in row-major format with row length = radix.
        // For each column calculate its DFT.
        int radix   = radices.back();
        int rest    = size / radix;

        radices.pop_back();

        // solve DFT of each column recursively
        std::vector<StridedVector<complex_t>> columns;
        for (int l = 0; l < radix; ++l)
        {

            columns.emplace_back(in, radix, l, rest);
            solve_dft_recursive(columns[l], radices);
        
        }

        // butterfly: DFT(in, i * radix + j)
        complex_t y;
        complex_t twiddle_factor_step;
        std::vector<complex_t> buffer(size, 0.0);
        
        for (int k0 = 0; k0 < rest; ++k0)
        {
            
            twiddle_factor_step = static_cast<complex_t>(std::polar(1.0L, -2 * PI / size * k0));

            for (int k1 = 0; k1 < radix; ++k1)
            {

                phase_step  = static_cast<complex_t>(std::polar(1.0L, -2 * PI / radix * k1));
                y           = 0.0;

                for (int j0 = radix - 1; j0 >= 0; --j0)
                {

                    y   = y * twiddle_factor_step * phase_step + columns[j0][k0];

                }

                buffer[k1 * rest + k0] = y;

            }

        }

        for (int k = 0; k < size; ++k)
            in[k] = buffer[k];
   
    }
}

// Cooley-Tuckey type implementation of the DFT by decimation in time, depth-first, mixed-radix
template <typename complex_t>
std::vector<complex_t> fft_recursive_depth_first(std::vector<complex_t>& x, std::vector<int>& radices)
{
    assert( std::accumulate(radices.begin(), radices.end(), 1, std::multiplies<int>()) == static_cast<int>(x.size()) );

    StridedVector<complex_t> strided_in{x}; 

    solve_dft_recursive<complex_t>(strided_in, radices);

    return x;
}


// Cooley-Tuckey type implementation of the DFT by decimation in time, breadth-first, mixed-radix
template <typename complex_t>
std::vector<complex_t> fft_iterative_breadth_first(std::vector<complex_t> &x, std::vector<int> &radices)
{
    assert( std::accumulate(radices.begin(), radices.end(), 1, std::multiplies<int>()) == static_cast<int>(x.size()) );

    std::vector<int> prev_radices{};
    // base cases
    int radix   = radices[0];
    int size    = radix;
    int nrows   = static_cast<int>(x.size()) / radix;
    
    complex_t phase_step;

    std::vector<complex_t> buffer(radix, 0.0);
    for (int low = 0; low < nrows; ++low)
    {
        
        for (int k = 0; k < radix; ++k)
        {
        
            phase_step  = static_cast<complex_t>(std::polar(1.0L, -2 * PI / size * k));
        
            for (int j = radix - 1; j >= 0; --j)
            {
        
                buffer[k]   *= phase_step;
                buffer[k]   += x[j * nrows + low];
        
            }
        }
        
        for (int k = 0; k < radix; ++k)
        {

            x[k * nrows + low] = buffer[k];
            buffer[k] = 0.0;

        }
    }

    prev_radices.push_back(radix);

    // iteration
    complex_t                   y;
    complex_t                   twiddle_factor_step;

    for (size_t i = 1; i < radices.size(); ++i) {

        radix   = radices[i];
        size    = size * radix;
        nrows   = nrows / radix;

        buffer  = std::vector<complex_t>(radix, 0.0);
        int l;

        for (int high = 0; high < size / radix; ++high)
        {
        
            l                   = reverse_digits(high, prev_radices.rbegin(), prev_radices.rend());
            twiddle_factor_step = static_cast<complex_t>(std::polar(1.0L, -2 * PI / size * l));

            for (int low = 0; low < nrows; ++low)
            {
            
                for (int k = 0; k < radix; ++k)
                {
            
                    phase_step  = static_cast<complex_t>(std::polar(1.0L, -2 * PI / radix * k));
            
                    for (int j = radix - 1; j >= 0; --j)
                    {
            
                        buffer[k] *= twiddle_factor_step * phase_step;
                        buffer[k] += x[(high * radix + j) * nrows + low];
            
                    }
                }
            
                for (int k = 0; k < radix; ++k)
                {

                    x[(high * radix + k) * nrows + low] = buffer[k];
                    buffer[k] = 0.0;

                }
            }
        }

        prev_radices.push_back(radix);

    }
    permute_by_digit_reversal<complex_t>(x, radices);

    return x;
}

#endif
