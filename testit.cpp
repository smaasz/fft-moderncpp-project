#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <complex>
#include <cstdlib>
#include <limits>
#include <string>
#include <getopt.h>
#include <fftw3.h>

#include "ffts.hpp"
#include "utils.hpp"

#define REPETITIONS 10

using namespace std;

// simple struct for test instances
template <typename complex_t>
struct TestInstance {
    int size;
    vector<complex_t> in;
    vector<complex_t> out;

    TestInstance() : size{}, in{}, out{} {};
    TestInstance(int size, vector<complex_t> in, vector<complex_t> out) : size{size}, in{in}, out{out} {};
};

// utility function to generate random test instances of a given size
template <typename complex_t>
TestInstance<complex_t> generate_test_instance(int size, int seed) {
    
    srand(seed);
    TestInstance<complex_t> test_instance;
    test_instance.size = size;
    
    for (int i = 0; i < size; ++i)
        test_instance.in.emplace_back( static_cast<double>(rand()) / RAND_MAX, static_cast<double>(rand()) / RAND_MAX );
    
    vector<complex<long double>>    test_in_ld(test_instance.in.begin(), test_instance.in.end());
    vector<complex<long double>>    test_out_ld = dft_matrix_mult<complex<long double>>(test_in_ld);
    
    test_instance.out = vector<complex_t>(test_out_ld.begin(), test_out_ld.end());
    
    return test_instance;
}

// enum type for the different algorithms
enum Algorithm {recursive_depth_first, iterative_breadth_first, fftw_lib};

// simple struct to store setup information i.e. parameters of the radix computation
struct SetupInfo {
    int                     radix_option;
    int                     radix_threshold;
    static constexpr int    not_used = std::numeric_limits<int>::max();
};

// function that tests the accuracy of a specified FFT algorithm
template <typename complex_t>
void test_accuracy(string const& text, Algorithm a, SetupInfo const& setup_info) {
    using std::fixed;
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    // generate test instances
    vector<TestInstance<complex_t>> test_instances; 
    test_instances.push_back(generate_test_instance<complex_t>(2 * 2 * 2 * 3 * 5 * 7, 43));
    test_instances.push_back(generate_test_instance<complex_t>(3 * 5 * 11 * 13, 43));
    test_instances.push_back(generate_test_instance<complex_t>(2 * 2 * 2 * 2 * 3 * 37, 43));
    test_instances.push_back(generate_test_instance<complex_t>(27000, 43));
    int size = 2 * 2 * 2 * 2 * 3 * 37;
    test_instances.emplace_back(size, vector<complex_t>(size, 0.0), vector<complex_t>(size, 0.0));

    cout << text << endl;

    {   
        cout << " size   time (ms)  accuracy (max-norm)   accuracy (two-norm)" << endl;

        for (auto test_instance : test_instances) {
            
            vector<complex_t>               out;
            vector<int>                     radices;
            duration<double, std::milli>    duration_ms;
            auto                            start_time_ms = high_resolution_clock::now();

            switch (a)
            {
            case recursive_depth_first:
                radices         = compute_radices(test_instance.size, setup_info.radix_option, setup_info.radix_threshold);
                start_time_ms   = high_resolution_clock::now();
                out             = fft_recursive_depth_first(test_instance.in, radices);
                duration_ms     = high_resolution_clock::now() - start_time_ms;
                break;
            case iterative_breadth_first:
                radices         = compute_radices(test_instance.size, setup_info.radix_option, setup_info.radix_threshold);
                start_time_ms   = high_resolution_clock::now();
                out             = fft_iterative_breadth_first(test_instance.in, radices);
                duration_ms     = high_resolution_clock::now() - start_time_ms;
                break;
            case fftw_lib:
                fftw_complex    *in;
                fftw_complex    *out_fftw;
                fftw_plan       p;

                //fftw_init_threads();

                in          = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * test_instance.size);
                out_fftw    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * test_instance.size);

                //fftw_plan_with_nthreads(4);

                //int max_threads = fftw_planner_nthreads();

                p = fftw_plan_dft_1d(test_instance.size, in, out_fftw, FFTW_FORWARD, FFTW_ESTIMATE);

                for (int i = 0; i < test_instance.size; ++i) {
                    in[i][0] = test_instance.in[i].real();
                    in[i][1] = test_instance.in[i].imag();
                }
                
                start_time_ms = high_resolution_clock::now();

                fftw_execute(p);

                duration_ms = high_resolution_clock::now() - start_time_ms;

                for (int i = 0; i < test_instance.size; ++i) {
                    out.emplace_back(out_fftw[i][0], out_fftw[i][1]);
                }

                fftw_destroy_plan(p);
                fftw_free(in);
                fftw_free(out_fftw);
                break;
            }
            
            int const default_precision = static_cast<int>(std::cout.precision());
            cout << setw(5) << test_instance.size
                    << setw(10) << setprecision(2) << fixed << duration_ms.count()
                    << setw(18) << fixed << setprecision(12) << max_norm(test_instance.out - out)
                    << setw(22) << fixed << setprecision(12) << two_norm(test_instance.out - out)
                    << endl;
            cout << setprecision(default_precision);
        }

        cout << endl;
    }
}

template <typename complex_t>
void test_speed(string const& text, Algorithm a, SetupInfo const& setup_info, bool use_powers_of_2) {
    using std::fixed;
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    // powers-of-2
    vector<TestInstance<complex_t>> test_instances;
    if (use_powers_of_2) 
    {
        int size = 32;
        for (size_t i = 0; i < 10; ++i) {
            size *= 2;
            test_instances.emplace_back(size, vector<complex_t>(size, 0.0), vector<complex_t>(size, 0.0));
        }
    } else
    {
        vector<int> sizes = {6, 9, 12, 15, 18, 24, 36, 80, 108, 210, 504, 1000, 1960, 4725, 10368, 27000, 75600, 165375};
        for (int size : sizes)
            test_instances.emplace_back(size, vector<complex_t>(size, 0.0), vector<complex_t>(size, 0.0));
    }

    cout << text << endl;

    {   
        cout << "  size   time (ms)   repetitions   average (ms)" << endl;

        for (auto test_instance : test_instances) {
            
            vector<complex_t>               out;
            vector<int>                     radices;
            duration<double, std::milli>    duration_ms;
            auto start_time_ms = high_resolution_clock::now();

            switch (a)
            {
            case recursive_depth_first:
                radices = compute_radices(test_instance.size, setup_info.radix_option, setup_info.radix_threshold);
                start_time_ms = high_resolution_clock::now();
                for (int i = 0; i < REPETITIONS; ++i)
                    out = fft_recursive_depth_first(test_instance.in, radices);
                duration_ms    = high_resolution_clock::now() - start_time_ms;
                break;
            case iterative_breadth_first:
                radices = compute_radices(test_instance.size, setup_info.radix_option, setup_info.radix_threshold);
                start_time_ms  = high_resolution_clock::now();
                
                for (int i = 0; i < REPETITIONS; ++i)
                    out = fft_iterative_breadth_first(test_instance.in, radices);
                
                duration_ms = high_resolution_clock::now() - start_time_ms;
                break;
            case fftw_lib:
                fftw_complex    *in;
                fftw_complex    *out_fftw;
                fftw_plan       p;

                //fftw_init_threads();

                in          = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * test_instance.size);
                out_fftw    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * test_instance.size);

                //fftw_plan_with_nthreads(4);

                //int max_threads = fftw_planner_nthreads();

                p = fftw_plan_dft_1d(test_instance.size, in, out_fftw, FFTW_FORWARD, FFTW_ESTIMATE);

                for (int i = 0; i < test_instance.size; ++i) {
                    in[i][0] = test_instance.in[i].real();
                    in[i][1] = test_instance.in[i].imag();
                }
                
                start_time_ms = high_resolution_clock::now();

                for (int i = 0; i < REPETITIONS; ++i)
                    fftw_execute(p);

                duration_ms = high_resolution_clock::now() - start_time_ms;

                fftw_destroy_plan(p);
                fftw_free(in);
                fftw_free(out_fftw);
                break;
            }
            
            int const default_precision = static_cast<int>(std::cout.precision());
            cout << setw(6) << test_instance.size
                    << setw(10) << setprecision(2) << fixed << duration_ms.count()
                    << setw(12) << fixed << setprecision(0) << REPETITIONS
                    << setw(14) << fixed << setprecision(4) << duration_ms.count() / static_cast<double>(REPETITIONS)
                    << endl;
            cout << setprecision(default_precision);
        }
    }
}



int main(int argc, char ** argv){
    
    constexpr char const* const options = "a:g:hnp:r:st:";
    constexpr char const* const usage = " [options]\n" \
        " -a n         Choose algorithm: 1 = iterative, 2 = recursive, 3 = FFTW (3)\n" \
        " -g n         Choose algorithm for radix generation: 1 = factors, 2 = factors reversed, 3 = thresholded (1)\n" \
        " -r n         Threshold for radix generation (not used)\n" \
        " -t n         Choose test: 1 = performance, 2 = accuracy (1)\n" \
        " -n           Use non-powers-of-2\n" \
        " -s           Use single precision\n" \
        " -p text      Print text before the test\n" \
        " -h           Show this help\n" \
        "\n";
    
    using std::fixed;
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    try 
    {
        int     algo                    = 3;
        int     algo_radix              = SetupInfo::not_used;
        int     radix_threshold         = SetupInfo::not_used;
        int     test_type               = 1;
        bool    use_powers_of_two       = true;
        bool    use_single_precision    = false;
        string  preamble                = "";
        int     c;

        while((c = getopt(argc, argv, options)) != -1)
        {
            switch(c)
            {
            case 'a' :
                algo = stoi(optarg);
                break;
            case 'g' :
                algo_radix = stoi(optarg);
                break;
            case 'r' :
                radix_threshold = stoi(optarg);
                break;
            case 't' :
                test_type = stoi(optarg);
                break;
            case 'n' :
                use_powers_of_two = false;
                break;
            case 's' :
                use_single_precision = true;
                break;
            case 'p' :
                preamble = string(optarg);
                break;
            case 'h' :
                cout << "usage: " << argv[0] << usage;
                return 0;
            case '?' :
            default :
                cerr << "unknown option [" << c << "]\n";
                cerr << "usage: " << argv[0] << usage;
                return -1;
            }
        }
        if (algo < 1 or algo > 3 or (algo_radix != SetupInfo::not_used and (algo_radix < 1 or algo_radix > 3)) or test_type < 1 or test_type > 2)
        {
            cerr << algo << " " << algo_radix << " " << test_type << endl;
            cerr << "usage: " << argv[0] << usage;
            return -2;
        }
        
        if ((algo == 1 or algo == 2) and algo_radix == SetupInfo::not_used)
        {
            cerr << "for algorithm " << algo << " a radix generation algorithm must be specified" << endl;
            cerr << "usage: " << argv[0] << usage;
            return -3;
        }

        Algorithm a;
        switch (algo)
        {
        case 1:
            a = iterative_breadth_first;
            break;
        case 2:
            a = recursive_depth_first;
            break;
        case 3:
            a = fftw_lib;
            break;
        }

        SetupInfo setup_info{algo_radix, radix_threshold};

        if (test_type == 1) {

            if (use_single_precision)
                test_speed<complex<float>>(preamble, a, setup_info, use_powers_of_two);
            else
                test_speed<complex<double>>(preamble, a, setup_info, use_powers_of_two);

        } else
        {

            if (use_single_precision)
                test_accuracy<complex<float>>(preamble, a, setup_info);
            else
                test_accuracy<complex<double>>(preamble, a, setup_info);

        }

    }
    catch(exception const& e)
    {
        cerr << argv[0] << ": Exception " << e.what() << " -- aborting\n";
        return -4;
    }


    // SetupInfo setup_info{1, SetupInfo::not_used};
    // test_accuracy<complex<float>>("iterative, float, radices=factors", iterative_breadth_first, setup_info);
    // test_accuracy<complex<double>>("iterative, double, radices=factors", iterative_breadth_first, setup_info);

    // test_accuracy<complex<float>>("recursive, float, radices=factors", recursive_depth_first, setup_info);
    // test_accuracy<complex<double>>("recursive, double, radices=factors", recursive_depth_first, setup_info);

    // setup_info = {SetupInfo::not_used, SetupInfo::not_used};
    // test_accuracy<complex<double>>("FFTW, double, radices=factors", fftw_lib, setup_info);

    // bool use_powers_of_2 = false;
    // setup_info = {3, 32};
    // test_speed<complex<float>>("iterative, float, radices=factors", iterative_breadth_first, setup_info, use_powers_of_2);
    // test_speed<complex<double>>("iterative, double, radices=factors", iterative_breadth_first, setup_info, use_powers_of_2);

    // test_speed<complex<float>>("recursive, float, radices=factors", recursive_depth_first, setup_info, use_powers_of_2);
    // test_speed<complex<double>>("recursive, double, radices=factors", recursive_depth_first, setup_info, use_powers_of_2);

    // setup_info = {SetupInfo::not_used, SetupInfo::not_used};
    // test_speed<complex<double>>("FFTW, double, radices=factors", fftw_lib, setup_info, use_powers_of_2);




    // int                 size    = 1 << 10;
    
    // {// test own implementation: recursive, depth-first
    //     cout << "Own implementation (recursive, depth-first): " << endl;

    //     vector<int>         radices = compute_radices(size, 3, 32);
    //     vector<complex_t>   x(size, 0.0);

    //     auto start_time_ms = high_resolution_clock::now();

    //     vector<complex_t>   y = fft_recursive_depth_first(x, radices);

    //     duration<double, std::milli> duration_ms = high_resolution_clock::now() - start_time_ms;
        
    //     const auto default_precision {std::cout.precision()};
    //     cout << "time= " << setprecision(0) << fixed << duration_ms.count() << " ms" << endl;
    //     cout << setprecision(default_precision);

    // }

    // {// test own implementation: iterative, breadth-first
    //     cout << "Own implementation (iterative, breadth-first): " << endl;

    //     vector<int>         radices = compute_radices(size, 3, 32);
    //     vector<complex_t>   x(size, 0.0);

    //     auto start_time_ms = high_resolution_clock::now();

    //     vector<complex_t>   y = fft_iterative_breadth_first(x, radices);

    //     duration<double, std::milli> duration_ms = high_resolution_clock::now() - start_time_ms;

    //     const auto default_precision {std::cout.precision()};
    //     cout << "time= " << setprecision(0) << fixed << duration_ms.count() << " ms" << endl;
    //     cout << setprecision(default_precision);

    // }

    // {// matrix multiplication
    //     cout << "Matrix multiplication: " << endl;

    //     vector<complex_t> x(size, 0.0);

    //     auto start_time_ms = high_resolution_clock::now();

    //     vector<complex_t> y = dft_matrix_mult<complex_t>(x);

    //     duration<double, std::milli> duration_ms = high_resolution_clock::now() - start_time_ms;

    //     const auto default_precision {std::cout.precision()};
    //     cout << "time= " << setprecision(0) << fixed << duration_ms.count() << " ms" << endl;
    //     cout << setprecision(default_precision);
    // }


    // {// test fftw
    //     cout << "FFTW3: " << endl;

    //     fftw_complex    *in;
    //     fftw_complex    *out;
    //     fftw_plan       p;

    //     //fftw_init_threads();

    //     in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size);
    //     out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size);

    //     //fftw_plan_with_nthreads(4);

    //     //int max_threads = fftw_planner_nthreads();

    //     p = fftw_plan_dft_1d(size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    //     for (int i = 0; i < size; ++i) {
    //         in[i][0] = 0.0;
    //         in[i][1] = 0.0;
    //     }
    //     in[0][0] = 1.0;
        
    //     auto start_time_ms = high_resolution_clock::now();

    //     fftw_execute(p);

    //     duration<double, std::milli> duration_ms = high_resolution_clock::now() - start_time_ms;

    //     const auto default_precision {std::cout.precision()};
    //     cout << "time= " << setprecision(0) << fixed << duration_ms.count() << " ms" << endl;
    //     cout << setprecision(default_precision);

    //     //cout << max_threads << endl;

    //     // for (int i = 0; i < size; ++i)
    //     //     cout << setprecision(4) << out[i][0] << " + i * " << out[i][1] << endl;

    //     fftw_destroy_plan(p);
    //     fftw_free(in);
    //     fftw_free(out);

    //     //fftw_cleanup_threads();
    // }

    // { // powers-of-2
    //     cout << "Powers of 2:" << endl;

    //     vector<int> power2 = {128};
    //     for (int i = 0; i < 8; ++i)
    //         power2.push_back(2 * power2.back());
    //     vector<duration<double, std::milli>> durations_ms(6, duration<double, std::milli>{});
    //     for (int length : power2) {
    //         vector<int>         radices = compute_radices(length, 3, 32);
    //         vector<complex_t>   x(length, 0.0);

    //         auto start_time_ms = high_resolution_clock::now();

    //         vector<complex_t>   y = fft_recursive_depth_first(x, radices);

    //         duration<double, std::milli> duration_ms = high_resolution_clock::now() - start_time_ms;

    //         const auto default_precision {std::cout.precision()};
    //         cout << "time= " << setprecision(0) << fixed << duration_ms.count() << " ms" << endl;
    //         cout << setprecision(default_precision);
    //     }
    // }

    // {
    //     cout << "Own implementation (iterative, breadth-first): " << endl;

    //     vector<int>         radices = compute_radices(8, 1, 32);
    //     vector<complex_t>   x(8, 0.0);
    //     x[0] = 3.0;
    //     x[1] = 2.0;

    //     vector<complex_t>   y = fft_iterative_breadth_first(x, radices);

    //     const auto default_precision{std::cout.precision()};
    //     cout << setprecision(4);
    //     for (auto element : y)
    //         cout << element << endl;
    //     cout << endl;
    //     cout << setprecision(default_precision);
    // }

    // {// test fftw
    //     cout << "FFTW3: " << endl;

    //     fftw_complex    *in;
    //     fftw_complex    *out;
    //     fftw_plan       p;

    //     //fftw_init_threads();

    //     in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 8);
    //     out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 8);

    //     //fftw_plan_with_nthreads(4);

    //     //int max_threads = fftw_planner_nthreads();

    //     p = fftw_plan_dft_1d(8, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    //     for (int i = 0; i < 8; ++i) {
    //         in[i][0] = 0.0;
    //         in[i][1] = 0.0;
    //     }
    //     in[0][0] = 3.0;
    //     in[1][0] = 2.0;
        
    //     auto start_time_ms = high_resolution_clock::now();

    //     fftw_execute(p);

    //     duration<double, std::milli> duration_ms = high_resolution_clock::now() - start_time_ms;
    //     cout << "time= " << setprecision(0) << fixed << duration_ms.count() << " ms" << endl; 

    //     // cout << max_threads << endl;

    //     for (int i = 0; i < 8; ++i)
    //         cout << setprecision(4) << out[i][0] << " + i * " << out[i][1] << endl;

    //     fftw_destroy_plan(p);
    //     fftw_free(in);
    //     fftw_free(out);

    //     //fftw_cleanup_threads();
    // }

    // { // test digit-reversal
    //     vector<int> v = {0,1,2,3};
    //     vector<int> radices{2,2};
    //     digit_reversal<int>(v, radices);
    //     for (auto element : v)
    //         std::cout << element << " ";
    //     std::cout << std::endl;
    // }

}