#include <vector>
#include <iostream>


// class for strided access with an offset of the underlying vector
template <typename T>
class StridedVector
{
private:
    std::vector<T>& v;
    int const stride;
    int const offset;
    size_t const strided_size;

public:
    T& operator[](int i)
    {
        assert( 0 <= i and i < static_cast<int>(strided_size) );

        return v[offset + i * stride];
    };
    T const& operator[](int i) const
    {
        assert( 0 <= i and i < static_cast<int>(strided_size) );

        return v[offset + i * stride];
    };

    size_t size() const
    {
        return strided_size;
    };

    explicit StridedVector(std::vector<T>& v) : v{v}, stride{1}, offset{0}, strided_size{v.size()} {};
    explicit StridedVector(std::vector<T>& v, int const stride, int const offset, int const size) : v{v}, stride{stride}, offset{offset}, strided_size(size)
    {
        assert(0 <= offset and offset < v.size() and stride > 0);
    };
    explicit StridedVector(StridedVector<T>& sv, int const stride, int const offset, int const size) : v{sv.v}, stride{sv.stride * stride}, offset{sv.offset + offset * sv.stride}, strided_size(size)
    {
        assert( 0 <= offset and offset < static_cast<int>(v.size()) and stride > 0 );
    };
    StridedVector(StridedVector const &sv) = default;
    StridedVector(StridedVector &&) = default;
    StridedVector &operator=(StridedVector const &) = default;
    StridedVector &operator=(StridedVector &&) = default;
    ~StridedVector() = default;
};