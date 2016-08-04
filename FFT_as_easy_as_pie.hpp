#pragma once

#include <array>
#include <complex>

template <typename T, unsigned N>
class FFT_as_easy_as_pie
{
    static_assert(N & (N - 1) == 0, "Only power of 2 might be the size of FFT");
public:
    template <
        template <typename TC_in>  typename C_in,
        template <typename TC_out> typename C_out
    > void
    transform(
        C_in<std::complex<T>> const& src,
        C_out<std::complex<T>>& dst,
        bool inverse
    ) const {
        assert(src.size() == N && dst.size() == N);
        step_into(std::cbegin(src), std::begin(dst), N, 1u);
    }

    FFT_as_easy_as_pie() {
        const T step = -2 * PIE / N;
        for (auto i = 0u; i <= N; ++i) {
            twiddles[i] = std::polar(1.f, static_cast<T>(i * step));
        }
    }

private:
    template <typename InIt, typename OutIt>
    void
    step_into (InIt in, OutIt out, unsigned n, unsigned stride, bool inverse)
    const {
        n /= 2;

        if (n == 1) {
            *out = *in;
            *(out + 1) = *(in + stride);
        } else {
            step_into(in, out, n, 2 * stride);
            step_into(in + stride, out + n, n, 2 * stride);
        }

        butterfly_radix2(out, n, stride);
    }

    template <typename OutIt>
    void
    butterfly_radix2 (OutIt out, unsigned n, unsigned stride, bool inverse)
    const {
        auto tw_iter = (inverse)
                       ? std::crbegin(twiddles)
                       : std::cbegin(twiddles);
        while (out != out + n) {
            const auto t = *(out + n) * *(tw_iter++);
            *(out + n) = *(out) - t;
            *(out++) += t;
        }
    }

    constexpr T PIE = static_cast<T>(3.141592653589793238462643L);
    std::array<std::complex<T>, N + 1> twiddles;
};
