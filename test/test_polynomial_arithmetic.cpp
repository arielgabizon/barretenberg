#include <gtest/gtest.h>

#include <barretenberg/curves/bn254/fr.hpp>
#include <barretenberg/polynomials/evaluation_domain.hpp>
#include <barretenberg/polynomials/polynomial_arithmetic.hpp>

using namespace barretenberg;

TEST(polynomials, evaluation_domain)
{
    size_t n = 256;
    evaluation_domain domain = evaluation_domain(n);

    EXPECT_EQ(domain.size, 256UL);
    EXPECT_EQ(domain.log2_size, 8UL);
}

TEST(polynomials, domain_roots)
{
    size_t n = 256;
    evaluation_domain domain = evaluation_domain(n);

    fr::field_t result;
    fr::field_t expected;
    expected = fr::one;
    fr::__pow_small(domain.root, n, result);

    EXPECT_EQ(fr::eq(result, expected), true);
}

TEST(polynomials, fft_with_small_degree)
{
    size_t n = 16;
    fr::field_t fft_transform[n];
    fr::field_t poly[n];

    for (size_t i = 0; i < n; ++i)
    {
        poly[i] = fr::random_element();
        fr::__copy(poly[i], fft_transform[i]);
    }

    evaluation_domain domain = evaluation_domain(n);
    domain.compute_lookup_table();
    polynomial_arithmetic::fft(fft_transform, domain);

    fr::field_t work_root;
    work_root = fr::one;
    fr::field_t expected;
    for (size_t i = 0; i < n; ++i)
    {
        expected = polynomial_arithmetic::evaluate(poly, work_root, n);
        EXPECT_EQ(fr::eq(fft_transform[i], expected), true);
        fr::__mul(work_root, domain.root, work_root);
    }
}

TEST(polynomials, basic_fft)
{
    size_t n = 1 << 14;
    fr::field_t* data = (fr::field_t*)aligned_alloc(32, sizeof(fr::field_t) * n * 2);
    fr::field_t* result = &data[0];
    fr::field_t* expected = &data[n];
    for (size_t i = 0; i < n; ++i)
    {
        result[i] = fr::random_element();
        fr::__copy(result[i], expected[i]);
    }

    evaluation_domain domain = evaluation_domain(n);
    domain.compute_lookup_table();
    polynomial_arithmetic::fft(result, domain);
    polynomial_arithmetic::ifft(result, domain);

    for (size_t i = 0; i < n; ++i)
    {
        EXPECT_EQ(fr::eq(result[i], expected[i]), true);
    }
    aligned_free(data);
}

TEST(polynomials, fft_ifft_consistency)
{
    size_t n = 256;
    fr::field_t result[n];
    fr::field_t expected[n];
    for (size_t i = 0; i < n; ++i)
    {
        result[i] = fr::random_element();
        fr::__copy(result[i], expected[i]);
    }

    evaluation_domain domain = evaluation_domain(n);
    domain.compute_lookup_table();
    polynomial_arithmetic::fft(result, domain);
    polynomial_arithmetic::ifft(result, domain);

    for (size_t i = 0; i < n; ++i)
    {
        EXPECT_EQ(fr::eq(result[i], expected[i]), true);
    }
}

TEST(polynomials, fft_coset_ifft_consistency)
{
    size_t n = 256;
    fr::field_t result[n];
    fr::field_t expected[n];
    for (size_t i = 0; i < n; ++i)
    {
        result[i] = fr::random_element();
        fr::__copy(result[i], expected[i]);
    }

    evaluation_domain domain = evaluation_domain(n);
    domain.compute_lookup_table();
    fr::field_t T0;
    fr::__mul(domain.generator, domain.generator_inverse, T0);
    EXPECT_EQ(fr::eq(T0, fr::one), true);

    polynomial_arithmetic::coset_fft(result, domain);
    polynomial_arithmetic::coset_ifft(result, domain);

    for (size_t i = 0; i < n; ++i)
    {
        EXPECT_EQ(fr::eq(result[i], expected[i]), true);
    }
}

TEST(polynomials, fft_coset_ifft_cross_consistency)
{
    size_t n = 2;
    fr::field_t expected[n];
    fr::field_t poly_a[4 * n];
    fr::field_t poly_b[4 * n];
    fr::field_t poly_c[4 * n];

    for (size_t i = 0; i < n; ++i)
    {
        poly_a[i] = fr::random_element();
        fr::__copy(poly_a[i], poly_b[i]);
        fr::__copy(poly_a[i], poly_c[i]);
        fr::__add(poly_a[i], poly_c[i], expected[i]);
        fr::__add(expected[i], poly_b[i], expected[i]);
    }

    for (size_t i = n; i < 4 * n; ++i)
    {
        poly_a[i] = fr::zero;
        poly_b[i] = fr::zero;
        poly_c[i] = fr::zero;
    }
    evaluation_domain small_domain = evaluation_domain(n);
    evaluation_domain mid_domain = evaluation_domain(2 * n);
    evaluation_domain large_domain = evaluation_domain(4 * n);
    small_domain.compute_lookup_table();
    mid_domain.compute_lookup_table();
    large_domain.compute_lookup_table();
    polynomial_arithmetic::coset_fft(poly_a, small_domain);
    polynomial_arithmetic::coset_fft(poly_b, mid_domain);
    polynomial_arithmetic::coset_fft(poly_c, large_domain);

    for (size_t i = 0; i < n; ++i)
    {
        fr::__add(poly_a[i], poly_c[4 * i], poly_a[i]);
        fr::__add(poly_a[i], poly_b[2 * i], poly_a[i]);
    }

    polynomial_arithmetic::coset_ifft(poly_a, small_domain);

    for (size_t i = 0; i < n; ++i)
    {
        EXPECT_EQ(fr::eq(poly_a[i], expected[i]), true);
    }
}

TEST(polynomials, compute_lagrange_polynomial_fft)
{
    size_t n = 256;
    evaluation_domain small_domain = evaluation_domain(n);
    evaluation_domain mid_domain = evaluation_domain(2 * n);
    small_domain.compute_lookup_table();
    mid_domain.compute_lookup_table();
    fr::field_t l_1_coefficients[2 * n];
    fr::field_t scratch_memory[2 * n + 4];
    for (size_t i = 0; i < 2 * n; ++i)
    {
    l_1_coefficients[i] = fr::zero;
    scratch_memory[i] = fr::zero;
    }
    polynomial_arithmetic::compute_lagrange_polynomial_fft(l_1_coefficients, small_domain, mid_domain);

    polynomial_arithmetic::copy_polynomial(l_1_coefficients, scratch_memory, 2 * n, 2 * n);

    polynomial_arithmetic::coset_ifft(l_1_coefficients, mid_domain);

    fr::field_t z = fr::random_element();
    fr::field_t shifted_z;
    fr::__mul(z, small_domain.root, shifted_z);
    fr::__mul(shifted_z, small_domain.root, shifted_z);

    fr::field_t eval;
    fr::field_t shifted_eval;

    eval = polynomial_arithmetic::evaluate(l_1_coefficients, shifted_z, small_domain.size);
    polynomial_arithmetic::fft(l_1_coefficients, small_domain);

    fr::__copy(scratch_memory[0], scratch_memory[2 * n]);
    fr::__copy(scratch_memory[1], scratch_memory[2 * n + 1]);
    fr::__copy(scratch_memory[2], scratch_memory[2 * n + 2]);
    fr::__copy(scratch_memory[3], scratch_memory[2 * n + 3]);
    fr::field_t* l_n_minus_one_coefficients = &scratch_memory[4];
    polynomial_arithmetic::coset_ifft(l_n_minus_one_coefficients, mid_domain);

    shifted_eval = polynomial_arithmetic::evaluate(l_n_minus_one_coefficients, z, small_domain.size);
    EXPECT_EQ(fr::eq(eval, shifted_eval), true);

    polynomial_arithmetic::fft(l_n_minus_one_coefficients, small_domain);

    EXPECT_EQ(fr::eq(l_1_coefficients[0], fr::one), true);

    for (size_t i = 1; i < n; ++i)
    {
        EXPECT_EQ(fr::eq(l_1_coefficients[i], fr::zero), true);
    }

    EXPECT_EQ(fr::eq(l_n_minus_one_coefficients[n - 2], fr::one), true);

    for (size_t i = 0; i < n; ++i)
    {
        if (i == (n - 2))
        {
            continue;
        }
        EXPECT_EQ(fr::eq(l_n_minus_one_coefficients[i], fr::zero), true);
    }
}

TEST(polynomials, divide_by_pseudo_vanishing_polynomial)
{
    size_t n = 256;
    fr::field_t a[4 * n];
    fr::field_t b[4 * n];
    fr::field_t c[4 * n];

    fr::field_t T0;
    for (size_t i = 0; i < n; ++i)
    {
        a[i] = fr::random_element();
        b[i] = fr::random_element();
        fr::__mul(a[i], b[i], c[i]);
        fr::__neg(c[i], c[i]);
        fr::__mul(a[i], b[i], T0);
        fr::__add(T0, c[i], T0);
    }
    for (size_t i = n; i < 4 * n; ++i)
    {
        a[i] = fr::zero;
        b[i] = fr::zero;
        c[i] = fr::zero;
    }

    // make the final evaluation not vanish
    // fr::one(c[n-1]);
    evaluation_domain small_domain = evaluation_domain(n);
    evaluation_domain mid_domain = evaluation_domain(4 * n);
    small_domain.compute_lookup_table();
    mid_domain.compute_lookup_table();

    polynomial_arithmetic::ifft(a, small_domain);
    polynomial_arithmetic::ifft(b, small_domain);
    polynomial_arithmetic::ifft(c, small_domain);

    polynomial_arithmetic::coset_fft(a, mid_domain);
    polynomial_arithmetic::coset_fft(b, mid_domain);
    polynomial_arithmetic::coset_fft(c, mid_domain);

    fr::field_t result[mid_domain.size];
    for (size_t i = 0; i < mid_domain.size; ++i)
    {
        fr::__mul(a[i], b[i], result[i]);
        fr::__add(result[i], c[i], result[i]);
    }

    polynomial_arithmetic::divide_by_pseudo_vanishing_polynomial(&result[0], small_domain, mid_domain);

    polynomial_arithmetic::coset_ifft(result, mid_domain);

    for (size_t i = n + 1; i < mid_domain.size; ++i)
    {

        EXPECT_EQ(fr::eq(result[i], fr::zero), true);
    }
}

TEST(polynomials, compute_kate_opening_coefficients)
{
    // generate random polynomial F(X) = coeffs
    size_t n = 256;
    fr::field_t coeffs[2 * n];
    fr::field_t W[2 * n];
    for (size_t i = 0; i < n; ++i)
    {
        coeffs[i] = fr::random_element();
        coeffs[i + n] = fr::zero;
    }
    polynomial_arithmetic::copy_polynomial(coeffs, W, 2 * n, 2 * n);

    // generate random evaluation point z
    fr::field_t z = fr::random_element();

    // compute opening polynomial W(X), and evaluation f = F(z)
    fr::field_t f = polynomial_arithmetic::compute_kate_opening_coefficients(W, W, z, n);

    // validate that W(X)(X - z) = F(X) - F(z)
    // compute (X - z) in coefficient form
    fr::field_t multiplicand[2 * n];
    fr::__neg(z, multiplicand[0]);
    multiplicand[1] = fr::one;
    for (size_t i = 2; i < 2 * n; ++i)
    {
        multiplicand[i] = fr::zero;
    }

    // set F(X) = F(X) - F(z)
    fr::__sub(coeffs[0], f, coeffs[0]);

    // compute fft of polynomials
    evaluation_domain domain = evaluation_domain(2 * n);
    domain.compute_lookup_table();
    polynomial_arithmetic::coset_fft(coeffs, domain);
    polynomial_arithmetic::coset_fft(W, domain);
    polynomial_arithmetic::coset_fft(multiplicand, domain);

    // validate that, at each evaluation, W(X)(X - z) = F(X) - F(z)
    fr::field_t result;
    for (size_t i = 0; i < domain.size; ++i)
    {
        fr::__mul(W[i], multiplicand[i], result);
        EXPECT_EQ(fr::eq(result, coeffs[i]), true);
    }
}

TEST(polynomials, get_lagrange_evaluations)
{
    size_t n = 16;

    evaluation_domain domain = evaluation_domain(n);
    domain.compute_lookup_table();
    fr::field_t z = fr::random_element();

    polynomial_arithmetic::lagrange_evaluations evals = polynomial_arithmetic::get_lagrange_evaluations(z, domain);

    fr::field_t vanishing_poly[2 * n];
    fr::field_t l_1_poly[n];
    fr::field_t l_n_minus_1_poly[n];

    for (size_t i = 0; i < n; ++i)
    {
        l_1_poly[i] = fr::zero;
        l_n_minus_1_poly[i] = fr::zero;
        vanishing_poly[i] = fr::zero;
    }
    l_1_poly[0] = fr::one;
    l_n_minus_1_poly[n - 2] = fr::one;

    fr::field_t n_mont{ { n, 0, 0, 0 } };
    fr::__to_montgomery_form(n_mont, n_mont);
    fr::__mul(n_mont, domain.root, vanishing_poly[n - 1]);

    polynomial_arithmetic::ifft(l_1_poly, domain);
    polynomial_arithmetic::ifft(l_n_minus_1_poly, domain);
    polynomial_arithmetic::ifft(vanishing_poly, domain);

    fr::field_t l_1_expected;
    fr::field_t l_n_minus_1_expected;
    fr::field_t vanishing_poly_expected;
    l_1_expected = polynomial_arithmetic::evaluate(l_1_poly, z, n);
    l_n_minus_1_expected = polynomial_arithmetic::evaluate(l_n_minus_1_poly, z, n);
    vanishing_poly_expected = polynomial_arithmetic::evaluate(vanishing_poly, z, n);
    EXPECT_EQ(fr::eq(evals.l_1, l_1_expected), true);
    EXPECT_EQ(fr::eq(evals.l_n_minus_1, l_n_minus_1_expected), true);
    EXPECT_EQ(fr::eq(evals.vanishing_poly, vanishing_poly_expected), true);
}