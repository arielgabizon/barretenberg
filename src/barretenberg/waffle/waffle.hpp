#pragma once

#include "stddef.h"
#include "stdint.h"
#include "malloc.h"

#include "../groups/pairing.hpp"
#include "../groups/scalar_multiplication.hpp"
#include "../polynomials/fft.hpp"

namespace waffle
{

struct plonk_srs
{
    g1::affine_element *monomials;
    g2::affine_element t2;
    size_t degree;
};

// contains the state of a PLONK proof, including witness values, instance values
// and Kate polynomial commitments
struct circuit_state
{
    // Kate polynomial commitments required for a proof of knowledge
    g1::affine_element W_L;
    g1::affine_element W_R;
    g1::affine_element W_O;
    g1::affine_element Z_1;
    g1::affine_element Z_2;
    g1::affine_element T;
    g1::affine_element PI_Z;
    g1::affine_element PI_Z_OMEGA;

    // random challenges
    fr::field_t gamma;
    fr::field_t beta;
    fr::field_t alpha;
    fr::field_t alpha_squared;
    fr::field_t alpha_cubed;
    fr::field_t alpha_quad;
    fr::field_t alpha_quintic;
    fr::field_t z;

    // pointers to witness vectors. Originally these are in Lagrange-base form,
    // during the course of proof construction, are replaced by their coefficient form
    fr::field_t *w_l;
    fr::field_t *w_r;
    fr::field_t *w_o;
    fr::field_t *z_1;
    fr::field_t *z_2;
    fr::field_t *t;

    // pointers to instance vectors. Originally in Lagrange-base form,
    // will be converted into coefficient form
    fr::field_t *q_c;
    fr::field_t *q_m;
    fr::field_t *q_l;
    fr::field_t *q_r;
    fr::field_t *q_o;
    fr::field_t *sigma_1;
    fr::field_t *sigma_2;
    fr::field_t *sigma_3;
    fr::field_t *s_id;

    fr::field_t *product_1;
    fr::field_t *product_2;
    fr::field_t *product_3;
    size_t n;
};

// Stores pointers to various polynomials required during proof construction
// Note: these pointers can overlap! We want to efficiently use available memory,
// and only a handful of these polynomials are required at any one time
struct fft_pointers
{
    fr::field_t* w_l_poly;
    fr::field_t* w_r_poly;
    fr::field_t* w_o_poly;
    fr::field_t* z_1_poly;
    fr::field_t* z_2_poly;

    fr::field_t* identity_poly;
    fr::field_t* gate_poly_mid;
    fr::field_t* gate_poly_long;

    fr::field_t* q_c_poly;
    fr::field_t* q_r_poly;
    fr::field_t* q_l_poly;
    fr::field_t* scratch_memory;
};

inline void init_srs(size_t size, plonk_srs& res)
{
    // res->monomials = (g1::affine_element *)aligned_alloc(32, sizeof(g1::affine_element) * 4 * size);
    printf("pre scalars\n");
    fr::field_t *scalars = (fr::field_t *)aligned_alloc(32, sizeof(fr::field_t) * size);
    printf("post scalars\n");
    // TODO: load this from disk. Dummy SRS for now
    fr::field_t x;
    fr::random_element(x);

    g1::affine_element X = g1::affine_one();

    printf("copying\n");
    g1::copy(&X, &res.monomials[0]);
    fr::one(scalars[0]);
    printf("about to do more copying\n");
    for (size_t i = 1; i < size; ++i)
    {
        g1::copy(&X, &res.monomials[i]);
        fr::mul(scalars[i - 1], x, scalars[i]);
    }
    printf("g2 group exponentiation\n");
    g2::affine_element g2_input = g2::affine_one();
    g2_input = g2::group_exponentiation(g2_input, x);
    // res.t2 = g2::affine_one();
    // res.t2 = g2::group_exponentiation(res.t2, x);
    printf("about to copy into res.t2\n");
    g2::copy_affine(g2_input, res.t2);
    printf("copied\n");
    res.degree = size;
    printf("post g2 group exponentiation\n");
    printf("pre point table\n");
    scalar_multiplication::generate_pippenger_point_table(res.monomials, res.monomials, size);
    printf("post point table\n");
    free(scalars);
    // blah blah blah
}


inline void compute_wire_coefficients(circuit_state &state, polynomials::evaluation_domain &domain, fft_pointers& ffts)
{
    const size_t n = state.n;

    polynomials::copy_polynomial(state.w_l, state.z_1, n, n);
    polynomials::copy_polynomial(state.w_r, state.z_2, n, n);
    polynomials::copy_polynomial(state.w_o, state.t, n, n);

    polynomials::ifft(state.w_l, domain.short_root_inverse, n);
    polynomials::ifft(state.w_r, domain.short_root_inverse, n);
    polynomials::ifft(state.w_o, domain.short_root_inverse, n);

    // compute [w_l], [w_r], [w_o]
    // scalar multiplication algorithm modifies scalars, which we want to preserve, so copy first
    polynomials::copy_polynomial(state.w_l, &ffts.w_l_poly[0], n, n);
    polynomials::copy_polynomial(state.w_r, &ffts.w_r_poly[0], n, n);
    polynomials::copy_polynomial(state.w_o, &ffts.w_o_poly[0], n, n);
}


inline void compute_z_coefficients(circuit_state& state, polynomials::evaluation_domain &domain, fft_pointers& ffts)
{
    const size_t n = state.n;
    // compute Z1, Z2
    fr::field_t T0;
    fr::field_t beta_n = { .data = { n, 0, 0, 0 } };
    fr::to_montgomery_form(beta_n, beta_n);
    fr::mul(beta_n, state.beta, beta_n);
    fr::field_t beta_n_2;
    fr::add(beta_n, beta_n, beta_n_2);

    // TODO: multithread this part!
    fr::field_t beta_identity = { .data = { 0, 0, 0, 0 } } ;
    // for the sigma permutation, as we compute each product term, store the intermediates in `product_1/2/3`.
    for (size_t i = 0; i < n; ++i)
    {
        fr::mul(state.sigma_1[i], state.beta, state.product_1[i]);
        fr::add(state.product_1[i], state.gamma, state.product_1[i]);
        fr::add(state.product_1[i], state.w_l[i], state.product_1[i]);

        fr::mul(state.sigma_2[i], state.beta, state.product_2[i]);
        fr::add(state.product_2[i], state.gamma, state.product_2[i]);
        fr::add(state.product_2[i], state.w_r[i], state.product_2[i]);

        fr::mul(state.sigma_3[i], state.beta, state.product_3[i]);
        fr::add(state.product_3[i], state.gamma, state.product_3[i]);
        fr::add(state.product_3[i], state.w_o[i], state.product_3[i]);

        fr::mul(state.product_1[i], state.product_2[i], state.z_2[i+1]);
        fr::mul(state.z_2[i+1], state.product_3[i], state.z_2[i+1]);
    
        fr::add(beta_identity, state.beta, beta_identity);
    
        fr::add(beta_identity, state.gamma, state.z_1[i+1]);
        fr::add(state.z_1[i+1], state.w_l[i], state.z_1[i+1]);

        fr::add(beta_identity, state.gamma, T0);
        fr::add(beta_n, T0, T0);
        fr::add(T0, state.w_r[i], T0);

        fr::mul(state.z_1[i+1], T0, state.z_1[i+1]);

        fr::add(beta_identity, state.gamma, T0);
        fr::add(beta_n_2, T0, T0);
        fr::add(T0, state.w_o[i], T0);

        fr::mul(state.z_1[i+1], T0, state.z_1[i+1]);
    }
    fr::one(state.z_1[0]);
    fr::one(state.z_2[0]);

    for (size_t i = 1; i < n; ++i)
    {
        fr::mul(state.z_1[i], state.z_1[i-1], state.z_1[i]);
        fr::mul(state.z_2[i], state.z_2[i-1], state.z_2[i]);
    }

    polynomials::ifft(state.z_1, domain.short_root_inverse, n);
    polynomials::ifft(state.z_2, domain.short_root_inverse, n);

    // compute [z_1], [z_2]
    // scalar multiplication algorithm modifies scalars, which we want to preserve, so copy first
    polynomials::copy_polynomial(state.z_1, &ffts.z_1_poly[0], n, n);
    polynomials::copy_polynomial(state.z_2, &ffts.z_2_poly[0], n, n);
}

inline void compute_wire_commitments(circuit_state &state, fr::field_t* wire_coefficients, plonk_srs &srs)
{
    size_t n = state.n;
    scalar_multiplication::multiplication_state mul_state[3];
    mul_state[0].num_elements = n;
    mul_state[0].scalars = &wire_coefficients[0];
    mul_state[0].points = srs.monomials;
    mul_state[1].num_elements = n;
    mul_state[1].scalars = &wire_coefficients[n];
    mul_state[1].points = srs.monomials;
    mul_state[2].num_elements = n;
    mul_state[2].scalars = &wire_coefficients[n];
    mul_state[2].points = srs.monomials;

    // scalar_multiplication::batched_scalar_multiplications(mul_state, 3);

    // TODO: make a method for normal-to-affine :/
    fq::copy(mul_state[0].output.x, state.W_L.x);
    fq::copy(mul_state[1].output.x, state.W_R.x);
    fq::copy(mul_state[2].output.x, state.W_O.x);
    fq::copy(mul_state[0].output.y, state.W_L.y);
    fq::copy(mul_state[1].output.y, state.W_R.y);
    fq::copy(mul_state[2].output.y, state.W_O.y);

    // compute beta, gamma
    // TODO: use keccak256
    fr::random_element(state.beta);
    fr::random_element(state.gamma);
}

inline void compute_z_commitments(circuit_state& state, fr::field_t* z_coefficients, plonk_srs& srs)
{
    size_t n = state.n;
    scalar_multiplication::multiplication_state mul_state[3];

    mul_state[0].num_elements = n;
    mul_state[0].scalars = &z_coefficients[0];
    mul_state[0].points = srs.monomials;
    mul_state[1].num_elements = n;
    mul_state[1].scalars = &z_coefficients[n];
    mul_state[1].points = srs.monomials;

    scalar_multiplication::batched_scalar_multiplications(mul_state, 2);

    // compute alpha
    // TODO: use keccak256, this is just for testing
    // precompute some powers of alpha for later on
    fr::random_element(state.alpha);
    fr::mul(state.alpha, state.alpha, state.alpha_squared);
    fr::mul(state.alpha_squared, state.alpha, state.alpha_cubed);
    fr::mul(state.alpha_cubed, state.alpha, state.alpha_quad);
    fr::mul(state.alpha_quad, state.alpha, state.alpha_quintic);
}

inline void compute_identity_grand_product_coefficients(circuit_state& state, polynomials::evaluation_domain &domain, fft_pointers& ffts)
{
    size_t n = state.n;
    // 4 (4n) 'registers'
    // fr::field_t* w_l_poly = &scratch_space[0];
    // fr::field_t* w_r_poly = &scratch_space[4 * n];
    // fr::field_t* w_o_poly = &scratch_space[8 * n];
    // fr::field_t* identity_poly = &scratch_space[12 * n];

    // compute 4n fft transforms for w_l, w_r, w_o and the identity permutation
    polynomials::copy_polynomial(state.w_l, ffts.w_l_poly, n, 4 * n);
    polynomials::copy_polynomial(state.w_r, ffts.w_r_poly, n, 4 * n);
    polynomials::copy_polynomial(state.w_o, ffts.w_o_poly, n, 4 * n);

    // TODO: optimize this!
    fr::one(state.s_id[0]);
    fr::field_t one;
    fr::one(one);
    for (size_t i = 1; i < n; ++i)
    {
        fr::add(state.s_id[i-1], one, state.s_id[i]);
    }
    polynomials::ifft(state.s_id, domain.short_root_inverse, domain.short_domain);
    polynomials::copy_polynomial(state.s_id, ffts.identity_poly, domain.short_domain, domain.long_domain);
    
    polynomials::fft_with_coset_and_constant(ffts.identity_poly, domain.long_root, domain.generator, state.beta, domain.long_domain);
    polynomials::fft_with_coset(ffts.w_l_poly, domain.long_root, domain.generator, domain.long_domain);
    polynomials::fft_with_coset(ffts.w_r_poly, domain.long_root, domain.generator, domain.long_domain);
    polynomials::fft_with_coset(ffts.w_o_poly, domain.long_root, domain.generator, domain.long_domain);

    // compute partial identity grand product
    fr::field_t T0;
    fr::field_t T1;
    fr::field_t T2;
    fr::field_t T3;
    fr::field_t beta_n = { .data = { n, 0, 0, 0 } };
    fr::to_montgomery_form(beta_n, beta_n);
    fr::mul(beta_n, state.beta, beta_n);
    fr::field_t beta_n_2;
    fr::add(beta_n, beta_n, beta_n_2);
    for (size_t i = 0; i < domain.long_domain; ++i)
    {
        fr::add(ffts.identity_poly[i], state.gamma, T0);
        fr::add(T0, ffts.w_l_poly[i], T1);
        fr::add(T0, ffts.w_r_poly[i], T2);
        fr::add(T0, ffts.w_o_poly[i], T3);
        fr::add(T2, beta_n, T2);
        fr::add(T3, beta_n_2, T3);
        fr::mul(T1, T2, ffts.identity_poly[i]);
        fr::mul(ffts.identity_poly[i], T3, ffts.identity_poly[i]);
    }
}

inline void compute_arithmetisation_coefficients(circuit_state& state, polynomials::evaluation_domain &domain, fft_pointers& ffts)
{
    size_t n = state.n;
    // compute q.o * w.o
    polynomials::ifft(state.q_o, domain.short_root_inverse, domain.short_domain);
    polynomials::copy_polynomial(state.q_o, ffts.gate_poly_mid, n, 2 * n);
    polynomials::fft_with_coset(ffts.gate_poly_mid, domain.mid_root, domain.generator, domain.mid_domain);

    // the fft transform on q.o is half that of w.o - access every other index of w.o
    for (size_t i = 0; i < domain.mid_domain; ++i)
    {
        fr::mul(ffts.w_o_poly[i * 2], ffts.gate_poly_mid[i], ffts.gate_poly_mid[i]);
    }

    // great! we've freed up w_o now
    // we can use that scratch space to compute q.l*w.l and q.r*w.r and q.c

    // compute q_c_poly and add into accumulator
    polynomials::ifft(state.q_c, domain.short_root_inverse, domain.short_domain);
    polynomials::copy_polynomial(state.q_c, ffts.q_c_poly, n, 2 * n);
    polynomials::fft_with_coset(ffts.q_c_poly, domain.mid_root, domain.generator, domain.mid_domain);
    for (size_t i = 0; i < domain.mid_domain; ++i)
    {
        fr::add(ffts.gate_poly_mid[i], ffts.q_c_poly[i], ffts.gate_poly_mid[i]);
    }
    // compute q.r * w.r
    polynomials::ifft(state.q_r, domain.short_root_inverse, domain.short_domain);
    polynomials::ifft(state.q_l, domain.short_root_inverse, domain.short_domain);

    polynomials::copy_polynomial(state.q_r, ffts.q_r_poly, n, 2 * n);
    polynomials::copy_polynomial(state.q_l, ffts.q_l_poly, n, 2 * n);
    polynomials::fft_with_coset(ffts.q_r_poly, domain.mid_root, domain.generator, domain.mid_domain);
    polynomials::fft_with_coset(ffts.q_l_poly, domain.mid_root, domain.generator, domain.mid_domain);

    // the fft transform on q.o is half that of w.o - access every other index of w.o
    fr::field_t T0;
    fr::field_t T1;
    for (size_t i = 0; i < domain.mid_domain; ++i)
    {
        fr::mul(ffts.w_r_poly[i * 2], ffts.q_r_poly[i], T0);
        fr::mul(ffts.w_l_poly[i * 2], ffts.q_l_poly[i], T1);
        fr::add(T1, T0, T1);
        fr::add(ffts.gate_poly_mid[i], T1, ffts.gate_poly_mid[i]);
        fr::mul(ffts.gate_poly_mid[i], state.alpha, ffts.gate_poly_mid[i]);
    }

    // The next step is to compute q_m.w_l.w_r - we need a 4n fft for this
    // requisition the memory that w_o was using
    polynomials::ifft(state.q_m, domain.short_root_inverse, domain.short_domain);
    polynomials::copy_polynomial(state.q_m, ffts.gate_poly_long, n, 4 * n);
    polynomials::fft_with_coset_and_constant(ffts.gate_poly_long, domain.long_root, domain.generator, state.alpha, domain.long_domain);

    for (size_t i = 0; i < domain.long_domain; ++i)
    {
        fr::mul(ffts.w_l_poly[i], ffts.w_r_poly[i], T0);
        fr::mul(ffts.gate_poly_long[i], T0, ffts.gate_poly_long[i]);
    }
}

inline void concatenate_arithmetic_and_identity_coefficients(circuit_state& state, polynomials::evaluation_domain &domain, fft_pointers& ffts)
{
    // we've now freed up use of w_l and w_r, but w_o is occupied with q_m.q_l.q_r

    // we can free up this memory by combining the grand-product term with q_m.q_l.q_r
    size_t n = state.n;

    // We assume that state.z_1 has had an inverse fourier transform performed on it, and is in coefficient form
    polynomials::copy_polynomial(state.z_1, ffts.z_1_poly, n, 4 * n);

    // when we transform z_1 into point-evaluation form, scale up by `alpha_squared` - saves us a mul later on
    polynomials::fft_with_coset_and_constant(ffts.z_1_poly, domain.long_root, domain.generator, state.alpha_squared, domain.long_domain);

    // we have the 4n component of the arithmetistion polynomial in `gate_poly_long`
    // and the identity grand product polynomial in `identity_poly` = I(X)
    // to finish the identity grand product, we need I(X).Z_1(X) - Z_1(X \omega^{-1})
    // the FFT of Z_1(X \omega^{-1}) is the same as Z_1(X), just shifted by 4 indices

    // instead of mapping indices, we just graft the first 4 array elements of Z_1, onto the end of Z_1,
    // and create a pointer to the shifted poly
    for (size_t i = 0; i < 4; ++i)
    {
        fr::copy(ffts.z_1_poly[i], ffts.z_1_poly[n+i]);
    }
    fr::field_t* shifted_z_1_poly = ffts.z_1_poly + 4; // tadaa
    for (size_t i = 0; i < domain.long_domain; ++i)
    {
        // multiply identity evaluation by normal z_1 evaluation
        fr::mul(ffts.identity_poly[i], ffts.z_1_poly[i], ffts.identity_poly[i]);
    
        // subtract the shifted evaluation from result
        fr::sub(ffts.identity_poly[i], shifted_z_1_poly[i], ffts.identity_poly[i]);

        // and add the gate poly evaluation into the result
        fr::add(ffts.identity_poly[i], ffts.gate_poly_long[i], ffts.identity_poly[i]);
    }
}

inline void compute_permutation_grand_product_coefficients(circuit_state& state, polynomials::evaluation_domain& domain, fr::field_t* gate_poly_mid, fr::field_t* quotient_poly, fr::field_t* z_1_poly, fr::field_t* scratch_memory)
{
    // The final steps are:
    // 1: Compute the permutation grand product
    // 2: Compute permutation check coefficients
    size_t n = state.n;
    fr::field_t* sigma_1_poly = &scratch_memory[0];
    fr::field_t* sigma_2_poly = &scratch_memory[domain.long_domain];
    fr::field_t* l_1_poly = &scratch_memory[2 * domain.long_domain];
    // // free memory: w_r, w_o
    // fr::field_t* l_1_poly = w_r_poly;
    // fr::field_t* l_n_poly = w_r_poly + domain.mid_domain;

    polynomials::compute_lagrange_polynomial_fft(l_1_poly, domain, &scratch_memory[2 * domain.long_domain + domain.mid_domain + 4]);

    fr::copy(l_1_poly[0], l_1_poly[domain.mid_domain]);
    fr::copy(l_1_poly[1], l_1_poly[domain.mid_domain + 1]);
    fr::copy(l_1_poly[2], l_1_poly[domain.mid_domain + 2]);
    fr::copy(l_1_poly[3], l_1_poly[domain.mid_domain + 3]);
    fr::copy(z_1_poly[0], z_1_poly[domain.long_domain]);
    fr::copy(z_1_poly[1], z_1_poly[domain.long_domain + 1]);
    fr::copy(z_1_poly[2], z_1_poly[domain.long_domain + 2]);
    fr::copy(z_1_poly[3], z_1_poly[domain.long_domain + 3]);
    fr::field_t* l_n_minus_1_poly = &l_1_poly[4];
    fr::field_t* shifted_z_1_poly = &z_1_poly[4];
    for (size_t i = 0; i < domain.mid_domain; ++i)
    {
        fr::mul(l_1_poly[i], z_1_poly[i * 2], l_1_poly[i]);
        fr::mul(l_n_minus_1_poly[i], shifted_z_1_poly[i * 2], l_n_minus_1_poly[i]);
        fr::mul(l_1_poly[i], state.alpha_squared, l_1_poly[i]);
        fr::mul(l_n_minus_1_poly[i], state.alpha_cubed, l_n_minus_1_poly[i]);
    }

    // we've computed the fft evaluation of z_1(X) * l_1(X)
    // and can use a shift to get to z_1(X.w^{-1}) * l_{n-1}(X)
    // z_1 is now free
    // => w_o, w_l are now free

    // fr::field_t* sigma_1_poly = w_l_poly;
    // fr::field_t* sigma_2_poly = w_o_poly;
    polynomials::copy_polynomial(state.sigma_1, sigma_1_poly, n, 4 * n);
    polynomials::copy_polynomial(state.sigma_2, sigma_2_poly, n, 4 * n);

    polynomials::fft_with_coset(sigma_1_poly, domain.long_root, domain.generator, domain.long_domain);
    polynomials::fft_with_coset(sigma_2_poly, domain.long_root, domain.generator, domain.long_domain);

    for (size_t i = 0; i < domain.long_domain; ++i)
    {
        fr::mul(sigma_1_poly[i], sigma_2_poly[i], sigma_1_poly[i]);
    }

    fr::field_t* sigma_3_poly = z_1_poly;
    polynomials::copy_polynomial(state.sigma_3, sigma_3_poly, n, 4 * n);
    polynomials::fft_with_coset(sigma_3_poly, domain.long_root, domain.generator, domain.long_domain);

    for (size_t i = 0; i < domain.long_domain; ++i)
    {
        fr::mul(sigma_1_poly[i], sigma_3_poly[i], sigma_1_poly[i]);
    }

    fr::field_t* z_2_poly = z_1_poly;
    polynomials::copy_polynomial(state.z_2, z_2_poly, n, 4 * n);
    polynomials::fft_with_coset_and_constant(z_2_poly, domain.long_root, domain.generator, state.alpha_cubed, domain.long_domain);
    fr::copy(z_2_poly[0], z_2_poly[domain.long_domain]);
    fr::copy(z_2_poly[1], z_2_poly[domain.long_domain + 1]);
    fr::copy(z_2_poly[2], z_2_poly[domain.long_domain + 2]);
    fr::copy(z_2_poly[3], z_2_poly[domain.long_domain + 3]);
    fr::field_t* shifted_z_2_poly = &z_2_poly[4];
    for (size_t i = 0; i < domain.long_domain; ++i)
    {
        size_t shifted_index = i < (domain.long_domain - 4) ? i + 4 : 4 - (domain.long_domain - i);
        fr::mul(sigma_1_poly[i], z_2_poly[i], sigma_1_poly[i]);
        fr::sub(sigma_1_poly[i], z_2_poly[shifted_index], sigma_1_poly[i]);

        // combine product term into quotient poly
        fr::mul(quotient_poly[i], sigma_1_poly[i], quotient_poly[i]);
    }

    fr::field_t T0;
    // accumulate degree-2n terms into gate_poly_mid
    for (size_t i = 0; i < domain.mid_domain; ++i)
    {
        fr::mul(z_2_poly[i], state.alpha, T0);
        fr::sub(l_1_poly[i], T0, l_1_poly[i]);

        fr::mul(shifted_z_2_poly[i], state.alpha_squared, T0);
        fr::sub(l_n_minus_1_poly[i], T0, l_n_minus_1_poly[i]);

        fr::add(gate_poly_mid[i], l_1_poly[i], gate_poly_mid[i]);
        fr::add(gate_poly_mid[i], l_n_minus_1_poly[i], gate_poly_mid[i]);
    }

    polynomials::ifft_with_coset(gate_poly_mid, domain.mid_root_inverse, domain.generator_inverse, domain.mid_domain);
    memset((void *)(gate_poly_mid + domain.mid_domain), 0, (domain.mid_domain) * sizeof(fr::field_t));

    polynomials::fft_with_coset(gate_poly_mid, domain.long_root, domain.generator, domain.long_domain);

    for (size_t i = 0; i < domain.long_domain; ++i)
    {
        fr::add(quotient_poly[i], gate_poly_mid[i], quotient_poly[i]);
    }

}

inline void compute_quotient_polynomial(circuit_state& state, polynomials::evaluation_domain &domain, const plonk_srs &)
{
    // q_c
    // w_o.q_o
    // w_l.q_l
    // w_r.q_r
    // (w_l + b.s_id + g) = sigma_1
    // (w_r + b.s_id + (g + n.b)) = sigma_2
    // (w_o + b.s_id + (g + 2n.b)) = sigma_3
    // (w_l + b.s_s1 + g)
    // (w_r + b.s_s2 + g)
    // (w_o + b.s_s3 + g)

    // By this point, we have already computed the z_1, z_2 coefficients, which are now represented in Lagrange-base form
    // Reserve 20n bytes of memory for our fft transforms
    size_t n = state.n;
    fr::field_t* scratch_space = (fr::field_t*)aligned_alloc(32, sizeof(fr::field_t) * n * 20);
    // 4 (4n) 'registers'
    fr::field_t* w_l_poly = &scratch_space[0];
    fr::field_t* w_r_poly = &scratch_space[4 * n];
    fr::field_t* w_o_poly = &scratch_space[8 * n];
    fr::field_t* identity_poly = &scratch_space[12 * n];
    // 2 (2n) 'registers'
    fr::field_t* gate_poly_mid = &scratch_space[16 * n];
    // fr::field_t* mid_work_register = &scratch_space[18 * n];

    // compute 4n fft transforms for w_l, w_r, w_o and the identity permutation
    polynomials::copy_polynomial(state.w_l, w_l_poly, n, 4 * n);
    polynomials::copy_polynomial(state.w_r, w_r_poly, n, 4 * n);
    polynomials::copy_polynomial(state.w_o, w_o_poly, n, 4 * n);

    polynomials::fft_with_coset(w_l_poly, domain.long_root, domain.generator, domain.long_domain);
    polynomials::fft_with_coset(w_r_poly, domain.long_root, domain.generator, domain.long_domain);
    polynomials::fft_with_coset(w_o_poly, domain.long_root, domain.generator, domain.long_domain);

    // TODO: optimize this out!
    fr::one(identity_poly[0]);
    for (size_t i = 1; i < n; ++i)
    {
        fr::add(identity_poly[i], identity_poly[i-1], identity_poly[i]);
    }
    memset((void *)(identity_poly + domain.short_domain), 0, (domain.long_domain - domain.short_domain) * sizeof(fr::field_t));
    polynomials::fft_with_coset_and_constant(identity_poly, domain.long_root, domain.generator, state.beta, domain.long_domain);
    
    // compute identity grand product
    fr::field_t T0;
    fr::field_t T1;
    fr::field_t T2;
    fr::field_t T3;
    fr::field_t beta_n = { .data = { n, 0, 0, 0 } };
    fr::to_montgomery_form(beta_n, beta_n);
    fr::mul(beta_n, state.beta, beta_n);
    fr::field_t beta_n_2;
    fr::add(beta_n, beta_n, beta_n_2);
    for (size_t i = 0; i < domain.long_domain; ++i)
    {
        fr::add(identity_poly[i], state.gamma, T0);
        fr::add(T0, w_l_poly[i], T1);
        fr::add(T0, w_r_poly[i], T2);
        fr::add(T0, w_o_poly[i], T3);
        fr::add(T2, beta_n, T2);
        fr::add(T3, beta_n_2, T3);
        fr::mul(T1, T2, identity_poly[i]);
        fr::mul(identity_poly[i], T3, identity_poly[i]);
    }

    // compute q.o * w.o
    polynomials::ifft(state.q_o, domain.short_root_inverse, domain.short_domain);
    polynomials::copy_polynomial(state.q_o, gate_poly_mid, n, 2 * n);
    polynomials::fft_with_coset(gate_poly_mid, domain.mid_root, domain.generator, domain.mid_domain);

    // the fft transform on q.o is half that of w.o - access every other index of w.o
    for (size_t i = 0; i < domain.mid_domain; ++i)
    {
        fr::mul(w_o_poly[i * 2], gate_poly_mid[i], gate_poly_mid[i]);
    }

    // great! we've freed up w_o now
    // we can use that scratch space to compute q.l*w.l and q.r*w.r

    // compute q.r * w.r
    fr::field_t* q_r_poly = w_o_poly; // hehe
    fr::field_t* q_l_poly = w_o_poly + domain.mid_domain;
    polynomials::ifft(state.q_r, domain.short_root_inverse, domain.short_domain);
    polynomials::ifft(state.q_l, domain.short_root_inverse, domain.short_domain);

    polynomials::copy_polynomial(state.q_r, q_r_poly, n, 2 * n);
    polynomials::copy_polynomial(state.q_l, q_l_poly, n, 2 * n);
    polynomials::fft_with_coset(q_r_poly, domain.mid_root, domain.generator, domain.mid_domain);
    polynomials::fft_with_coset(q_l_poly, domain.mid_root, domain.generator, domain.mid_domain);

    // the fft transform on q.o is half that of w.o - access every other index of w.o
    for (size_t i = 0; i < domain.mid_domain; ++i)
    {
        fr::mul(w_r_poly[i * 2], q_r_poly[i], T0);
        fr::mul(w_l_poly[i * 2], q_l_poly[i], T1);
        fr::add(T1, T0, T1);
        fr::add(gate_poly_mid[i], T1, gate_poly_mid[i]);
    }

    // The next step is to compute q_m.w_l.w_r - we need a 4n fft for this
    // requisition the memory that w_o was using
    fr::field_t* q_m_poly = w_o_poly;
    polynomials::ifft(state.q_m, domain.short_root_inverse, domain.short_domain);
    polynomials::copy_polynomial(state.q_m, q_m_poly, n, 4 * n);
    polynomials::fft_with_coset_and_constant(q_m_poly, domain.long_root, domain.generator, state.alpha, domain.long_domain);

    for (size_t i = 0; i < domain.long_domain; ++i)
    {
        fr::mul(w_l_poly[i], w_r_poly[i], T0);
        fr::mul(q_m_poly[i], T0, q_m_poly[i]);
    }

    // we've now freed up use of w_l and w_r, but w_o is occupied with q_m.q_l.q_r

    // we can free up this memory by combining the grand-product term with q_m.q_l.q_r
    fr::field_t* z_1_poly = w_l_poly; // <-- already performed an ifft transform for this one
    polynomials::copy_polynomial(state.z_1, z_1_poly, n, 4 * n);
    polynomials::fft_with_coset_and_constant(z_1_poly, domain.long_root, domain.generator, state.alpha_squared, domain.long_domain);
    // we want both z_1(X) and z_1(X \omega^-1) - repeat the final 4 indices of the fft to get this
    for (size_t i = 0; i < 4; ++i)
    {
        fr::copy(z_1_poly[4 * n - 4 + i], z_1_poly[4 * n + i]);
    }
    // fr::field_t* shifted_z_1_poly = z_1_poly + 4;

    for (size_t i = 0; i < domain.long_domain; ++i)
    {
        size_t shifted_index = i < (domain.long_domain - 4) ? i + 4 : 4 - (domain.long_domain - i);
        fr::mul(identity_poly[i], z_1_poly[i], identity_poly[i]);
        fr::sub(identity_poly[i], z_1_poly[shifted_index], identity_poly[i]);
        fr::add(identity_poly[i], q_m_poly[i], identity_poly[i]);
    }

    // create placeholder local variable for the quotient polynomial
    fr::field_t* quotient_poly = identity_poly;
    // free memory: w_r, w_o
    fr::field_t* l_1_poly = w_r_poly;
    fr::field_t* l_n_poly = w_r_poly + domain.mid_domain;

    // TODO: implement an O(n) version of this
    polynomials::compute_split_lagrange_polynomial_fft(l_1_poly, l_n_poly, domain);

    for (size_t i = 0; i < domain.mid_domain; ++i)
    {
        size_t shifted_index = i < (domain.mid_domain - 2) ? (i + 2) * 2 : (2 - (domain.mid_domain - i)) * 2;
        fr::mul(l_1_poly[i], z_1_poly[i * 2], l_1_poly[i]);
        fr::mul(l_n_poly[i], z_1_poly[shifted_index], l_n_poly[i]);
        fr::mul(l_1_poly[i], state.alpha_squared, l_1_poly[i]);
        fr::mul(l_n_poly[i], state.alpha_cubed, l_n_poly[i]);
    }

    // z_1 is now free
    // => w_o, w_l are now free

    fr::field_t* sigma_1_poly = w_l_poly;
    fr::field_t* sigma_2_poly = w_o_poly;
    polynomials::copy_polynomial(state.sigma_1, sigma_1_poly, n, 4 * n);
    polynomials::copy_polynomial(state.sigma_2, sigma_2_poly, n, 4 * n);

    polynomials::fft_with_coset(sigma_1_poly, domain.long_root, domain.generator, domain.long_domain);
    polynomials::fft_with_coset(sigma_2_poly, domain.long_root, domain.generator, domain.long_domain);

    for (size_t i = 0; i < domain.long_domain; ++i)
    {
        fr::mul(sigma_1_poly[i], sigma_2_poly[i], sigma_1_poly[i]);
    }

    fr::field_t* sigma_3_poly = w_o_poly;
    polynomials::copy_polynomial(state.sigma_3, sigma_3_poly, n, 4 * n);
    polynomials::fft_with_coset(sigma_3_poly, domain.long_root, domain.generator, domain.long_domain);

    for (size_t i = 0; i < domain.long_domain; ++i)
    {
        fr::mul(sigma_1_poly[i], sigma_3_poly[i], sigma_3_poly[i]);
    }

    fr::field_t* z_2_poly = w_o_poly;
    polynomials::copy_polynomial(state.sigma_2, sigma_2_poly, n, 4 * n);
    polynomials::fft_with_coset_and_constant(z_2_poly, domain.long_root, domain.generator, state.alpha_cubed, domain.long_domain);

    for (size_t i = 0; i < domain.long_domain; ++i)
    {
        size_t shifted_index = i < (domain.long_domain - 4) ? i + 4 : 4 - (domain.long_domain - i);
        fr::mul(sigma_1_poly[i], z_2_poly[i], sigma_1_poly[i]);
        fr::sub(sigma_1_poly[i], z_2_poly[shifted_index], sigma_1_poly[i]);

        // combine product term into quotient poly
        fr::mul(quotient_poly[i], sigma_1_poly[i], quotient_poly[i]);
    }

    // accumulate degree-2n terms into gate_poly_mid
    for (size_t i = 0; i < domain.mid_domain; ++i)
    {
        fr::mul(z_2_poly[i], state.alpha, T0);
        fr::sub(l_1_poly[i], T0, l_1_poly[i]);

        size_t shifted_index = i < (domain.mid_domain - 2) ? (i + 2) * 2 : (2 - (domain.mid_domain - i)) * 2;
        fr::mul(z_2_poly[shifted_index], state.alpha_squared, T0);
        fr::sub(l_n_poly[i], T0, l_n_poly[i]);

        fr::add(gate_poly_mid[i], l_1_poly[i], gate_poly_mid[i]);
        fr::add(gate_poly_mid[i], l_n_poly[i], gate_poly_mid[i]);
    }

    polynomials::ifft_with_coset(gate_poly_mid, domain.mid_root_inverse, domain.generator_inverse, domain.mid_domain);
    memset((void *)(gate_poly_mid + domain.mid_domain), 0, (domain.mid_domain) * sizeof(fr::field_t));

    polynomials::fft_with_coset(gate_poly_mid, domain.long_root, domain.generator, domain.long_domain);

    for (size_t i = 0; i < domain.long_domain; ++i)
    {
        fr::add(quotient_poly[i], gate_poly_mid[i], quotient_poly[i]);
    }

    // final step is to divide by the vanishing polynomial...
}
} // namespace waffle