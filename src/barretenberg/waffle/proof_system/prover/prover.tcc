#include "./prover.hpp"
#include <d2d1helper.h>
#include <sspsidl.h>

#include "../../../curves/bn254/fr.hpp"
#include "../../../curves/bn254/g1.hpp"
#include "../../../curves/bn254/g2.hpp"
#include "../../../curves/bn254/scalar_multiplication/scalar_multiplication.hpp"
#include "../../../io/io.hpp"
#include "../../../polynomials/polynomial_arithmetic.hpp"

#include "../../reference_string/reference_string.hpp"

#include "../widgets/base_widget.hpp"
#include "../linearizer.hpp"
#include "../challenge.hpp"
#include "../permutation.hpp"

using namespace barretenberg;

namespace waffle
{
template<size_t program_width>
Prover<program_width>::Prover(const size_t __n) :
n(__n),
fft_state(n),
reference_string(n)
{
}

template<size_t program_width>
Prover<program_width>::Prover(Prover<program_width> &&other) :
n(other.n),
w(std::move(other.w)),
fft_state(std::move(other.fft_state)),
sigma(std::move(other.sigma))
//
//sigma_1_mapping(std::move(other.sigma_1_mapping)),
//sigma_2_mapping(std::move(other.sigma_2_mapping)),
//sigma_3_mapping(std::move(other.sigma_3_mapping))
{
    for (size_t i = 0; i < other.widgets.size(); ++i)
    {
        widgets.emplace_back(std::move(other.widgets[i]));
    }
    reference_string = std::move(other.reference_string);
    update_needs_w_shifted();
}

template<size_t program_width>
Prover<program_width>& Prover<program_width>::operator=(Prover<program_width> &&other)
{
    n = other.n;
    w = std::move(other.w);
    fft_state = waffle::CircuitFFTState(std::move(other.fft_state));
    sigma = std::move(other.sigma);
    widgets.resize(0);
    for (size_t i = 0; i < other.widgets.size(); ++i)
    {
        widgets.emplace_back(std::move(other.widgets[i]));
    }
    reference_string = std::move(other.reference_string);
    return *this;
}

template<size_t program_width>
void Prover<program_width>::compute_wire_commitments()
{
    scalar_multiplication::multiplication_state mul_state[program_width];

    for (size_t i = 0; i < program_width; ++i) {
        mul_state[i]= {reference_string.monomials, w[i].get_coefficients(), n, {}};
    }

    scalar_multiplication::batched_scalar_multiplications(mul_state, program_width);

//    // TODO: make a method for normal-to-affine copies :/
//    fq::__copy(mul_state[0].output.x, proof.W_L.x);
//    fq::__copy(mul_state[1].output.x, proof.W_R.x);
//    fq::__copy(mul_state[2].output.x, proof.W_O.x);
//    fq::__copy(mul_state[0].output.y, proof.W_L.y);
//    fq::__copy(mul_state[1].output.y, proof.W_R.y);
//    fq::__copy(mul_state[2].output.y, proof.W_O.y);

    for (size_t i = 0; i < program_width; ++i) {
        fq::__copy(mul_state[i].output.x, proof.W[i].x);
        fq::__copy(mul_state[i].output.y, proof.W[i].y);
    }

    // compute beta, gamma
    challenges.gamma = compute_gamma(proof);
    challenges.beta = compute_beta(proof, challenges.gamma);
}

template<size_t program_width>
void Prover<program_width>::compute_z_commitment()
{
    scalar_multiplication::multiplication_state mul_state{
        reference_string.monomials,
        z.get_coefficients(),
        n,
        g1::element()};

    scalar_multiplication::batched_scalar_multiplications(&mul_state, 1);

    // TODO: make a method for normal-to-affine copies :/
    fq::__copy(mul_state.output.x, proof.Z_1.x);
    fq::__copy(mul_state.output.y, proof.Z_1.y);

    // compute alpha
    // TODO: does this really belong here?
    challenges.alpha = compute_alpha(proof);
}

template<size_t program_width>
void Prover<program_width>::compute_quotient_commitment()
{
    std::array<scalar_multiplication::multiplication_state, program_width> mul_state;
    for (size_t i = 0; i < program_width; ++i) {
        mul_state[i]= {reference_string.monomials, &fft_state.quotient_large.get_coefficients()[i*n],n {} };
    }



    for (size_t i = 0; i < program_width; ++i) {
        g1::jacobian_to_affine(mul_state[i].output, proof.T[i]);
    }

//g1::jacobian_to_affine(mul_state[0].output, proof.T_LO);
//    g1::jacobian_to_affine(mul_state[1].output, proof.T_MID);
//    g1::jacobian_to_affine(mul_state[2].output, proof.T_HI);

 //   challenges.z = compute_evaluation_challenge(proof);
}

template<size_t program_width>
void Prover<program_width>::compute_wire_coefficients()
{
    for (size_t i = 0; i < program_width; ++i) {
        fft_state.w[i] = polynomial(w[i], n);
    }

    for (size_t i = 0; i < program_width; ++i) {
        w[i].ifft(fft_state.small_domain);
    }

}

template<size_t program_width>
void Prover<program_width>::compute_z_coefficients()
{
    std::array<polynomial,program_width*2> accumulators;
    for (auto i=0; i< accumulators.length(); ++i)
    {
        accumulators[i] = polynomial(n + 1, n + 1);
    }
#ifndef NO_MULTITHREADING
#pragma omp parallel for
#endif
    for (size_t j = 0; j < fft_state.small_domain.num_threads; ++j)
    {
        fr::field_t work_root;
        fr::field_t thread_root;
        fr::__pow_small(fft_state.small_domain.root, j * fft_state.small_domain.thread_size, thread_root);
        fr::__mul(thread_root, challenges.beta, work_root);
        std::array<fr::field_t, program_width-1> k;
        fr::compute_coset_generators(program_width-1, 1 << 30, k);
//        fr::field_t k1 = fr::multiplicative_generator;
//        fr::field_t k2 = fr::alternate_multiplicative_generator;

        for (size_t i = (j * fft_state.small_domain.thread_size); i < ((j + 1) * fft_state.small_domain.thread_size); ++i)
        {
            std::array<fr::field_t, program_width> T;
            T[0] = work_root;
            for(size_t l = 0; l < program_width; ++l){
                if (l>0)
                    fr::__mul(work_root, k[l], T[l]);

                fr::__add(T[l], challenges.gamma, T[l]);
                fr::__add(T[l], fft_state.w_fft[l][i], accumulators[l][i + 1]);
                fr::__mul(sigma[l][i], challenges.beta, T[l]);
                fr::__add(T[l], challenges.gamma, T[l]);
                fr::__add(T[l], fft_state.w_fft[l][i], accumulators[program_width+l][i + 1]);
            }
            fr::__mul(work_root, fft_state.small_domain.root, work_root);


        }
    }

    // step 2: compute the constituent components of Z1(X), Z2(X). This is a small bottleneck, as we have
    // 6 non-parallelizable processes
#ifndef NO_MULTITHREADING
#pragma omp parallel for
#endif
    for (size_t i = 0; i < 2 * program_width; ++i)
    {
        fr::field_t *coeffs = accumulators[i].get_coefficients();
        coeffs[0] = fr::one;
        for (size_t j = 1; j < fft_state.small_domain.size - 1; ++j)
        {
            fr::__mul(coeffs[j + 1], coeffs[j], coeffs[j + 1]);
        }
    }

    polynomial z2(n, n);
    z = polynomial(n, n);
    // step 3: concatenate together the accumulator elements into Z1(X), Z2(X)
    // THIS CODE ASSUMES program_width>=2!
    ITERATE_OVER_DOMAIN_START(fft_state.small_domain);
        fr::__mul(accumulators[0][i], accumulators[1][i], z[i]);
        fr::__mul(accumulators[3][i], accumulators[4][i], z2[i]);
        for(size_t l =2; l< program_width; ++l) {
            fr::__mul(z[i], accumulators[l][i], z[i]);
            fr::__mul(z2[i], accumulators[l+program_width][i], z2[i]);
        }
    ITERATE_OVER_DOMAIN_END;

    fr::batch_invert(z2.get_coefficients(), fft_state.small_domain.size);

    ITERATE_OVER_DOMAIN_START(fft_state.small_domain);
        fr::__mul(z[i], z2[i], z[i]);
    ITERATE_OVER_DOMAIN_END;

    z.ifft(fft_state.small_domain);
}

template<size_t program_width>
void Prover<program_width>::compute_permutation_grand_product_coefficients(polynomial& z_fft)
{
    // Our permutation check boils down to two 'grand product' arguments,
    // that we represent with a single polynomial Z(X).
    // We want to test that Z(X) has been constructed correctly.
    // When evaluated at elements of w \in H, the numerator of Z(w) will equal the
    // identity permutation grand product, and the denominator will equal the copy permutation grand product.

    // The identity that we need to evaluate is: Z(X.w).(permutation grand product) = Z(X).(identity grand product)
    // i.e. The next element of Z is equal to the current element of Z, multiplied by (identity grand product) / (permutation grand product)

    // This method computes `Z(X).(identity grand product).{alpha}`.
    // The random `alpha` is there to ensure our grand product polynomial identity is linearly independent from the
    // other polynomial identities that we are going to roll into the quotient polynomial T(X).

    // Specifically, we want to compute (in the original 3 wire case):
    // (w_l(X) + \beta.sigma1(X) + \gamma).(w_r(X) + \beta.sigma2(X) + \gamma).(w_o(X) + \beta.sigma3(X) + \gamma).Z(X).alpha
    // Once we divide by the vanishing polynomial, this will be a degree 3n polynomial.

    // Step 1: convert sigma1(X), sigma2(X), sigma3(X) from point-evaluation form into coefficient form.
    // When we do this, scale the coefficients up by `beta` - we can get this for free by rolling it into the ifft transform
    for (int i = 0; i < program_width; ++i) {
        sigma[i].ifft_with_constant(fft_state.small_domain, challenges.beta);
    }

    // Step 2: convert sigma1(X), sigma2(X), sigma3(X), Z(X) back into point-evaluation form, but this time evaluated
    // at the 4n'th roots of unity.
    
    // Step 2a: Make copies when doing this we'll need the coefficient form polynomials later
    std::array<polynomial, program_width> sigma_fft;
    for (int i = 0; i < program_width; ++i) {
        sigma_fft[i] = polynomial(sigma[i], fft_state.large_domain.size);
    }
    // z_fft = polynomial(z, fft_state.large_domain.size + 4);

    // add `gamma` to sigma_1(X), sigma2(X), sigma3(X), so that we don't have to add it into each evaluation
    for (int i = 0; i < program_width; ++i) {
        fr::__add(sigma_fft[i][0], challenges.gamma, sigma_fft[i][0]); // sigma1_fft = \beta.sigma_1(X) + \gamma
    }

    // before performing our fft, add w_l(X), w_r(X), w_o(X) into sigma1_fft, sigma2_fft, sigma3_fft,
    // (cheaper to add n terms in coefficient form, than 4n terms over our extended evaluation domain)
    ITERATE_OVER_DOMAIN_START(fft_state.small_domain);
    for (int j = 0; j < program_width; ++j) {
        fr::__add(sigma_fft[j][i], w[j][i], sigma_fft[j][i]); // sigma1_fft = w_l(X) + \beta.sigma_1(X) + \gamma
    }
    ITERATE_OVER_DOMAIN_END;

    // Step 2c: perform fft transforms to map into point-evaluation form.
    // (use coset fft so that we get evaluations over the roots of unity * multiplicative generator)
    // (if we evaluated at the raw roots of unity, dividing by the vanishing polynomial would require dividing by zero)
    for (int i = 0; i < program_width; ++i) {
        sigma_fft[i].coset_fft(fft_state.large_domain);
    }
    // Multiply Z(X) by \alpha^2 when performing fft transform - we get this for free if we roll \alpha^2 into the multiplicative generator
    z_fft.coset_fft_with_constant(fft_state.large_domain, challenges.alpha);

    // We actually want Z(X.w), not Z(X)! But that's easy to get. z_fft contains Z(X) evaluated at the 4n'th roots of unity.
    // So z_fft(i) = Z(w^{i/4})
    // i.e. z_fft(i + 4) = Z(w^{i/4}.w)
    // => if virtual term 'foo' contains a 4n fft of Z(X.w), then z_fft(i + 4) = foo(i)
    // So all we need to do, to get Z(X.w) is to offset indexes to z_fft by 4.
    // If `i >= 4n  4`, we need to wrap around to the start - so just append the 4 starting elements to the end of z_fft
    z_fft.add_lagrange_base_coefficient(z_fft[0]);
    z_fft.add_lagrange_base_coefficient(z_fft[1]);
    z_fft.add_lagrange_base_coefficient(z_fft[2]);
    z_fft.add_lagrange_base_coefficient(z_fft[3]);


    // Step 4: Set the quotient polynomial to be equal to
    // (w_l(X) + \beta.sigma1(X) + \gamma).(w_r(X) + \beta.sigma2(X) + \gamma).(w_o(X) + \beta.sigma3(X) + \gamma).Z(X).alpha
    ITERATE_OVER_DOMAIN_START(fft_state.large_domain);
    for (int j = 1; j < program_width; ++j) {
        fr::__mul(sigma_fft[0][i],
                  sigma_fft[j][i],
                  sigma_fft[j][i]); // sigma1_fft = (w_l(X) + B.sigma_1(X) + \gamma).(w_r(X) + B.sigma_2(X) + \gamma)
    }
    fr::__mul(sigma_fft[0][i], z_fft[i + 4], sigma_fft[0][i]); // sigma1_fft = (w_l(X) + B.sigma_1(X) + \gamma).(w_r(X) + B.sigma_2(X) + \gamma).(w_o(X) + B.sigma_3(X) + \gamma).Z(X.omega)
    fr::__neg(sigma_fft[0][i], fft_state.quotient_large[i]); // Q(X) -= (w_l(X) + B.sigma_1(X) + \gamma).(w_r(X) + B.sigma_2(X) + \gamma).(w_o(X) + B.sigma_3(X) + \gamma).Z(X.omega)
    ITERATE_OVER_DOMAIN_END;
}

template<size_t program_width>
void Prover<program_width>::compute_identity_grand_product_coefficients(polynomial &z_fft)
{
//    fr::field_t right_shift = fr::multiplicative_generator;
//    fr::field_t output_shift = fr::alternate_multiplicative_generator;
    std::array<fr::field_t, program_width-1> k;
    fr::compute_coset_generators(program_width-1, 1 << 30, k);

#ifndef NO_MULTITHREADING
    #pragma omp parallel for
#endif
    for (size_t j = 0; j < fft_state.large_domain.num_threads; ++j)
    {
//        fr::field_t T0;
//        fr::field_t T1;
//        fr::field_t T2;
        std::array<fr::field_t, program_width> T;
        fr::field_t beta_id;

        fr::field_t work_root;
        fr::__pow_small(fft_state.large_domain.root, j * fft_state.large_domain.thread_size, work_root);
        fr::__mul(work_root, fr::multiplicative_generator, work_root);
        for (size_t i = (j * fft_state.large_domain.thread_size); i < ((j + 1) * fft_state.large_domain.thread_size); ++i)
        {
            fr::__mul(work_root, challenges.beta, beta_id);
            fr::__add(beta_id, challenges.gamma, T[0]);
            fr::__add(T[0], fft_state.w_fft[0][i], T[0]);
            for (int l = 1; l < program_width; ++l) {
                fr::__mul(beta_id, k[l - 1], T[1]);
                fr::__add(T[l], challenges.gamma, T[l]);
                fr::__add(T[l], fft_state.w_fft[l][i], T[l]);

                // combine three identity product terms, with z_1_poly evaluation

                fr::__mul(T[0], T[l], T[0]);
            }
            fr::__mul(T[0], z_fft[i], T[0]);
            fr::__add(fft_state.quotient_large[i], T[0], fft_state.quotient_large[i]);
            fr::__mul(work_root, fft_state.large_domain.root, work_root);
        }
    }

    // // We can shrink the evaluation domain by 2 for the wire polynomials and Z(X), to save on memory
    // fft_state.w_l_fft.shrink_evaluation_domain(2);
    // fft_state.w_r_fft.shrink_evaluation_domain(2);
    // fft_state.w_o_fft.shrink_evaluation_domain(2);
    // z_fft.shrink_evaluation_domain(2);

    // size = 2n, max size = 2n + 4 (appending 4 coefficients after poly arithmetic call)
    polynomial l_1(n + n, n + n + 4);
    polynomial_arithmetic::compute_lagrange_polynomial_fft(l_1.get_coefficients(), fft_state.small_domain, fft_state.mid_domain);
    l_1.add_lagrange_base_coefficient(l_1[0]);
    l_1.add_lagrange_base_coefficient(l_1.at(1));
    l_1.add_lagrange_base_coefficient(l_1.at(2));
    l_1.add_lagrange_base_coefficient(l_1.at(3));

    // accumulate degree-2n terms into gate_poly_mid
    fr::field_t alpha_squared = fr::sqr(challenges.alpha);

    ITERATE_OVER_DOMAIN_START(fft_state.mid_domain);
        fr::field_t T4;
        fr::field_t T6;

        // Step 1: Compute (Z(X.w) - 1).(\alpha^3).L{n-1}(X)
        // The \alpha^3 term is so that we can subsume this polynomial into the quotient polynomial,
        // whilst ensuring the term is linearly independent form the other terms in the quotient polynomial

        // We want to verify that Z(X) equals `1` when evaluated at `w_n`, the 'last' element of our multiplicative subgroup H.
        // But PLONK's 'vanishing polynomial', Z*_H(X), isn't the true vanishing polynomial of subgroup H.
        // We need to cut a root of unity out of Z*_H(X), specifically `w_n`, for our grand product argument.
        // When evaluating Z(X) has been constructed correctly, we verify that Z(X.w).(identity permutation product) = Z(X).(sigma permutation product),
        // for all X \in H. But this relationship breaks down for X = w_n, because Z(X.w) will evaluate to the *first* element of our grand product argument.
        // The last element of Z(X) has a dependency on the first element, so the first element cannot have a dependency on the last element.

        // TODO: With the reduction from 2 Z polynomials to a single Z(X), the above no longer applies
        // TODO: Fix this to remove the (Z(X.w) - 1).L_{n-1}(X) check
    
        // To summarise, we can't verify claims about Z(X) when evaluated at `w_n`.
        // But we can verify claims about Z(X.w) when evaluated at `w_{n-1}`, which is the same thing
    
        // To summarise the summary: If Z(w_n) = 1, then (Z(X.w) - 1).L_{n-1}(X) will be divisible by Z_H*(X)
        // => add linearly independent term (Z(X.w) - 1).(\alpha^3).L{n-1}(X) into the quotient polynomial to check this

        // z_fft already contains evaluations of Z(X).(\alpha^2)
        // at the (2n)'th roots of unity
        // => to get Z(X.w) instead of Z(X), index element (i+2) instead of i
        fr::__sub(z_fft[2 * i + 4], challenges.alpha, T6); // T6 = (Z(X.w) - 1).(\alpha^2)
        fr::__mul(T6, challenges.alpha, T6); // T6 = (Z(X.w) - 1).(\alpha^3)
        fr::__mul(T6, l_1[i + 4], T6); // T6 = (Z(X.w) - 1).(\alpha^3).L{n-1}(X)

        // Step 2: Compute (Z(X) - 1).(\alpha^4).L1(X)
        // We need to verify that Z(X) equals `1` when evaluated at the first element of our subgroup H
        // i.e. Z(X) starts at 1 and ends at 1
        // The `alpha^4` term is so that we can add this as a linearly independent term in our quotient polynomial

        fr::__sub(z_fft[2 * i], challenges.alpha, T4); // T4 = (Z(X) - 1).(\alpha^2)
        fr::__mul(T4, alpha_squared, T4); // T4 = (Z(X) - 1).(\alpha^4)
        fr::__mul(T4, l_1[i], T4); // T4 = (Z(X) - 1).(\alpha^2).L1(X)
    
        // Add T4 and T6 into the degree 2n component of the quotient polynomial
        fr::__add(T4, T6, fft_state.quotient_mid[i]);
    ITERATE_OVER_DOMAIN_END;
}

template<size_t program_width>
void Prover<program_width>::compute_quotient_polynomial()
{
    fft_state.quotient_large.resize(4 * n);
    fft_state.quotient_mid.resize(2 * n);

    compute_wire_coefficients();

    compute_wire_commitments();

    compute_z_coefficients();

    compute_z_commitment();

    for (int i = 0; i < ; ++i) {
        fft_state.w_fft[i] = polynomial(w[i], 4 * n + 4);
        fft_state.w_fft[i].coset_fft(fft_state.large_domain);
        fft_state.w_fft[i].add_lagrange_base_coefficient(fft_state.w_fft[i][0]);
        fft_state.w_fft[i].add_lagrange_base_coefficient(fft_state.w_fft[i][1]);
        fft_state.w_fft[i].add_lagrange_base_coefficient(fft_state.w_fft[i][2]);

    }

    polynomial z_fft(z, 4 * n + 4);

    compute_permutation_grand_product_coefficients(z_fft);

    compute_identity_grand_product_coefficients(z_fft);

    fr::field_t alpha_base = fr::sqr(fr::sqr(challenges.alpha));
    fr::mul(challenges.alpha, alpha_base);
    for (size_t i = 0; i < widgets.size(); ++i)
    {
        alpha_base = widgets[i]->compute_quotient_contribution(alpha_base, challenges.alpha, fft_state);
    }

    polynomial_arithmetic::divide_by_pseudo_vanishing_polynomial(fft_state.quotient_mid.get_coefficients(), fft_state.small_domain, fft_state.mid_domain);
    polynomial_arithmetic::divide_by_pseudo_vanishing_polynomial(fft_state.quotient_large.get_coefficients(), fft_state.small_domain, fft_state.large_domain);


    fft_state.quotient_mid.coset_ifft(fft_state.mid_domain);
    fft_state.quotient_large.coset_ifft(fft_state.large_domain);


    ITERATE_OVER_DOMAIN_START(fft_state.mid_domain);
        fr::__add(fft_state.quotient_large[i], fft_state.quotient_mid[i], fft_state.quotient_large[i]);
    ITERATE_OVER_DOMAIN_END;

}

template<size_t program_width>
fr::field_t Prover<program_width>::compute_linearisation_coefficients()
{
    r.resize_unsafe(n);
    // ok... now we need to evaluate polynomials. Jeepers
    fr::field_t beta_inv;
    fr::__invert(challenges.beta, beta_inv);
    fr::field_t shifted_z;
    fr::__mul(challenges.z, fft_state.small_domain.root, shifted_z);

    // evaluate the prover and instance polynomials.
    // (we don't need to evaluate the quotient polynomial, that can be derived by the verifier)
    for (size_t k = 0; k < program_width; ++k) {
        proof.w_eval[k] = w[k].evaluate(challenges.z, n);

    }



    for (int j = 0; j < program_width; ++j) {
        if (needs_w_shifted[j]) {
            proof.w_shifted_eval[j] = w[j].evaluate(shifted_z, n);
        }
    }

    for (size_t k = 0; k < program_width-1; ++k) {
        proof.sigma_eval[k] = sigma[k].evaluate(challenges.z, n);
    }
    proof.z_1_shifted_eval = z.evaluate(shifted_z, n);

    for (size_t i = 0; i < widgets.size(); ++i)
    {
        widgets[i]->compute_proof_elements(proof, challenges.z);
    }
    fr::field_t t_eval = fft_state.quotient_large.evaluate(challenges.z, 3 * n);
    // we scaled the sigma polynomials up by beta, so scale back down
    for (size_t k = 0; k < program_width-1; ++k) {
        fr::__mul(proof.sigma_eval[k], beta_inv, proof.sigma_eval[k]);
    }

    polynomial_arithmetic::lagrange_evaluations lagrange_evals = polynomial_arithmetic::get_lagrange_evaluations(challenges.z, fft_state.small_domain);
    plonk_linear_terms linear_terms = compute_linear_terms(proof, challenges, lagrange_evals.l_1, n);

    ITERATE_OVER_DOMAIN_START(fft_state.small_domain);
        fr::field_t T0;
        fr::field_t T1;
        fr::__mul(z[i], linear_terms.z_1, T0);
        fr::__mul(sigma[program_width-1][i], linear_terms.sigma_last , T1);
        // we scaled sigma_3[i] by beta, need to correct for that...
        fr::__mul(T1, beta_inv, T1);
        fr::__add(T0, T1, r[i]);
    ITERATE_OVER_DOMAIN_END;

    fr::field_t alpha_base = fr::sqr(fr::sqr(challenges.alpha));
    for (size_t i = 0; i < widgets.size(); ++i)
    {
        alpha_base = widgets[i]->compute_linear_contribution(alpha_base, challenges.alpha, proof, fft_state.small_domain, r);
    }

    proof.linear_eval = r.evaluate(challenges.z, n);
    return t_eval;
}

template<size_t program_width>
void Prover<program_width>::compute_opening_elements()
{

    fr::field_t t_eval = compute_linearisation_coefficients();

    challenges.nu = compute_linearisation_challenge(proof, t_eval);

    size_t nu_degree = 2 * program_width + 2;
    fr::field_t nu_powers[nu_degree];
    fr::__copy(challenges.nu, nu_powers[0]);
    for (size_t i = 1; i < nu_degree; ++i)
    {
        fr::__mul(nu_powers[i - 1], nu_powers[0], nu_powers[i]);
    }

    fr::field_t beta_inv;
    fr::__invert(challenges.beta, beta_inv);

    // Next step: compute the two Kate polynomial commitments, and associated opening proofs
    // We start with the opening proofs. We have two evaluation points: z and z.omega
    // We need to create random linear combinations of each individual polynomial and combine them
    // The poly opened at z is a combination of : the quotient poly T, the linearized poly r, the w witness polynomials,
    // and the w-1 first sigma polynomials describing the permutation.
    // The poly opened at z.omega is the permutation argument polynomial Z
    polynomial opening_poly(n, n);
    polynomial shifted_opening_poly(n, n);
    fr::field_t z_pow_n;
    fr::field_t z_pow_2_n;
    fr::__pow_small(challenges.z, n, z_pow_n);
    fr::__pow_small(challenges.z, 2 * n, z_pow_2_n);

    ITERATE_OVER_DOMAIN_START(fft_state.small_domain);
        fr::field_t T;
        fr::field_t T1;
        fr::field_t T2;
        fr::field_t T3;
        fr::field_t T4;
        fr::field_t T5;
        fr::field_t T8;
        fr::field_t T9;
        std::array<fr::field_t, program_width> W;
        fr::field_t acc;

        // contributions from the quotient polynomial (first contribution added in end)
        fr::__mul(fft_state.quotient_large[i+n], z_pow_n, acc);
        fr::__mul(fft_state.quotient_large[i+n+n], z_pow_2_n, T);
        fr::__add(acc, T, acc);
        fr::__mul(r[i], nu_powers[0], T);
        for (int k = 0; k < program_width; ++k) {
            fr::__mul(w[k][i], nu_powers[k+1], W[k]);
            fr::__add(acc, W[k], acc);
        }

        std::array<fr::field_t, program_width-1> S; //temp variables to accumulate sigmas' contribution
        fr::field_t S_sum;
        fr::__mul(sigma[0][i], nu_powers[program_width+1], S_sum);

        //add random coefficients to the sigma's (except the last that is included in the linearizer r
        for (int l = 1; l < program_width-1; ++l) {
            fr::__mul(sigma[l][i], nu_powers[program_width+l+1], S[1]);
            fr::__add(S_sum,S[l],S_sum);
        }
        // we added a \beta multiplier to sigma_1(X), sigma_2(X), (and also sigma_3(X), s_id(X)) - need to undo that here
        fr::__mul(S_sum, beta_inv, S_sum);

        fr::__add(acc, S_sum, acc);
        fr::__mul(z[i], nu_powers[nu_degree-2], shifted_opening_poly[i]);
        fr::__add(fft_state.quotient_large[i], acc, opening_poly[i]);
    ITERATE_OVER_DOMAIN_END;

    fr::field_t nu_base = nu_powers[nu_degree-1];

    for (int k = 0; k < program_width ; ++k) {

        if (needs_w_shifted[k]) {
            ITERATE_OVER_DOMAIN_START(fft_state.small_domain);
            fr::field_t T0;
            fr::__mul(nu_base, w[k][i], T0);
            fr::__add(shifted_opening_poly[i], T0, shifted_opening_poly[i]);
            ITERATE_OVER_DOMAIN_END;
            nu_base = fr::mul(nu_base, challenges.nu);
        }
    }


    for (size_t i = 0; i < widgets.size(); ++i) {
        nu_base = widgets[i]->compute_opening_poly_contribution(
            &opening_poly[0],fft_state.small_domain, nu_base, nu_powers[0]);
    }

    fr::field_t shifted_z;
    fr::__mul(challenges.z, fft_state.small_domain.root, shifted_z);

    opening_poly.compute_kate_opening_coefficients(challenges.z);

    shifted_opening_poly.compute_kate_opening_coefficients(shifted_z);

    g1::element PI_Z = scalar_multiplication::pippenger(opening_poly.get_coefficients(), reference_string.monomials, n);
    g1::element PI_Z_OMEGA =
        scalar_multiplication::pippenger(shifted_opening_poly.get_coefficients(), reference_string.monomials, n);

    g1::jacobian_to_affine(PI_Z, proof.PI_Z);
    g1::jacobian_to_affine(PI_Z_OMEGA, proof.PI_Z_OMEGA);
}

template <size_t program_width>
plonk_proof<program_width> Prover<program_width>::construct_proof()
{
    for (int i = 0; i < program_width; ++i) {
        compute_permutation_lagrange_base_single(sigma[i], sigma_mapping[i], fft_state.small_domain);

    }
    compute_quotient_polynomial();
    compute_quotient_commitment();
    compute_opening_elements();
    return proof;
}

template <size_t program_width>
void Prover<program_width>::update_needs_w_shifted() {
    for (size_t i = 0; i < widgets.size(); ++i)
    {
        for (int j = 0; j < program_width; ++j) {
            needs_w_shifted[j] |= widgets[i]->version.has_dependency(j);

        }
    }
}

template <size_t program_width>
void Prover<program_width>::reset()
{
    for (int j = program_width; j < ; ++j) {

        w[j].fft(fft_state.small_domain);
        sigma[j] = polynomial(0, 0);
        fft_state.w_fft[j] = polynomial(0, 0);
    }    fft_state.quotient_mid = polynomial(0, 0);
    fft_state.quotient_large = polynomial(0, 0);
    r = polynomial(0, 0);
    for (size_t i = 0; i < widgets.size(); ++i)
    {
        widgets[i]->reset(fft_state.small_domain);
    }
}
} // namespace waffle