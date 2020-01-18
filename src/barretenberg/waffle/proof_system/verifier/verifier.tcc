#include "./verifier.hpp"
#include <sspsidl.h>

#include "../../../curves/bn254/fq12.hpp"
#include "../../../curves/bn254/g1.hpp"
#include "../../../curves/bn254/g2.hpp"
#include "../../../curves/bn254/pairing.hpp"
#include "../../../curves/bn254/scalar_multiplication/scalar_multiplication.hpp"
#include "../../../polynomials/evaluation_domain.hpp"
#include "../../../polynomials/polynomial_arithmetic.hpp"

#include "../../../types.hpp"

#include "../../reference_string/reference_string.hpp"

#include "../challenge.hpp"
#include "../linearizer.hpp"
#include "../widgets/base_widget.hpp"

using namespace barretenberg;

namespace waffle {
template<size_t program_width>
Verifier<program_width>::Verifier(const size_t subgroup_size)
    : n(subgroup_size)
{}

template<size_t program_width>
Verifier<program_width>::Verifier(Verifier&& other)
    : n(other.n)
{
    reference_string = std::move(other.reference_string);
    for (int j = 0; j < program_width; ++j) {
        g1::copy_affine(other.SIGMA[j], SIGMA[j]);
    }
    for (size_t i = 0; i < other.verifier_widgets.size(); ++i) {
        verifier_widgets.emplace_back(std::move(other.verifier_widgets[i]));
    }
    update_needs_w_shifted();
}

template<size_t program_width>
Verifier<program_width>& Verifier<program_width>::operator=(Verifier<program_width>&& other)
{
    n = other.n;
    reference_string = std::move(other.reference_string);

    for (int j = 0; j < program_width; ++j) {
        g1::copy_affine(other.SIGMA[j], SIGMA[j]);
    }
    verifier_widgets.resize(0);
    for (size_t i = 0; i < other.verifier_widgets.size(); ++i) {
        verifier_widgets.emplace_back(std::move(other.verifier_widgets[i]));
    }
    return *this;
}

template<size_t program_width>
bool Verifier<program_width>::verify_proof(const waffle::plonk_proof<program_width>& proof)
{
    evaluation_domain domain = evaluation_domain(n);

    bool inputs_valid = g1::on_curve(proof.T_LO)
                        // && g1::on_curve(proof.T_MID)
                        // && g1::on_curve(proof.T_HI)
                        // && g1::on_curve(proof.W_L)
                        // && g1::on_curve(proof.W_R)
                        // && g1::on_curve(proof.W_O)
                        && g1::on_curve(proof.Z_1) && g1::on_curve(proof.PI_Z);
    // && g1::on_curve(proof.PI_Z_OMEGA);

    if (!inputs_valid) {
        printf("inputs not valid!\n");
        return false;
    }


    for (int j = 0; j < program_width; ++j) {
        if (!g1::on_curve(SIGMA[j])) {
            printf("instance not valid!\n");
            return false;
        }
    }


    bool widget_instance_valid = true;
    for (size_t i = 0; i < verifier_widgets.size(); ++i) {
        widget_instance_valid = widget_instance_valid && verifier_widgets[i]->verify_instance_commitments();
    }
    if (!widget_instance_valid) {
        printf("widget instance not valid!\n");
        return false;
    }

    for (int j = 0; j < program_width-1; ++j) {
        if (!fr::eq(proof.sigma_eval[j], fr::zero)) {
            printf("proof field elements not valid!\n");
            return false;
        }
    }
    if (!fr::eq(proof.linear_eval, fr::zero)) {
        printf("proof field elements not valid!\n");
        return false;
    }

    // reconstruct challenges
    plonk_challenges challenges;
    fr::field_t alpha_pow[4];
    fr::field_t nu_pow[10];
    challenges.alpha = compute_alpha(proof);
    challenges.gamma = compute_gamma(proof);
    challenges.beta = compute_beta(proof, challenges.gamma);
    challenges.z = compute_evaluation_challenge(proof);

    polynomial_arithmetic::lagrange_evaluations lagrange_evals =
        polynomial_arithmetic::get_lagrange_evaluations(challenges.z, domain);

    // compute the terms we need to derive R(X)
    plonk_linear_terms linear_terms = compute_linear_terms(proof, challenges, lagrange_evals.l_1, n);

    // reconstruct evaluation of quotient polynomial from prover messages
    fr::field_t t_eval;
    fr::field_t T0; //T0 acts as an accumulator for all terms
    fr::field_t T1;
    fr::field_t T2;
    fr::__copy(challenges.alpha, alpha_pow[0]);
    for (size_t i = 1; i < 4; ++i) {
        fr::__mul(alpha_pow[i - 1], alpha_pow[0], alpha_pow[i]);
    }


    // computing the term Z(g\cdot x)\cdot \prod_{w_i(x)+ beta\cdot sigma_i(x)+ gamma}  (last witness doesn't have beta term as optimization)
    fr::__mul(proof.sigma_eval[0], challenges.beta, T0);
    fr::__add(proof.w_eval[0], challenges.gamma, T1);
    fr::__add(T0, T1, T0);
    for (size_t i = 1; i < program_width-2; ++i) {
        fr::__mul(proof.sigma_eval[i], challenges.beta, T2);
        fr::__add(proof.w_r_eval, challenges.gamma, T1);
        fr::__add(T2, T1, T2);
        fr::__mul(T0,T2,T0);
    }
    fr::__add(proof.w_eval[program_width-1], challenges.gamma, T1);
    fr::__mul(T0, T1, T0);
    fr::__mul(T0, proof.z_1_shifted_eval, T0);
    fr::__mul(T0, alpha_pow[0], T0);


    fr::__sub(proof.z_1_shifted_eval, fr::one, T1);
    fr::__mul(T1, lagrange_evals.l_n_minus_1, T1);
    fr::__mul(T1, alpha_pow[1], T1);

    fr::__mul(lagrange_evals.l_1, alpha_pow[2], T2);

    fr::__sub(T1, T2, T1);
    fr::__sub(T1, T0, T1);

    fr::__add(T1, proof.linear_eval, t_eval);

    fr::__invert(lagrange_evals.vanishing_poly, T0);
    fr::__mul(t_eval, T0, t_eval);

    std::array<fr::field_t, program_width-1> z_pow;
    make_z_pow(challenges.z, program_width-1, n, z_pow);

    challenges.nu = compute_linearisation_challenge(proof, t_eval);

    fr::field_t u = compute_kate_separation_challenge(proof, t_eval);
    fr::__copy(challenges.nu, nu_pow[0]);
    for (size_t i = 1; i < 9; ++i) {
        fr::__mul(nu_pow[i - 1], nu_pow[0], nu_pow[i]);
    }

    // reconstruct Kate opening commitments from committed values
    fr::__mul(linear_terms.q_m, nu_pow[0], linear_terms.q_m);
    fr::__mul(linear_terms.q_l, nu_pow[0], linear_terms.q_l);
    fr::__mul(linear_terms.q_r, nu_pow[0], linear_terms.q_r);
    fr::__mul(linear_terms.q_o, nu_pow[0], linear_terms.q_o);
    fr::__mul(linear_terms.q_c, nu_pow[0], linear_terms.q_c);
    fr::__mul(linear_terms.z_1, nu_pow[0], linear_terms.z_1);
    fr::__mul(linear_terms.sigma_last, nu_pow[0], linear_terms.sigma_last);

    fr::__mul(nu_pow[6], u, T0);
    fr::__add(linear_terms.z_1, T0, linear_terms.z_1);

    fr::field_t batch_evaluation;
    fr::__copy(t_eval, batch_evaluation);
    fr::__mul(nu_pow[0], proof.linear_eval, T0);
    fr::__add(batch_evaluation, T0, batch_evaluation);

    for (size_t k = 0; k < program_width; ++k) {

    fr::__mul(nu_pow[k+1], proof.w_eval[k], T0);
    fr::__add(batch_evaluation, T0, batch_evaluation);
    }


    for (size_t k = 0; k < program_width-1; ++k) {
        fr::__mul(nu_pow[k+program_width+1], proof.sigma_eval[k], T0);
        fr::__add(batch_evaluation, T0, batch_evaluation);
    }

    fr::__mul(nu_pow[2*program_width], u, T0);
    fr::__mul(T0, proof.z_1_shifted_eval, T0);
    fr::__add(batch_evaluation, T0, batch_evaluation);

    fr::field_t nu_base = nu_pow[2*program_width+1];

    for (size_t k = 0; k < program_width; ++k) {
        if (needs_w_shifted[k]) {
            fr::__mul(proof.w_l_shifted_eval, nu_base, T0);
            fr::__mul(T0, u, T0);
            fr::__add(batch_evaluation, T0, batch_evaluation);
            fr::__mul(nu_base, nu_pow[0], nu_base);
        }
    }
    for (size_t i = 0; i < verifier_widgets.size(); ++i) {
        nu_base =
            verifier_widgets[i]->compute_batch_evaluation_contribution(batch_evaluation, nu_base, nu_pow[0], proof);
    }

    fr::__neg(batch_evaluation, batch_evaluation);

    fr::field_t z_omega_scalar;
    fr::__mul(challenges.z, domain.root, z_omega_scalar);
    fr::__mul(z_omega_scalar, u, z_omega_scalar);

    std::vector<fr::field_t> scalars;
    std::vector<g1::affine_element> elements;

    elements.emplace_back(proof.Z_1);
    scalars.emplace_back(linear_terms.z_1);

    fr::__copy(nu_pow[2*program_width+1], nu_base);

    for (size_t k = 0; k < program_width; ++k) {
        if (g1::on_curve(proof.W_L)) {
            elements.emplace_back(proof.W_L);
            if (needs_w_shifted[k]) {
                fr::__mul(nu_base, u, T0);
                fr::__add(T0, nu_pow[k+1], T0);
                scalars.emplace_back(T0);
                // scalars.emplace_back(fr::add(nu_pow[1], nu_base));
                fr::__mul(nu_base, nu_pow[0], nu_base);
            } else {
                scalars.emplace_back(nu_pow[k+1]);
            }
        }
    }

    for (size_t k = 0; k < program_width-1; ++k) {
        elements.emplace_back(SIGMA[k]);
        scalars.emplace_back(nu_pow[program_width+k+1]);
    }

    elements.emplace_back(SIGMA[program_width]);
    scalars.emplace_back(linear_terms.sigma_last);

    elements.emplace_back(g1::affine_one);
    scalars.emplace_back(batch_evaluation);

    if (g1::on_curve(proof.PI_Z_OMEGA)) {
        elements.emplace_back(proof.PI_Z_OMEGA);
        scalars.emplace_back(z_omega_scalar);
    }

    elements.emplace_back(proof.PI_Z);
    scalars.emplace_back(challenges.z);

    for (size_t k = 1; k < program_width; ++k) {
        if (g1::on_curve(proof.T[k])) {
            elements.emplace_back(proof.T[k]);
            scalars.emplace_back(z_pow[k-1]);
        }
    }

    VerifierBaseWidget::challenge_coefficients coeffs{
        fr::sqr(fr::sqr(challenges.alpha)), challenges.alpha, nu_base, challenges.nu, challenges.nu
    };

    for (size_t i = 0; i < verifier_widgets.size(); ++i) {
        coeffs = verifier_widgets[i]->append_scalar_multiplication_inputs(coeffs, proof, elements, scalars);
    }

    size_t num_elements = elements.size();
    elements.resize(num_elements * 2);
    scalar_multiplication::generate_pippenger_point_table(&elements[0], &elements[0], num_elements);
    g1::element P[2];

    P[0] = g1::group_exponentiation_inner(proof.PI_Z_OMEGA, u);
    P[1] = scalar_multiplication::pippenger(&scalars[0], &elements[0], num_elements);

    g1::mixed_add(P[1], proof.T_LO, P[1]);
    g1::mixed_add(P[0], proof.PI_Z, P[0]);
    g1::__neg(P[0], P[0]);
    g1::batch_normalize(P, 2);

    g1::affine_element P_affine[2];
    fq::__copy(P[0].x, P_affine[1].x);
    fq::__copy(P[0].y, P_affine[1].y);
    fq::__copy(P[1].x, P_affine[0].x);
    fq::__copy(P[1].y, P_affine[0].y);

    fq12::field_t result =
        pairing::reduced_ate_pairing_batch_precomputed(P_affine, reference_string.precomputed_g2_lines, 2);

    return fq12::eq(result, fq12::one);
}

template <size_t program_width>
void Verifier<program_width>::update_needs_w_shifted() {
    for (size_t i = 0; i < verifier_widgets.size(); ++i)
    {
        for (int j = 0; j < program_width; ++j) {
            needs_w_shifted[j] |= verifier_widgets[i]->version.has_dependency(j);

        }
    }
}

} // namespace waffle