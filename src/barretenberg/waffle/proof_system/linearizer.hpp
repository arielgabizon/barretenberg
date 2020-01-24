#pragma once

#include "../../types.hpp"
#include "../../curves/bn254/fr.hpp"

#include "../waffle_types.hpp"

namespace waffle
{

template<size_t program_width>
    struct plonk_linear_terms
    {
        std::array<barretenberg::fr::field_t, program_width> w_l;
        barretenberg::fr::field_t z_1;
        barretenberg::fr::field_t sigma_last;
    };

    // This linearisation trick was originated from Mary Maller and the SONIC paper. When computing Kate commitments to the PLONK polynomials, we wish to find the minimum number of polynomial evaluations that the
    // prover must send to the verifier. I.e. we want to find the minimum number of polynomial evaluations that are needed, so that the remaining
    // polynomial evaluations can be expressed as a linear sum of polynomials. The verifier can derive the prover's commitment to this linear polynomial
    // from the original commitments - the prover can provide an evaluation of this linear polynomial, instead of the evaluations of its consitutent polynomials.
    // This shaves 6 field elements off of the proof size!
    template <size_t program_width>
    inline plonk_linear_terms<program_width> compute_linear_terms(const plonk_proof& proof, const g1_challenges& challenges,
        const barretenberg::fr::field_t& l_1, const size_t)
    {
        plonk_linear_terms<program_width> result;
        barretenberg::fr::field_t T0;
        barretenberg::fr::field_t T1;
        barretenberg::fr::field_t T2;


        barretenberg::fr::field_t right_shift = barretenberg::fr::multiplicative_generator;
        barretenberg::fr::field_t output_shift = barretenberg::fr::alternate_multiplicative_generator;

        barretenberg::fr::field_t alpha_pow[6];
        barretenberg::fr::__copy(challenges.alpha, alpha_pow[0]);
        for (size_t i = 1; i < 6; ++i)
        {
            barretenberg::fr::__mul(alpha_pow[i-1], alpha_pow[0], alpha_pow[i]);
        }

        barretenberg::fr::__mul(challenges.z, challenges.beta, T0);
        barretenberg::fr::__add(T0, proof.w_eval[0], T0);
        barretenberg::fr::__add(T0, challenges.gamma, T0);

        barretenberg::fr::__mul(challenges.z, challenges.beta, T1);
        barretenberg::fr::__mul(T1, right_shift, T1);
        barretenberg::fr::__add(T1, proof.w_eval[1], T1);
        barretenberg::fr::__add(T1, challenges.gamma, T1);

        barretenberg::fr::__mul(challenges.z, challenges.beta, T2);
        barretenberg::fr::__mul(T2, output_shift, T2);
        barretenberg::fr::__add(T2, proof.w_eval[2], T2);
        barretenberg::fr::__add(T2, challenges.gamma, T2);

        barretenberg::fr::__mul(T2, T1, T1);
        barretenberg::fr::__mul(T1, T0, T0);
        barretenberg::fr::__mul(T0, alpha_pow[0], result.z_1);

        barretenberg::fr::__mul(proof.sigma_1_eval, challenges.beta, T0);
        barretenberg::fr::__add(T0, proof.w_eval[0], T0);
        barretenberg::fr::__add(T0, challenges.gamma, T0);

        barretenberg::fr::__mul(proof.sigma_2_eval, challenges.beta, T1);
        barretenberg::fr::__add(T1, proof.w_eval[1], T1);
        barretenberg::fr::__add(T1, challenges.gamma, T1);


        barretenberg::fr::__mul(T1, T0, T0);
        barretenberg::fr::__mul(T0, proof.z_1_shifted_eval, T0);
        barretenberg::fr::__mul(T0, alpha_pow[0], result.sigma_last);
        barretenberg::fr::__neg(result.sigma_last, result.sigma_last);
        barretenberg::fr::__mul(result.sigma_last, challenges.beta, result.sigma_last);


        barretenberg::fr::__mul(l_1, alpha_pow[2], T0);
        barretenberg::fr::__add(result.z_1, T0, result.z_1);

        return result;
    }
}