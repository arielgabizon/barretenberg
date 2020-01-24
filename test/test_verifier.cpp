#include <gtest/gtest.h>

#include <barretenberg/polynomials/polynomial_arithmetic.hpp>
#include <barretenberg/waffle/proof_system/preprocess.hpp>
#include <barretenberg/waffle/proof_system/prover/prover.hpp>
#include <barretenberg/waffle/proof_system/verifier/verifier.hpp>
#include <barretenberg/waffle/proof_system/widgets/arithmetic_widget.hpp>
#include <memory>

namespace
{

using namespace barretenberg;

void generate_test_data(waffle::Prover<3>& state)
{
    size_t n = state.n;
    std::unique_ptr<waffle::ProverArithmeticWidget> widget = std::make_unique<waffle::ProverArithmeticWidget>(n);
    // state.widgets.emplace_back(std::make_unique<waffle::ProverArithmeticWidget>(n));


    fr::field_t T0;
    // even indices = mul gates, odd incides = add gates

    state.w[0].resize(n);
    state.w[1].resize(n);
    state.w[2].resize(n);

    for (size_t i = 0; i < n / 4; ++i)
    {
        state.w[0].at(2 * i) = fr::random_element();
        state.w[1].at(2 * i) = fr::random_element();
        fr::__mul(state.w[0].at(2 * i), state.w[1].at(2 * i), state.w[2].at(2 * i));
        fr::__add(state.w[2][2 * i], state.w[0][2 * i], state.w[2][2 * i]);
        fr::__add(state.w[2][2 * i], state.w[1][2 * i], state.w[2][2 * i]);
        fr::__add(state.w[2][2 * i], fr::one, state.w[2][2 * i]);
        fr::__copy(fr::one, widget->q_l.at(2 * i));
        fr::__copy(fr::one, widget->q_r.at(2 * i));
        fr::__copy(fr::neg_one(), widget->q_o.at(2 * i));
        fr::__copy(fr::one, widget->q_c.at(2 * i));
        fr::__copy(fr::one, widget->q_m.at(2 * i));

        state.w[0].at(2 * i + 1) = fr::random_element();
        state.w[1].at(2 * i + 1) = fr::random_element();
        state.w[2].at(2 * i + 1) = fr::random_element();

        fr::__add(state.w[0].at(2 * i + 1), state.w[1].at(2 * i + 1), T0);
        fr::__add(T0, state.w[2].at(2 * i + 1), widget->q_c.at(2 * i + 1));
        fr::__neg(widget->q_c.at(2 * i + 1), widget->q_c.at(2 * i + 1));
        widget->q_l.at(2 * i + 1) = fr::one;
        widget->q_r.at(2 * i + 1) = fr::one;
        widget->q_o.at(2 * i + 1) = fr::one;
        widget->q_m.at(2 * i + 1) = fr::zero;
    }
    size_t shift = n / 2;
    polynomial_arithmetic::copy_polynomial(&state.w[0].at(0), &state.w[0].at(shift), shift, shift);
    polynomial_arithmetic::copy_polynomial(&state.w[1].at(0), &state.w[1].at(shift), shift, shift);
    polynomial_arithmetic::copy_polynomial(&state.w[2].at(0), &state.w[2].at(shift), shift, shift);
    polynomial_arithmetic::copy_polynomial(&widget->q_m.at(0), &widget->q_m.at(shift), shift, shift);
    polynomial_arithmetic::copy_polynomial(&widget->q_l.at(0), &widget->q_l.at(shift), shift, shift);
    polynomial_arithmetic::copy_polynomial(&widget->q_r.at(0), &widget->q_r.at(shift), shift, shift);
    polynomial_arithmetic::copy_polynomial(&widget->q_o.at(0), &widget->q_o.at(shift), shift, shift);
    polynomial_arithmetic::copy_polynomial(&widget->q_c.at(0), &widget->q_c.at(shift), shift, shift);

    // create basic permutation - second half of witness vector is a copy of the first half
    state.sigma_map[0].resize(n);
    state.sigma_map[1].resize(n);
    state.sigma_map[2].resize(n);

    for (size_t i = 0; i < n / 2; ++i)
    {
        state.sigma_map[0][shift + i] = (uint32_t)i;
        state.sigma_map[1][shift + i] = (uint32_t)i + (1U << 30U);
        state.sigma_map[2][shift + i] = (uint32_t)i + (1U << 31U);
        state.sigma_map[0][i] = (uint32_t)(i + shift);
        state.sigma_map[1][i] = (uint32_t)(i + shift) + (1U << 30U);
        state.sigma_map[2][i] = (uint32_t)(i + shift) + (1U << 31U);
    }
    // make last permutation the same as identity permutation
    state.sigma_map[0][shift - 1] = (uint32_t)shift - 1;
    state.sigma_map[1][shift - 1] = (uint32_t)shift - 1 + (1U << 30U);
    state.sigma_map[2][shift - 1] = (uint32_t)shift - 1 + (1U << 31U);
    state.sigma_map[0][n - 1] = (uint32_t)n - 1;
    state.sigma_map[1][n - 1] = (uint32_t)n - 1 + (1U << 30U);
    state.sigma_map[2][n - 1] = (uint32_t)n - 1 + (1U << 31U);

    state.w[0].at(n - 1) = fr::zero;
    state.w[1].at(n - 1) = fr::zero;
    state.w[2].at(n - 1) = fr::zero;
    widget->q_c.at(n - 1) = fr::zero;
    widget->q_l.at(n - 1) = fr::zero;
    widget->q_r.at(n - 1) = fr::zero;
    widget->q_o.at(n - 1) = fr::zero;
    widget->q_m.at(n - 1) = fr::zero;

    state.w[0].at(shift - 1) = fr::zero;
    state.w[1].at(shift - 1) = fr::zero;
    state.w[2].at(shift - 1) = fr::zero;
    widget->q_c.at(shift - 1) = fr::zero;

    state.widgets.emplace_back(std::move(widget));
}
} // namespace

TEST(verifier, verify_arithmetic_proof_small)
{
    size_t n = 4;

    waffle::Prover<3> state(n);

    generate_test_data(state);

    waffle::Verifier verifier = waffle::preprocess(state);

    // construct proof
    waffle::plonk_proof proof = state.construct_proof();

    // verify proof
    bool result = verifier.verify_proof(proof);

    EXPECT_EQ(result, true);
}

TEST(verifier, verify_arithmetic_proof)
{
    size_t n = 1 << 14;

    waffle::Prover<3> state(n);

    generate_test_data(state);

    waffle::Verifier verifier = waffle::preprocess(state);

    // construct proof
    waffle::plonk_proof proof = state.construct_proof();

    // verify proof
    bool result = verifier.verify_proof(proof);

    EXPECT_EQ(result, true);
}