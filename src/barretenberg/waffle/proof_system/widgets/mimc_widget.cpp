#include "./mimc_widget.hpp"

#include "../../../types.hpp"

#include "../../../curves/bn254/fr.hpp"
#include "../../../curves/bn254/g1.hpp"
#include "../../../curves/bn254/scalar_multiplication.hpp"

using namespace barretenberg;

namespace waffle
{
ProverMiMCWidget::ProverMiMCWidget(const size_t _n) :
    ProverBaseWidget(
        static_cast<size_t>(WidgetVersionControl::Dependencies::REQUIRES_W_O_SHIFTED),
        static_cast<size_t>(WidgetVersionControl::Features::HAS_MIMC_SELECTORS)
    ),
    n(_n)
{
    q_mimc_selector.resize(_n);
    q_mimc_coefficient.resize(_n);
}

ProverMiMCWidget::ProverMiMCWidget(const ProverMiMCWidget &other) :
    ProverBaseWidget(other),
    n(other.n)
{
    q_mimc_selector = polynomial(other.q_mimc_selector);
    q_mimc_coefficient = polynomial(other.q_mimc_coefficient);
}

ProverMiMCWidget::ProverMiMCWidget(ProverMiMCWidget &&other) :
    ProverBaseWidget(other),
    n(other.n)
{
    q_mimc_selector = polynomial(other.q_mimc_selector);
    q_mimc_coefficient = polynomial(other.q_mimc_coefficient);
}

ProverMiMCWidget& ProverMiMCWidget::operator=(const ProverMiMCWidget &other)
{
    q_mimc_selector = polynomial(other.q_mimc_selector);
    q_mimc_coefficient = polynomial(other.q_mimc_coefficient);
    n = other.n;
    version = WidgetVersionControl(other.version);
    return *this;
}

ProverMiMCWidget& ProverMiMCWidget::operator=(ProverMiMCWidget &&other)
{
    q_mimc_selector = polynomial(other.q_mimc_selector);
    q_mimc_coefficient = polynomial(other.q_mimc_coefficient);
    n = other.n;
    version = WidgetVersionControl(other.version);
    return *this;
}

fr::field_t ProverMiMCWidget::compute_quotient_contribution(const barretenberg::fr::field_t& alpha_base, const barretenberg::fr::field_t &alpha_step, CircuitFFTState& circuit_state)
{
    q_mimc_selector.ifft(circuit_state.small_domain);
    q_mimc_coefficient.ifft(circuit_state.small_domain);

    polynomial q_mimc_selector_fft = polynomial(q_mimc_selector, circuit_state.large_domain.size);
    polynomial q_mimc_coefficient_fft = polynomial(q_mimc_coefficient, circuit_state.large_domain.size);

    q_mimc_selector_fft.coset_fft_with_constant(circuit_state.large_domain, alpha_base);
    q_mimc_coefficient_fft.coset_fft(circuit_state.large_domain);
    ITERATE_OVER_DOMAIN_START(circuit_state.large_domain);
        fr::field_t T0;
        fr::field_t T1;
        fr::field_t T2;
        fr::__add_with_coarse_reduction(circuit_state.w_o_fft[i], circuit_state.w_l_fft[i], T0); // T0 = w_o + w_l
        fr::__add_with_coarse_reduction(T0, q_mimc_coefficient_fft[i], T0);  // T0 = (w_o + w_l + q_c)
        fr::__sqr_with_coarse_reduction(T0, T1);              // T1 = (w_o + w_l + q_c)^2
        fr::__mul_with_coarse_reduction(T1, T0, T1);          // T1 = (w_o + w_l + q_c)^3
        fr::__sub_with_coarse_reduction(T1, circuit_state.w_r_fft[i], T1); // T1 = (w_o + w_l + q_c)^3 - w_r
        fr::__sqr_with_coarse_reduction(circuit_state.w_r_fft[i], T2); // T2 = w_r^2
        fr::__mul_with_coarse_reduction(T2, T0, T2);  // T2 = (w_o + w_l + q_c).w_r^2 
        fr::__sub_with_coarse_reduction(T2, circuit_state.w_o_fft[i + 4], T2); // T2 = (w_o + w_l + q_c).w_r^2 - w_{o.next}
        fr::__mul_with_coarse_reduction(T2, alpha_step, T2);  // T2 = (w_o + w_l + q_c).w_r^2 - w_{o.next}).alpha
        fr::__add_with_coarse_reduction(T1, T2, T1);  // T1 = ((w_o + w_l + q_c)^3 - w_r) + (w_o + w_l + q_c).w_r^2 - w_{o.next}).alpha

        fr::__mul(T1, q_mimc_selector_fft[i], T1); // T1 = (((w_o + w_l + q_c)^3 - w_r) + (w_o + w_l + q_c).w_r^2 - w_{o.next}).alpha).q_mimc
        fr::__add(circuit_state.quotient_large[i], T1, circuit_state.quotient_large[i]);
    ITERATE_OVER_DOMAIN_END;

    return fr::mul(alpha_base, fr::sqr(alpha_step));
}

void ProverMiMCWidget::compute_proof_elements(plonk_proof &proof, const fr::field_t &z)
{
    proof.q_mimc_coefficient_eval = q_mimc_coefficient.evaluate(z, n);
}

fr::field_t ProverMiMCWidget::compute_linear_contribution(const fr::field_t &alpha_base, const fr::field_t &alpha_step, const waffle::plonk_proof &proof, const evaluation_domain& domain, polynomial &r)
{
    fr::field_t mimc_T0 = fr::add(fr::add(proof.w_o_eval, proof.w_l_eval), proof.q_mimc_coefficient_eval);
    fr::field_t mimc_a = fr::sqr(mimc_T0);
    mimc_a = fr::mul(mimc_a, mimc_T0);
    mimc_a = fr::sub(mimc_a, proof.w_r_eval);
    fr::field_t mimc_term = fr::mul(fr::sub(fr::mul(fr::sqr(proof.w_r_eval), mimc_T0), proof.w_o_shifted_eval), alpha_step);
    mimc_term = fr::mul(fr::add(mimc_term, mimc_a), alpha_base);

    ITERATE_OVER_DOMAIN_START(domain);
        fr::field_t T0;
        fr::__mul(mimc_term, q_mimc_selector[i], T0);
        fr::__add(r[i], T0, r[i]);
    ITERATE_OVER_DOMAIN_END;
    return fr::mul(alpha_base, fr::sqr(alpha_step));
}

fr::field_t ProverMiMCWidget::compute_opening_poly_contribution(barretenberg::fr::field_t* poly, const evaluation_domain& domain, const fr::field_t &nu_base, const fr::field_t &nu_step)
{
    ITERATE_OVER_DOMAIN_START(domain);
        fr::field_t T0;
        fr::__mul(q_mimc_coefficient[i], nu_base, T0);
        fr::__add(poly[i], T0, poly[i]);
    ITERATE_OVER_DOMAIN_END;
    return fr::mul(nu_base, nu_step);
}

std::unique_ptr<VerifierBaseWidget> ProverMiMCWidget::compute_preprocessed_commitments(const evaluation_domain& domain, const ReferenceString &reference_string) const
{
    polynomial polys[2]{
        polynomial(q_mimc_coefficient, domain.size),
        polynomial(q_mimc_selector, domain.size)
    };

    for (size_t i = 0; i < 2; ++i)
    {
        polys[i].ifft(domain);
    }

    scalar_multiplication::multiplication_state mul_state[2]{
        { reference_string.monomials, polys[0].get_coefficients(), domain.size, {}},
        { reference_string.monomials, polys[1].get_coefficients(), domain.size, {}},
    };

    scalar_multiplication::batched_scalar_multiplications(mul_state, 2);
    std::vector<barretenberg::g1::affine_element> commitments;
    commitments.resize(2);

    for (size_t i = 0; i < 2; ++i)
    {
        g1::jacobian_to_affine(mul_state[i].output, commitments[i]);
    }
    std::unique_ptr<VerifierBaseWidget> result = std::make_unique<VerifierMiMCWidget>(commitments);
    return result;
}

void ProverMiMCWidget::reset(const evaluation_domain& domain)
{
    q_mimc_coefficient.fft(domain);
    q_mimc_selector.fft(domain);
}

// ###

VerifierMiMCWidget::VerifierMiMCWidget(std::vector<barretenberg::g1::affine_element> &instance_commitments) :
    VerifierBaseWidget(
        static_cast<size_t>(WidgetVersionControl::Dependencies::REQUIRES_W_O_SHIFTED),
        static_cast<size_t>(WidgetVersionControl::Features::HAS_MIMC_SELECTORS)
    )
{
    ASSERT(instance_commitments.size() == 2);
    instance = std::vector<g1::affine_element>{
        instance_commitments[0],
        instance_commitments[1]
    };
}

barretenberg::fr::field_t VerifierMiMCWidget::compute_batch_evaluation_contribution(barretenberg::fr::field_t &batch_eval, barretenberg::fr::field_t &nu_base, barretenberg::fr::field_t &nu_step, const plonk_proof &proof)
{
    fr::field_t T0;
    fr::__mul(proof.q_mimc_coefficient_eval, nu_base, T0);
    fr::__add(batch_eval, T0, batch_eval);
    return fr::mul(nu_base, nu_step);
}

VerifierBaseWidget::challenge_coefficients VerifierMiMCWidget::append_scalar_multiplication_inputs(
    const VerifierBaseWidget::challenge_coefficients &challenge,
    const waffle::plonk_proof &proof,
    std::vector<barretenberg::g1::affine_element> &points,
    std::vector<barretenberg::fr::field_t> &scalars)
{
    if (g1::on_curve(instance[0]))
    {
        points.push_back(instance[0]);
        scalars.push_back(challenge.nu_base);
    }

    fr::field_t mimc_T0 = fr::add(fr::add(proof.w_o_eval, proof.w_l_eval), proof.q_mimc_coefficient_eval);
    fr::field_t mimc_a = fr::sqr(mimc_T0);
    mimc_a = fr::mul(mimc_a, mimc_T0);
    mimc_a = fr::sub(mimc_a, proof.w_r_eval);
    fr::field_t q_mimc_term = fr::mul(fr::sub(fr::mul(fr::sqr(proof.w_r_eval), mimc_T0), proof.w_o_shifted_eval), challenge.alpha_step);
    q_mimc_term = fr::mul(fr::add(q_mimc_term, mimc_a), challenge.alpha_base);
    q_mimc_term = fr::mul(q_mimc_term, challenge.linear_nu);

    if (g1::on_curve(instance[1]))
    {
        points.push_back(instance[1]);
        scalars.push_back(q_mimc_term);
    }

    return VerifierBaseWidget::challenge_coefficients{
        fr::mul(challenge.alpha_base, fr::sqr(challenge.alpha_step)),
        challenge.alpha_step,
        fr::mul(challenge.nu_base, challenge.nu_step),
        challenge.nu_step,
        challenge.linear_nu
    };
}
}