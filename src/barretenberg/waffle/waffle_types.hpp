#pragma once

#include "../curves/bn254/fq.hpp"
#include "../curves/bn254/fr.hpp"
#include "../curves/bn254/g1.hpp"

namespace waffle
{
struct plonk_challenges
{
    barretenberg::fr::field_t beta;
    barretenberg::fr::field_t gamma;
    barretenberg::fr::field_t alpha;
    barretenberg::fr::field_t z;
    barretenberg::fr::field_t nu;
};

void make_z_pow(const barretenberg::fr::field_t z, const size_t length, const uint64_t n, barretenberg::fr::field_t* z_pow)
{
    for (size_t i = 0; i < length; ++i) {
        barretenberg::fr::__pow_small(z, n * (i + 1), z_pow[i]);
    }
}

template<size_t program_width>
struct plonk_proof
{

    // Kate polynomial commitments required for a proof of knowledge
    std::array<barretenberg::g1::affine_element, program_width> W;
//
//    barretenberg::g1::affine_element W_L;
//    barretenberg::g1::affine_element W_R;
//    barretenberg::g1::affine_element W_O;
    barretenberg::g1::affine_element Z;

    std::array<barretenberg::g1::affine_element, program_width> T;
//
//    barretenberg::g1::affine_element T_LO;
//    barretenberg::g1::affine_element T_MID;
//    barretenberg::g1::affine_element T_HI;
    barretenberg::g1::affine_element PI_Z;
    barretenberg::g1::affine_element PI_Z_OMEGA;


    std::array<barretenberg::fr::field_t, program_width> w_eval;


    std::array<barretenberg::fr::field_t, program_width-1> sigma_eval;


//            barretenberg::fr::field_t sigma_1_eval;
//    barretenberg::fr::field_t sigma_2_eval;
    barretenberg::fr::field_t z_1_shifted_eval;
    barretenberg::fr::field_t linear_eval;


    std::array<barretenberg::fr::field_t, program_width> w_shifted_eval;
//
//    barretenberg::fr::field_t w_l_shifted_eval;
//    barretenberg::fr::field_t w_r_shifted_eval;
//    barretenberg::fr::field_t w_o_shifted_eval;
    barretenberg::fr::field_t q_c_eval;
    barretenberg::fr::field_t q_mimc_coefficient_eval;
    std::vector<barretenberg::fr::field_t> custom_gate_evaluations;
};
}