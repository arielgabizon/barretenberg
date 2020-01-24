#pragma once

#include <map>

#include "../curves/bn254/fq.hpp"
#include "../curves/bn254/fr.hpp"
#include "../curves/bn254/g1.hpp"

namespace waffle
{
struct g1_challenges
{
    barretenberg::fr::field_t beta;
    barretenberg::fr::field_t gamma;
    barretenberg::fr::field_t alpha;
    barretenberg::fr::field_t z;
    barretenberg::fr::field_t nu;
};

struct plonk_proof
{

    std::map<std::string, std::pair<barretenberg::g1::affine_element, size_t> > g1_commitments;
    // Kate polynomial commitments required for a proof of knowledge
    std::vector<barretenberg::g1::affine_element> W;
//
//    barretenberg::g1::affine_element W_L;
//    barretenberg::g1::affine_element W_R;
//    barretenberg::g1::affine_element W_O;
    barretenberg::g1::affine_element Z;

    std::vector<barretenberg::g1::affine_element> T;
//
//    barretenberg::g1::affine_element T[0];
//    barretenberg::g1::affine_element T[1];
//    barretenberg::g1::affine_element T[2];
    barretenberg::g1::affine_element PI_Z;
    barretenberg::g1::affine_element PI_Z_OMEGA;


    std::vector<barretenberg::fr::field_t> w;

    barretenberg::fr::field_t sigma_1_eval;
    barretenberg::fr::field_t sigma_2_eval;
    barretenberg::fr::field_t z_1_shifted_eval;
    barretenberg::fr::field_t linear_eval;

    std::vector<barretenberg::fr::field_t> w_eval;
    std::vector<barretenberg::fr::field_t> w_shifted_eval;
//
//    barretenberg::fr::field_t w_l_shifted_eval;
//    barretenberg::fr::field_t w_r_shifted_eval;
//    barretenberg::fr::field_t w_o_shifted_eval;
    barretenberg::fr::field_t q_c_eval;
    barretenberg::fr::field_t q_mimc_coefficient_eval;
    std::vector<barretenberg::fr::field_t> custom_gate_evaluations;
};
}