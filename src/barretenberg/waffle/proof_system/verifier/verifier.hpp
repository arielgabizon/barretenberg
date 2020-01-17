#pragma once

#include "../../../types.hpp"

#include "../../reference_string/reference_string.hpp"
#include "../../waffle_types.hpp"

#include "../widgets/base_widget.hpp"

namespace waffle
{
template<size_t program_width>
class Verifier
{
  public:
    Verifier(const size_t subgroup_size = 0);
    Verifier(Verifier&& other);
    Verifier(const Verifier& other) = delete;
    Verifier& operator=(const Verifier& other) = delete;
    Verifier& operator=(Verifier&& other);

    bool verify_proof(const plonk_proof<program_width>& proof);

    ReferenceString reference_string;
    std::array<barretenberg::g1::affine_element, program_width> SIGMA;
    std::vector<std::unique_ptr<VerifierBaseWidget>> verifier_widgets;
    std::array<bool, program_width> needs_w_shifted = { 0 };
    size_t n;
    void update_needs_w_shifted();
};
} // namespace waffle
