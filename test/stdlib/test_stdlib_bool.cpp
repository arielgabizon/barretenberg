#include <gtest/gtest.h>

#include <barretenberg/waffle/composer/standard_composer.hpp>
#include <barretenberg/waffle/proof_system/preprocess.hpp>
#include <barretenberg/waffle/proof_system/widgets/arithmetic_widget.hpp>
#include <barretenberg/waffle/proof_system/prover/prover.hpp>
#include <barretenberg/waffle/proof_system/verifier/verifier.hpp>

#include <barretenberg/polynomials/polynomial_arithmetic.hpp>

#include <barretenberg/waffle/stdlib/common.hpp>
#include <barretenberg/waffle/stdlib/bool/bool.hpp>
#include <barretenberg/waffle/stdlib/field/field.hpp>

#include <memory>


using namespace barretenberg;
using namespace plonk;

typedef stdlib::field_t<waffle::StandardComposer> field_t;
typedef stdlib::bool_t<waffle::StandardComposer> bool_t;
typedef stdlib::witness_t<waffle::StandardComposer> witness_t;


bool get_value(bool_t &input)
{
    return static_cast<bool>(barretenberg::fr::from_montgomery_form(field_t(input).get_value()).data[0]);
}

TEST(stdlib_bool, test_basic_operations)
{
    waffle::StandardComposer composer = waffle::StandardComposer();
    bool_t a(&composer);
    bool_t b(&composer);
    a = stdlib::witness_t(&composer, barretenberg::fr::one);
    b = stdlib::witness_t(&composer, barretenberg::fr::zero);
    a = a ^ b; // a = 1
    b = !b;   // b = 1 (witness 0)
    bool_t d = (a == b); //
    d = false;           // d = 0
    bool_t e = a | d;   // e = 1 = a
    bool_t f = e ^ b;  // f = 0
    d = (!f) & a; // d = 1
    waffle::Prover prover = composer.preprocess();

    EXPECT_EQ(fr::eq(fr::from_montgomery_form(prover.w_l[0]), {{ 1, 0, 0, 0 }}), true);
    EXPECT_EQ(fr::eq(fr::from_montgomery_form(prover.w_r[0]), {{ 1, 0, 0, 0 }}), true);
    EXPECT_EQ(fr::eq(fr::from_montgomery_form(prover.w_o[0]), {{ 1, 0, 0, 0 }}), true);

    EXPECT_EQ(fr::eq(fr::from_montgomery_form(prover.w_l[1]), {{ 0, 0, 0, 0 }}), true);
    EXPECT_EQ(fr::eq(fr::from_montgomery_form(prover.w_r[1]), {{ 0, 0, 0, 0 }}), true);
    EXPECT_EQ(fr::eq(fr::from_montgomery_form(prover.w_o[1]), {{ 0, 0, 0, 0 }}), true);

    EXPECT_EQ(fr::eq(fr::from_montgomery_form(prover.w_l[2]), {{ 1, 0, 0, 0 }}), true);
    EXPECT_EQ(fr::eq(fr::from_montgomery_form(prover.w_r[2]), {{ 0, 0, 0, 0 }}), true);
    EXPECT_EQ(fr::eq(fr::from_montgomery_form(prover.w_o[2]), {{ 1, 0, 0, 0 }}), true);

    EXPECT_EQ(fr::eq(fr::from_montgomery_form(prover.w_l[3]), {{ 1, 0, 0, 0 }}), true);
    EXPECT_EQ(fr::eq(fr::from_montgomery_form(prover.w_r[3]), {{ 0, 0, 0, 0 }}), true);
    EXPECT_EQ(fr::eq(fr::from_montgomery_form(prover.w_o[3]), {{ 1, 0, 0, 0 }}), true);

    EXPECT_EQ(fr::eq(fr::from_montgomery_form(prover.w_l[4]), {{ 1, 0, 0, 0 }}), true);
    EXPECT_EQ(fr::eq(fr::from_montgomery_form(prover.w_r[4]), {{ 0, 0, 0, 0 }}), true);
    EXPECT_EQ(fr::eq(fr::from_montgomery_form(prover.w_o[4]), {{ 0, 0, 0, 0 }}), true);

    EXPECT_EQ(fr::eq(fr::from_montgomery_form(prover.w_l[5]), {{ 0, 0, 0, 0 }}), true);
    EXPECT_EQ(fr::eq(fr::from_montgomery_form(prover.w_r[5]), {{ 1, 0, 0, 0 }}), true);
    EXPECT_EQ(fr::eq(fr::from_montgomery_form(prover.w_o[5]), {{ 1, 0, 0, 0 }}), true);

    EXPECT_EQ(prover.n, 8UL);
}


TEST(stdlib_bool, xor)
{
    waffle::StandardComposer composer = waffle::StandardComposer();
    for (size_t i = 0; i < 32; ++i)
    {
        bool_t a = witness_t(&composer, (bool)(i % 2));
        bool_t b = witness_t(&composer, (bool)(i % 3 == 1));
        a ^ b;
    }
    waffle::Prover prover = composer.preprocess();
    waffle::Verifier verifier = waffle::preprocess(prover);

    waffle::plonk_proof proof = prover.construct_proof();

    bool result = verifier.verify_proof(proof);
    EXPECT_EQ(result, true);
}

TEST(stdlib_bool, xor_constants)
{
    waffle::StandardComposer composer = waffle::StandardComposer();
    for (size_t i = 0; i < 32; ++i)
    {
        bool_t a = witness_t(&composer, (bool)(i % 2));
        bool_t b = witness_t(&composer, (bool)(i % 3 == 1));
        a ^ b;
    }
    for (size_t i = 0; i < 32; ++i)
    {
        if (i % 2 == 0)
        {
            bool_t a = witness_t(&composer, (bool)(i % 2));
            bool_t b(&composer, (bool)(i % 3 == 1));
            a ^ b;
        }
        else
        {
            bool_t a(&composer, (bool)(i % 2));
            bool_t b = witness_t(&composer, (bool)(i % 3 == 1));
            a ^ b;
        }
    }
    waffle::Prover prover = composer.preprocess();
    waffle::Verifier verifier = waffle::preprocess(prover);

    waffle::plonk_proof proof = prover.construct_proof();

    bool result = verifier.verify_proof(proof);
    EXPECT_EQ(result, true);
}


TEST(stdlib_bool, xor_twin_constants)
{
    waffle::StandardComposer composer = waffle::StandardComposer();
    bool_t c;
    for (size_t i = 0; i < 32; ++i)
    {
        bool_t a(&composer, (i % 1) == 0);
        bool_t b(&composer, (i % 1) == 1);
        c = c ^ a ^ b;
    }
    c = c ^ bool_t(witness_t(&composer, true));
    for (size_t i = 0; i < 32; ++i)
    {
        if (i % 2 == 0)
        {
            bool_t a = witness_t(&composer, (bool)(i % 2));
            bool_t b(&composer, (bool)(i % 3 == 1));
            c = c ^ a ^ b;
        }
        else
        {
            bool_t a(&composer, (bool)(i % 2));
            bool_t b = witness_t(&composer, (bool)(i % 3 == 1));
            c = c ^ a ^ b;
        }
    }
    waffle::Prover prover = composer.preprocess();
    waffle::Verifier verifier = waffle::preprocess(prover);

    waffle::plonk_proof proof = prover.construct_proof();

    bool result = verifier.verify_proof(proof);
    EXPECT_EQ(result, true);
}

TEST(stdlib_bool, and)
{
    waffle::StandardComposer composer = waffle::StandardComposer();
    for (size_t i = 0; i < 32; ++i)
    {
        bool_t a = witness_t(&composer, (bool)(i % 1));
        bool_t b = witness_t(&composer, (bool)(i % 2 == 1));
        a & b;
    }
    waffle::Prover prover = composer.preprocess();
    waffle::Verifier verifier = waffle::preprocess(prover);

    waffle::plonk_proof proof = prover.construct_proof();

    bool result = verifier.verify_proof(proof);
    EXPECT_EQ(result, true);
}

TEST(stdlib_bool, and_constants)
{
    waffle::StandardComposer composer = waffle::StandardComposer();
    for (size_t i = 0; i < 32; ++i)
    {
        bool_t a = witness_t(&composer, (bool)(i % 2));
        bool_t b = witness_t(&composer, (bool)(i % 3 == 1));
        a & b;
    }
    for (size_t i = 0; i < 32; ++i)
    {
        if (i % 2 == 0)
        {
            bool_t a = witness_t(&composer, (bool)(i % 2));
            bool_t b(&composer, (bool)(i % 3 == 1));
            a & b;
        }
        else
        {
            bool_t a(&composer, (bool)(i % 2));
            bool_t b = witness_t(&composer, (bool)(i % 3 == 1));
            a & b;
        }
    }
    waffle::Prover prover = composer.preprocess();
    waffle::Verifier verifier = waffle::preprocess(prover);

    waffle::plonk_proof proof = prover.construct_proof();

    bool result = verifier.verify_proof(proof);
    EXPECT_EQ(result, true);
}

TEST(stdlib_bool, or)
{
    waffle::StandardComposer composer = waffle::StandardComposer();
    for (size_t i = 0; i < 32; ++i)
    {
        bool_t a = witness_t(&composer, (bool)(i % 2));
        bool_t b = witness_t(&composer, (bool)(i % 3 == 1));
        a | b;
    }
    waffle::Prover prover = composer.preprocess();
    waffle::Verifier verifier = waffle::preprocess(prover);

    waffle::plonk_proof proof = prover.construct_proof();

    bool result = verifier.verify_proof(proof);
    EXPECT_EQ(result, true);
}


TEST(stdlib_bool, or_constants)
{
    waffle::StandardComposer composer = waffle::StandardComposer();
    for (size_t i = 0; i < 32; ++i)
    {
        bool_t a = witness_t(&composer, (bool)(i % 2));
        bool_t b = witness_t(&composer, (bool)(i % 3 == 1));
        a | b;
    }
    for (size_t i = 0; i < 32; ++i)
    {
        if (i % 2 == 0)
        {
            bool_t a = witness_t(&composer, (bool)(i % 2));
            bool_t b(&composer, (bool)(i % 3 == 1));
            a | b;
        }
        else
        {
            bool_t a(&composer, (bool)(i % 2));
            bool_t b = witness_t(&composer, (bool)(i % 3 == 1));
            a | b;
        }
    }
    waffle::Prover prover = composer.preprocess();
    waffle::Verifier verifier = waffle::preprocess(prover);

    waffle::plonk_proof proof = prover.construct_proof();

    bool result = verifier.verify_proof(proof);
    EXPECT_EQ(result,true);
}

TEST(stdlib_bool, eq)
{
    waffle::StandardComposer composer = waffle::StandardComposer();
    bool a_alt[32];
    bool b_alt[32];
    bool c_alt[32];
    bool d_alt[32];
    for (size_t i = 0; i < 32; ++i)
    {
        if (i % 2 == 0)
        {
            a_alt[i] = bool(i % 2);
            b_alt[i] = false;
            c_alt[i] = a_alt[i] ^ b_alt[i];
            d_alt[i] = a_alt[i] == c_alt[i];
        }
        else
        {
            a_alt[i] = true;
            b_alt[i] = false;
            c_alt[i] = false;
            d_alt[i] = false;
        }
    }
    bool_t a[32];
    bool_t b[32];
    bool_t c[32];
    bool_t d[32];
    for (size_t i = 0; i < 32; ++i)
    {
        if (i % 2 == 0)
        {
            a[i] = witness_t(&composer, (bool)(i % 2));
            b[i] = witness_t(&composer, (bool)(0));
            c[i] = a[i] ^ b[i];
            d[i] = a[i] == c[i];
        }
        else
        {
            a[i] = witness_t(&composer, (bool)(1));
            b[i] = witness_t(&composer, (bool)(0));
            c[i] = a[i] & b[i];
            d[i] = a[i] == c[i];
        }
    }
    for (size_t i = 0; i < 32; ++i)
    {
        EXPECT_EQ(get_value(a[i]), a_alt[i]);
        EXPECT_EQ(get_value(b[i]), b_alt[i]);
        EXPECT_EQ(get_value(c[i]), c_alt[i]);
        EXPECT_EQ(get_value(d[i]), d_alt[i]);
    }
    waffle::Prover prover = composer.preprocess();
    waffle::Verifier verifier = waffle::preprocess(prover);

    waffle::plonk_proof proof = prover.construct_proof();

    bool result = verifier.verify_proof(proof);
    EXPECT_EQ(result, true);
}

TEST(stdlib_bool, test_simple_proof)
{
    waffle::StandardComposer composer = waffle::StandardComposer();
    bool_t a(&composer);
    bool_t b(&composer);
    a = stdlib::witness_t(&composer, barretenberg::fr::one);
    b = stdlib::witness_t(&composer, barretenberg::fr::zero);
    // bool_t c(&composer);
    a = a ^ b; // a = 1
    b = !b;   // b = 1 (witness 0)
    bool_t c = (a == b); // c = 1
    bool_t d(&composer); // d = ?
    d = false;           // d = 0
    bool_t e = a | d;   // e = 1 = a
    bool_t f = e ^ b;  // f = 0
    d = (!f) & a; // d = 1
    for (size_t i = 0; i < 64; ++i)
    {
        a = witness_t(&composer, (bool)(i % 2));
        b = witness_t(&composer, (bool)(i % 3 == 1));
        c = a ^ b;
        a = b ^ c;
        c = a;
        a = b;
        f = b;
    }
    waffle::Prover prover = composer.preprocess();
    waffle::Verifier verifier = waffle::preprocess(prover);

    waffle::plonk_proof proof = prover.construct_proof();

    bool result = verifier.verify_proof(proof);
    EXPECT_EQ(result, true);
}
