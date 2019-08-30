#pragma once

namespace fq
{
    struct field_t
    {
        alignas(32) uint64_t data[4];
    };
}

namespace fr
{
    struct field_t
    {
        alignas(32) uint64_t data[4];
    };

    struct field_wide_t
    {
        alignas(64) uint64_t data[8];
    };
}

namespace fq2
{
    struct fq2_t
    {
        fq::field_t c0;
        fq::field_t c1;
    };
}