#include "./evaluation_domain.hpp"

#include "../fields/fr.hpp"
#include "../assert.hpp"

#include "math.h"
#include "memory.h"

#ifndef NO_MULTITHREADING
#include "omp.h"
#endif

using namespace barretenberg;

namespace
{
// constexpr size_t MIN_THREADED_SIZE = 32;

size_t compute_num_threads(const size_t size)
{
#ifndef NO_MULTITHREADING
    size_t num_threads = static_cast<size_t>(omp_get_max_threads());
#else
    size_t num_threads = 1;
#endif
    if (size <= num_threads)// || size <= MIN_THREADED_SIZE)
    {
        num_threads = 1;
    }
    return num_threads;
}

void compute_lookup_table_single(const fr::field_t& input_root, const size_t size, fr::field_t* const roots, std::vector<fr::field_t*>& round_roots)
{
    const size_t num_rounds = static_cast<size_t>(log2(size));

    round_roots.emplace_back(&roots[0]);
    for (size_t i = 1; i < num_rounds - 1; ++i)
    {
        round_roots.emplace_back(round_roots.back() + (1UL << i));
    }

    for (size_t i = 0; i < num_rounds - 1; ++i)
    {
        const size_t m = 1UL << (i + 1);
        const fr::field_t round_root = fr::pow_small(input_root, (size / (2 * m)));
        fr::field_t* const current_round_roots = round_roots[i];
        fr::one(current_round_roots[0]);
        for (size_t j = 1; j < m; ++j)
        {
            fr::__mul_without_reduction(current_round_roots[j - 1], round_root, current_round_roots[j]);
        }
    }
}
}

evaluation_domain::evaluation_domain(const size_t domain_size):
    size(domain_size),
    num_threads(compute_num_threads(domain_size)),
    thread_size(domain_size / num_threads),
    log2_size(static_cast<size_t>(log2(size))),
    log2_thread_size(static_cast<size_t>(log2(thread_size))),
    log2_num_threads(static_cast<size_t>(log2(num_threads))),
    root(fr::get_root_of_unity(log2_size)),
    root_inverse(fr::invert(root)),
    domain(fr::to_montgomery_form({{size,0,0,0}})),
    domain_inverse(fr::invert(domain)),
    generator(fr::multiplicative_generator()),
    generator_inverse(fr::multiplicative_generator_inverse()),
    roots(nullptr)
{
    ASSERT((1UL << log2_size) == size);
    ASSERT((1UL << log2_thread_size) == thread_size);
    ASSERT((1UL << log2_num_threads) == num_threads);
}

evaluation_domain::evaluation_domain(const evaluation_domain& other):
    size(other.size),
    num_threads(compute_num_threads(other.size)),
    thread_size(other.size / thread_size),
    log2_size(static_cast<size_t>(log2(size))),
    log2_thread_size(static_cast<size_t>(log2(thread_size))),
    log2_num_threads(static_cast<size_t>(log2(num_threads))),
    root(fr::get_root_of_unity(log2_size)),
    root_inverse(fr::invert(root)),
    domain(fr::to_montgomery_form({{size,0,0,0}})),
    domain_inverse(fr::invert(domain)),
    generator(fr::multiplicative_generator()),
    generator_inverse(fr::multiplicative_generator_inverse())
{
    ASSERT((1UL << log2_size) == size);
    ASSERT((1UL << log2_thread_size) == thread_size);
    ASSERT((1UL << log2_num_threads) == num_threads);
    if (other.roots != nullptr)
    {
        size_t mem_size = sizeof(fr::field_t) * size * 2;
        roots = static_cast<fr::field_t*>(aligned_alloc(32, mem_size));
        memcpy(static_cast<void*>(roots), static_cast<void*>(other.roots), mem_size);
        round_roots.resize(log2_size - 1);
        inverse_round_roots.resize(log2_size - 1);
        round_roots[0] = &roots[0];
        inverse_round_roots[0] = &roots[size];
        for (size_t i = 1; i < log2_size - 1; ++i)
        {
            round_roots[i] = round_roots[i - 1] + (1UL << i);
            inverse_round_roots[i] = inverse_round_roots[i - 1] + (1UL << i);
        }
    }
    else
    {
        roots = nullptr;
    }
}

evaluation_domain::evaluation_domain(evaluation_domain&& other):
    size(other.size),
    num_threads(compute_num_threads(other.size)),
    thread_size(other.size / thread_size),
    log2_size(static_cast<size_t>(log2(size))),
    log2_thread_size(static_cast<size_t>(log2(thread_size))),
    log2_num_threads(static_cast<size_t>(log2(num_threads))),
    root(fr::get_root_of_unity(log2_size)),
    root_inverse(fr::invert(root)),
    domain(fr::to_montgomery_form({{size,0,0,0}})),
    domain_inverse(fr::invert(domain)),
    generator(fr::multiplicative_generator()),
    generator_inverse(fr::multiplicative_generator_inverse())
{
    if (roots != nullptr)
    {
        free(roots);
    }
    roots = other.roots;
    round_roots = std::move(other.round_roots);
    inverse_round_roots = std::move(other.inverse_round_roots);
    // std::copy(other.round_roots.begin(), other.round_roots.end(), round_roots.begin());
    // std::copy(other.inverse_round_roots.begin(), other.inverse_round_roots.end(), round_roots.begin());
}

evaluation_domain::~evaluation_domain()
{
    if (roots != nullptr)
    {
        free(roots);
    }
}

void evaluation_domain::compute_lookup_table()
{
    ASSERT(roots == nullptr);
    roots = (fr::field_t*)(aligned_alloc(32, sizeof(fr::field_t) * size * 2));
    compute_lookup_table_single(root, size, roots, round_roots);
    compute_lookup_table_single(root_inverse, size, &roots[size], inverse_round_roots);
}