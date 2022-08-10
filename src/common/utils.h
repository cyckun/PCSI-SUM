#pragma once

#ifndef _UTILS_H
#define _UTILS_H

#include <vector>
#include <cinttypes>
#include <algorithm>
#include <random>
#include <unordered_set>

#include "HashingTables/common/hashing.h"
#include "constants.h"

namespace PCSI {

// generate n elements of length bitlen under the specific seed
std::vector<uint64_t> random_ele_generator(const uint32_t n, const uint32_t bitlen, const uint32_t seed) {
    std::vector<uint64_t> eles;
    eles.reserve(n);

    // initial prng
    std::mt19937 engine(seed);
    while (eles.size() != n) {
        // select element uniformly
        std::uniform_int_distribution<uint64_t> dis(0, (1ull << bitlen) - 1);
        while (eles.size() != n) {
            eles.push_back(dis(engine));
        }

        // using unordered_set to distinct
        std::unordered_set<uint64_t> tmp;
        for (auto ele : eles) {
            tmp.insert(ele);
        }
        eles.assign(tmp.begin(), tmp.end());
    }
    std::sort(eles.begin(), eles.end());
    // process cuckoo hashing
    for (auto &ele : eles) {
        ele = ENCRYPTO::HashingTable::ElementToHash(ele) & _61_mask;
    }
    return eles;
}

}

#endif