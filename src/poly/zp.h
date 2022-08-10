#pragma once

#ifndef _ZP_H
#define _ZP_H

#include <x86intrin.h> // simd, support 128-bit register
#include "NTL/ZZ.h"
#include "NTL/ZZ_p.h"
#include "gmp.h"

#include <sstream>
#include <string>
#include <vector>

using namespace NTL;

class ZpLongEle {

public:
    static const unsigned long long p = 2305843009213693951;
    unsigned long long ele;

    ZpLongEle() { ele = 0; };
    ZpLongEle(unsigned long long e) { ele = e; };
    
    // check equal
    bool operator!=(const ZpLongEle& other) { return other.ele != ele; };

    // assign
    inline ZpLongEle& operator=(const ZpLongEle& other) { ele = other.ele; return *this; };
    
    // add
    ZpLongEle operator+(const ZpLongEle& other) { 
        ZpLongEle ans;
        ans.ele = (ele + other.ele) % p;
        return ans;
    };

    // minus
    ZpLongEle operator-(const ZpLongEle& other) { 
        ZpLongEle ans;
        ans.ele = (ele - other.ele) % p;
        return ans;
    };

    ZpLongEle operator*(const ZpLongEle& other) {
        ZpLongEle ans;

        unsigned long long high;
        // store the low 64-bits of the result in low, and store the high 64-bits in hi
        unsigned long long low = _mulx_u64(ele, other.ele, &high);
        
        unsigned long long low61 = (low & p);
        unsigned long long res = (low61 + (low >> 61) + (high << 3)) % p;

        ans.ele = res;

        return ans;
    }   

    // div (using exgcd)
    // \frac{this}{other} = this \cdot ( other^{-1} mod p )
    ZpLongEle operator/(const ZpLongEle& other) {
        ZpLongEle ans;

        // c = p
        // a = this
        // b = other
        // inv = b^{-1} mod p
        mpz_t c, inv, a, b;
        mpz_init_set_str(c, "2305843009213693951", 10);
        mpz_init(inv);
        mpz_init(a);
        mpz_init(b);

        mpz_set_ui(a, ele);
        mpz_set_ui(b, other.ele);

        mpz_invert(inv, b, c);

        mpz_mul(inv, inv, a);
        mpz_mod(inv, inv, c);
        
        ans.ele = mpz_get_ui(inv);

        return ans;
    }

    ZpLongEle& operator+=(const ZpLongEle& other) {
        ele = (ele + other.ele) % p;
        return *this;
    }

    ZpLongEle& operator*=(const ZpLongEle& other) {

        unsigned long long high;
        // store the low 64-bits of the result in low, and store the high 64-bits in hi
        unsigned long long low = _mulx_u64(ele, other.ele, &high);
        
        unsigned long long low61 = (low & p);
        unsigned long long res = (low61 + (low >> 61) + (high << 3)) % p;

        ele = res;

        return *this;
    }
};

#endif