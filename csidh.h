#ifndef CSIDH_H
#define CSIDH_H

#include "u512.h"
#include "fp.h"
#include "mont.h"

/* specific to p, should perhaps be somewhere else */
#define num_primes 74
fp invs_[9];

const unsigned primes[num_primes];

typedef struct private_key {
    int8_t e[num_primes];
} private_key;

typedef struct public_key {
    fp A; /* Montgomery coefficient: represents y^2 = x^3 + Ax^2 + x */
} public_key;

extern const public_key base;

void csidh_private(private_key *priv, const int8_t *max_exponent);
void action(public_key *out, public_key const *in, private_key const *priv,
		uint8_t num_intervals, int8_t const *max_exponent, unsigned int const num_isogenies, uint8_t const my);
bool csidh(public_key *out, public_key const *in, private_key const *priv,
		uint8_t const num_intervals, int8_t const *max_exponent, unsigned int const num_isogenies, uint8_t const my);
void elligator(proj *P, proj *Pd, const fp *A);
bool validate(public_key const *in);


int32_t lookup(size_t pos, int8_t const *priv);
uint32_t isequal(uint32_t a, uint32_t b);
void cmov(int8_t *r, const int8_t *a, uint32_t b);


#endif
