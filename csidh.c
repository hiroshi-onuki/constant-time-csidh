#include <string.h>
#include <assert.h>

#include "csidh.h"
#include "rng.h"

const unsigned primes[num_primes] = {    359, 353, 349, 347, 337, 331, 317, 313, 311,
307, 293, 283, 281, 277, 271, 269, 263, 257, 251, 241, 239, 233, 229,
227, 223, 211, 199, 197, 193, 191, 181, 179, 173, 167, 163, 157, 151,
149, 139, 137, 131, 127, 113, 109, 107, 103, 101, 97, 89, 83, 79, 73,
71, 67, 61, 59, 53, 47, 43, 41, 37, 31, 29, 23, 19, 17, 13, 11, 7, 5, 3,
587, 373, 367 };

const u512 four_sqrt_p = { { 0x85e2579c786882cf, 0x4e3433657e18da95,
		0x850ae5507965a0b3, 0xa15bc4e676475964, } };

const public_key base = { 0 }; /* A = 0 */

/* get priv[pos] in constant time  */
int32_t lookup(size_t pos, int8_t const *priv)
{
	int b;
	int8_t r = priv[0];
	for(size_t i=1;i<num_primes;i++)
	{
		b = isequal(i, pos);
		//ISEQUAL(i, pos, b);
		//b = (uint8_t)(1-((-(i ^ pos)) >> 31));
		cmov(&r, &priv[i], b);
		//CMOV(&r, &priv[i], b);
	}
	return r;
}

/* check if a and b are equal in constant time  */
uint32_t isequal(uint32_t a, uint32_t b)
{
	//size_t i;
	uint32_t r = 0;
	unsigned char *ta = (unsigned char *)&a;
	unsigned char *tb = (unsigned char *)&b;
	r = (ta[0] ^ tb[0]) | (ta[1] ^ tb[1]) | (ta[2] ^ tb[2]) |  (ta[3] ^ tb[3]);
	r = (-r);
	r = r >> 31;
	return (int)(1-r);
}


/* decision bit b has to be either 0 or 1 */
void cmov(int8_t *r, const int8_t *a, uint32_t b)
{
	uint32_t t;
	b = -b; /* Now b is either 0 or 0xffffffff */
	t = (*r ^ *a) & b;
	*r ^= t;
}


void csidh_private(private_key *priv, const int8_t *max_exponent) {
	memset(&priv->e, 0, sizeof(priv->e));
	for (size_t i = 0; i < num_primes;) {
		int8_t buf[64];
		randombytes(buf, sizeof(buf));
		for (size_t j = 0; j < sizeof(buf); ++j) {
			if (buf[j] <= max_exponent[i] && buf[j] >= -max_exponent[i]) {
				priv->e[i] = lookup(j, buf);
				if (++i >= num_primes)
					break;
			}
		}
	}
}

/* compute [(p+1)/l] P for all l in our list of primes. */
/* divide and conquer is much faster than doing it naively,
 * but uses more memory. */
static void cofactor_multiples(proj *P, const proj *A, size_t lower,
		size_t upper) {
	assert(lower < upper);

	if (upper - lower == 1)
		return;

	size_t mid = lower + (upper - lower + 1) / 2;

	u512 cl = u512_1, cu = u512_1;
	for (size_t i = lower; i < mid; ++i)
		u512_mul3_64(&cu, &cu, primes[i]);
	for (size_t i = mid; i < upper; ++i)
		u512_mul3_64(&cl, &cl, primes[i]);

	xMUL(&P[mid], A, &P[lower], &cu);
	xMUL(&P[lower], A, &P[lower], &cl);

	cofactor_multiples(P, A, lower, mid);
	cofactor_multiples(P, A, mid, upper);
}

/* never accepts invalid keys. */
bool validate(public_key const *in) {
	const proj A = { in->A, fp_1 };

	do {

		proj P[num_primes];
		fp_random(&P->x);
		P->z = fp_1;

		/* maximal 2-power in p+1 */
		xDBL(P, &A, P);
		xDBL(P, &A, P);

		cofactor_multiples(P, &A, 0, num_primes);

		u512 order = u512_1;

		for (size_t i = num_primes - 1; i < num_primes; --i) {

			/* we only gain information if [(p+1)/l] P is non-zero */
			if (memcmp(&P[i].z, &fp_0, sizeof(fp))) {

				u512 tmp;
				u512_set(&tmp, primes[i]);
				xMUL(&P[i], &A, &P[i], &tmp);

				if (memcmp(&P[i].z, &fp_0, sizeof(fp)))
					/* P does not have order dividing p+1. */
					return false;

				u512_mul3_64(&order, &order, primes[i]);

				if (u512_sub3(&tmp, &four_sqrt_p, &order)) /* returns borrow */
					/* order > 4 sqrt(p), hence definitely supersingular */
					return true;
			}
		}

		/* P didn't have big enough order to prove supersingularity. */
	} while (1);
}

/* compute x^3 + Ax^2 + x */
static void montgomery_rhs(fp *rhs, fp const *A, fp const *x) {
	fp tmp;
	*rhs = *x;
	fp_sq1(rhs);
	fp_mul3(&tmp, A, x);
	fp_add2(rhs, &tmp);
	fp_add2(rhs, &fp_1);
	fp_mul2(rhs, x);
}


/* generates curve points */
void elligator(proj *P, proj *Pd, const fp *A) {

	fp u2, u2m1, tmp, rhs;
	bool issquare;

	fp_random(&u2);
	fp_sq1(&u2);				// u^2
	fp_sub3(&u2m1, &u2, &fp_1); // u^2 - 1

	fp_sq2(&tmp, &u2m1);	// (u^2 - 1)^2
	fp_sq2(&rhs, A);		// A^2
	fp_mul2(&rhs, &u2);		// A^2u^2
	fp_add2(&rhs, &tmp);	// A^2u^2 + u(u^2 - 1)^2
	fp_mul2(&rhs, A);		// (A^2u^2 + u(u^2 - 1)^2)A
	fp_mul2(&rhs, &u2m1);	// (A^2u^2 + u(u^2 - 1)^2)A(u^2 - 1)

	fp_set(&P->x, 0);
	fp_add2(&P->x, A);
	fp_set(&P->z, 0);
	fp_add2(&P->z, &u2m1);
	fp_set(&Pd->x, 0);
	fp_sub2(&Pd->x, A);
	fp_mul2(&Pd->x, &u2);
	fp_set(&Pd->z, 0);
	fp_add2(&Pd->z, &u2m1);

	// swap (x:z) and (xd:zd)
	issquare = fp_issquare(&rhs);
	fp_cswap(&P->x, &Pd->x, !issquare);
	fp_cswap(&P->z, &Pd->z, !issquare);
}

/* constant-time. */
void action(public_key *out, public_key const *in, private_key const *priv,
		uint8_t num_batches, int8_t const *max_exponent, unsigned int const num_isogenies, uint8_t const my) {

	//factors k for different batches
	u512 k[3] = {{ .c = {0x1b5933af628d005c, 0x9d4af02b1d7b7f56, 0x8977a8435092262a, 0xb86302ff54a37ca2, 0xd6e09db2af04d095, 0x5c73f, 0x0, 0x0} },
				{.c = {0xd97b8b6bc6f6be1c, 0x315872c44ea6e448, 0x1aae7c54fd380c86, 0x237ec4cf2da454a2, 0x3733f9e3d9fea1b4, 0x1fdc0e, 0x0, 0x0} },
				{.c = {0x629ea97b02169a84, 0xc4b9616a12d48d22, 0x492a10278ad7b45a, 0xc44ac4dce55b87f8, 0x9e12876886632d6e, 0xe0c0c5, 0x0, 0x0} }};


	u512 p_order = { .c = {0x24403b2c196b9323, 0x8a8759a31723c208, 0xb4a93a543937992b, 0xcdd1f791dc7eb773, 0xff470bd36fd7823b, 0xfbcf1fc39d553409, 0x9478a78dd697be5c, 0x0ed9b5fb0f251816}};

	int8_t ec = 0, m = 0;
	uint8_t count = 0;
	uint8_t elligator_index = 0;
	uint8_t last_iso[3], bc, ss;
	proj P, Pd, K;
	u512 cof, l;
	bool finished[num_primes] = {0};
	int8_t e[num_primes] = {0};
	int8_t counter[num_primes] = {0};
	int8_t s, ps;
	unsigned int isog_counter = 0;

	//index for skipping point evaluations
	last_iso[0] = 72;
	last_iso[1] = 73;
	last_iso[2] = 71;

	memcpy(e, priv->e, sizeof(priv->e));

	memcpy(counter, max_exponent, sizeof(counter));

	proj A = { in->A, fp_1 };

	while (isog_counter < num_isogenies) {
		m = (m + 1) % num_batches;
		
		if(count == my*num_batches) {  //merge the batches after my rounds
			m = 0;
			last_iso[0] = 73;    //doesn't skip point evaluations anymore after merging batches
			u512_set(&k[m], 4);  //recompute factor k
			num_batches = 1;

			// no need for constant-time, depends only on randomness
			for (uint8_t i = 0; i < num_primes; i++) {
				if(counter[i]==0) {
					u512_mul3_64(&k[m], &k[m], primes[i]);
				}
			}
		}

		assert(!memcmp(&A.z, &fp_1, sizeof(fp)));

		if(memcmp(&A.x, &fp_0, sizeof(fp))) {
			elligator(&P, &Pd, &A.x);
		} else {
			fp_enc(&P.x, &p_order); // point of full order on E_a with a=0
			fp_sub3(&Pd.x, &fp_0, &P.x);
			P.z = fp_1;
			Pd.z = fp_1;
		}

		xMUL(&P, &A, &P, &k[m]);
		xMUL(&Pd, &A, &Pd, &k[m]);
		ps = 1;

		for (uint8_t i = m; i < num_primes; i = i + num_batches) {
			if(finished[i] == true) {  //depends only on randomness
				continue;
			} else {
				cof = u512_1;
				for (uint8_t j = i + num_batches; j < num_primes; j = j + num_batches) {
					if (finished[j] == false)  //depends only on randomness
						u512_mul3_64(&cof, &cof, primes[j]);
				}

				ec = lookup(i, e);  //check in constant-time if normal or dummy isogeny must be computed
				bc = isequal(ec, 0);
				s = (uint8_t)ec >> 7;
				ss = !isequal(s, ps);
				ps = s;

				fp_cswap(&P.x, &Pd.x, ss);
				fp_cswap(&P.z, &Pd.z, ss);
				xMUL(&K, &A, &P, &cof);
				u512_set(&l, primes[i]);
				xMUL(&Pd, &A, &Pd, &l);

				if (memcmp(&K.z, &fp_0, sizeof(fp))) {  //depends only on randomness

						if (i == last_iso[m])
						{
							lastxISOG(&A, &K, primes[i], bc);	// doesn't compute the images of points
						}
						else
						{
							xISOG(&A, &P, &Pd, &K, primes[i], bc);
						}

						e[i] = ec - (1 ^ bc) + (s << 1);
						counter[i] = counter[i] - 1;
						isog_counter = isog_counter + 1;



				}

			}


			if(counter[i]==0) {   //depends only on randomness
				finished[i] = true;
				u512_mul3_64(&k[m], &k[m], primes[i]);
			}
		}




		fp_inv(&A.z);
		fp_mul2(&A.x, &A.z);
		A.z = fp_1;
		count = count + 1;

	}

	out->A = A.x;

}


/* includes public-key validation. */
bool csidh(public_key *out, public_key const *in, private_key const *priv,
		uint8_t const num_batches, int8_t const *max_exponent, unsigned int const num_isogenies, uint8_t const my) {
	if (!validate(in)) {
		fp_random(&out->A);
		return false;
	}
	action(out, in, priv, num_batches, max_exponent, num_isogenies, my);

	return true;
}


