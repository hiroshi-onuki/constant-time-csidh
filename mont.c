
#include <assert.h>

#include "mont.h"
#include "u512.h"

void xDBLADD(proj *R, proj *S, proj const *P, proj const *Q, proj const *PQ, proj const *A24)
{
    fp tmp0, tmp1, tmp2;        //requires precomputation of A24=(A+2C:4C)

    fp_add3(&tmp0, &P->x, &P->z);
    fp_sub3(&tmp1, &P->x, &P->z);
    fp_sq2(&R->x, &tmp0);
    fp_sub3(&tmp2, &Q->x, &Q->z);
    fp_add3(&S->x, &Q->x, &Q->z);
    fp_mul2(&tmp0, &tmp2);
    fp_sq2(&R->z, &tmp1);
    fp_mul2(&tmp1, &S->x);
    fp_sub3(&tmp2, &R->x, &R->z);
    fp_mul2(&R->z, &A24->z);
    fp_mul2(&R->x, &R->z);
    fp_mul3(&S->x, &A24->x, &tmp2);
    fp_sub3(&S->z, &tmp0, &tmp1);
    fp_add2(&R->z, &S->x);
    fp_add3(&S->x, &tmp0, &tmp1);
    fp_mul2(&R->z, &tmp2);
    fp_sq1(&S->z);
    fp_sq1(&S->x);
    fp_mul2(&S->z, &PQ->x);
    fp_mul2(&S->x, &PQ->z);
}

void xDBL(proj *Q, proj const *A, proj const *P)
{
    fp a, b, c;
    fp_add3(&a, &P->x, &P->z);
    fp_sq1(&a);
    fp_sub3(&b, &P->x, &P->z);
    fp_sq1(&b);
    fp_sub3(&c, &a, &b);
    fp_add2(&b, &b); fp_add2(&b, &b); /* multiplication by 4 */
    fp_mul2(&b, &A->z);
    fp_mul3(&Q->x, &a, &b);
    fp_add3(&a, &A->z, &A->z); /* multiplication by 2 */
    fp_add2(&a, &A->x);
    fp_mul2(&a, &c);
    fp_add2(&a, &b);
    fp_mul3(&Q->z, &a, &c);
}

void xADD(proj *S, proj const *P, proj const *Q, proj const *PQ)
{
    fp a, b, c, d;
    fp_add3(&a, &P->x, &P->z);
    fp_sub3(&b, &P->x, &P->z);
    fp_add3(&c, &Q->x, &Q->z);
    fp_sub3(&d, &Q->x, &Q->z);
    fp_mul2(&a, &d);
    fp_mul2(&b, &c);
    fp_add3(&c, &a, &b);
    fp_sub3(&d, &a, &b);
    fp_sq1(&c);
    fp_sq1(&d);
    fp_mul3(&S->x, &PQ->z, &c);
    fp_mul3(&S->z, &PQ->x, &d);
}

/* Montgomery ladder. */
/* P must not be the unique point of order 2. */
/* not constant-time! */
/* factors are independent from the secret -> no constant-time ladder */
void xMUL(proj *Q, proj const *A, proj const *P, u512 const *k)
{
    proj R = *P;
    proj A24;
    const proj Pcopy = *P; /* in case Q = P */

    Q->x = fp_1;
    Q->z = fp_0;

    fp_add3(&A24.x, &A->z, &A->z);    //precomputation of A24=(A+2C:4C)
    fp_add3(&A24.z, &A24.x, &A24.x);
    fp_add2(&A24.x, &A->x);

    unsigned long i = 512;
    while (--i && !u512_bit(k, i));

    do {

        bool bit = u512_bit(k, i);

        if (bit) { proj T = *Q; *Q = R; R = T; } /* not constant-time */
        //fp_cswap(&Q->x, &R.x, bit);
        //fp_cswap(&Q->z, &R.z, bit);

        xDBLADD(Q, &R, Q, &R, &Pcopy, &A24);

        if (bit) { proj T = *Q; *Q = R; R = T; } /* not constant-time */
        //fp_cswap(&Q->x, &R.x, bit);
        //fp_cswap(&Q->z, &R.z, bit);

    } while (i--);
}

//simultaneous square-and-multiply, computes x^exp and y^exp 
void exp_by_squaring_(fp* x, fp* y, uint64_t exp)
{
	fp result1, result2;
	fp_set(&result1, 1);
	fp_set(&result2, 1);

    while (exp)
    {
        if (exp & 1){
          fp_mul2(&result1, x);
          fp_mul2(&result2, y);
        }
	
        fp_sq1(x);
        fp_sq1(y);
        exp >>= 1;
    }

    fp_set(x, 0);
    fp_add2(x, &result1);
    fp_set(y, 0);
    fp_add2(y, &result2);


}


/* computes the isogeny or dummy isogeny with kernel point K of order k */
/* returns the new curve coefficient A and the image of P for real isogenies*/
/* returns the old curve coefficient A and [k]P for dummy isogenies */
void xISOG(proj *A, proj *P, proj *Pd, proj *K, uint64_t k, int mask)
{
    assert (k >= 3);
    assert (k % 2 == 1);

    fp tmp0, tmp1, tmp2, tmp3, tmp4, Psum, Pdif, Pdsum, Pddif;
    proj Q, Qd, Aed, prod;
    proj Acopy = *A;
    proj Pdcopy = *Pd;

    fp_add3(&Aed.z, &A->z, &A->z);  //compute twisted Edwards curve coefficients
    fp_add3(&Aed.x, &A->x, &Aed.z);
    fp_sub3(&Aed.z, &A->x, &Aed.z);
   
    fp_add3(&Psum, &P->x, &P->z);   //precomputations
    fp_sub3(&Pdif, &P->x, &P->z);
    fp_add3(&Pdsum, &Pd->x, &Pd->z);   //precomputations
    fp_sub3(&Pddif, &Pd->x, &Pd->z);

    fp_sub3(&prod.x, &K->x, &K->z);
    fp_add3(&prod.z, &K->x, &K->z);
    
    fp_mul3(&tmp1, &prod.x, &Psum);
    fp_mul3(&tmp0, &prod.z, &Pdif);
    fp_add3(&Q.x, &tmp0, &tmp1);
    fp_sub3(&Q.z, &tmp0, &tmp1);
    fp_mul3(&tmp1, &prod.x, &Pdsum);   // for P'
    fp_mul3(&tmp0, &prod.z, &Pddif);
    fp_add3(&Qd.x, &tmp0, &tmp1);
    fp_sub3(&Qd.z, &tmp0, &tmp1);

    // CONSTANT TIME :
    proj *R = K;
    proj *S = P;

    // CONSTANT TIME :
    fp_cswap(&R->x, &S->x, mask);
    fp_cswap(&R->z, &S->z, mask);

    proj M[3] = {*R};  //K for real iso, P for dum iso
    xDBL(&M[1], A, R);

    for (uint64_t i = 1; i < k / 2; ++i) {

        if (i >= 2)
            xADD(&M[i % 3], &M[(i - 1) % 3], R, &M[(i - 2) % 3]);

	fp_sub3(&tmp1, &M[i % 3].x, &M[i % 3].z);
	fp_add3(&tmp0, &M[i % 3].x, &M[i % 3].z);
	fp_mul2(&prod.x, &tmp1);
	fp_mul2(&prod.z, &tmp0);
	fp_mul3(&tmp3, &tmp1, &Psum);  // for P
	fp_mul3(&tmp4, &tmp0, &Pdif);
	fp_add3(&tmp2, &tmp3, &tmp4);
	fp_mul2(&Q.x, &tmp2);
	fp_sub3(&tmp2, &tmp3, &tmp4);
	fp_mul2(&Q.z, &tmp2);
	fp_mul3(&tmp3, &tmp1, &Pdsum);  // for P'
	fp_mul3(&tmp4, &tmp0, &Pddif);
	fp_add3(&tmp2, &tmp3, &tmp4);
	fp_mul2(&Qd.x, &tmp2);
	fp_sub3(&tmp2, &tmp3, &tmp4);
	fp_mul2(&Qd.z, &tmp2);

    }

    if (k>3)
	xADD(&M[((k-1) / 2) % 3], &M[(((k-1) / 2)-1) % 3], R, &M[(((k-1) / 2)-2) % 3]);
    proj Pdummy = *R, Pcopy = *R;

    xADD(&Pdummy, &M[((k-1) / 2) % 3],  &M[(((k-1) / 2)-1) % 3], &Pcopy);

    // point evaluation
    fp_sq1(&Q.x);
    fp_sq1(&Q.z);
    fp_mul2(&P->x, &Q.x);
    fp_mul2(&P->z, &Q.z);
    fp_sq1(&Qd.x);
    fp_sq1(&Qd.z);
    fp_mul2(&Pd->x, &Qd.x);
    fp_mul2(&Pd->z, &Qd.z);


    //compute Aed.x^k, Aed.z^k
    exp_by_squaring_(&Aed.x, &Aed.z, k);

    //compute prod.x^8, prod.z^8
    fp_sq1(&prod.x);
    fp_sq1(&prod.x);
    fp_sq1(&prod.x);
    fp_sq1(&prod.z);
    fp_sq1(&prod.z);
    fp_sq1(&prod.z);

    //compute image curve parameters
    fp_mul2(&Aed.z, &prod.x);
    fp_mul2(&Aed.x, &prod.z);

    //compute Montgomery params
    fp_add3(&A->x, &Aed.x, &Aed.z);
    fp_sub3(&A->z, &Aed.x, &Aed.z);
    fp_add2(&A->x, &A->x);

    // CONSTANT TIME : swap back
    fp_cswap(&A->x, &Acopy.x, mask);
    fp_cswap(&A->z, &Acopy.z, mask);

    // CONSTANT TIME :
    fp_cswap(&P->x, &Pdummy.x, mask);
    fp_cswap(&P->z, &Pdummy.z, mask);
    fp_cswap(&Pd->x, &Pdcopy.x, mask);
    fp_cswap(&Pd->z, &Pdcopy.z, mask);

}

/* computes the last real/dummy isogeny per batch with kernel point K of order k */
/* real isogeny: returns the new curve coefficient A, no point evaluation */
/* dummy isogeny: returns the old curve coefficient A, no point evaluation */
void lastxISOG(proj *A, proj const *K, uint64_t k, int mask)
{
    assert (k >= 3);
    assert (k % 2 == 1);

    fp tmp0, tmp1;
    proj Aed, prod;
    proj Acopy = *A;

    fp_add3(&Aed.z, &A->z, &A->z);  //compute twisted Edwards curve coefficients
    fp_add3(&Aed.x, &A->x, &Aed.z);
    fp_sub3(&Aed.z, &A->x, &Aed.z);
   
    fp_sub3(&prod.x, &K->x, &K->z);
    fp_add3(&prod.z, &K->x, &K->z);

    proj M[3] = {*K};
    xDBL(&M[1], A, K);

    for (uint64_t i = 1; i < k / 2; ++i) {

        if (i >= 2)
            xADD(&M[i % 3], &M[(i - 1) % 3], K, &M[(i - 2) % 3]);

	fp_sub3(&tmp1, &M[i % 3].x, &M[i % 3].z);
    	fp_add3(&tmp0, &M[i % 3].x, &M[i % 3].z);
	fp_mul2(&prod.x, &tmp1);
        fp_mul2(&prod.z, &tmp0);

    }

    //compute Aed.x^k, Aed.z^k
    exp_by_squaring_(&Aed.x, &Aed.z, k);

    //compute prod.x^8, prod.z^8
    fp_sq1(&prod.x);
    fp_sq1(&prod.x);
    fp_sq1(&prod.x);
    fp_sq1(&prod.z);
    fp_sq1(&prod.z);
    fp_sq1(&prod.z);

    //compute image curve parameters
    fp_mul2(&Aed.z, &prod.x);
    fp_mul2(&Aed.x, &prod.z);

    //compute Montgomery params
    fp_add3(&A->x, &Aed.x, &Aed.z);
    fp_sub3(&A->z, &Aed.x, &Aed.z);
    fp_add2(&A->x, &A->x);

    // CONSTANT TIME : swap back
    fp_cswap(&A->x, &Acopy.x, mask);
    fp_cswap(&A->z, &Acopy.z, mask);

}


