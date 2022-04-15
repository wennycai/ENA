/*
 * RELIC is an Efficient LIbrary for Cryptography
 * Copyright (C) 2007-2019 RELIC Authors
 *
 * This file is part of RELIC. RELIC is legal property of its developers,
 * whose names are not listed here. Please refer to the COPYRIGHT file
 * for contact information.
 *
 * RELIC is free software; you can redistribute it and/or modify it under the
 * terms of the version 2.1 (or later) of the GNU Lesser General Public License
 * as published by the Free Software Foundation; or version 2.0 of the Apache
 * License as published by the Apache Software Foundation. See the LICENSE files
 * for more details.
 *
 * RELIC is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the LICENSE files for more details.
 *
 * You should have received a copy of the GNU Lesser General Public or the
 * Apache License along with RELIC. If not, see <https://www.gnu.org/licenses/>
 * or <https://www.apache.org/licenses/>.
 */

/**
 * @file
 *
 * Implementation of utilities for prime elliptic curves over cubic
 * extensions.
 *
 * @ingroup epx
 */

#include "relic_core.h"

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

int ep3_is_infty(ep3_t p) {
	return (fp3_is_zero(p->z) == 1);
}

void ep3_set_infty(ep3_t p) {
	fp3_zero(p->x);
	fp3_zero(p->y);
	fp3_zero(p->z);
	p->norm = 1;
}

void ep3_copy(ep3_t r, ep3_t p) {
	fp3_copy(r->x, p->x);
	fp3_copy(r->y, p->y);
	fp3_copy(r->z, p->z);
	r->norm = p->norm;
}

int ep3_cmp(ep3_t p, ep3_t q) {
    ep3_t r, s;
    int result = RLC_EQ;

    ep3_null(r);
    ep3_null(s);

    TRY {
        ep3_new(r);
        ep3_new(s);

        if ((!p->norm) && (!q->norm)) {
            /* If the two points are not normalized, it is faster to compare
             * x1 * z2^2 == x2 * z1^2 and y1 * z2^3 == y2 * z1^3. */
            fp3_sqr(r->z, p->z);
            fp3_sqr(s->z, q->z);
            fp3_mul(r->x, p->x, s->z);
            fp3_mul(s->x, q->x, r->z);
            fp3_mul(r->z, r->z, p->z);
            fp3_mul(s->z, s->z, q->z);
            fp3_mul(r->y, p->y, s->z);
            fp3_mul(s->y, q->y, r->z);
        } else {
            if (!p->norm) {
                ep3_norm(r, p);
            } else {
                ep3_copy(r, p);
            }

            if (!q->norm) {
                ep3_norm(s, q);
            } else {
                ep3_copy(s, q);
            }
        }

        if (fp3_cmp(r->x, s->x) != RLC_EQ) {
            result = RLC_NE;
        }

        if (fp3_cmp(r->y, s->y) != RLC_EQ) {
            result = RLC_NE;
        }
    } CATCH_ANY {
        THROW(ERR_CAUGHT);
    } FINALLY {
        ep3_free(r);
        ep3_free(s);
    }

    return result;
}

void ep3_rand(ep3_t p) {
	bn_t n, k;

	bn_null(k);
	bn_null(n);

	TRY {
		bn_new(k);
		bn_new(n);

		ep3_curve_get_ord(n);
		bn_rand_mod(k, n);

		ep3_mul_gen(p, k);
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
	FINALLY {
		bn_free(k);
		bn_free(n);
	}
}

void ep3_rhs(fp3_t rhs, ep3_t p) {
	fp3_t t0;
	fp3_t t1;

	fp3_null(t0);
	fp3_null(t1);

	TRY {
		fp3_new(t0);
		fp3_new(t1);

		/* t0 = x1^2. */
		fp3_sqr(t0, p->x);
		/* t1 = x1^3. */
		fp3_mul(t1, t0, p->x);

		ep3_curve_get_a(t0);
		fp3_mul(t0, p->x, t0);
		fp3_add(t1, t1, t0);

		ep3_curve_get_b(t0);
		fp3_add(t1, t1, t0);

		fp3_copy(rhs, t1);

	} CATCH_ANY {
		THROW(ERR_CAUGHT);
	} FINALLY {
		fp3_free(t0);
		fp3_free(t1);
	}
}


int ep3_is_valid(ep3_t p) {
	ep3_t t;
	int r = 0;

	ep3_null(t);

	TRY {
		ep3_new(t);

		ep3_norm(t, p);

		ep3_rhs(t->x, t);
		fp3_sqr(t->y, t->y);

		r = (fp3_cmp(t->x, t->y) == RLC_EQ) || ep3_is_infty(p);
	} CATCH_ANY {
		THROW(ERR_CAUGHT);
	} FINALLY {
		ep3_free(t);
	}
	return r;
}

void ep3_tab(ep3_t *t, ep3_t p, int w) {
	if (w > 2) {
		ep3_dbl(t[0], p);
#if defined(EP_MIXED)
		ep3_norm(t[0], t[0]);
#endif
		ep3_add(t[1], t[0], p);
		for (int i = 2; i < (1 << (w - 2)); i++) {
			ep3_add(t[i], t[i - 1], t[0]);
		}
#if defined(EP_MIXED)
		ep3_norm_sim(t + 1, t + 1, (1 << (w - 2)) - 1);
#endif
	}
	ep3_copy(t[0], p);
}

void ep3_print(ep3_t p) {
	fp3_print(p->x);
	fp3_print(p->y);
	fp3_print(p->z);
}

int ep3_size_bin(ep3_t a, int pack) {
	ep3_t t;
	int size = 0;

	ep3_null(t);

	if (ep3_is_infty(a)) {
		return 1;
	}

	TRY {
		ep3_new(t);

		ep3_norm(t, a);

		size = 1 + 2 * RLC_FP_BYTES;
		if (!pack) {
			size += 2 * RLC_FP_BYTES;
		}
	} CATCH_ANY {
		THROW(ERR_CAUGHT);
	} FINALLY {
		ep3_free(t);
	}

	return size;
}

void ep3_read_bin(ep3_t a, const uint8_t *bin, int len) {
	if (len == 1) {
		if (bin[0] == 0) {
			ep3_set_infty(a);
			return;
		} else {
			THROW(ERR_NO_BUFFER);
			return;
		}
	}

	if (len != (2 * RLC_FP_BYTES + 1) && len != (4 * RLC_FP_BYTES + 1)) {
		THROW(ERR_NO_BUFFER);
		return;
	}

	a->norm = 1;
	fp_set_dig(a->z[0], 1);
	fp_zero(a->z[1]);
	fp_zero(a->z[2]);
	fp3_read_bin(a->x, bin + 1, 2 * RLC_FP_BYTES);
	if (len == 2 * RLC_FP_BYTES + 1) {
		switch(bin[0]) {
			case 2:
				fp3_zero(a->y);
				break;
			case 3:
				fp3_zero(a->y);
				fp_set_bit(a->y[0], 0, 1);
				fp_zero(a->y[1]);
				fp_zero(a->y[2]);
				break;
			default:
				THROW(ERR_NO_VALID);
				break;
		}
		ep3_upk(a, a);
	}

	if (len == 4 * RLC_FP_BYTES + 1) {
		if (bin[0] == 4) {
			fp3_read_bin(a->y, bin + 2 * RLC_FP_BYTES + 1, 2 * RLC_FP_BYTES);
		} else {
			THROW(ERR_NO_VALID);
		}
	}
}

void ep3_write_bin(uint8_t *bin, int len, ep3_t a, int pack) {
	ep3_t t;

	ep3_null(t);

	if (ep3_is_infty(a)) {
		if (len < 1) {
			THROW(ERR_NO_BUFFER);
		} else {
			bin[0] = 0;
			return;
		}
	}

	TRY {
		ep3_new(t);

		ep3_norm(t, a);

		if (pack) {
			if (len < 2 * RLC_FP_BYTES + 1) {
				THROW(ERR_NO_BUFFER);
			} else {
				ep3_pck(t, t);
				bin[0] = 2 | fp_get_bit(t->y[0], 0);
				fp3_write_bin(bin + 1, 2 * RLC_FP_BYTES, t->x);
			}
		} else {
			if (len < 4 * RLC_FP_BYTES + 1) {
				THROW(ERR_NO_BUFFER);
			} else {
				bin[0] = 4;
				fp3_write_bin(bin + 1, 2 * RLC_FP_BYTES, t->x);
				fp3_write_bin(bin + 2 * RLC_FP_BYTES + 1, 2 * RLC_FP_BYTES, t->y);
			}
		}
	} CATCH_ANY {
		THROW(ERR_CAUGHT);
	} FINALLY {
		ep3_free(t);
	}
}
