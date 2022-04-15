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
 * Implementation of hashing to a prime elliptic curve over a quadratic
 * extension.
 *
 * @ingroup epx
 */

#include "relic_core.h"
#include "relic_md.h"

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

void ep3_map(ep3_t p, const uint8_t *msg, int len) {
	bn_t x;
	fp3_t t0;
	uint8_t digest[RLC_MD_LEN];

	bn_null(x);
	fp3_null(t0);

	TRY {
		bn_new(x);
		fp3_new(t0);

		md_map(digest, msg, len);
		bn_read_bin(x, digest, RLC_MIN(RLC_FP_BYTES, RLC_MD_LEN));

		fp_prime_conv(p->x[0], x);
		fp_zero(p->x[1]);
		fp_set_dig(p->z[0], 1);
		fp_zero(p->z[1]);

		while (1) {
			ep3_rhs(t0, p);

			if (fp3_srt(p->y, t0)) {
				p->norm = 1;
				break;
			}

			fp_add_dig(p->x[0], p->x[0], 1);
		}

		switch (ep_curve_is_pairf()) {
			default:
				/* Now, multiply by cofactor to get the correct group. */
				ep3_curve_get_cof(x);
				if (bn_bits(x) < RLC_DIG) {
					ep3_mul_dig(p, p, x->dp[0]);
					if (bn_sign(x) == RLC_NEG) {
						ep3_neg(p, p);
					}
				} else {
					ep3_mul(p, p, x);
				}
				break;
		}
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
	FINALLY {
		bn_free(x);
		fp3_free(t0);
	}
}
