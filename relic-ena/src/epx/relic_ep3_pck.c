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
 * Implementation of point compression on prime elliptic curves over cubic
 * extensions.
 *
 * @ingroup ep
 */

#include "relic_core.h"

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/
//should fix,为什么可以用最低位判断
void ep3_pck(ep3_t r, ep3_t p) {
	int b = fp_get_bit(p->y[0], 0);
	fp3_copy(r->x, p->x);
	fp3_zero(r->y);
	fp_set_bit(r->y[0], 0, b);
	fp_zero(r->y[1]);
	fp_zero(r->y[2]);
	fp_set_dig(r->z[0], 1);
	fp_zero(r->z[1]);
	fp_zero(r->z[2]);
	r->norm = 1;
}

int ep3_upk(ep3_t r, ep3_t p) {
	fp3_t t;
	int result = 0;

	fp3_null(t);

	TRY {
		fp3_new(t);

		ep3_rhs(t, p);//t=x1^3 + a * x1 + b

		/* t0 = sqrt(x1^3 + a * x1 + b). */
		result = fp3_srt(t, t);

		if (result) {
			/* Verify if least significant bit of the result matches the
			 * compressed y-coordinate. */
			if (fp_get_bit(t[0], 0) != fp_get_bit(p->y[0], 0)) {
				fp3_neg(t, t);
			}
			fp3_copy(r->x, p->x);
			fp3_copy(r->y, t);
			fp_set_dig(r->z[0], 1);
			fp_zero(r->z[1]);
			fp_zero(r->z[2]);
			r->norm = 1;
		}
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
	FINALLY {
		fp3_free(t);
	}
	return result;
}
