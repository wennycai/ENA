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
 * Implementation of point normalization for points used in pairing computation.
 *
 * @ingroup pp
 */

#include "relic_core.h"
#include "relic_md.h"
#include "relic_pp.h"
#include "relic_conf.h"
#include "relic_fp_low.h"

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

#if EP_ADD == PROJC || !defined(STRIP)

/**
 * Normalizes a point represented in projective coordinates.
 *
 * @param r			- the result.
 * @param p			- the point to normalize.
 */
void pp_norm_imp(ep2_t r, ep2_t p) {
	fp2_inv(r->z, p->z);
	fp2_mul(r->x, p->x, r->z);
	fp2_mul(r->y, p->y, r->z);
	fp_zero(r->z[0]);
	fp_zero(r->z[1]);
	fp_set_dig(r->z[0], 1);

	r->norm = 1;
}

#endif /* EP_ADD == PROJC */

#if EP_ADD == PROJC || !defined(STRIP)

/**
 * Normalizes a point represented in projective coordinates.
 *
 * @param r			- the result.
 * @param p			- the point to normalize.
 */
void pp_norm_imp3(ep3_t r, ep3_t p) {
	fp3_inv(r->z, p->z);
	fp3_mul(r->x, p->x, r->z);
	fp3_mul(r->y, p->y, r->z);
	fp_zero(r->z[1]);
	fp_zero(r->z[2]);
	fp_set_dig(r->z[0], 1);

	r->norm = 1;
}

#endif /* EP_ADD == PROJC */

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

void pp_norm_k2(ep_t r, ep_t p) {
	ep_norm(r, p);
}

void pp_norm_k12(ep2_t r, ep2_t p) {
	if (ep2_is_infty(p)) {
		ep2_set_infty(r);
		return;
	}

	if (p->norm) {
		/* If the point is represented in affine coordinates, we just copy it. */
		ep2_copy(r, p);
	}
#if EP_ADD == PROJC || !defined(STRIP)
	pp_norm_imp(r, p);
#endif
}
void pp_norm_k18(ep3_t r, ep3_t p) {
	if (ep3_is_infty(p)) {
		ep3_set_infty(r);
		return;
	}

	if (p->norm) {
		/* If the point is represented in affine coordinates, we just copy it. */
		ep3_copy(r, p);
	}
#if EP_ADD == PROJC || !defined(STRIP)
	pp_norm_imp3(r, p);
#endif
}
