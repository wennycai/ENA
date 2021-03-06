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
 * Implementation of point negation on elliptic prime curves over cubic
 * extensions.
 *
 * @ingroup epx
 */

#include "relic_core.h"

/*============================================================================*/
	/* Public definitions                                                         */
/*============================================================================*/

#if EP_ADD == BASIC || !defined(STRIP)

void ep3_neg_basic(ep3_t r, ep3_t p) {
	if (ep3_is_infty(p)) {
		ep3_set_infty(r);
		return;
	}

	if (r != p) {
		fp3_copy(r->x, p->x);
		fp3_copy(r->z, p->z);
	}

	fp3_neg(r->y, p->y);

	r->norm = 1;
}

#endif

#if EP_ADD == PROJC || !defined(STRIP)

void ep3_neg_projc(ep3_t r, ep3_t p) {
	if (ep3_is_infty(p)) {
		ep3_set_infty(r);
		return;
	}

	if (r != p) {
		fp3_copy(r->x, p->x);
		fp3_copy(r->z, p->z);
	}

	fp3_neg(r->y, p->y);

	r->norm = p->norm;
}

#endif
