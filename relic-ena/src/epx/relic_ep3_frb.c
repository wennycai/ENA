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
 * Implementation of frobenius action on prime elliptic curves over
 * cubic extensions.
 *
 * @ingroup epx
 */

#include "relic_core.h"

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/
//case 4 has question
void ep3_frb(ep3_t r, ep3_t p, int i) {
	ep3_copy(r, p);

	switch (i) {
		case 1:
			fp3_frb(r->x, r->x, 1);
			fp3_frb(r->y, r->y, 1);
			fp3_frb(r->z, r->z, 1);
			fp3_mul_frb(r->x, r->x, 1,1, 2);
			fp3_mul_frb(r->y, r->y, 1,1, 3);
			break;
		case 2:
			fp3_frb(r->x, r->x, 2);
			fp3_frb(r->y, r->y, 2);
			fp3_frb(r->z, r->z, 2);
			fp3_mul_frb(r->x, r->x, 1, 2, 2);
			fp3_mul_frb(r->y, r->y, 1, 2, 3);
			break;
		case 3:
			fp3_mul_frb(r->x, r->x, 1, 3, 2);
			fp3_neg(r->y, r->y);
			break;
		case 4:
				fp3_frb(r->x, r->x, 1);
				fp3_frb(r->y, r->y, 1);
				fp3_mul_frb(r->x, r->x, 1,4, 2);
				fp3_mul_frb(r->y, r->y, 1,4, 3);
			break;
	}
}
