/**
 * Copyright (C) 2015. Mario Rincon-Nigro.
 *
 * This file is a part of Flowie.
 *
 * Flowie is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Flowie is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Flowie.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <math/random.hpp>

#include <stdlib.h>

namespace flowie
{

std::vector<int> rand_range(int lo, int hi, unsigned int n) {
    std::vector<int> elements;
    
    for(unsigned int i = 0; i < n; i++) {
	int random = rand() % (hi - lo) + lo;
	elements.push_back(random);
    }

    return elements;
}

}

