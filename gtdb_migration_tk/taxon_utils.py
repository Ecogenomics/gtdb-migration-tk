###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import os
import re

def canonical_strain_id(strain_id):
    """Generate canonical strain identifier."""

    strain_id = strain_id.strip()
    strain_id = re.sub(r'\(.+\)', ' ', strain_id)
    strain_id = ' '.join(strain_id.split())
    strain_id = re.sub('[\W_]+', '', strain_id).upper()
    
    return strain_id
    