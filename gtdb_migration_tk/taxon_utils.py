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
    
    if strain_id.startswith('strain '):
        strain_id = strain_id.replace('strain ', '', 1)
    if strain_id.startswith('Strain '):
        strain_id = strain_id.replace('Strain ', '', 1)
        
    strain_id = re.sub(r'\(.+\)', ' ', strain_id)
    strain_id = ' '.join(strain_id.split())
    strain_id = re.sub('[\W_]+', '', strain_id).upper()
    
    return strain_id
    
    
def check_format_strain(strain_id):
    """Check if strain ID has an acceptable format."""

    # run old format check
    old_format = False
    if old_format:
        # print in red
        print(f'\033[91mWE ARE USING THE OLD PARSING.\033[0m')
        print("MAKE SURE TO UPDATE THE STRAIN ID FORMAT.")

    # skip IDs without a number
    if old_format:
        if not any(char.isdigit() for char in strain_id):
            return False

        if all(c.isdigit() or c.isupper() for c in strain_id):
            return True
    else:

        if all(c.isdigit() or c.isupper() for c in strain_id) and len(strain_id) >=2 :
            return True

    special_characters = ['-', '.', ' ']
    processed_strain = str(strain_id)
    if len(processed_strain) < 1:
        return False
    for spechar in special_characters:
        processed_strain = processed_strain.replace(spechar, '')

    if all(c.isdigit() or c.isupper() for c in processed_strain):
        return True
        
    if strain_id.count(' ') == 0 and all(c.isdigit() or c.isupper() or c.lower() for c in processed_strain):
        return True


    return False
    