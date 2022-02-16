#!/bin/bash

# Converts .hic files to .mcool files, preserving all the resolutions stored
# in the source file
# 
# Usage: convert_hic_to_cool.sh input.hic output.mcool
#
# Copyright 2022 Luca Nanni
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

hic_path=${1}
output_path=${2}

hic2cool convert ${1} ${2} -r 0