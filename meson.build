# ***********************************************************************
#
#   Copyright (C) 2025 -- The 4D-STAR Collaboration
#   File Author: Emily Boudreaux
#   Last Modified: June 21, 2025
#
#   GridFire is free software; you can use it and/or modify
#   it under the terms and restrictions the GNU General Library Public
#   License version 3 (GPLv3) as published by the Free Software Foundation.
#
#   GridFire is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#   See the GNU Library General Public License for more details.
#
#   You should have received a copy of the GNU Library General Public License
#   along with this software; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
#
# *********************************************************************** #
project('GridFire', ['c', 'cpp'], version: 'v1.0.3', default_options: ['cpp_std=c++23'], meson_version: '>=1.5.0')

gridfire_args = []

# Start by running the code which validates the build environment
subdir('build-check')

# Configure the logging level
subdir('build-extra/log-level')

# Then build the external dependencies
subdir('build-config')

# Build the main source code
subdir('src')

# Build the Python bindings
subdir('build-python')

# Build the test suite
subdir('tests')

# Build the tool suite
subdir('tools')

# Build the benchmark suite
subdir('benchmarks')

# Build the pkg-config file
subdir('build-extra/pkg-config')


