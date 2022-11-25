###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
#       Date: 2012-12-16
#
#  Copyright (C) 2012 Feel++ Consortium
#
# Distributed under the GPL(GNU Public License):
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#

cmake_minimum_required(VERSION 3.0)
project(mor)

enable_testing()

# find Feel++ package
find_package(Feel++
  PATHS /usr/share/feelpp/feel/cmake/modules
  REQUIRED)

find_package(Feel++-MOR
  PATHS /usr/share/feelpp/mor/cmake/
  REQUIRED)

find_package(Feel++-Toolboxes
  PATHS /usr/share/feelpp/toolboxes/cmake/
  REQUIRED)

if(NOT FEELPP_FOUND)
    message(WARNING "Feel++ was not found on your system. Make sure to install it and specify the FEELPP_DIR to reference the installation directory.")
else()
  message(STATUS "[feelpp] found Feel++ library in ${Feel++_DIR}")
endif()

if( NOT Feel++-Toolboxes_FOUND )
    message(WARNING "Feel++-Toolboxes was not found on your system. Make sure to install it and specify the FEELPP_TOOLBOXES_DIR to reference the installation directory.")
else()
  message(STATUS "[feelpp-toolboxes] found Feel++-Toolboxes library in ${Feel++-Toolboxes_DIR}")
endif()

if ( NOT Feel++-MOR_FOUND )
    message(WARNING "Feel++-MOR was not found on your system. Make sure to install it and specify the FEELPP_MOR_DIR to reference the installation directory.")
else()
  message(STATUS "[feelpp-mor] found Feel++-MOR library in ${Feel++-MOR_DIR}")
endif()

find_package(OpenTURNS 1.20
  REQUIRED)
message(STATUS "Found OpenTURNS: ${OpenTURNS_DIR}  (found version \"${OpenTURNS_VERSION}\")")

# find_package(OpenMP
#   REQUIRED)
# message(STATUS "Found OpenMP")


add_subdirectory(src)
add_subdirectory(test)