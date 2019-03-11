# -*- coding: utf-8 -*-
# GRID is a numerical integration library for quantum chemistry.
#
# Copyright (C) 2011-2019 The GRID Development Team
#
# This file is part of GRID.
#
# GRID is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# GRID is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --
""" Computes de Becke weighting function for every point in the grid """

import numpy as np

__all__ = ['becke_helper_atom']


def distanceBecke( point0, point1 ):
    """ Computes the distance between 2 points in the space xyz 
        Parameters
        ----------

        point0 :
            The coordinates of the first point

        point1 :
            The coordinates of the second point
    """
    tmp  = (point1[0] - point0[0])**2 
    tmp += (point1[1] - point0[1])**2 
    tmp += (point1[2] - point0[2])**2 

    return np.sqrt(tmp)
    


def becke_helper_atom( npoint, points, weights, natom, radii, centers, select, order)
    """ Computes the Becke weighting function for every point in the grid 

        Parameters
        ----------

        npoint :
            The number of grid points.

        points :
            The grid points. Row-mayor storage of (npoint,3) array.

        weights :
            The output, i.e. the Becke weights for the grid points. Note that the
            becke weights is **multiplied** with the original contents of the array!
        
        natom  :
            The number of atoms.

        radii :
            The radii of the atoms.

        centers :
            The posicion of the atoms.

        select :
            The select atom for wich the Becke weights are to be computed.

        order :
            The order of the switching function in the Becke scheme.

        See Becke's paper for the details:
        A. D. Becke, The Journal of Chemical Physics 88, 2547 (1988)
        URL http://dx.doi.org/10.1063/1.454033
    """
    alphas = np.zeros( int( natom*(natom+1)/2 ) )
    offset = 0
    
    for iatom in range( natom ) :
        for jatom in range ( iatom+1 ) :
            alpha = (radii[iatom] - radii[jatom])/(radii[iatom] + radii[jatom])
            alpha = alpha/(alpha**2 -1)
            if alpha > 0.45 :
                 alpha = 0.45
            elif alpha < -0.45 :
                alpha = -0.45
            alphas[offset] = alpha
            offset += 1

#   Precompute interatomic distances
    atomic_dists = np.zeros( int(natom*(natom+1)/2 ) )
    offset = 0

    for iatom in range( natom ) :
        for jatom in range ( iatom+1) :
            atomic_dists[offset] = distance( centers[3*iatom], centers[3*jatom] )
            offset += 1

#   Calculate the Becke Weights
    for ipoint in range ( npoint-1,-1,-1 ) :
        num = 0
        den = 0 
        for iatom in range ( natom ) :
            p = 1
            for jatom in range ( natom ) :
                if iatom == jatom :
                    continue
            # compute offset for alpha and interatomic distance
                if iatom < jatom :
                    offset = int( jatom*(jatom+1)/2 ) + iatom
                    term   = 1
                else:
                    offset = int( iatom*(iatom+1)/2 ) + jatom
                    tem    = 0
            # diatomic switching function
            
                s = (distance(points,centers[3*iatom]) - 
                     distance(points,centers[3*jatom]))/atomic_dists[offset]
                s = s + alphas[offset]*( 1 - 2*term )*( 1 - s**2)

                for k in range (1,order+1)
                    s = 0.5*s*(3-s**2)

                s = 0.5*(1-s)

                p *= s

            if iatom == selec :
                nom = p
            
            denom += p

        weights[ipoint] *= num/den;

        points += 3

        weights++

    


