# -*- coding: utf-8 -*-
# %% [markdown]
#  # Implementation of the GRAS algorithm used for balancing consumer expenditure surveys and input-output data.
#
# %% [markdown]
# Code taken from:
# https://github.com/rich-wood/pygras/blob/main/gras.py

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# %%
import numpy as np

def gras(tabin, coltot, rowtot, iter_in=None):
    """
    Created on Mon Mar 27
    Richard Wood
    translated from Matlab code
    
    Balance tabin to column total (coltot) and row total (rowtot)
    
    Inputs:
    tabin: initial estimate of table to be balanced
    coltot: column totals to be reached of the table (must be same dimensions of table)
    rowtot: row totals to be reached of the table (must be same dimensions of table)
    iter_in: number of iterations (optional)
    
    Outputs:
    tabout: balanced table
    r1: row scaling vectors
    s1: column scaling vectors
    """
    
    # Check that the totals match, otherwise it will never converge!
    if abs(np.sum(rowtot) - np.sum(coltot)) > 1e-7:
        print('Warning: row and column totals do not match')
        print(abs(np.sum(rowtot)))
        print(abs(np.sum(coltot)))
    
    # Split the input table into a positive and negative 
    postab = np.maximum(tabin, 0)
    negtab = -np.minimum(tabin, 0)
    
    # Set up some dimension variables
    tabdim1, tabdim2 = postab.shape
    rdim, sdim = len(rowtot), len(coltot)
    
    # If the maximum number of iterations is externally defined, use it, otherwise use 100 iterations
    if iter_in is not None:
        MAXITER = iter_in
    else:
        MAXITER = 100
    print(f'MAXITER: {MAXITER}')
    
    # Initialise the row and column scaling vectors (r1 and s1) to unity
    if rdim > 0:
        r1 = np.ones((rdim, MAXITER+1))
    s1 = np.ones((MAXITER+1, sdim))
    
    for k in range(MAXITER):
        # s1[k, :] = scalcer(postab, negtab, r1[:, k], coltot, sdim, tabdim1)
        # r1[:, k+1] = rcalcer(postab, negtab, s1[k, :], rowtot, rdim, tabdim2)
        tmp=scalcer(postab,negtab,r1[0:rdim,k],coltot,sdim,tabdim1)
        s1[k,0:sdim] = tmp.transpose()
        tmp=rcalcer(postab,negtab,s1[k,0:sdim],rowtot,rdim,tabdim2)
        r1[0:rdim,k+1] = tmp.transpose()
        
        if np.sum(np.abs(r1[:, k+1] - r1[:, k])) < 1e-8:
            if k > 1:
                if np.sum(np.abs(s1[k, :] - s1[k-1, :])) < 1e-8:
                    print(f'Threshold reached: {k}')
                    break
    
    if k == MAXITER:
        print(f'Max runs reached: {k}')
    
    r2 = np.where(np.abs(r1) < 1e-10, 1, r1)
    s2 = np.where(np.abs(s1) < 1e-10, 1, s1)
    r2 = np.minimum(r2, 1e10)
    s2 = np.minimum(s2, 1e10)
    
    postab1 = np.hstack((np.matmul(postab[:, :sdim] , np.diag(s1[k, :])), postab[:, sdim:]))
    postab2 = np.vstack((np.diag(r1[:, k+1]) @ postab1[:rdim, :], postab1[rdim:, :]))


    # negtab = np.concatenate((np.dot(np.linalg.inv(np.diag(r2[:,k+1])), np.concatenate((negtab[0:r1.shape[0], 0:s1.shape[1]], np.linalg.inv(np.diag(s2[k,:]))), axis=1)),
    #                      negtab[0:r1.shape[0], s1.shape[1]:negtab.shape[1]]), axis=1)

    negtab1 = np.concatenate((np.dot(negtab[:, 0:s1.shape[1]], np.linalg.inv(np.diag(np.maximum(1e-5, s2[k,:])))),
                          negtab[0:negtab.shape[0], s1.shape[1]:negtab.shape[1]]), axis=1)

    negtab2 = np.concatenate((np.dot(np.linalg.inv(np.diag(np.maximum(1e-5, r2[:,k+1]))), negtab1[0:r1.shape[0], :]),
                          negtab1[r1.shape[0]:negtab1.shape[0], 0:negtab.shape[1]]), axis=0)

    tabout = postab2 - negtab2

    
    return tabout



def rcalcer(p, n, s, u, dim1, dim2):
    # Calculate the scaling factor on rows
    # p is positive matrix
    # n is negative matrix
    # s is scalar
    # dim1 is no of rows of p/n
    # dim2 is no of cols of p/n

    r = np.ones((dim1,1))

    for i in range(dim1):
        pscale = pcalcer(p, s, i, dim2)
        nscale = ncalcer(n, s, i, dim2)
        if pscale == 0:
            if nscale == 0:
                r[i] = 1
            else:
                r[i] = (u[i] + np.sqrt((u[i])**2+4*nscale))/(2)
        else:
            r[i] = (u[i] + np.sqrt((u[i])**2+4*pscale*nscale))/(2*pscale)
            # r[i] = (np.exp(1)*u[i] + np.sqrt((np.exp(1)*u[i])**2+4*pscale*nscale))/(2*pscale)

    return r

def pcalcer(p, s, i, dim2):
    pofr = 0
    for j in range(dim2):
        if j > len(s)-1:
            pofr += p[i,j]
        else:
            pofr += p[i,j]*s[j]

    return pofr

def ncalcer(n, s, i, dim2):
    nofr = 0
    for j in range(dim2):
        # if s[j] != 0:
        if j > len(s)-1:
            nofr += n[i,j]
        else:
            nofr += n[i,j]/(s[j]+1e-10)
        # end

    return nofr

    
def scalcer(p, n, r, v, dim1, dim2):
    """
    Calculate the scaling factor on columns
    p is positive matrix
    n is negative matrix
    r is scalar
    dim1 is no of rows of p/n
    dim2 is no of cols of p/n
    """
    s = np.ones(dim1)
    for j in range(dim1):
        pscale = spcalcer(p, r, j, dim2)
        nscale = sncalcer(n, r, j, dim2)
        if pscale == 0:
            if nscale == 0:
                s[j] = 1
            else:
                s[j] = (v[j] + np.sqrt(v[j]**2 + 4 * nscale)) / 2
        else:
            s[j] = (v[j] + np.sqrt(v[j]**2 + 4 * pscale * nscale)) / (2 * pscale)
            # s[j] = (np.exp(1) * v[j] + np.sqrt((np.exp(1) * v[j])**2 + 4 * pscale * nscale)) / (2 * pscale)
    return s


def spcalcer(p, r, j, dim2):
    """
    Calculate the scaling factor for positive matrix
    p is positive matrix
    r is scalar
    j is column index
    dim1 is no of rows of p
    """
    pofr = 0
    for i in range(dim2):
        if i >= r.size:
            pofr += p[i, j]
        else:
            pofr += p[i, j] * r[i]
    return pofr


def sncalcer(n, r, j, dim2):
    """
    Calculate the scaling factor for negative matrix
    n is negative matrix
    r is scalar
    j is column index
    dim1 is no of rows of n
    """
    nofr = 0
    for i in range(dim2):
        # if r[i] != 0:
        if i >= r.size:
            nofr += n[i, j]
        else:
            nofr += n[i, j] / (r[i] + 1e-10)
    return nofr
