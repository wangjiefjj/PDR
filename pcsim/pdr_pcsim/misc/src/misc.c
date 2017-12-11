#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include "misc.h"




/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
DBL** mallocArray2D_DBL(U32 row, U32 col)
{
    DBL **pmatrix = NULL;
    U32 i = 0;
    U32 sizeofDBL = sizeof(DBL);
    U32 sizeofPtrDBL = sizeof(DBL *);

    pmatrix = (DBL **)malloc(sizeofPtrDBL * row);

    if (pmatrix == NULL)
    {
        return NULL;
    }

    for (i = 0; i < row; i++)
    {
        pmatrix[i] = (DBL *)malloc(sizeofDBL * col);

        if (pmatrix[i] == NULL)
        {
            U32 j;

            for (j = 0; j < i; j++)
            {
                free(pmatrix[j]);
            }
            free(pmatrix);
            
            return NULL;
        }
        memset(pmatrix[i], 0, sizeofDBL * col);
    }

    return pmatrix;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
U32 freeArray2D_DBL(DBL **pmatrix, U32 row, U32 col)
{
    U32 i = 0;
    PARAMETER_NOT_USED(col);

    for (i = 0; i < row; i++)
    {
        free(pmatrix[i]);
    }

    return 0;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
U32 computeMeanStd(FLT* const mean, FLT* const std, const FLT array[][CHN], U32 count)
{
    U32 i = 0;
    FLT sum = 0;
    
    if (count < 2)
    {
        return -1;
    }

    for (i = 0; i < count; i++)
    {
        FLT det = 0;

        det = sqrtf(array[i][CHX]*array[i][CHX] + array[i][CHY]*array[i][CHY] + array[i][CHZ]*array[i][CHZ]);
        sum += det;
    }
    *mean = sum / count;

    sum = 0;
    for (i = 0; i < count; i++)
    {
        FLT det = 0;

        det = sqrtf(array[i][CHX]*array[i][CHX] + array[i][CHY]*array[i][CHY] + array[i][CHZ]*array[i][CHZ]);
        sum += (det - *mean) * (det - *mean);
    }
    *std = sqrtf(sum / (count - 1.0F));

    return 0;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
void dcm2euler(const FLT cbn[3][3], FLT* const pyaw, FLT* const ppitch, FLT* const proll)
{
    *pyaw = atan2f(cbn[CHY][CHX], cbn[CHX][CHX]);
    *ppitch = asinf(-cbn[CHZ][CHX]);
    *proll = atan2f(cbn[CHZ][CHY], cbn[CHZ][CHZ]);
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
void euler2dcm(FLT cbn[3][3], FLT fyaw, FLT fpitch, FLT froll)
{
    cbn[0][0] = cosf(fpitch) * cosf(fyaw);
    cbn[0][1] = sinf(froll) * sinf(fpitch) * cosf(fyaw) - cosf(froll) * sinf(fyaw);
    cbn[0][2] = cosf(froll) * sinf(fpitch) * cosf(fyaw) + sinf(froll) * sinf(fyaw);

    cbn[1][0] = cosf(fpitch) * sinf(fyaw);
    cbn[1][1] = sinf(froll) * sinf(fpitch) * sinf(fyaw) + cosf(froll) * cosf(fyaw);
    cbn[1][2] = cosf(froll) * sinf(fpitch) * sinf(fyaw) - sinf(froll) * cosf(fyaw);

    cbn[2][0] = -sinf(fpitch);
    cbn[2][1] = sinf(froll) * cosf(fpitch);
    cbn[2][2] = cosf(froll) * cosf(fpitch);
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
void euler2q(FLT q[4], FLT fyaw, FLT fpitch, FLT froll)
{
    q[0] = cosf(froll / 2) * cosf(fpitch / 2) * cosf(fyaw / 2) + sinf(froll / 2) * sinf(fpitch / 2) * sinf(fyaw / 2);
    q[1] = sinf(froll / 2) * cosf(fpitch / 2) * cosf(fyaw / 2) - cosf(froll / 2) * sinf(fpitch / 2) * sinf(fyaw / 2);
    q[2] = cosf(froll / 2) * sinf(fpitch / 2) * cosf(fyaw / 2) + sinf(froll / 2) * cosf(fpitch / 2) * sinf(fyaw / 2);
    q[3] = cosf(froll / 2) * cosf(fpitch / 2) * sinf(fyaw / 2) - sinf(froll / 2) * sinf(fpitch / 2) * cosf(fyaw / 2);
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
void q2dcm(FLT q[4], FLT cbn[3][3])
{
    cbn[0][0] = q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3];
    cbn[0][1] = 2 * (q[1] * q[2] - q[0] * q[3]);
    cbn[0][2] = 2 * (q[1] * q[3] + q[0] * q[2]);

    cbn[1][0] = 2 * (q[1] * q[2] + q[0] * q[3]);
    cbn[1][1] = q[0] * q[0] - q[1] * q[1] + q[2] * q[2] - q[3] * q[3];
    cbn[1][2] = 2 * (q[2] * q[3] - q[0] * q[1]);

    cbn[2][0] = 2 * (q[1] * q[3] - q[0] * q[2]);
    cbn[2][1] = 2 * (q[2] * q[3] + q[0] * q[1]);
    cbn[2][2] = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3];
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
void qNorm(FLT q[4])
{
    FLT fnorm;

    fnorm = sqrtf(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
    if (fnorm > 0.001F)
    {
        fnorm = 1.0F / fnorm;
        q[0] *= fnorm;
        q[1] *= fnorm;
        q[2] *= fnorm;
        q[3] *= fnorm;
    }
    else
    {
        // return with identity quaternion since the quaternion is corrupted
        q[0] = 1.0F;
        q[1] = 0.0F;
        q[2] = 0.0F;
        q[3] = 0.0F;
    }

    // correct a negative scalar component if the function was called with negative q0
    if (q[0] < 0.0F)
    {
        q[0] = -q[0];
        q[1] = -q[1];
        q[2] = -q[2];
        q[3] = -q[3];
    }


}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
void f3x3matrixTranspose(FLT matrix[3][3])
{
    U32 i = 0;
    U32 j = 0;
    FLT temp;

    for (i = 1; i < 3; i++)
    {
        for (j = 0; j < i; j++)
        {
            temp = matrix[i][j];
            matrix[i][j] = matrix[j][i];
            matrix[j][i] = temp;
        }
    }
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
void fmatrixAeqI(FLT *A[], U32 rc)
{
    // rc = rows and columns in A

    FLT *pAij;	    // pointer to A[i][j]
    U32 i, j;		// loop counters

    for (i = 0; i < rc; i++)
    {
        // set pAij to &A[i][j=0]
        pAij = A[i];
        for (j = 0; j < rc; j++)
        {
            *(pAij++) = 0.0F;
        }
        A[i][i] = 1.0F;
    }
    return;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
void f3x3matrixEqI(FLT matrix[3][3])
{
    U32 i = 0;
    U32 j = 0;

    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            matrix[i][j] = 0.0F;
        }
        matrix[i][i] = 1.0F;
    }
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
void f3x3matrixEqScalar(FLT matrix[3][3], FLT scalar)
{
    U32 i = 0;
    U32 j = 0;

    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            matrix[i][j] = scalar;
        }
    }
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
void f3x3matrixEqAxScalar(FLT A[][3], FLT scalar)
{
    FLT *pAij;
    U32 i, j;

    for (i = 0; i < 3; i++)
    {
        // set pAij to &A[i][j=0]
        pAij = A[i];
        for (j = 0; j < 3; j++)
        {
            *(pAij++) *= scalar;
        }
    }
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
void f3x3matrixEqMinusA(FLT A[][3])
{
    FLT *pAij;
    U32 i, j;

    for (i = 0; i < 3; i++)
    {
        // set pAij to &A[i][j=0]
        pAij = A[i];
        for (j = 0; j < 3; j++)
        {
            *pAij = -*pAij;
            pAij++;
        }
    }
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
U32 matrixMult(const DBL** const matrixA, const DBL** const matrixB, U32 rowA, U32 colA, U32 rowB, U32 colB, DBL** const matrixC)
{
    U32 i = 0;
    U32 j = 0;
    U32 k = 0;
    DBL temp = 0.0F;

    if (colA != rowB)
    {
        return -1;
    }

    for (i = 0; i < rowA; i++)
    {
        for (j = 0; j < colB; j++)
        {
            temp = 0.0F;
            for (k = 0; k < colA; k++)
            {
                temp += matrixA[i][k] * matrixB[k][j];
            }
            matrixC[i][j] = temp;
        }
    }

    return 0;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
void fmatrixAeqInvA(FLT *A[], S32 iColInd[], S32 iRowInd[], S32 iPivot[], S32 isize, U32 *pierror)
{
    FLT largest;
    FLT scaling;
    FLT recippiv;
    FLT ftmp;
    S32 i, j, k, l, m;
    U32 iPivotRow, iPivotCol;

    iPivotRow = iPivotCol = 0;
    *pierror = 0;
    for (j = 0; j < isize; j++)
    {
        iPivot[j] = 0;
    }
    for (i = 0; i < isize; i++)
    {
        largest = 0.0F;
        for (j = 0; j < isize; j++)
        {
            if (iPivot[j] != 1)
            {
                for (k = 0; k < isize; k++)
                {
                    if (iPivot[k] == 0)
                    {
                        if (fabsf(A[j][k]) >= largest)
                        {
                            iPivotRow = j;
                            iPivotCol = k;
                            largest = (FLT)fabsf(A[iPivotRow][iPivotCol]);
                        }
                    }
                    else if (iPivot[k] > 1)
                    {
                        fmatrixAeqI(A, isize);
                        *pierror = 1;
                        return;
                    }
                }
            }
        }
        iPivot[iPivotCol]++;
        if (iPivotRow != iPivotCol)
        {
            for (l = 0; l < isize; l++)
            {
                ftmp = A[iPivotRow][l];
                A[iPivotRow][l] = A[iPivotCol][l];
                A[iPivotCol][l] = ftmp;
            }
        }
        iRowInd[i] = iPivotRow;
        iColInd[i] = iPivotCol;
        if (A[iPivotCol][iPivotCol] == 0.0F)
        {
            fmatrixAeqI(A, isize);
            *pierror = 1;
            return;
        }
        recippiv = 1.0F / A[iPivotCol][iPivotCol];
        A[iPivotCol][iPivotCol] = 1.0F;
        for (l = 0; l < isize; l++)
        {
            if (A[iPivotCol][l] != 0.0F)
                A[iPivotCol][l] *= recippiv;
        }
        for (m = 0; m < isize; m++)
        {
            if (m != iPivotCol)
            {
                scaling = A[m][iPivotCol];
                A[m][iPivotCol] = 0.0F;
                for (l = 0; l < isize; l++)
                {
                    if ((A[iPivotCol][l] != 0.0F) && (scaling != 0.0F))
                        A[m][l] -= A[iPivotCol][l] * scaling;
                }
            }
        }
    }
    for (l = isize - 1; l >= 0; l--)
    {
        i = iRowInd[l];
        j = iColInd[l];

        if (i != j)
        {
            for (k = 0; k < isize; k++)
            {
                ftmp = A[k][i];
                A[k][i] = A[k][j];
                A[k][j] = ftmp;
            }
        }
    }
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
void eigencompute10(FLT A[][10], FLT eigval[], FLT eigvec[][10], S32 n)
{
#define NITERATIONS 15

    FLT cot2phi, tanhalfphi, tanphi, sinphi, cosphi;
    FLT ftmp;
    FLT residue;
    S32 ir, ic;
    S32 j;
    S32 ctr;

    for (ir = 0; ir < n; ir++)
    {
        for (ic = 0; ic < n; ic++)
        {
            eigvec[ir][ic] = 0.0F;
        }
        eigvec[ir][ir] = 1.0F;
        eigval[ir] = A[ir][ir];
    }
    ctr = 0;
    do
    {
        residue = 0.0F;
        for (ir = 0; ir < n - 1; ir++)
        {
            for (ic = ir + 1; ic < n; ic++)
            {
                residue += fabsf(A[ir][ic]);
            }
        }
        if (residue > 0.0F)
        {
            for (ir = 0; ir < n - 1; ir++)
            {
                for (ic = ir + 1; ic < n; ic++)
                {
                    if (fabsf(A[ir][ic]) > 0.0F)
                    {
                        cot2phi = 0.5F * (eigval[ic] - eigval[ir]) / (A[ir][ic]);
                        tanphi = 1.0F / (fabsf(cot2phi) + sqrtf(1.0F + cot2phi * cot2phi));
                        if (cot2phi < 0.0F)
                        {
                            tanphi = -tanphi;
                        }
                        cosphi = 1.0F / sqrtf(1.0F + tanphi * tanphi);
                        sinphi = tanphi * cosphi;
                        tanhalfphi = sinphi / (1.0F + cosphi);
                        ftmp = tanphi * A[ir][ic];
                        eigval[ir] -= ftmp;
                        eigval[ic] += ftmp;
                        A[ir][ic] = 0.0F;
                        for (j = 0; j < n; j++)
                        {
                            ftmp = eigvec[j][ir];
                            eigvec[j][ir] = ftmp - sinphi * (eigvec[j][ic] + tanhalfphi * ftmp);
                            eigvec[j][ic] = eigvec[j][ic] + sinphi * (ftmp - tanhalfphi * eigvec[j][ic]);
                        }
                        for (j = 0; j <= ir - 1; j++)
                        {
                            ftmp = A[j][ir];
                            A[j][ir] = ftmp - sinphi * (A[j][ic] + tanhalfphi * ftmp);
                            A[j][ic] = A[j][ic] + sinphi * (ftmp - tanhalfphi * A[j][ic]);
                        }
                        for (j = ir + 1; j <= ic - 1; j++)
                        {
                            ftmp = A[ir][j];
                            A[ir][j] = ftmp - sinphi * (A[j][ic] + tanhalfphi * ftmp);
                            A[j][ic] = A[j][ic] + sinphi * (ftmp - tanhalfphi * A[j][ic]);
                        }
                        for (j = ic + 1; j < n; j++)
                        {
                            ftmp = A[ir][j];
                            A[ir][j] = ftmp - sinphi * (A[ic][j] + tanhalfphi * ftmp);
                            A[ic][j] = A[ic][j] + sinphi * (ftmp - tanhalfphi * A[ic][j]);
                        }
                    }
                }
            }
        }
    } while ((residue > 0.0F) && (ctr++ < NITERATIONS));
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
void udDecompose(DBL** const matrix, U32 n)
{
    const DBL EPS = 1.0e-8;
    U32 i = 0;
    U32 j = 0;
    U32 k = 0;
    DBL amax = 0.0;
    DBL alpha = 0.0;
    DBL beta = 0.0;

    // Search for maximum of diagonal elements and scale down by 1.0e+8
    // for threshold for setting alpha to zero.
    amax = 0.0;

    for (i = 0; i <= n - 1; i++)
    {
        if (amax < matrix[i][i])
        {
            amax = matrix[i][i];
        }
    }

    amax = amax * EPS;

    for (j = n; j >= 2; j--)
    {
        if (matrix[j - 1][j - 1] > amax)
        {
            alpha = (DBL)1.0 / matrix[j - 1][j - 1] ;

            for (k = 1; k <= j - 1; k++)
            {
                beta   =  matrix[k - 1][j - 1];
                matrix[k - 1][j - 1] =  alpha * beta;

                for (i = 1; i <= k; i++)
                {
                    matrix[i - 1][k - 1] -= beta * matrix[i - 1][j - 1];
                }
            }   // k in 1..j-1
        }
        else
        {
            for (k = 1; k <= j - 1; k++)
            {
                matrix[k - 1][j - 1] = 0.0;
            }
        }
    }    // j in reverse 2..n
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
void multPhimUp(const DBL** const phim, const DBL* const u, U32 n, DBL** const w)
{
    U32 i = 0;
    U32 j = 0;
    U32 k = 0;
    U32 l = 0;
    U32 jm1 = 0;
    U32 np2 = n + 2;
    U32 nn1 = (n * (n + 1)) / 2;
    DBL sum = 0.0;

    for (i = 0; i < n; i++)
    {
        w[i][0] = phim[i][0];
    }

    for (l = 2; l <= n; l++)
    {
        j = np2 - l;
        nn1 = nn1 - j;
        jm1 = j - 1;
        for (i = 1; i <= n; i++)
        {
            if (j > n)
            {
                sum = 0.0;
                jm1 = n;
            }
            else
            {
                sum = phim[i - 1][j - 1];
            }
            
            for (k = 1; k <= jm1; k++)
            {
                sum = sum + phim[i - 1][k - 1] * u[nn1 + k - 1];
            }

            w[i - 1][j - 1] = sum;
        } // for (i = 1; i <= n; i++)

    } // for (l = 2; l <= n; l++)

}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
void storeUq(const U32 l, const U32 m, const DBL** const u, const U32 istart, const U32 jstart, DBL** const w)
{
    U32 i = 0;
    U32 j = 0;
    U32 row_r = 0;
    U32 col_r = 0;
    U32 row_u = 0;
    U32 col_u = 0;
    U32 istm1 = 0;
    U32 jstm1 = 0;
    U32 jm1 = 0;
    U32 jp1 = 0;

    // Begin
    // Copy u into r(istart..istart+m-1,jstart..jstart+m-1)
    istm1 = istart - 1;
    jstm1 = jstart - 1;

    for (j = 1; j <= m; j++)
    {
        jm1   = j - 1;
        col_r = jstm1 + j;
        col_u = l + j;

        for (i = 1; i <= jm1; i++)
        {
            row_r = istm1 + i;
            row_u = l + i;
            w[row_r - 1][col_r - 1] = u[row_u - 1][col_u - 1];
        }

        row_r = istm1 + j;
        w[row_r - 1][col_r - 1] = 1.0;
        jp1 = j + 1;

        for (i = jp1; i <= m; i++)
        {
            row_r = istm1 + i;
            w[row_r - 1][col_r - 1] = 0.0;
        }
    }
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
void udTimeUpdate(const U32 iw, const U32 jw, DBL** const w, const DBL* const dw, DBL* const u)
{
    U32 jwt2 = 2 * jw;
    U32 jm1 = 0;
    U32 jwpk = 0;
    U32 i = 0;
    U32 ij = 0;
    U32 j = 0;
    U32 j_index = 0;
    U32 k = 0;
    DBL w1k = 0.0;
    DBL dinv = 0.0;
    DBL sum = 0.0;
    DBL *v = NULL;

    v = (DBL*)malloc(sizeof(DBL) * iw * 4);

    if (iw > 1)
    {
        for (j = iw; j >= 2; j--) //j=n,...,2
        {
            sum = 0.0;

            for (k = 1; k <= jw; k++)
            {
                v[k - 1] = w[j - 1][k - 1];
                jwpk = jw + k;
                v[jwpk - 1] = dw[k - 1] * v[k - 1]; //D*w(j)//D---dw
                sum  = v[k - 1] * v[jwpk - 1] + sum;
            }//d(j)=w(j)*D*w(j)

            w[j - 1][j - 1] = sum;
            dinv   = sum;///dinv=d(j)
            jm1    = j - 1;

            if (sum <= (DBL)0.0)
            {
                // w(j,.) := 0 when dinv = 0 (dinv = norm(w(j,.))**2);
                for (k = 1; k <= jm1; k++)
                {
                    w[j - 1][k - 1] = 0.0;
                }
            }
            else
            {
                for (k = 1; k <= jm1; k++) //i=1,2,...j-1
                {
                    sum  = 0.0;

                    for (i = 1; i <= jw; i++)
                    {
                        sum = w[k - 1][i - 1] * v[jw + i - 1] + sum;
                    }

                    sum = sum / dinv;//u(i,j)=w(i)*D*w(j)/d(j)

                    for (i = 1; i <= jw; i++)
                    {
                        w[k - 1][i - 1] = w[k - 1][i - 1] - sum * v[i - 1];    //w(i)=w(i)-u(i,j)*w(j)
                    }

                    w[j - 1][k - 1] = sum;
                }
            }   // if (sum <=  0.0)
        }   // for (j=iw; j>=2; j--) => (j in reverse 2..iw)

        // The lower part of w is u' (u transpose)
    }   // if (iw > 1) */

    sum = 0.0;

    for (k = 1; k <= jw; k++)
    {
        w1k = w[1 - 1][k - 1];
        sum = sum + dw[k - 1] * w1k * w1k; //d(1)=w(1)*D*w(1)
    }

    u[0] = sum;//u(1,1)=d(1)

    if (iw > 1)
    {
        ij = 1;

        for (j_index = 2; j_index <= iw; j_index++)
        {
            for (i = 1; i <= j_index; i++)
            {
                ij = ij + 1;
                u[ij - 1] = w[j_index - 1][i - 1];
            }
        }  // for (j_index=2; j_index<=iw; j_index++)
    }  // if (iw > 1)

    free(v);
}

/*-------------------------------------------------------------------------*/
/**
  @brief    
  @param    
  @return   
  

 */
/*--------------------------------------------------------------------------*/
void udMeasUpdate(DBL* const u, DBL* const x, const U32 n, const DBL r, const DBL* const a, const DBL z, DBL* const alpha, DBL* const res)
{
    U32 np1 = 0;
    U32 np2 = 0;
    U32 ntot = 0;
    U32 j = 0;
    U32 jj = 0;
    U32 jm1 = 0;
    U32 jjn = 0;
    U32 kj = 0;
    U32 j_index = 0;
    U32 k = 0;
    U32 l = 0;
    DBL beta = 0.0;
    DBL gamma = 0.0;
    DBL sum = 0.0;
    DBL temp = 0.0;
    DBL p = 0.0;
    DBL res_local = 0.0;
    DBL s = 0.0;
    DBL *gx = (DBL*)malloc(sizeof(DBL) * n);
    DBL *utax = (DBL*)malloc(sizeof(DBL) * n);
    DBL *g = (DBL*)malloc(sizeof(DBL) * n);
    DBL *uta = (DBL*)malloc(sizeof(DBL) * n);

    memset(gx, 0, sizeof(DBL) * n);
    memset(utax, 0, sizeof(DBL) * n);
    memset(g, 0, sizeof(DBL) * n);
    memset(uta, 0, sizeof(DBL) * n);

    // Begin
    np1  = n + 1;
    np2  = n + 2;
    ntot = (n * np1) / 2;
    // Compute residual z := z-ax
    temp = 0;

    for (j = 0; j < n; j++) //pay attention here
    {
        temp += a[j] * x[j];
    }

    res_local = z - temp;
    // Compute uta = u'a;
    //           g = d(u'a);
    jjn = ntot;

    for (l = 2; l <= n; l++)
    {
        j   = np2 - l;
        jj  = jjn - j;
        sum = a[j - 1];
        jm1 = j - 1;

        for (k = 1; k <= jm1; k++)
        {
            sum = sum + u[jj + k - 1] * a[k - 1];
        }

        uta[j - 1]     = sum;
        g[j - 1] = sum * u[jjn - 1];
        jjn = jj;
    } // for (l=2; l<=n; l++)

    uta[0] = a[0];
    g[0] = u[0] * uta[0];
    // uta = u'a;
    // g   = d(u'a);
    sum = r + g[0] * uta[0];
    gamma = 0.0;

    if (sum > (DBL)0.0)
    {
        gamma = (DBL)1.0 / sum;
    }

    if (uta[0] != (DBL)0.0)
    {
        u[0] = u[0] * r * gamma;
    }

    kj = 2;

    for (j_index = 2; j_index <= n; j_index++)
    {
        beta =  sum;
        temp =  g[j_index - 1];
        sum  =  sum + temp * uta[j_index - 1];
        p    = -uta[j_index - 1] * gamma;
        jm1  = j_index - 1;

        for (k = 1; k <= jm1; k++)
        {
            s            = u[kj - 1];
            u[kj - 1]     = s + p * g[k - 1];
            g[k - 1] = g[k - 1] + temp * s;
            kj           = kj + 1;
        }

        if (fabs(temp) >= (DBL)1.0e-20)
        {
            gamma = (DBL)1.0 / sum;
            u[kj - 1] = u[kj - 1] * beta * gamma;
        }

        kj = kj + 1;
    }    // for (j_index=2; j_index<=n; j_index++)

    *alpha = sum;
    *res   = res_local;

    for (j = 1; j <= n; j++)
    {
        g[j - 1] *= gamma;          // g = gamma*g;
        x[j - 1] += res_local * g[j - 1]; // x = x + res_local*g;
    }

    free(gx);
    free(utax);
    free(g);
    free(uta);
}
