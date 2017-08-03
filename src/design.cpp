/*
 * design.cpp
 *
 *  Created on: 09.11.2009
 *      Author: daniel
 */

#include <design.h>
#include <fpUcHandling.h>
#include <dataStructure.h>

#include <rcppExport.h>
#include <types.h>

#include <algorithm>


// get different concatenated columns of matrix
AMatrix
getMultipleCols(const AMatrix& m,
                const PosIntVector& s) // expect that the column numbers s are 1-based!
{
    AMatrix ret(m.n_rows, s.size());

    PosIntVector::size_type cols = 0; // invariant: about to process column nr. "cols"
    for (PosIntVector::const_iterator i = s.begin(); i != s.end(); i++)
    {
        ret.col(cols++) = m.col(*i - 1); // so we must subtract 1 here
    }

    return ret;
}


// ***************************************************************************************************//

// build Fp basis matrix from vector, power indices and power set for one covariate
AMatrix
getFpMatrix(const std::vector<AVector>& tcols,
            const Powers& powerinds,
            const DataValues& data)
{
    const int logInd = 3; // this index corresponds to power = 0, i.e. log.
    const int nrow = tcols.at(1).n_elem;

    // get the log column
    const AVector& logColumn = tcols.at(logInd);

    // this will be the returned matrix
    AMatrix ret(nrow, powerinds.size());

    // start recursion
    int lastInd = logInd;
    AVector lastCol = arma::ones<AVector>(nrow);

    // there is at least one power present
    Powers::size_type cols = 0; // invariant: about to process column number "cols"

    for (Powers::const_iterator
         now = powerinds.begin();
         now != powerinds.end();
         now++)
    {
        if (*now == lastInd)
        { // repeated powers case:

            // elementwise multiplication (Schur product):
            lastCol = lastCol % logColumn;
        }
        else
        { // normal case

            lastInd = *now;
            lastCol = tcols.at(lastInd);
        }
        // center the column
        ret.col(cols++) = lastCol - arma::mean(lastCol);
    }

    return ret;
}

// ***************************************************************************************************//

// construct centered design matrix including intercept for the model
// optionally the intercept column is not included
AMatrix
getDesignMatrix(const ModelPar &mod,
                const DataValues &data,
                const FpInfo &fpInfo,
                const UcInfo& ucInfo,
                const FixInfo& fixInfo,
                bool includeIntercept)
{
    // total number of columns
    int nColumns = mod.size(ucInfo, fixInfo); 
    if(includeIntercept)
    {
        ++nColumns;
    }

    // initialize the return matrix
    AMatrix ret(data.nObs, nColumns);

    // invariant: nextColumn is the next column to be written
    PosInt nextColumn = 0;

    // start with the intercept column?
    if(includeIntercept)
    {
        ret.col(0) = data.onesVector;
        ++nextColumn;
    }

    // go on with centered fp matrices
    for (PosInt i = 0; i != fpInfo.nFps; i++)
    {
        Powers powersi = mod.fpPars.at(i);

        if (! powersi.empty())
        {
            // what is the end column in the return matrix?
            PosInt endColumn = nextColumn + powersi.size() - 1;

            // insert the centered FP matrix into the return matrix
            ret.cols(nextColumn, endColumn) = getFpMatrix(fpInfo.tcols.at(i), powersi, data);

            // correct invariant
            nextColumn = endColumn + 1;
        }
    }

    // centered uc matrices
    for(IntSet::const_iterator
            g = mod.ucPars.begin();
            g != mod.ucPars.end();
            ++g)
    {
        // the C index is
        Int i = *g - 1;

        // what is the column list for this covariate?
        PosIntVector thisColList = ucInfo.ucColList.at(i);

        // what is the end column in the return matrix?
        PosInt endColumn = nextColumn + thisColList.size() - 1;

        // insert the centered UC matrix into the return matrix
        ret.cols(nextColumn, endColumn) = getMultipleCols(data.centeredDesign, thisColList);

        // correct invariant
        nextColumn = endColumn + 1;
    }
    
    
    // centered fix matrices
    for(IntSet::const_iterator
          g = mod.fixPars.begin();
        g != mod.fixPars.end();
        ++g)
    {
      // the C index is
      Int i = *g - 1;
      
      // what is the column list for this covariate?
      PosIntVector thisColList = fixInfo.fixColList.at(i);
      
      // what is the end column in the return matrix?
      PosInt endColumn = nextColumn + thisColList.size() - 1;
      
      // insert the centered UC matrix into the return matrix
      ret.cols(nextColumn, endColumn) = getMultipleCols(data.centeredDesign, thisColList);
      
      // correct invariant
      nextColumn = endColumn + 1;
    }
    

    return ret;
}

// ***************************************************************************************************//



// End of file.
