/*
 * fpUcHandling.cpp
 *
 *  Created on: 09.11.2009
 *      Author: daniel
 *      
 *  13/07/2015 Replace assert() with Rccp:Stop()
 */

#include <fpUcHandling.h>

//#include <cassert>
#include <algorithm>

#include <rcppExport.h>

// ***************************************************************************************************//

// small helper function to structure code:
// get the maximum power set for a range of FP maximum degrees
MyDoubleVector
getMaxPowerSet(const PosIntVector& fpmaxs)
{
    // always in the power set:
    static const double fixedpowers[] = { -2, -1, -0.5, 0, 0.5, 1, 2, 3 };
    //        corresponding indices        0   1     2  3    4  5  6  7

    // convert to double vector
    MyDoubleVector maxPowerset(fixedpowers, fixedpowers + 8);

    if (! fpmaxs.empty())
    {
        // consider if the power set is even bigger
        const PosInt biggestMaxDegree = * std::max_element(fpmaxs.begin(),
                                                           fpmaxs.end());
        for (PosInt more = 4; more <= biggestMaxDegree; ++more)
        {
            maxPowerset.push_back(more);
        }
    }

    // return the maximum power "set"
    return (maxPowerset);
}


// ***************************************************************************************************//

// fpInfo //

MyDoubleVector
FpInfo::inds2powers(const Powers& m) const // convert inds m (a "Powers" object) into powers vector p (a MyDoubleVector)
{
    MyDoubleVector ret;

    for (Powers::const_iterator j = m.begin();
            j != m.end();
            j++)
    {
        ret.push_back(powerset[*j]);
    }

    return ret;
}


// convert std::vector of powers into Powers object (i.e. a multiset) of ints of the power indexes.
Powers
FpInfo::vec2inds(const MyDoubleVector& p) const
{
    Powers ret;

    for(MyDoubleVector::const_iterator j = p.begin();
            j != p.end();
            j++)
    {
        // the index which we look for is the difference of the start of the powerset and the
        // iterator pointing to the location where the power was found.
        int index = find(powerset.begin(), powerset.end(), *j) - powerset.begin();
        ret.insert(index);
    }

    return ret;
}

// ***************************************************************************************************//


// Box-Tidwell transform
static inline double
boxtidwell(double x, double power)
{
    return (power != 0) ? pow(x, power) : log(x);
}

// ***************************************************************************************************//

// build array of vectors of AVectors holding the required transformed values for the design matrices
// do not! center the column. This is done inside getFpMatrix, because the repeated powers case cannot
// be treated here!!
AVectorArray
getTransformedCols(const PosIntVector& fpcards,
                   const PosIntVector& fppos,
                   const PosIntVector& fpmaxs,
                   const AMatrix& x)
{
    // initialize the return value
    AVectorArray transformedCols;

    // get maximum powerset
    MyDoubleVector maxPowerset = getMaxPowerSet(fpmaxs);

    // process each FP term
    PosIntVector::const_iterator card = fpcards.begin();
    for (PosIntVector::const_iterator pos = fppos.begin(); pos != fppos.end(); ++card, ++pos)
    {
        // the original column
        const AVector thisCol = x.col(*pos - 1);

        // start vector of columns for this FP term
        std::vector<AVector> thisFp;

        // for every possible power
        for (PosInt j = 0; j != *card; j++)
        {
            // start with original column
            AVector thisTransform = thisCol;

            // then transform it according to the power:
            for (PosInt k = 0; k < thisTransform.n_rows; ++k)
            {
                // assert(thisTransform(k) > 0);
              if(!(thisTransform(k) > 0)) Rcpp::stop("fpUcHandling.cpp:getTransformedCols: thisTransform(k) not greater than 0");

                thisTransform(k) = boxtidwell(thisTransform(k),
                                              maxPowerset[j]);

                // assert(! ISNAN(thisTransform(k)));
                if(ISNAN(thisTransform(k))) Rcpp::stop("fpUcHandling.cpp:getTransformedCols: thisTransform(k) is NAN");
            }

            // and put it into vector of columns
            thisFp.push_back(thisTransform);
        }

        // push vector of all powers for this FP term into the array
        transformedCols.push_back(thisFp);
    }

    // return the array
    return (transformedCols);
}

// ***************************************************************************************************//


// convert frequency vector into multiset
Powers
freqvec2Powers(IntVector& vec, const int &vecLength)
{
    Powers ret;
    for (int power = 0; power != vecLength; power++)
    {
        for (int times = 0; times != vec[power]; times++)
            ret.insert(power);
    }
    return ret;
}

// ***************************************************************************************************//

// End of file.
