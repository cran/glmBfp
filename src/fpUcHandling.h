/*
 * fpUcHandling.h
 *
 *  Created on: 09.11.2009
 *      Author: daniel
 */

#ifndef FPUCHANDLING_H_
#define FPUCHANDLING_H_

#include <types.h>
#include <numeric>
#include <rcppExport.h>

// ***************************************************************************************************//

// small helper function to structure code:
// get the maximum power set for a range of FP maximum degrees
MyDoubleVector
getMaxPowerSet(const PosIntVector& fpmaxs);

// ***************************************************************************************************//

// build array of vectors of ColumnVectors holding the required transformed values for the design matrices
// do not! center the column. This is done inside getFpMatrix, because the repeated powers case cannot
// be treated here.
AVectorArray
getTransformedCols(const PosIntVector& fpcards,
                   const PosIntVector& fppos,
                   const PosIntVector& fpmaxs,
                   const AMatrix& x);


// ***************************************************************************************************//


struct FpInfo
{ // collects all information on fractional polynomials needed to be passed down
    PosInt nFps;
    MyDoubleVector powerset;
    PosIntVector fpcards;
    PosIntVector fppos;
    PosIntVector fpmaxs;
    StrVector fpnames;
    AVectorArray tcols;
    PosInt maxFpDim;

    // number of possible univariate fps for each FP?
    IntVector numberPossibleFps;

    // what is the multiset expressing a linear inclusion of a covariate?
    Powers linearPowers;


    // ctr
    FpInfo(const PosIntVector& fpcards,
           const PosIntVector& fppos,
           const PosIntVector& fpmaxs,
           const StrVector& fpnames,
           const AMatrix& x) :
        nFps(fpmaxs.size()), powerset(getMaxPowerSet(fpmaxs)), fpcards(fpcards),
                fppos(fppos), fpmaxs(fpmaxs), fpnames(fpnames), tcols(getTransformedCols(fpcards, fppos, fpmaxs, x)),
                maxFpDim(std::accumulate(fpmaxs.begin(), fpmaxs.end(), 0)),
                numberPossibleFps(),
                linearPowers()
    {
        // numbers of possible univariate fps?
        for(PosInt i=0; i != nFps; ++i)
        {
            int thisNumber = 0;
            for(PosInt deg = 0; deg <= fpmaxs[i]; ++deg)
            {
                thisNumber += Rf_choose(fpcards[i] - 1 + deg, deg);
            }
            numberPossibleFps.push_back(thisNumber);
        }

        // insert the index 5 for linear power 1
        linearPowers.insert(5);

    }

    // convert inds m into powers vector
    MyDoubleVector
    inds2powers(const Powers& m) const;

    // convert std::vector of powers into Powers object (i.e. a multiset)
    Powers
    vec2inds(const MyDoubleVector& p) const;

};

// ***************************************************************************************************//

// collects all information on uncertain fixed form covariates groups
struct UcInfo
{
    const PosIntVector ucSizes;
    const PosInt maxUcDim;
    const PosIntVector ucIndices;
    const std::vector <PosIntVector> ucColList;
    const PosInt nUcGroups;

    UcInfo(const PosIntVector& ucSizes,
           const PosInt maxUcDim,
           const PosIntVector& ucIndices,
           const std::vector<PosIntVector>& ucColList) :
        ucSizes(ucSizes), maxUcDim(maxUcDim), ucIndices(ucIndices),
                ucColList(ucColList), nUcGroups(ucColList.size())
    {
    }
};

// ***************************************************************************************************//

// collects all information on fixed form covariates groups
struct FixInfo
{
	const PosIntVector fixSizes;
	const PosInt maxFixDim;
	const PosIntVector fixIndices;
	const std::vector <PosIntVector> fixColList;
	const PosInt nFixGroups;

	FixInfo(const PosIntVector& fixSizes,
		const PosInt maxFixDim,
		const PosIntVector& fixIndices,
		const std::vector<PosIntVector>& fixColList) :
		fixSizes(fixSizes), maxFixDim(maxFixDim), fixIndices(fixIndices),
		fixColList(fixColList), nFixGroups(fixColList.size())
	{
	}
};

// ***************************************************************************************************//


// return iterator of random element of myset; should be enclosed in getRNGstate() etc.
template<class T>
    typename T::iterator
    discreteUniform(const T& container)
    {
        if (container.empty())
        {
            Rf_error("\ncontainer in call to discreteUniform is empty!\n");
        }

        double u = unif_rand();

        typename T::size_type size = container.size();
        typename T::const_iterator i = container.begin();
        typename T::size_type j = 1;

        while (u > 1.0 / size * j)
        {
            i++;
            j++;
        }

        return i;
    }

// ***************************************************************************************************//

// get random int x with lower <= x < upper; should be enclosed in getRNGstate() etc.
template<class INT>
INT
discreteUniform(const INT& lower, const INT& upper)
{
    if (lower >= upper)
    {
        Rf_error("\nlower = %d >= %d = upper in discreteUniform call\n", lower,
                 upper);
    }

    double u = unif_rand();

    INT size = upper - lower;
    INT ret = lower;

    while (u > 1.0 / size * (ret - lower + 1))
    {
        ret++;
    }

    return ret;
}

// ***************************************************************************************************//


// delete a number from a set
template<class T>
    typename std::set<T>
    removeElement(std::set<T> input, T element)
    {
        typename std::set<T>::iterator iter = input.begin();
        while (iter != input.end())
        {
            if (*iter == element)
                // A copy of iter is passed into erase(), ++ is executed after erase().
                // Thus iter remains valid
                input.erase(iter++);
            else
                ++iter;
        }

        return input;
    }

// ***************************************************************************************************//

// construct a sequence 1:maximum
template<class T>
    typename std::set<T>
    constructSequence(T maximum)
    {
        std::set<T> ret;

        for (T i = 1; i <= maximum; ++i)
        {
            ret.insert(ret.end(), i);
        }

        return ret;
    }

// ***************************************************************************************************//

// convert frequency vector into multiset
Powers
freqvec2Powers(IntVector& vec, const int &vecLength);

// ***************************************************************************************************//


#endif /* FPUCHANDLING_H_ */
