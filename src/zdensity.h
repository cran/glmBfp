/*
 * zdensity.h
 *
 *  Created on: 16.11.2009
 *      Author: daniel
 */

#ifndef ZDENSITY_H_
#define ZDENSITY_H_

#include <iwls.h>
#include <coxfit.h>
#include <types.h>

// 21/11/2012: major update due to inclusion of TBF methodology

class NegLogUnnormZDens {
public:

    // call the function object:
    // z is the argument,
    double
    operator()(double z);

    // constructor
    NegLogUnnormZDens(const ModelPar &mod,
                      const DataValues& data,
                      const FpInfo& fpInfo,
                      const UcInfo& ucInfo,
                      const FixInfo& fixInfo,
                      const GlmModelConfig& config,
                      // return the approximate *conditional* density f(y | z, mod) by operator()?
                      // otherwise return the approximate unnormalized *joint* density f(y, z | mod).
                      const Book& bookkeep,
                      PosInt nIter=40);

    // try to get the TBF log marginal likelihood
    double
    getTBFLogMargLik() const;

    // get the maximum log conditional marginal likelihood
    // and put the local EB estimate into zMode
    double
    getTBFMaxLogCondMargLik(double& zMode) const;

    // get residual deviance
    double
    getResidualDeviance() const
    {
        return modResidualDeviance;
    }

    // destructor
    ~NegLogUnnormZDens()
    {
        delete iwlsObject;
        delete coxfitObject;
    }

private:
    // save the model reference and fp info, so we can write a nice warning message if the
    // IWLS fails for some z
    const ModelPar& mod;
    const FpInfo& fpInfo;
    const GlmModelConfig& config;
    const Book& bookkeep;

    // also save the original start linear predictor
    const AVector linPredStart;

    // pointer to an IWLS object
    Iwls * iwlsObject;

    // pointer to a Coxfit object
    Coxfit * coxfitObject;

    // number of IWLS iterations
    PosInt nIter;

    // the size of the model
    const int modSize;

    // the residual deviance of the model (only filled with correct value if TBF approach is used)
    double modResidualDeviance;
};


#endif /* ZDENSITY_H_ */
