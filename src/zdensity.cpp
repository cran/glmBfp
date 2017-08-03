#include <zdensity.h>
#include <stdexcept>

#include <rcppExport.h>
#include <linalgInterface.h>
#include <coxfit.h>

// 03/07/2013: use offsets

// constructor
NegLogUnnormZDens::NegLogUnnormZDens(const ModelPar &mod,
                                     const DataValues& data,
                                     const FpInfo& fpInfo,
                                     const UcInfo& ucInfo,
                                     const FixInfo& fixInfo,
                                     const GlmModelConfig& config,
                                     // return the approximate *conditional* density f(y | z, mod) by operator()?
                                     // otherwise return the approximate unnormalized *joint* density f(y, z | mod).
                                     const Book& bookkeep,
                                     PosInt nIter) :
                                     mod(mod),
                                     fpInfo(fpInfo),
                                     config(config),
                                     bookkeep(bookkeep),
                                     linPredStart(config.linPredStart),
                                     iwlsObject(0),
                                     coxfitObject(0),
                                     nIter(nIter),
                                     modSize(mod.size(ucInfo, fixInfo)), 
                                     modResidualDeviance(R_NaReal)
{
    if(bookkeep.doGlm)
    {
        iwlsObject = new Iwls(mod, data, fpInfo, ucInfo, fixInfo, config,
                              config.linPredStart,
                              // take the same original start value for each model, but then update
                              // it inside the iwls object when new calls to the functor are made.
                              (bookkeep.useFixedg || bookkeep.empiricalBayes),
                              EPS, // take EPS as the convergence epsilon
                              bookkeep.debug,
                              bookkeep.tbf);

        if(bookkeep.tbf)
        {
            // compute the residual deviance for this model
            modResidualDeviance = config.nullModelDeviance - iwlsObject->computeDeviance(100);

            // echo detailed progress in debug mode
            if(bookkeep.debug)
            {
                Rprintf("\nNegLogUnnormZDens: Residual deviance of this model is %f", modResidualDeviance);
            }
        }

    }
    else
    {
        coxfitObject = new Coxfit(data.response,
                                  data.censInd,
                                  getDesignMatrix(mod, data, fpInfo, ucInfo, fixInfo, false),
                                  config.weights,
                                  config.offsets,
                                  1);

        // compute the residual deviance for this model
        modResidualDeviance = coxfitObject->computeResidualDeviance();

        // echo detailed progress in debug mode
        if(bookkeep.debug)
        {
            Rprintf("\nNegLogUnnormZDens: Residual deviance of this Cox model is %f", modResidualDeviance);
        }
    }
}

// try to get the TBF log marginal likelihood
double
NegLogUnnormZDens::getTBFLogMargLik() const
{
    double ret = R_NaReal;

    if(bookkeep.tbf)
    {
        // ask g-prior object to give log marginal likelihood
        // if a closed form is not available (as for most classes) then
        // NA is returned.
        ret = config.gPrior->getTBFLogMargLik(modResidualDeviance, modSize);
    } else {
        std::ostringstream stream;
        stream << "getTBFLogMargLik asked from NegLogUnnormZDens, but TBF methodology is not used!";
        throw std::domain_error(stream.str().c_str());
    }

    return ret;
}

// get the maximum log conditional marginal likelihood
// and put the local EB estimate into zMode
double
NegLogUnnormZDens::getTBFMaxLogCondMargLik(double& zMode) const
{
    double ret = R_NaReal;

    if(bookkeep.tbf)
    {
        if(modResidualDeviance > modSize)
        {
            double gMode = modResidualDeviance / modSize - 1.0;
            zMode = log(gMode);

            ret = - modSize / 2.0 * log1p(gMode) + gMode / (gMode + 1.0) * modResidualDeviance / 2.0;
        }
        else
        {
            zMode = R_NegInf;
            ret = 0.0;
        }
    } else {
        std::ostringstream stream;
        stream << "getTBFMaxLogCondMargLik asked from NegLogUnnormZDens, but TBF methodology is not used!";
        throw std::domain_error(stream.str().c_str());
    }

    return ret;
}

// call the function object
double
NegLogUnnormZDens::operator()(double z)
{
    // map back to the original covariance factor scale
    const double g = exp(z);

    // optional status message
    if(bookkeep.debug)
    {
        Rprintf("\nNegLogUnnormZDens: Computing function call NegLogUnnormZDens(%f) ...", z);
    }

    // if g is numerically zero, we cannot proceed
    if(g == 0.0)
    {
        std::ostringstream stream;
        stream << "g numerically zero";
        throw std::domain_error(stream.str().c_str());
    }

    // the return value will be placed in here:
    double ret = 0.0;

    // first consider the TBF case
    if(bookkeep.tbf)
    {
        if(bookkeep.debug)
        {
            Rprintf("\nNegLogUnnormZDens: Using TBF approach");
        }

        double logBF = - modSize / 2.0 * log1p(g) + g / (g + 1.0) * modResidualDeviance / 2.0;

        // Note that the following distinction is done *inside*
        // Iwls::computeLogUnPosteriorDens for the generalized hyper-g prior approach.
        if(bookkeep.useFixedg || bookkeep.empiricalBayes)
        {
            // return the negative conditional log marginal likelihood approximation, - log(f(y | z))
            ret = - logBF;
        }
        else
        {
            // return - log(f(y | z) f(z)) = - log(f(y, z)), the negative unnormalized log posterior of z
            double logGprior = config.gPrior->logDens(g) + z;
            ret = - logBF - logGprior;
        }
    }
    else // hyper-g prior
    {
        if(bookkeep.debug)
        {
            Rprintf("\nNegLogUnnormZDens: Using generalised hyper-g prior approach");
        }

        // protect against errors in IWLS and non-finite result
        try
        {
            // compute Gaussian approximation:
            PosInt requiredIter = iwlsObject->startWithLastLinPred(nIter, // iterate at most nIter times (because we want convergence)
                                                                  g);    // for this covariance factor, for the current linear predictor start.

            // if we did not converge in the maximum number of iterations, then start again from the original linear predictor,
            // and allow even more iterations.
            if(requiredIter > nIter)
            {
                PosInt higherIter = 2 * nIter;
                requiredIter = iwlsObject->startWithNewLinPred(higherIter,
                                                              g,
                                                              linPredStart);

                // if we still do not have converged, print a warning message if we are verbose.
                if(requiredIter > higherIter && bookkeep.debug)
                {
                    Rprintf("\nnNegLogUnnormZDens: IWLS did not converge in %d iterations (even after restart)!", higherIter);
                }

                // do not warn the user, because this need not be a serious problem.
            }

            if(bookkeep.debug)
            {
                Rprintf("\nNegLogUnnormZDens: finished after %d IWLS iterations", requiredIter);
            }

            // get iwls results
            const IwlsResults iwlsResults = iwlsObject->getResults();

            // then the return value is:
            ret = 0.5 * iwlsResults.logPrecisionDeterminant -
                    iwlsObject->nCoefs * M_LN_SQRT_2PI -    // This is dependent on the model dimension! (we had to add it here because
                    // computeLogUnPosteriorDens also contains now the analogous normalization constant.
                    iwlsObject->computeLogUnPosteriorDens(Parameter(iwlsResults.coefs, z));


            // the Raudenbush & Yang & Yosef correction to the Laplace approximation
            // in the canonical link case

            // Warning: this may not work yet for some reason for some data sets,
            // strange negative correctionFactor values can occur!
            // (both for the removed Fortran and for the actual C++ code...)
            if(bookkeep.higherOrderCorrection && config.canonicalLink)
            {

                // initialize E(T_4), E(T_6):
                // (excluding constant factors which are multiplied at the end)
                double exp_T4 = 0.0;
                double exp_T6 = 0.0;

                // for the sum over m_i^(3) * x_i * B_i:
                AVector k = arma::zeros(iwlsObject->nCoefs);

                // design matrix is also needed
                const AMatrix& design = iwlsObject->design;

                // iterate through the observations:
                for(PosInt i = 0; i < iwlsObject->nObs; ++i)
                {
                    // calculate the fitted probability and resulting derivatives m_i^(x)
                    // along the way.
                    const double mu = config.link->linkinv(iwlsResults.linPred(i));
                    double m3, m4, m6;

                    // see page 147 in Raudenbush et al. for the formulas of m3, m4, m6
                    if(config.familyString == "binomial")
                    {
                        const double m2 = mu * (1.0 - mu);
                        m3 = m2 * (1.0 - 2.0 * mu);
                        m4 = m2 * (1.0 - 6.0 * m2);
                        m6 = m4 * (1.0 - 12.0 * m2) - 12.0 * m3 * m3;
                    }
                    else if (config.familyString == "poisson")
                    {
                        m3 = mu;
                        m4 = mu;
                        m6 = mu;
                    }
                    else
                    {
                        Rf_error("Higher order correction not implemented for this family.");
                    }

                    // add to the sums:
                    AVector tmp = arma::trans(design.row(i));
                    trs(false, false, iwlsResults.qFactor, tmp);
                    const double B = arma::dot(tmp, tmp);
                    const double B2 = B * B;

                    exp_T4 += m4 * B2;
                    exp_T6 += m6 * B2 * B;
                    k += (m3 * B) * arma::trans(design.row(i));
                }

                // calculate k^T * Cov * k, and overwrite k with the intermediate result:
                trs(false, false, iwlsResults.qFactor, k);
                double exp_T3T3 = arma::dot(k, k);

                // So the correction factor for the conditional marginal likelihood is (with the new terms):
                double correctionFactor = - 1.0 / 8.0 * exp_T4 -
                        1.0 / 48.0 * exp_T6 +
                        5.0 / 24.0 * exp_T3T3;

                if(bookkeep.debug)
                {
                    Rprintf("\nNegLogUnnormZDens: Higher-order correction factor is %f", 1.0 + correctionFactor);
                }

                if(correctionFactor > -1.0)
                {
                    // since we return here the negative log of the conditional marg lik:
                    ret -= log1p(correctionFactor);
                } else {
                    Rf_warning("negative value for correction factor! We are not using the higher order correction here.");
                }

            } // end if(higherOrderCorrection)

            // check finiteness of return value
            if(! R_finite(ret))
            {
                std::ostringstream stream;
                stream << "NegLogUnnormZDens() got non-finite result " << ret << " for z=" << z;
                throw std::domain_error(stream.str().c_str());
            }
        }
        // catch errors in IWLS and non-finite result
        catch (std::domain_error& error)
        {
            if(bookkeep.debug)
            {
                Rprintf("For z=%f, the density value could not be computed because\n%s",
                        z, error.what());
            }
            Rf_warning("for z=%f, the density value could not be computed for the following model:\n%s\nCheck for near-collinearity of covariates.",
                       z, mod.print(fpInfo).c_str());

            // return NaN. This can be handled by the Brent optimize routine! It apparently replaces it (implicitely) by
            // the highest positive number.
            return R_NaN;
        }
    }

    // return value
    return ret;

}
