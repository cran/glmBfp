/*
 * evalZdensity.cpp
 *
 *  Created on: 16.05.2010
 *      Author: daniel
 */


#include <rcppExport.h>
#include <dataStructure.h>
#include <types.h>
#include <zdensity.h>
#include <fpUcHandling.h>
#include <stdexcept>

using namespace Rcpp;

// ***************************************************************************************************//


// R call is:
//
//    samples <- .External(cpp_evalZdensity,
//                         config,
//                         attrs$data,
//                         attrs$fpInfos,
//                         attrs$ucInfos,
//                         attrs$distribution,
//                         options)

SEXP
cpp_evalZdensity(SEXP r_interface)
{
    // ----------------------------------------------------------------------------------
    // extract arguments
    // ----------------------------------------------------------------------------------

    r_interface = CDR(r_interface);
    List rcpp_config(CAR(r_interface));

    r_interface = CDR(r_interface);
    List rcpp_data(CAR(r_interface));

    r_interface = CDR(r_interface);
    List rcpp_fpInfos(CAR(r_interface));

    r_interface = CDR(r_interface);
    List rcpp_ucInfos(CAR(r_interface));

    r_interface = CDR(r_interface);
    List rcpp_distribution(CAR(r_interface));

    r_interface = CDR(r_interface);
    List rcpp_options(CAR(r_interface));

    // ----------------------------------------------------------------------------------
    // unpack the R objects
    // ----------------------------------------------------------------------------------

    // data:

    const NumericMatrix n_x = rcpp_data["x"];
    const AMatrix x(n_x.begin(), n_x.nrow(),
                   n_x.ncol());

    const NumericMatrix n_xCentered = rcpp_data["xCentered"];
    const AMatrix xCentered(n_xCentered.begin(), n_xCentered.nrow(),
                           n_xCentered.ncol());

    const NumericVector n_y = rcpp_data["y"];
    const AVector y(n_y.begin(), n_y.size());

    const IntVector censInd = as<IntVector>(rcpp_data["censInd"]);

    // FP configuration:

    // vector of maximum fp degrees
    const PosIntVector fpmaxs = as<PosIntVector>(rcpp_fpInfos["fpmaxs"]);
    // corresponding vector of fp column indices
    const PosIntVector fppos = rcpp_fpInfos["fppos"];
    // corresponding vector of power set cardinalities
    const PosIntVector fpcards = rcpp_fpInfos["fpcards"];
    // names of fp terms
    const StrVector fpnames = rcpp_fpInfos["fpnames"];


    // UC configuration:

    const PosIntVector ucIndices = rcpp_ucInfos["ucIndices"];
    List rcpp_ucColList = rcpp_ucInfos["ucColList"];

    std::vector<PosIntVector> ucColList;
    for (R_len_t i = 0; i != rcpp_ucColList.length(); ++i)
    {
        ucColList.push_back(as<PosIntVector>(rcpp_ucColList[i]));
    }


    // distributions info:

    const double nullModelLogMargLik = as<double>(rcpp_distribution["nullModelLogMargLik"]);
    const double nullModelDeviance = as<double>(rcpp_distribution["nullModelDeviance"]);
    const double empiricalMean = as<double>(rcpp_distribution["yMean"]);
    S4 rcpp_gPrior = rcpp_distribution["gPrior"];
    List rcpp_family = rcpp_distribution["family"];



    // options:

    const MyDoubleVector zValues = rcpp_options["zValues"];
//    const bool conditional = as<bool>(rcpp_options["conditional"]);
//    const bool debug = as<bool>(rcpp_options["debug"]);
//    const bool higherOrderCorrection = as<bool>(rcpp_options["higherOrderCorrection"]);

    const Book bookkeep(as<bool>(rcpp_distribution["tbf"]),
                        as<bool>(rcpp_distribution["doGlm"]),
                        false,
                        as<bool>(rcpp_options["conditional"]),
                        false, //useFixedc
                        100,
                        false,
                        as<bool>(rcpp_options["debug"]),
                        as<std::string>(rcpp_distribution["modelPrior"]),
                        10,
                        10,
                        10,
                        false,
                        as<bool>(rcpp_options["debug"]),
                        as<bool>(rcpp_options["higherOrderCorrection"]));


    // ----------------------------------------------------------------------------------
    // further process arguments
    // ----------------------------------------------------------------------------------

    // data:

     // only the intercept is always included, that is fixed, in the model
     IntSet fixedCols;
     fixedCols.insert(1);

     // totalnumber is set to 0 because we do not care about it.
     const DataValues data(x, xCentered, y, censInd, 0, fixedCols);

     // FP configuration:
     const FpInfo fpInfo(fpcards, fppos, fpmaxs, fpnames, x);

     // UC configuration:

     // determine sizes of the UC groups, and the total size == maximum size reached together by all
     // UC groups.
     PosIntVector ucSizes;
     PosInt maxUcDim = 0;
     for (std::vector<PosIntVector>::const_iterator cols = ucColList.begin(); cols != ucColList.end(); ++cols)
     {
         PosInt thisSize = cols->size();

         maxUcDim += thisSize;
         ucSizes.push_back(thisSize);
     }
     const UcInfo ucInfo(ucSizes, maxUcDim, ucIndices, ucColList);

     // search configuration:
     const GlmModelConfig config(rcpp_family, nullModelLogMargLik, nullModelDeviance, as<double>(rcpp_distribution["fixedg"]), rcpp_gPrior,
                                 data.response, bookkeep.debug, bookkeep.useFixedc, empiricalMean);
     // config of this model:
     const ModelPar thisModelConfig(rcpp_config, fpInfo);


     // ----------------------------------------------------------------------------------
     // evaluate the z density
     // ----------------------------------------------------------------------------------

     // get negative log unnormalized z density: a function object.
     NegLogUnnormZDens negLogUnnormZDens(thisModelConfig,
                                         data,
                                         fpInfo,
                                         ucInfo,
                                         config,
                                         bookkeep);

     // evaluate it at the given z values.
     NumericVector results;

     for (MyDoubleVector::const_iterator z = zValues.begin(); z != zValues.end(); ++z)
     {
         results.push_back(negLogUnnormZDens(* z));
     }

     // return the results vector
     return results;
}
