/*
 * design.h
 *
 *  Created on: 09.11.2009
 *      Author: daniel
 */

#ifndef DESIGN_H_
#define DESIGN_H_

#include <dataStructure.h>
#include <types.h>

// construct centered design matrix including intercept for the model
// optionally the intercept column is not included
AMatrix
getDesignMatrix(const ModelPar &mod,
                const DataValues &data,
                const FpInfo &fpInfo,
                const UcInfo& ucInfo,
                const FixInfo& fixInfo,
                bool includeIntercept = true);



#endif /* DESIGN_H_ */
