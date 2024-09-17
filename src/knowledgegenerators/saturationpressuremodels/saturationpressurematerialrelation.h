#ifndef SATURATIONPRESSUREMATERIALRELATION_H
#define SATURATIONPRESSUREMATERIALRELATION_H

//#include "knowledgegenerators.h"
#include "../../base/thing.h"

class SaturationPressureMaterialRelation : public MaterialRelation {

public:
    std::string getClassName() const { return "SaturationPressureMaterialRelation"; }

    /// Finds and run first related Software object model.
    void run()
    { findNearest<SoftwareModel>()->run(); }
};

#endif //SATURATIONPRESSUREMATERIALRELATION_H
