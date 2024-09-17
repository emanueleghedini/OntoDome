#ifndef SURFACETENSIONMATERIALRELATION_H
#define SURFACETENSIONMATERIALRELATION_H

//#include "knowledgegenerators.h"
#include "../../base/thing.h"

class SurfaceTensionMaterialRelation : public MaterialRelation {

public:
    std::string getClassName() const { return "SurfaceTensionMaterialRelation"; }

    /// Find and run the first related Software object model.
    void run()
    { findNearest<SoftwareModel>()->run(); }
};

#endif // SURFACETENSIONMATERIALRELATION_H
