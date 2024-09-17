#ifndef KNOWLEDGEGENERATORS_H
#define KNOWLEDGEGENERATORS_H

#include "../base/thing.h"

//#include "surfacetensionmodels/surfacetensionpolynomialmodel.h"
//#include "saturationpressuremodels/sapolynomialmodel.h"
//#include "gasmodels/gasmodel.h"
//#include "gasmodels/gasmodelcv.h"
//#include "nanomodels/nucleation/cnt.h"
//#include "nanomodels/moments/momentmodelpratsinis.h"
//#include "statemodels/stateinterpolator.h"

//class KnowledgeGeneratorsDB {
//protected:
//  std::vector<Thing*> list;

//public:
//  // Constructor
//  KnowledgeGeneratorsDB() {
//    // KnowledgeGenerators must be added manually by the developer
////    list.push_back(new STPolynomialModel);
////    list.push_back(new SAPolynomialModel);
//  }

//  // search for knowledge generators which are model for the given entity
//  auto ismodelfor (Thing* criteria) {
////  auto ismodelfor (std::string criteria) {
//    std::vector<std::string> champions;

//    for (auto i : list) {
//      auto mod_for = i->getRelations<hasSoftwareModel>();
//      for (auto ii : mod_for) {
//        auto modeled = ii->getRange()->getClassName();
//        if (modeled == criteria->getClassName())  {
////        if (modeled == criteria)  {
//          champions.push_back(i->getClassName());
//        }
//      }
//    }
//    return champions;
//  }

//};

#endif // KNOWLEDGEGENERATORS_H
