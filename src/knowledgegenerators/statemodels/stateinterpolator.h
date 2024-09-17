#ifndef STATEINTERPOLATOR_H
#define STATEINTERPOLATOR_H

#include "../../base/thing.h"

class StateInterpolator : public SoftwareModel {
public:
  std::string getClassName() const { return "State temporal intepolator"; }

  template <class TT> auto intepolate_state(TT* obj, double time) {
      // Get the pair embracing given time step from states map
      std::pair<double,Matter*> down,up;
      for (auto st : obj->states) {
          if (st.first > time) {
              up = st;
              break;
          }
          down = st;
      }

      // Abort if requested time is out of available range
      if (up.first == 0) { abort(); }

      // Create the interpolated state to be filled with new datas
      auto interpolated = new Matter;

        for (auto o1 : down.second->getRelatedObject<Thing>()) {
          for (auto o2 : up.second->getRelatedObject<Thing>()) {
            if (o1->getClassName() == o2->getClassName()) {
              if (dynamic_cast<ScalarQuantity*>(o1)) {
                o1 = dynamic_cast<ScalarQuantity*>(o1);
                auto v1 = o1->getRelatedObject<Scalar>()[0]->data;
                auto v2 = o2->getRelatedObject<Scalar>()[0]->data;
                auto vf = v1 + (time - down.first)/(up.first - down.first)*(v2-v1);
                auto dup = o1;
                dup->getRelatedObject<Scalar>()[0]->data = vf;
                interpolated->createRelationTo<hasProperty,Thing>(dup);
              }
              else if (dynamic_cast<VectorQuantity*>(o1)) {
                  auto v1 = o1->getRelatedObject<Vector>()[0]->data;
                  auto v2 = o2->getRelatedObject<Vector>()[0]->data;
                  std::vector<double> vec;
                  for (auto i : v1) {
                    vec.push_back(v1[i] + (time - down.first)/(up.first - down.first)*(v2[i]-v1[i]));
                  }
                  auto dup = o1;
                  dup->getRelatedObject<Vector>()[0]->data = vec;
                  interpolated->createRelationTo<hasProperty,Thing>(dup);
              }
              else if (dynamic_cast<Matter*>(o1)) {
                if (o1->getUuid() == o2->getUuid()) {
                  auto dup = o1;
                  interpolated->createRelationTo<hasPart,Thing>(dup);
                }
              }
            }
          }
        }

      return interpolated;
  }
};

#endif // STATEINTERPOLATOR_H
