// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"

#include <climits>

namespace Acts {

namespace Test {

class SurfaceMaterialStub : public ISurfaceMaterial {
  using ISurfaceMaterial::ISurfaceMaterial;

  ISurfaceMaterial& operator*=(double /*scale*/) override { return *this; };

  const MaterialSlab& materialSlab(const Vector2& /*lp*/) const override {
    return m_fullMaterial;
  }

  const MaterialSlab& materialSlab(const Vector3& /*gp*/) const override {
    return m_fullMaterial;
  }

  const MaterialSlab& materialSlab(size_t /*bin0*/,
                                   size_t /*bin1*/) const override {
    return m_fullMaterial;
  }

  std::ostream& toStream(std::ostream& sl) const override {
    sl << "SurfaceMaterialStub";
    return sl;
  };

  MaterialSlab m_fullMaterial{};
};

/// Test the constructors
BOOST_AUTO_TEST_CASE(ISurfaceMaterial_factor_test) {
  double splitFactor = 42.0;
  SurfaceMaterialStub stub{splitFactor};

  BOOST_CHECK_EQUAL(stub.factor(NavigationDirection::Forward,
                                MaterialUpdateStage::FullUpdate),
                    1.0);

  BOOST_CHECK_EQUAL(stub.factor(NavigationDirection::Backward,
                                MaterialUpdateStage::FullUpdate),
                    1.0);

  BOOST_CHECK_EQUAL(stub.factor(NavigationDirection::Forward,
                                MaterialUpdateStage::PostUpdate),
                    splitFactor);

  BOOST_CHECK_EQUAL(stub.factor(NavigationDirection::Backward,
                                MaterialUpdateStage::PreUpdate),
                    splitFactor);

  BOOST_CHECK_EQUAL(
      stub.factor(NavigationDirection::Forward, MaterialUpdateStage::PreUpdate),
      1 - splitFactor);

  BOOST_CHECK_EQUAL(stub.factor(NavigationDirection::Backward,
                                MaterialUpdateStage::PostUpdate),
                    1 - splitFactor);
}

}  // namespace Test
}  // namespace Acts
