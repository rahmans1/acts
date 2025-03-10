// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <any>
#include <vector>

/// @note this is foreseen for the 'Geometry' module

namespace Acts {

class Surface;

namespace Experimental {

class Portal;
class Detector;
class DetectorVolume;

/// @brief A navigation state struct that is
/// holding the current navigation information
/// about volume, surfaces, and portals
struct NavigationState {
  /// @brief  A surface candidate and its intersection
  ///
  /// candidates can either be surfaces or portals (which contain a surface)
  struct SurfaceCandidate {
    /// A candidate intersection, in Surface view
    ObjectIntersection<Surface> objectIntersection;
    /// A candidate is either a detector Surface
    const Surface* surface = nullptr;
    /// Or a portal
    const Portal* portal = nullptr;
    /// The boundary check used for the candidate, boundary checks
    /// can differ for sensitive surfaces and portals
    BoundaryCheck bCheck = true;
  };

  /// Surface candidate vector alias, this allows to use e.g. boost_small vector
  /// or other stl like containers
  using SurfaceCandidates = std::vector<SurfaceCandidate>;

  /// The current position
  Vector3 position = Vector3(0., 0., 0.);

  /// The current direction
  Vector3 direction = Vector3(0., 0., 0.);

  /// The current absolute momentum
  ActsScalar absMomentum = 0.;

  /// The current charge
  ActsScalar charge = 0.;

  /// The current magnetic field
  Vector3 magneticField = Vector3(0., 0., 0.);

  /// The current detector in processing
  const Detector* currentDetector = nullptr;

  /// The current volume in processing, i.e. the position is inside
  const DetectorVolume* currentVolume = nullptr;

  /// The current surface, i.e the position is on surface
  const Surface* currentSurface = nullptr;

  /// That are the candidate surfaces to process
  SurfaceCandidates surfaceCandidates = {};
  SurfaceCandidates::iterator surfaceCandidate = surfaceCandidates.end();

  /// Boundary directives for surfaces
  BoundaryCheck surfaceBoundaryCheck = true;

  /// An overstep tolerance
  ActsScalar overstepTolerance = -0.1;

  /// Auxilliary attached information
  std::any auxilliary;
};

/// Filler of the current volume
struct DetectorVolumeFiller {
  /// Helper struct that allows to fill a volume into the
  /// navigation state, it allows to use common navigation
  /// structs for volume, portal, surfaces
  ///
  /// @param nState the navigation state
  /// @param volume the volume that is filled
  inline static void fill(NavigationState& nState,
                          const DetectorVolume* volume) {
    nState.currentVolume = volume;
  }
};

/// Fillers and attachers for surfaces to act on the navigation state
struct SurfacesFiller {
  /// Helper struct that allows to fill surfaces into the candidate vector it
  /// allows to use common navigation structs for volume, portal, surfaces
  ///
  /// @param nState the navigation state
  /// @param surfaces the surfaces that are filled in
  inline static void fill(NavigationState& nState,
                          const std::vector<const Surface*>& surfaces) {
    std::for_each(surfaces.begin(), surfaces.end(), [&](const auto& s) {
      nState.surfaceCandidates.push_back(NavigationState::SurfaceCandidate{
          ObjectIntersection<Surface>{}, s, nullptr,
          nState.surfaceBoundaryCheck});
    });
  }
};

/// Fillers and attachers for portals to act on the navigation state
struct PortalsFiller {
  /// Helper struct that allows to fill surfaces into the candidate vector it
  /// allows to use common navigation structs for volume, portal, surfaces
  ///
  /// @param nState the navigation state
  /// @param portals the portals that are filled in
  inline static void fill(NavigationState& nState,
                          const std::vector<const Portal*>& portals) {
    std::for_each(portals.begin(), portals.end(), [&](const auto& p) {
      nState.surfaceCandidates.push_back(NavigationState::SurfaceCandidate{
          ObjectIntersection<Surface>{}, nullptr, p, true});
    });
  }
};

}  // namespace Experimental
}  // namespace Acts
