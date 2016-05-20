// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DD4hepLayerHelper.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DD4HEPGEOMETRYTOOLS_DD4HEPLAYERHELPER_H
#define ACTS_DD4HEPGEOMETRYTOOLS_DD4HEPLAYERHELPER_H 1

// Core module
#include "ACTS/Utilities/Definitions.hpp"
// Geometry module
#include "ACTS/Detector/TrackingVolume.hpp"
#include "ACTS/Utilities/BinnedArray.hpp"
#include "ACTS/Tools/ITrackingVolumeBuilder.hpp"
// DD4hep
#include "DD4hep/Detector.h"

namespace Acts {
    class Surface;
    class Volume;
}

namespace Acts {
    
    /** @ class DD4hepLayerHelper
        
     The DD4hepLayerHelper translates Layers with its containing modules (if given) from the DD4hep geometry into Layers and Surfaces from TrackingGeometry.
     @TODO find replacement for Gaudi exeption and message stream
     
    */
    
    class DD4hepGeometryHelper;
    
    typedef std::vector<const Acts::Surface*>                        SurfaceVector;
    typedef Acts::BinnedArray<const Acts::Surface*>                   SurfaceArray;
    
   class DD4hepLayerHelper  {
       
   public:
       /** constructor */
       DD4hepLayerHelper();
       
       /** destructor */
       ~DD4hepLayerHelper();
       
       /** Creates a triple of all the three possible Layer types of the given volume detector element*/
       const Acts::LayerTriple* createLayerTriple(DD4hep::Geometry::DetElement& motherDetElement);
       /** returns a pair of a possible barrel and endcap volume*/
       const Acts::VolumeTriple* volumeTriple();
       
   private:
       /**constructs all layers contained by this volume*/
       void constructLayers(DD4hep::Geometry::DetElement& detElement);
       /**creates the cylindrical shaped layers*/
       void createCylinderLayers(DD4hep::Geometry::DetElement& motherDetElement, std::shared_ptr<Acts::Transform3D> motherTransform);
       /**creates disc shaped layers*/
       void createDiscLayers(DD4hep::Geometry::DetElement& motherDetElement, std::shared_ptr<Acts::Transform3D> motherTransform);
       /**method to create a binned array out of  a vector of surfaces for cylinder layers (surfaces binned in z and phi)*/
       SurfaceArray* createCylinderBinnedSurfaceArray(DD4hep::Geometry::DetElement& motherDetElement, std::shared_ptr<const Acts::Transform3D> motherTransform, double Zmin, double Zmax) const;
       /**method to create a binned array out of a vector of surfaces for disc layers (surfaces binned in r and phi)*/
       SurfaceArray* createDiscBinnedSurfaceArray(DD4hep::Geometry::DetElement& motherDetElement, std::shared_ptr<const Acts::Transform3D> motherTransform) const;
       /**method which creates a vector of contained surfaces by a given DD4hep mother DetElement - used internally to extract the surfaces of a layer*/
       void createSurfaceVector(DD4hep::Geometry::DetElement& motherDetElement, std::shared_ptr<const Acts::Transform3D> motherTransform, SurfaceVector& surfaces) const;
       /**helper method to get the bin values for a binned array out of overlapping modules*/
       std::vector<float> orderRValues(std::vector<std::pair<float,float>> old) const;
        /**helper method to sort pairs of doubles*/
       static bool sortFloatPairs(std::pair<float,float> ap, std::pair<float,float> bp);
       
       /**possible layers of the negative end cap*/
       Acts::LayerVector*                                m_negativeLayers;
       /**possible layers of the central barrel*/
       Acts::LayerVector*                                m_centralLayers;
       /**possible layers of the positive end cap*/
       Acts::LayerVector*                                m_positiveLayers;
       /**the barrel volume of the current hierarchy*/
       VolumePtr                                         m_barrelVolume;
       /**the negative endcap volume of the current hierarchy*/
       VolumePtr                                         m_nEndcapVolume;
       /**the positive endcap volume of the current hierarchy*/
       VolumePtr                                         m_pEndcapVolume;
       
    };
    
} //end of namespace Acts

inline const Acts::VolumeTriple* Acts::DD4hepLayerHelper::volumeTriple() { return  new VolumeTriple(m_nEndcapVolume, m_barrelVolume,m_pEndcapVolume);}


#endif //ACTS_DD4HEPGEOMETRYTOOLS_DD4HEPLAYERHELPER_H
