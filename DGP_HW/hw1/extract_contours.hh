////////////////////////////////////////////////////////////////////////////////
// extract_contours.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Contour extraction framework code for ECS 289H HW 1.
*/
//  Author:  Julian Panetta (jpanetta), jpanetta@ucdavis.edu
////////////////////////////////////////////////////////////////////////////////
#ifndef EXTRACT_CONTOUR_HH
#define EXTRACT_CONTOUR_HH

#include <Eigen/Dense>
#include <functional>

using PtRef = Eigen::Ref<const Eigen::Vector3d>;
using ScalarField = std::function<double(PtRef)>;

// Extract the contours for the scalar field `sf` defined over the triangle
// mesh (V, F).
// The contour points are returned in `P` and the edges connecting them in `E`.
// If `linearInterpolation` is false, the contour points should be generated
// at the edge midpoint. If it is true, it should be linearly interpolated
// based on the contour value.
// Bonus:
// If `linearInterpolation` is true and `snapEpsilon > 0`, contour points whose
// relative distance along the edge from one of the edge endpoints is less than
// `snapEpsilon` will be snapped to that endpoint. Contour edges that
// degenerate to zero length due to this snapping should be neglected.
void extract_contours(Eigen::Ref<const Eigen::MatrixX3d> V, 
                      Eigen::Ref<const Eigen::MatrixX3i> F,
                      const ScalarField &sf,
                      Eigen::MatrixX3d &P,
                      Eigen::MatrixX2i &E,
                      bool linearInterpolation = true,
                      double snapEpsilon = 0.0);

#endif /* end of include guard: EXTRACT_CONTOUR_HH */
