#include "functions.h"
#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <vector>

namespace bm = boost::math;
namespace bg = boost::geometry;

using GeoPoint = bg::model::d2::point_xy<double>;
using GeoLine = bg::model::linestring<GeoPoint>;
using GeoPolygon = bg::model::polygon<GeoPoint>;

// Statistical values
const double chi95 = sqrt(boost::math::quantile(
  boost::math::chi_squared(2), 0.95));

const double delta0 = bm::quantile(bm::normal(), 0.975) + bm::quantile(bm::normal(), 0.80);


constexpr double PI = 3.1415926535897932384626433832795;
constexpr double SIGMA_DIST = 0.002;
constexpr double SIGMA_ANGLE = 2.0 * PI / 648000.0;
constexpr double TOLERANCE = 1e-89;

bool visible(const Eigen::Vector2d& a, const Eigen::Vector2d& b, const GeoPolygon building){
  GeoLine sight{
    {a.x(), a.y()},
    {b.x(), b.y()}
  };
  return bg::distance(tankCenter, sight) >= tankRadius && !bg::crosses(sight, building);
}


Eigen::RowVectorXd designRow(const std::vector<Station>& s, int i, int j, bool isAzimuth){
  Eigen::Vector2d v = s[j].position - s[i].position;
  Eigen::Vector2d d = v / v.norm();
  if (isAzimuth){
    d = Eigen::Vector2d(v.y(), -v.x()) / v.squaredNorm();
  };
  Eigen::RowVectorXd row = Eigen::RowVectorXd::Zero(2 * (s.size() - 2));

  if (i >= 2){
    row.segment<2>(2 * (i - 2)) -= d.transpose();
  };

  if (j >= 2){
    row.segment<2>(2 * (j - 2)) += d.transpose();
  }
  return row;
};