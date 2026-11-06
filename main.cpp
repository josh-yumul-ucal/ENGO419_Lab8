#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <iomanip>
#include "functions.h"

namespace bm = boost::math;
namespace bg = boost::geometry;

using GeoPoint = bg::model::d2::point_xy<double>;
using GeoLine = bg::model::linestring<GeoPoint>;
using GeoPolygon = bg::model::polygon<GeoPoint>;

int main(){

  // Control Points:
  const double aCP[]{800.0, 600.0};
  const double bCP[]{1200.0, 500.0};

  // New Points:
  const double oneNP[]{1220.00, 1100.00};
  const double twoNP[]{820.00, 880.00};

  // Treat control points as fixed 
  Station A_{"A", Eigen::Vector2d{aCP}, true};
  Station B_{"B", Eigen::Vector2d{bCP}, true};

  Station P1{"P1", Eigen::Vector2d{oneNP}, false};
  Station P2{"P2", Eigen::Vector2d{twoNP}, false};


  // List of pointers to points for use in makeStations function
  const Station* initial[]{&A_, &B_, &P1, &P2};
  // vector of building points

  // Approximation of building corners and tank center
  const std::vector<Eigen::Vector2d> buildingPoints ={
    {1300.0, 980.0}, // top left corner
    {1300.0, 740.0}, // bottom left corner
    {1630.0, 740.0}, // bottom right corner
    {1630.0, 820.0}, // far right corner
    {1500.0, 820.0}, // inside corner
    {1500.0, 980.0} // top right corner
  };


  GeoPolygon building;

  for (const Eigen::Vector2d& point : buildingPoints){
    bg::append(building.outer(), GeoPoint(point.x(), point.y()));
  };

  bg::correct(building);

  // Approximate tank radius: 1170 - 1050 = 120m

  // Iteration 1:
  // Task 1: Visual approximation of added points
  Station P3{"P3", Eigen::Vector2d{1600.0, 1100.0}, false};

  std::vector<Station> addedStationsV1;
  addedStationsV1.push_back(P3);

  std::vector<Station> preliminary = makeStations(initial);
  std::vector<Station> stationsV1 = makeStations(addedStationsV1, initial);

  std::cout << std::fixed;
  std::cout.precision(8);

  analyze("Preliminary: no new stations",
        preliminary, 1, buildingPoints, building);

  analyze("Iteration 1: P3 added",
  stationsV1, 1, buildingPoints, building);


  return 0;
};
