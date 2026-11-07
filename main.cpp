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

  std::vector<Station> addedStations;
  addedStations.push_back(P3);

  std::vector<Station> preliminary = makeStations(initial);
  std::vector<Station> stationsV1 = makeStations(addedStations, initial);

  std::cout << std::fixed;
  std::cout.precision(8);

  designResult preliminaryTrial = analyze("Preliminary: no new stations",
        preliminary, 1, buildingPoints, building);

  designResult iter1 = analyze("Iteration 1: P3 added",
  stationsV1, 1, buildingPoints, building);
  designResult iter2 = analyze("Iteration 2: P3 added, 3 rounds of observations", stationsV1, 3, buildingPoints, building);

  Station P4{"P4", Eigen::Vector2d{1700.0, 600.0}, false};
  addedStations.push_back(P4);
  std::vector<Station> stationsV2 = makeStations(addedStations, initial);
  designResult iter3 = analyze("Iteration 3: P4 added", stationsV2, 1, buildingPoints, building);

  designResult iter4 = analyze("Iteration 4: P4 added,  2 rounds of observations", stationsV2, 2, buildingPoints, building);

  designResult iter5 = analyze("Iteration 5: P3 added only, 4 rounds of observations", stationsV1, 4, buildingPoints, building);

  Station P5{"P5", Eigen::Vector2d{1300, 600}, false};
  addedStations.push_back(P5);
  std::vector<Station> stationsV3 = makeStations(addedStations, initial);
  designResult iter6 = analyze("Iteration 6: P5 added", stationsV3, 1, buildingPoints, building);

  designResult iter7 = analyze("Iteration 7: P5 added, 2 rounds of observations", stationsV3, 2, buildingPoints, building);

  designResult iter8 = analyze("Iteration 8: P3 added only, 7 rounds of observations", stationsV1, 7, buildingPoints, building);

  designResult iter9 = analyze("Iteration 9: P5 added, 3 rounds of measurements", stationsV3, 3, buildingPoints, building);

  std::vector<Station> addedStationsV4{P3, P5};
  std::vector<Station> stationsV4 = makeStations(addedStationsV4, initial);

  designResult iter10 = analyze("Iteration 10: P3 + P5, 4 rounds of observations", stationsV4, 4, buildingPoints, building);

  saveResults(iter8, building, "iter8");
  saveResults(iter4, building, "iter4");
  saveResults(iter3, building, "iter3");
  saveResults(preliminaryTrial, building, "preliminary");
  saveResults(iter9, building, "iter9");
  saveResults(iter10, building, "iter10");

  return 0;
};
