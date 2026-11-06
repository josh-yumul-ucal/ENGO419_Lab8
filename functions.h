#pragma once
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <iostream>
#include <Eigen/Dense>
#include <stdexcept>
#include <vector>
namespace bg = boost::geometry;

using GeoPoint = bg::model::d2::point_xy<double>;
using GeoLine = bg::model::linestring<GeoPoint>;
using GeoPolygon = bg::model::polygon<GeoPoint>;

extern const double chi95;
extern const double delta0;
inline const Eigen::Vector2d tankCenter{1050.0, 750.0};
inline const double tankRadius{120.0};

enum class ObservationKind{
  Distance,
  Angle
};

struct Observation{
  std::string label;
  ObservationKind kind;
  Eigen::VectorXd designRow;
  double sigma;
};

struct Station{
  std::string name;
  Eigen::Vector2d position;
  bool fixed;
};

struct Ellipse{
  double majorAxis;
  double minorAxis;
  double orientation;
};


// function that adds a station to the list of stations in the design 
// template required in order to function to deduce num of elements in stationPtrs
template <std::size_t N>
std::vector<Station> makeStations(const std::vector<Station>& addedStations, const Station* const (&stationPtrs)[N]){
  std::vector<Station> stations;

  for (const Station* existingStations : stationPtrs){
    stations.push_back(*existingStations);
  };

  for (const Station& added : addedStations){
    if (added.fixed){
      throw std::invalid_argument("New stations cannot be fixed: " + added.name);
    };

    for (const Station& existingStations : stations){
      if (added.name == existingStations.name){
        throw std::invalid_argument("Duplicate station name" + added.name);
      };

    };
    stations.push_back(added);
  };
  return stations;
};

// Overload to compute preliminary analysis
template <std::size_t N>
std::vector<Station> makeStations(const Station* const (&stationPtrs)[N]) {
    std::vector<Station> stations;

    for (const Station* existing : stationPtrs) {
        stations.push_back(*existing);
    }

    return stations;
}

bool visible(const Eigen::Vector2d& a, const Eigen::Vector2d& b, const GeoPolygon& building);


Eigen::RowVectorXd designRow(const std::vector<Station>& s, int i, int j, bool isAzimuth);

void analyze(const std::string& title, const std::vector<Station>& s, int measurementSets, const std::vector<Eigen::Vector2d>& buildingPoints, const GeoPolygon& building);

void printEllipse(const std::string& label, const Eigen::Matrix2d& c);