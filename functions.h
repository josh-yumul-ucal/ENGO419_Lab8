#pragma once
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <iostream>
#include <Eigen/Dense>
#include <vector>

extern const double chi95;
extern const double delta0;
extern const Eigen::Vector2d tankCenter{1050.0, 750.0};
extern const double tankRadius{120.0};


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

bool visible(const Eigen::Vector2d& a, const Eigen::Vector2d& b);

struct Network{
  std::vector<Station> stations;
  std::vector<std::pair<int, int>> visibleLines;
  std::vector<Observation> observations;
  Eigen::MatrixXd A;
}; 


struct Analysis{
  Eigen::MatrixXd Cx;
  Ellipse point1;
  Ellipse relative12;
  bool coversCorners = false;
  bool passesPrecision = false;
  double minimumRedundancy = 1.0;
  double maximumPoint1Shift = 0.0;
};

Eigen::RowVectorXd designRow(const std::vector<Station>& s, int i, int j, bool isAzimuth);

void analyze(const string& title, const std::vector<Station>& s, int iterations);