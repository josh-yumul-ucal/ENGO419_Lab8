#include "functions.h"
#include <iostream>
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <vector>
#include <string>
#include <fstream>


constexpr double PI = 3.1415926535897932384626433832795;
constexpr const double SIGMA_DIST = 0.002;
constexpr const double SIGMA_ANGLE = 2.0 * PI / 648000.0;
constexpr double TOLERANCE = 1e-89;
constexpr double P1_MAXIMUM_VALUE = 3.0;
constexpr double RELATIVE_MAXIMUM_VALUE = 5.0;

namespace bm = boost::math;
namespace bg = boost::geometry;

// Forward declaration of printEllipse cause I did it after I made analyze
void printEllipse(const std::string& label, const Eigen::Matrix2d& Cx);

using GeoPoint = bg::model::d2::point_xy<double>;
using GeoLine = bg::model::linestring<GeoPoint>;
using GeoPolygon = bg::model::polygon<GeoPoint>;

// Statistical values
const double chi95 = sqrt(boost::math::quantile(
  boost::math::chi_squared(2), 0.95));

const double delta0 = bm::quantile(bm::normal(), 0.975) + bm::quantile(bm::normal(), 0.80);

bool visible(const Eigen::Vector2d& a, const Eigen::Vector2d& b, const GeoPolygon& building){
  GeoLine sight{
    {a.x(), a.y()},
    {b.x(), b.y()}
  };

  GeoPoint tankPoint(tankCenter.x(), tankCenter.y());

  return bg::distance(tankPoint, sight) >= tankRadius && !bg::crosses(sight, building);
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


designResult analyze(const std::string& title, const std::vector<Station>& s, 
  int measurementSets, const std::vector<Eigen::Vector2d>& buildingPoints, const GeoPolygon& building){
  int n = s.size();
  std::vector<Eigen::RowVectorXd> rows;
  std::vector<double> sigmas;

  // distances A-B is skipped because they're control points
  for (int i = 0; i < n; i++){
    for (int j = std::max(i + 1, 2); j < n; j++){
      if (visible(s[i].position, s[j].position, building)){
        rows.push_back(designRow(s, i, j, false));
        sigmas.push_back(SIGMA_DIST);
      }
    }
  }

  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      for (int k = j + 1; k < n; k++){
        if (j != i && k != i && visible(s[i].position, s[j].position, building) && visible(s[i].position, s[k].position, building)){
          rows.push_back(designRow(s, i, k, true) - designRow(s, i, j, true));
          sigmas.push_back(SIGMA_ANGLE);
        }
      }
    }
  }

  Eigen::MatrixXd A(rows.size(), 2 * (n - 2));
  Eigen::VectorXd sigma(rows.size());
  for (size_t i = 0; i < rows.size(); i++){
    A.row(i) = rows[i];
    sigma(i) = sigmas[i] / sqrt(double(measurementSets));
  }

  // Standard covariance computations to find covariances of parameters
  Eigen::VectorXd weights = sigma.array().square().inverse();
  Eigen::MatrixXd Cx = (A.transpose() * weights.asDiagonal() * A).inverse();

  // Relative ellipse stuff
  Eigen::Matrix2d C11 = Cx.block<2, 2>(0, 0), C22 = Cx.block<2, 2>(2,2), C12 = Cx.block<2, 2>(0, 2);
  Eigen::MatrixXd Crel = C11 + C22 - C12 - C12.transpose();

  // Redundancy number computations
  double minR = 1, maxShift = 0;
  for (int i = 0; i < A.rows(); i++){
    Eigen::RowVectorXd a = A.row(i);
    // need (0, 0) to convert matrix multiplicationt to scalar value or else type error
    double r = 1.0 - (a * Cx * a.transpose())(0, 0) * weights(i); // a * Cx * a.T / std

    // Store the smallest R
    minR = std::min(minR, r);
    double mdb = delta0 * sigma(i) / sqrt(r);
    Eigen::VectorXd shift = Cx * a.transpose() * mdb * weights(i);
    maxShift = std::max(maxShift, shift.head<2>().norm());
  }

  // Corner = covered if any station can see it
  int cornersSeen = 0;
  for (const Eigen::Vector2d& corners : buildingPoints){
    bool seen = false;
    for (const Station& station : s){
      seen = seen || visible(station.position, corners, building);
    };
    cornersSeen += seen;
  }

  // Print statements
  std::cout << "\n" << title << "\n";
  std::cout << "Observations: " << A.rows() << ", unknowns: " << A.cols() << ", readings: " << A.rows() * measurementSets << "\n";
  std::cout << "Corners Visible: " << cornersSeen << " of " << buildingPoints.size() << "\n";
  std::cout << "Minimum redundancy: " << minR << "\n";
  std::cout << "Maximum P1 shift from an MDB [mm]: " << 1000 * maxShift << "\n";
  std::cout << "Cx [m^2]: \n" << Cx << "\n";

  // Ellipse parameters print statements
  printEllipse("95% P1 Ellipse", C11);
  printEllipse("95% P1-P2 relative ellipse", Crel);

  return {title, s, measurementSets, static_cast<int>(A.rows()), static_cast<int>(A.cols()), cornersSeen, 
    minR, maxShift, Cx, Crel};
};


void printEllipse(const std::string& label, const Eigen::Matrix2d& Cx){
  double sigmaX2 = Cx(0, 0);
  double sigmaY2 = Cx(1, 1);
  double sigmaXY = Cx(0, 1);

  // Error ellipse computations
  double a2 = 0.5 * (sigmaX2 + sigmaY2 + sqrt(pow(sigmaX2 - sigmaY2, 2) + 4 * pow(sigmaXY, 2)));
  double b2 = 0.5 * (sigmaX2 + sigmaY2 - sqrt(pow(sigmaX2 - sigmaY2, 2) + 4 * pow(sigmaXY, 2)));

  // 95% CI Scale
  double Smax = 1000 * chi95 * sqrt(a2);
  double Smin = 1000 * chi95 * sqrt(b2);

  double orientation = 0.5 * atan2(2 * sigmaXY, sigmaX2 - sigmaY2) * 180 / PI;
  
  // Print to terminal
  std::cout << label << " major/minor (mm): " << Smax << " / " << Smin << ", orientation (deg): " << orientation << "\n";
}

static const char* yesNo(bool value){
  return value ? "yes" : "no";
};

Ellipse computeEllipse(const Eigen::Matrix2d &Cx)
{
  double sigmaX2 = Cx(0, 0);
  double sigmaY2 = Cx(1, 1);
  double sigmaXY = Cx(0, 1);

  // Error ellipse computations
  double a2 = 0.5 * (sigmaX2 + sigmaY2 + sqrt(pow(sigmaX2 - sigmaY2, 2) + 4 * pow(sigmaXY, 2)));
  double b2 = 0.5 * (sigmaX2 + sigmaY2 - sqrt(pow(sigmaX2 - sigmaY2, 2) + 4 * pow(sigmaXY, 2)));

  // 95% CI Scale
  double Smax = 1000 * chi95 * sqrt(a2);
  double Smin = 1000 * chi95 * sqrt(b2);

  double orientation = 0.5 * atan2(2 * sigmaXY, sigmaX2 - sigmaY2) * 180 / PI;

  return {Smax, Smin, orientation};
}

std::ofstream openCSV(const std::string &path){
  std::ofstream file(path);
  if(!file){
    throw std::runtime_error("Could not write to " + path);
  };
  file << std::fixed << std::setprecision(8);
  return file;
  
};

void saveResults(const designResult &result, const GeoPolygon &building, const std::string name) {
  std::ofstream summary = openCSV(name + "_summary.csv");
  std::ofstream ellipses = openCSV(name + "_ellipses.csv");
  std::ofstream stations = openCSV(name + "_stations.csv");
  std::ofstream sightlines = openCSV(name + "_sightlines.csv");

  // CSV headers
  summary << "design,added_stations,measurement_sets,observations,unknowns,readings,"
        << "corners_seen,min_redundancy,max_p1_shift_mm,"
        << "p1_major_mm,p1_minor_mm,p1_orientation_deg,"
        << "rel_major_mm,rel_minor_mm,rel_orientation_deg,"
        << "meets_p1_spec,meets_rel_spec\n";

  ellipses << "design,name,center_x,center_y,major_mm,minor_mm,orientation_deg\n";
  stations << "design,name,x,y,fixed\n";
  sightlines << "design,station_a,station_b\n";

  // Design results
  const std::string design = "\"" + result.title + "\"";
  const std::vector<Station>& s = result.stations;


  std::string added;
  for (size_t i = 4; i < s.size(); i++){
    added += (added.empty() ? "" : " " ) + s[i].name;
  };

  Eigen::Matrix2d C11 = result.Cx.block<2, 2>(0,0);
  Ellipse P1 = computeEllipse(C11);
  Ellipse relative = computeEllipse(result.Crel);

  summary << design << ","
  << (added.empty() ? "none" : added) << ","
  << result.measurementSets << ","
  << result.observations << ","
  << result.unknowns << ","
  << result.observations * result.measurementSets << ","
  << result.cornersSeen << ","
  << result.minRedundancy << ","
  << 1000 * result.maxP1Shift << ","
  << P1.majorAxis << ","
  << P1.minorAxis << ","
  << P1.orientation << ","
  << relative.majorAxis << ","
  << relative.minorAxis << ","
  << relative.orientation << ","
  << yesNo(P1.majorAxis < P1_MAXIMUM_VALUE) << ","
  << yesNo(relative.majorAxis < RELATIVE_MAXIMUM_VALUE) << "\n";


  for (size_t i = 0; i < s.size(); i++){
    stations << design << ","
    << s[i].name << ","
    << s[i].position.x() << ","
    << s[i].position.y() << ","
    << yesNo(s[i].fixed) << "\n";

    for (size_t j = i + 1; j < s.size(); j++){
      if (visible(s[i].position, s[j].position, building)){
        sightlines << design << "," 
        << s[i].name << ","
        << s[j].name << "\n";
      }
    }
    // P1 error ellipse parameters
    if (i >= 2){
    Eigen::Matrix2d C = result.Cx.block<2, 2>(2 * (i - 2), 2 * (i - 2));

    Ellipse e = computeEllipse(C);
    ellipses << design << ","
    << s[i].name << ","
    << s[i].position.x() << ","
    << s[i].position.y() << ","
    << e.majorAxis << ","
    << e.minorAxis << ","
    << e.orientation << "\n";
    }
  }

  // Relative error ellipse parameters
  Eigen::Vector2d midpoint =
    0.5 * (s[2].position + s[3].position);

  ellipses << design << ",P1-P2,"
          << midpoint.x() << ","
          << midpoint.y() << ","
          << relative.majorAxis << ","
          << relative.minorAxis << ","
          << relative.orientation << "\n";

}
