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


constexpr double PI = 3.1415926535897932384626433832795;
constexpr const double SIGMA_DIST = 0.002;
constexpr const double SIGMA_ANGLE = 2.0 * PI / 648000.0;
constexpr double TOLERANCE = 1e-89;

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


void analyze(const std::string& title, const std::vector<Station>& s, 
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
};