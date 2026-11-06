# ENGO419 Lab 8

Written in C++
Completed by:
- Josh Yumul
- Tyler Birkbeck

## Goal(s) and Reasoning behind network design
Goal: Complete the network with minimum additional points, with lab specified accuracy requirements

We do not use inner constraints to solve the datum defect because there are already reliable control points in points A and B, which we can treat as points with variances of zero, this also improves the absolute error ellipses of new points 1 and 2 (not relative error ellipses between 1 and 2 as that is reliant on the relative uncertainty between the coordinate differences of points 1 and 2)

## Determining visibility of building points
Using Boost.Geometry library to simulate the geometry, and then using its objects to represent the 2D area and then tests for obstruction.

### GeoPoint
```cpp
using GeoPoint = bg::model::d2::point_xy<double>;
```
represents a single (x, y) coordinate

### GeoLine
```cpp
using GeoLine = bg::model::linestrin<GeoPoint>;
```
represents connected line segments between points

### GeoPolygon
```cpp
using GeoPolygon = bg::model::polygon<GeoPoint>;
```
represents a closed boundary around an area defined by corners (used for the building)

### Explaining bool visible()
It computes a GeoLine using two coordinates called <sight>, and then it uses the tankCenter and tank radius to represent the tank, and a GeoPolygon to represent the building. The line: 
```cpp
return bg::distance(tankPoint, sight) >= tankRadius && !bg::crosses(sight, building);
```
Checks whether the line of sight between two coordinates is obstructed by either the building or the tank. No obstruction from either object means that the two coordinates have a line of sight to one another. This function is used within analyze() to compute which building corners are visible to the point network. 

## Example of a non-optimal design 
The network with only a single added point at (1600.0, 1100.0) is insufficient

## Reasoning behind conclusions
<fill in later> 