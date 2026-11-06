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

## Iteration Analysis
*Iteration 1: P3 Added*
One round of observations, and a third point added at (1600.0, 1100.0). Significantly better minimum redundancy number, and MDBs max shift is 12.4 rather than 13.3. P1 95% absolute error ellipse major axis parameter did not decrease significantly. Neither did P1-P2 95% relative error ellipse major axis parameter. Continue with iterations.

*Iteration 2: P3 Added, 3 rounds of observations*
42 total measurements. P1 95% absolute error ellipse major axis parameter and P1-P2 95% relative error ellipse major axis parameter decreased significantly; 7.8 -> 4.5 / 6.3 -> 3.6. Absolute error of P1 still unacceptable, although maximum P1 shift from minimum detectable blunder has decreased significantly; 12.4 -> 7.2. Continue with iterations. 

*Iteration 3: P4 Added*
26 total measurements. Larger minimum redundancy number and smaller maximum P1 shift from MDB with less measurements would indicate a better designed network. Unfortunately, the major ellipse parameters are larger than the previous iteration, evident when inspecting covariance of parameters (all values seem to be larger by around 2x-3x). Continue with iterations.

*Iteration 4: P4 Added, 2 rounds of observations*
52 total measurements. Larger minimum redundancy and smaller maximum P1 shift from MDB like previous iteration, but significantly smaller elements in covariance matrix (comparable to iteration 2). Major axis parameters still not acceptable.

*Iteration 5: P3 added only, 4 rounds of observations"
Same as iteration 2 but with slightly smaller major axis elements. Not sure whether to commit to adding more rounds of observations or adding another station. Still not acceptable values.

*Iteration 6: P5 added*
48 total measurements. Major axis parameter of absolute P1 error ellipse still too large, but first iteration to have an acceptable relative P1-P2 error ellipse. Larger redundancy values and the smallest P1 shift from MDB indicates the strongest network so far. Added another round of observations for final iteration.

*Iteration 7: P5 added, 2 rounds of observations*
Same pattern as before, slightly smaller error ellipse parameters, but it seems to be impossible to create a network that can support an absolute error ellipse that is acceptable for P1.

*Iteration 8: P3 only, 7 rounds of observations*
Error ellipse specifications are finally met. Unholy 98 total measurements, but minimal points added to network. Will use as final network. 

## Example of a non-optimal design 
The network with only a single added point at (1600.0, 1100.0) is insufficient

## Reasoning behind conclusions
<fill in later> 