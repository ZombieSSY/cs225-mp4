/**
 * @file kdtree.cpp
 * Implementation of KDTree class.
 */
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

template <int Dim>
bool KDTree<Dim>::smallerDimVal(const Point<Dim>& first,
                                const Point<Dim>& second, int curDim) const
{
    /**
     * @todo Implement this function!
     */
  if (first[curDim] < second[curDim])
    return true;
  else if (first[curDim] > second[curDim])
    return false;
  else
    return first < second; // If tie, return by operator<
}

template <int Dim>
bool KDTree<Dim>::shouldReplace(const Point<Dim>& target,
                                const Point<Dim>& currentBest,
                                const Point<Dim>& potential) const
{
    /**
     * @todo Implement this function!
     */
  int dimension = Dim;
  return shouldReplace(target, currentBest, potential, dimension);
}

// Helper function for "shouldReplace"
template <int Dim>
bool KDTree<Dim>::shouldReplace(const Point<Dim>& target,
                                const Point<Dim>& currentBest,
                                const Point<Dim>& potential,
				int curDim) const
{
  // Calculation for distance
  int currDistance = 0;
  int potDistance = 0;
  for (int i = 0; i < curDim; i++){
    int d1 = pow(target[i] - currentBest[i], 2);
    int d2 = pow(target[i] - potential[i], 2);
    currDistance += d1;
    potDistance += d2;
  }

  // Check and return a boolean
  if (currDistance < potDistance)
    return false;
  else if (currDistance > potDistance)
    return true;
  else
    return potential < currentBest; // If tie, return by operator<
}

template <int Dim>
KDTree<Dim>::KDTree(const vector<Point<Dim>>& newPoints)
{
    /**
     * @todo Implement this function!
     */
  points = newPoints;
  if (points.size() != 0){
    int n = points.size();
    KDTreeHelper(points, 0, n-1, 0);
  }
}

// Helper function for KDTree
template <int Dim>
void KDTree<Dim>::KDTreeHelper(vector<Point<Dim>>& newPoints, int a, int b, int curDim)
{
  if (a > b)
    return;
  int mid = (a + b)/2;    
  MedianSelect(newPoints, a, b, mid, curDim);//Build the root for KDTree
  // Continue by recursion
  if (a < mid)
    KDTreeHelper(newPoints, a, mid - 1, (curDim + 1) % Dim);
  if (mid < b)
    KDTreeHelper(newPoints, mid + 1, b, (curDim + 1) % Dim);
}

// Implementation for the quickselect
template <int Dim>
Point<Dim> KDTree<Dim>::MedianSelect(vector<Point<Dim>>& newPoints, int left, int right, int pos, int curDim)
{
  int mid = (left + right)/2;
  if (left == right)
    return newPoints[left];

  while (left < right){
    int part = Partition(newPoints, left, right, mid, curDim);
    if (part == pos)
      return newPoints[part];
    else if (pos > part)
      left = part + 1;
    else
      right = part - 1;
  }
  return newPoints[left];
}

// Seperate the vector
// With the points in the left side all less than the mid, and the points in the right side all bigger than the mid
template <int Dim>
int KDTree<Dim>::Partition(vector<Point<Dim>>& newPoints, int left, int right, int mid, int curDim)
{
  Point<Dim> k = newPoints[mid];
  int leftMark = left;
  int rightMark = right;
  
  iter_swap(newPoints.begin() + mid, newPoints.begin() + rightMark);
  for (int i = left; i < right; i++){
    if (newPoints[i] == k || smallerDimVal(newPoints[i], k, curDim)){
      iter_swap(newPoints.begin() + i, newPoints.begin() + leftMark);
      leftMark += 1;
    }
  }
  iter_swap(newPoints.begin() + leftMark, newPoints.begin() + right);

  return leftMark;
}

// Implementation for findNearestNeighbor
// Start from the middle point, which is the root
template <int Dim>
Point<Dim> KDTree<Dim>::findNearestNeighbor(const Point<Dim>& query) const
{
    /**
     * @todo Implement this function!
     */
  int a = 0;
  int b = points.size() - 1;
  int mid = (a + b)/2;
  Point<Dim> middle = points[mid];
  return checkNearest(query, middle, points, a, b, 0);
}

// Helper function for the NNS
// Compare the current point with the provided point to find the closest point
template <int Dim>
Point<Dim> KDTree<Dim>::checkNearest(const Point<Dim>& query, const Point<Dim>& current, const vector<Point<Dim>> newPoints, int a, int b, int curDim) const
{
  // If a == b, return new one directly
  if (a == b){
    if (shouldReplace(query, current, newPoints[a]))
      return newPoints[a];
    else
      return current;
  }
  
  Point<Dim> best = current;
  int mid = (a + b)/2;
  // To check whether needs to replace the current point
  if (smallerDimVal(newPoints[mid], query, curDim) && b > mid)
    best = checkNearest(query, current, newPoints, mid + 1, b, (curDim + 1) % Dim);  
  if (smallerDimVal(query, newPoints[mid], curDim) && a < mid)
    best = checkNearest(query, current, newPoints, a, mid - 1, (curDim + 1) % Dim);
  
  if (shouldReplace(query, best, newPoints[mid]))
    best = newPoints[mid];
  
  int d = pow((newPoints[mid][curDim] - query[curDim]), 2);
  int r = 0;
  for (int i = 0; i < Dim; i++)
    r += pow(query[i] - best[i], 2);
  // Compare the distance
  if (d <= r){
    // Go to the right;
    if (mid < b)
      best = checkNearest(query, best, newPoints, mid + 1, b, (curDim + 1) % Dim);
    // Go to the left;
    if (mid > a)
      best = checkNearest(query, best, newPoints, a, mid - 1, (curDim + 1) % Dim);
  }
  return best;
}
