
#include <iostream>
#include <fstream>
#include <math.h>
#include <cstdlib>
#include <fstream>
#include <algorithm>
#include <vector>
#include <numeric>
#include <stdio.h>
#include <cstdlib>
using namespace std ;
#define pi 3.1415



#define thres_lateral 2.0 // 2.4 0r 2.0 for plant cells and 0.5 for wing disc
#define searchAreaForNghbr 20.0 * thres_lateral
#define thres_corners 2.0 *thres_lateral        //to find a vertices for cells who does not have a comon neighbor
#define thres_cyclic4  0.7 * thres_lateral      //when 4 cells are neighbors instead of 3 which is more common
#define thres_intersect 0.4 *thres_lateral      // vertices shared between 3 cells (1.1 for wingDisc 0.4 for plant)
#define thres_refineNodes 2 * thres_lateral   //used to delete close boundary nodes and xNew nodes
#define thres_noBoundary 20
#define thres_Production 20
#define thres_layer 3             //CLV3 layer 3+



double Dist2D (double x1 , double y1 , double x2 , double y2 ) ;
vector<double> Dist_PointToVec2D ( double x1, double y1, vector<double> X2, vector<double> Y2 ) ;
vector<vector <double> > Dist_VecToVec2D ( vector<double> , vector<double> , vector<double>, vector<double>  ) ;
vector<int> Indices_MinMatrix (vector<vector<double > >) ;
vector<int> Indices_MaxMatrix (vector<vector<double> >  ) ;
vector<double> Dist_pointToVec1D ( double x1, vector<double> X2) ;

double DotProduct (double x1 , double y1 , double x2 , double y2) ;
double AngleOfTwoVectors (double x1 , double y1 , double x2 , double y2) ;
double MagnitudeVec (double x , double y) ;
double TriangleArea (double x1, double y1, double x2 , double y2) ;
double sum_over_vec(const vector<vector<double> >& v, int a ) ;
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------

class plusNum
{
private:
    double num  ;
public:
    plusNum(double n) : num(n) {}
    
    
    // This operator overloading enables calling
    // operator function () on objects of increment
    double operator () (double arr_num) const
    {
        return num + arr_num;
    }
} ;
//---------------------------------------------------------------------------------------------
class productNum
{
private:
    double num ;
public:
    productNum(double n) : num(n) {}
    
    // This operator overloading enables calling
    // operator function () on objects of increment
    double operator () (double arr_num) const
    {
        return num * arr_num;
    }
} ;

class linearConfig
{
private:
    double a ;
    double b ;
public:
    linearConfig(double n , double m )
    {
        a = n ;
        b = m ;
    }
    
    // This operator overloading enables calling
    // operator function () on objects of increment
    double operator () (double arr_num1 , double arr_num2) const
    {
        return a * arr_num1 + b * arr_num2 ;
    }
} ;

class productVec
{
private:
    
public:
    // This operator overloading enables calling
    // operator function () on objects of increment
    double operator () (double arr_num1 , double arr_num2) const
    {
        return  arr_num1 * arr_num2 ;
    }
} ;

double NormalCDFInverse2(double p);
double RationalApproximation2(double t);


