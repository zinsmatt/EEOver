#include <array>
#include <cmath>
#include "solvers.h"


#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/foreach.hpp>
#include <boost/geometry/geometries/adapted/c_array.hpp>
typedef boost::geometry::model::d2::point_xy<double, boost::geometry::cs::cartesian> point_2d;
typedef boost::geometry::model::polygon<point_2d> polygon_2d;
using namespace boost::geometry;


#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;


using Ellipse = std::array<double, 5>;

int ellipse2poly(float PHI_1, float A1, float B1, float H1, float K1, polygon_2d * po)
{
    polygon_2d  poly;
    float xc = H1;
    float yc = K1;
    float a = A1;
    float b = B1;
    float w = PHI_1;
    
    const int N = 20; // Default discretization
    float t = 0;
    int i = 0;
    double coor[N*2+1][2];
  
    double x, y;
    float step = M_PI/N;
    float sinphi = sin(w);
    float cosphi = cos(w);
    for(i=0; i<2*N+1; i++)
    {   
        x = xc + a*cos(t)*cosphi - b*sin(t)*sinphi;
        y = yc + a*cos(t)*sinphi + b*sin(t)*cosphi;
        if(fabs(x) < 1e-4) x = 0;
        if(fabs(y) < 1e-4) y = 0;

        coor[i][0] = x;
        coor[i][1] = y;
        t += step;
    }

    assign_points(poly, coor);
    correct(poly);
    *po = poly;
    return 1;
}

// overlaping area of two polygons
float getOverlapingAreaPoly(polygon_2d poly, polygon_2d poly2)
{
    float overAreaPoly = 0.0;    
    std::deque<polygon_2d> output;
    bool ret = boost::geometry::intersection(poly, poly2, output);
    if(!ret) {
         std::cout << "Could not calculate the overlap of the polygons\n";
         exit(-1);
    }
    BOOST_FOREACH(polygon_2d const& p, output)
    {
        overAreaPoly = boost::geometry::area(p);
    }
    return overAreaPoly;
}

double compute_overlap_boost(const Ellipse& ell0, const Ellipse& ell1)
{
    polygon_2d poly0, poly1;
    ellipse2poly(ell0[2], ell0[0], ell0[1], ell0[3], ell0[4], &poly0);
    ellipse2poly(ell1[2], ell1[0], ell1[1], ell1[3], ell1[4], &poly1);
    double inter = getOverlapingAreaPoly(poly0, poly1);
    return inter;
}

double compute_iou_boost(const Ellipse& ell0, const Ellipse& ell1)
{
    double area0 = ell0[0]*ell0[1]*M_PI;
    double area1 = ell1[0]*ell1[1]*M_PI;
    double inter = compute_overlap_boost(ell0, ell1);
    return inter / (area0 + area1 - inter);
}


double compute_iou_toms(const Ellipse& ell0, const Ellipse& ell1)
{
    double x[4], y[4];
    int n_roots;
    int ret;
    // first try with TOMS
    double inter = ellipse_ellipse_overlap_netlibs(ell0[2], ell0[0], ell0[1], ell0[3], ell0[4], 
                                                   ell1[2], ell1[0], ell1[1], ell1[3], ell1[4],
                                                   x, y, &n_roots, &ret);
    double area0 = ell0[0]*ell0[1]*M_PI;
    double area1 = ell1[0]*ell1[1]*M_PI;
    if (inter < 0 || inter*2 > area0 + area1)
    {
        inter = compute_overlap_boost(ell0, ell1);
    }
    return inter / (area0 + area1 - inter);
}

double compute_iou_gsl(const Ellipse& ell0, const Ellipse& ell1, int choice)
{
    double x[4], y[4];
    int n_roots;
    int ret;
    // first try with GSL
    double inter = ellipse_ellipse_overlap_gsl(ell0[2], ell0[0], ell0[1], ell0[3], ell0[4], 
                                               ell1[2], ell1[0], ell1[1], ell1[3], ell1[4],
                                               x, y, &n_roots, &ret, choice);
    double area0 = ell0[0]*ell0[1]*M_PI;
    double area1 = ell1[0]*ell1[1]*M_PI;
    if (inter < 0 || inter*2 > area0 + area1)
    {
        inter = compute_overlap_boost(ell0, ell1);
    }
    return inter / (area0 + area1 - inter);
}



double compute_iou_gems(const Ellipse& ell0, const Ellipse& ell1)
{
    double x[4], y[4];
    int n_roots;
    int ret;
    // first try with GEMS
    double inter = ellipse_ellipse_overlap_gems(ell0[2], ell0[0], ell0[1], ell0[3], ell0[4], 
                                                ell1[2], ell1[0], ell1[1], ell1[3], ell1[4],
                                                x, y, &n_roots, &ret);
    double area0 = ell0[0]*ell0[1]*M_PI;
    double area1 = ell1[0]*ell1[1]*M_PI;
    if (inter < 0 || inter*2 > area0 + area1)
    {
        inter = compute_overlap_boost(ell0, ell1);
    }
    return inter / (area0 + area1 - inter);
}



PYBIND11_MODULE(ellipses_iou, m) {
    m.doc() = R"pbdoc(
        Ellipses IoU
        -----------------------

        .. currentmodule:: Ellipses_IoU

        .. autosummary::
           :toctree: _generate
			compute_iou
    )pbdoc";

    m.def("compute_iou_toms", &compute_iou_toms, R"pbdoc(
        Compute the IoU between two ellipses with TOMS. (probably the best solver)
    )pbdoc");
    m.def("compute_iou_boost", &compute_iou_boost, R"pbdoc(
        Compute the IoU between two ellipses with BOOST. The ellipses are turned into polygons
    )pbdoc");
    m.def("compute_iou_gsl", &compute_iou_gsl, R"pbdoc(
        Compute the IoU between two ellipses with GSL. The third parameter can be 1 or 2
        depending on the version you want to use. "1" is for the default and "2" is faster
        but less reliable.
    )pbdoc");
    m.def("compute_iou_gems", &compute_iou_gems, R"pbdoc(
        Compute the IoU between two ellipses with GEMS. (probably the worst one)
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
