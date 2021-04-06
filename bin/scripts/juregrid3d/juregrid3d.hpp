//
// Copyright 2015 by Forschungszentrum Juelich GmbH
// author: Joern Ungermann
//
#include <iostream>
#include <memory>
#include <array>
#include <vector>
#include <string>
#include <exception>
#include <limits>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Delaunay_triangulation_3.h>

#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include <CGAL/natural_neighbor_coordinates_3.h>
#include <CGAL/Location_policy.h>

typedef std::vector<std::pair<int, double>> IndiceWeights;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT CoordType;

typedef Kernel::Point_3 Point3;
typedef CGAL::Triangulation_vertex_base_with_info_3<int, Kernel> Vb3;
typedef CGAL::Triangulation_data_structure_3<Vb3> Tds3;
typedef CGAL::Delaunay_triangulation_3<Kernel, Tds3, CGAL::Fast_location> DelaunayTriangulation3;
//typedef CGAL::Delaunay_triangulation_3<Kernel, Tds3> DelaunayTriangulation3;
typedef DelaunayTriangulation3::Vertex_handle VertexHandle3;
typedef DelaunayTriangulation3::Cell_handle CellHandle3;

constexpr double Earth_Mean_Radius = 6371.23;

#define jassert(cond, msg)                       \
    if (!(cond))                                 \
    {                                            \
        std::stringstream message;               \
        message << msg;                          \
        throw std::runtime_error(message.str()); \
    }

/// Converts from degree to radians
inline constexpr double deg2rad(double x)
{
    return x * M_PI / 180.;
}

/// Converts from radians to degree
inline constexpr double rad2deg(double x)
{
    return x * M_1_PI * 180.;
}

inline double fixLongitude(double lon)
{
    if (lon >= 180)
        lon -= 360;
    if (lon < -180)
        lon += 360;
    jassert((-180 <= lon && lon < 180),
            "lon not legal?!" << lon);
    return lon;
}

inline double fixLatitude(double lat)
{
    if (lat > 90)
        lat = 180 - lat;
    if (lat < -90)
        lat = -180 - lat;
    jassert((-90 <= lat && lat <= 90),
            "lat not legal?!" << lat);
    return lat;
}

// *********************************************************************************
//

struct CoordTransform
{
    /// Converts from lon/lat to another grid
    virtual Point3 convert(const Point3 &point) const = 0;
    virtual bool wrap() const = 0;
    virtual ~CoordTransform(){};

    /// Converts longitude and latitude to a stereographic projection
    static std::pair<double, double>
    convertLonLatToStereographic(
        double lon, double lat, double lon0, double lat0)
    {
        lon = deg2rad(lon);
        lat = deg2rad(lat);
        lon0 = deg2rad(lon0);
        lat0 = deg2rad(lat0);
        double r = 2 * Earth_Mean_Radius / (1. + sin(lat0) * sin(lat) + cos(lat0) * cos(lat) * cos(lon - lon0));
        return std::make_pair(
            r * cos(lat) * sin(lon - lon0),
            r * (cos(lat0) * sin(lat) - sin(lat0) * cos(lat) * cos(lon - lon0)));
    }
};

struct CoordTransformStereographic : public CoordTransform
{
    double lon_center_;
    double lat_center_;
    bool wrapping_;
    double scale_;

    CoordTransformStereographic(double loce, double lace, bool wr, double scale) : lon_center_(loce), lat_center_(lace), wrapping_(wr), scale_(scale){};

    virtual Point3 convert(const Point3 &point) const
    {
        double lon = point[0];
        double lat = point[1];
        auto xy = CoordTransform::convertLonLatToStereographic(
            lon, lat, this->lon_center_, this->lat_center_);
        return Point3(xy.first, xy.second, point[2] * this->scale_);
    }

    virtual bool wrap() const
    {
        return this->wrapping_;
    }
};

struct CoordTransformStereographicLevels : public CoordTransform
{
    double lon_center_;
    double lat_center_;
    bool wrapping_;
    double distance_;
    std::vector<double> levels_;

    CoordTransformStereographicLevels(double loce, double lace, bool wr,
                                      const std::vector<double> &scales) : lon_center_(loce), lat_center_(lace), wrapping_(wr),
                                                                           distance_(scales[0]), levels_(scales.begin() + 1, scales.end()){};

    virtual Point3 convert(const Point3 &point) const
    {
        double lon = point[0];
        double lat = point[1];
        auto upper = std::lower_bound(this->levels_.begin(), this->levels_.end() - 1, point[2]);
        if (upper == this->levels_.begin() && *upper == point[2])
        {
            ++upper;
        }
        jassert(upper != this->levels_.begin(),
                "levels not encompassing " << *upper << " " << *(this->levels_.begin()) << " " << *(this->levels_.end() - 1) << " " << point[2]);
        auto lower = upper - 1;
        jassert(*lower <= point[2] && point[2] <= *upper,
                "levels not encompassing! " << this->levels_.front() << " " << *lower << " " << point[2] << " " << *upper << " " << this->levels_.back());
        double alt = (lower - this->levels_.begin()) + (point[2] - *lower) / (*upper - *lower);
        auto xy = CoordTransform::convertLonLatToStereographic(
            lon, lat, this->lon_center_, this->lat_center_);
        return Point3(xy.first, xy.second, alt * this->distance_);
    }

    virtual bool wrap() const
    {
        return this->wrapping_;
    }
};

// *********************************************************************************
//

/// Different strategies to leverage the underlying Delaunay triangulation for
/// interpolation
struct InterpolationStrategy
{
    virtual IndiceWeights interpolate(
        const DelaunayTriangulation3 &T, const Point3 &point) const = 0;
    /// Factory function returing pointer to selcted strategy
    static InterpolationStrategy *create(const std::string strategy);
    virtual ~InterpolationStrategy(){};
};

/// Facilitates CGAL functionality to implement a natural neighbour interpolation
struct InterpolationNaturalNeighbour : public InterpolationStrategy
{
    /// Pointer to last found cell to speed up interpolation of closely spaced points.
    mutable CellHandle3 start3_;

    IndiceWeights interpolate(const DelaunayTriangulation3 &T, const Point3 &point) const
    {
        IndiceWeights result;
        if (T.dimension() == 3)
        { // skip empty triangulations that crash in the
            std::vector<std::pair<VertexHandle3, CoordType>> vertices;
            CoordType norm;

            CGAL::laplace_natural_neighbor_coordinates_3(
                T, point, std::back_inserter(vertices), norm,
                this->start3_ != CellHandle3() ? this->start3_ : CellHandle3());

            double normalize = 0;
            double wgt = -1;
            for (const auto &vertex : vertices)
            {
                normalize += vertex.second;
                if (vertex.second > wgt)
                {
                    wgt = vertex.second;
                    this->start3_ = vertex.first->cell();
                }
            }

            for (const auto &vertex : vertices)
            {
                result.emplace_back(vertex.first->info(), vertex.second / normalize);
            }

            if (result.empty())
            {
                std::cerr << "WARNING - no points. switching to nearest neighbour\n";
                auto vertex = T.nearest_vertex(point);
                result.emplace_back(vertex->info(), 1.);
                this->start3_ = vertex->cell();
            }
        }
        return result;
    }
};

/// Uses a Delaunay triangulation for identifying tetrahedrons and interpolationg
/// between the vertices using barymetric weights.
struct InterpolationBarycentric : public InterpolationStrategy
{
    /// Pointer to last found cell to speed up interpolation of closely spaced points.
    mutable CellHandle3 start3_;

    template <typename Array>
    static double dot(const Array &a, const Array &b)
    {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    /// Solves 3x3 system by Gaussian Elimination and Pivoting.
    static std::array<double, 3> solveThreeByThreeSystem(
        const std::array<double, 9> &matrix_orig, std::array<double, 3> b_orig)
    {
        std::array<double, 9> matrix(matrix_orig);
        std::array<double, 3> b(b_orig);
        std::array<int, 3> perm({{0, 1, 2}});
        if (std::abs(matrix[3]) > std::abs(matrix[0]))
            std::swap(perm[0], perm[1]);
        if (std::abs(matrix[6]) > std::abs(matrix[perm[0] * 3]))
            std::swap(perm[0], perm[2]);

        for (int i = 1; i < 3; ++i)
        {
            double fac = matrix[perm[i] * 3] / matrix[perm[0] * 3];
            for (int j = 0; j < 3; ++j)
            {
                matrix[perm[i] * 3 + j] -= matrix[perm[0] * 3 + j] * fac;
            }
            b[perm[i]] -= b[perm[0]] * fac;
        }
        if (std::abs(matrix[perm[2] * 3 + 1]) > std::abs(matrix[perm[1] * 3 + 1]))
        {
            std::swap(perm[1], perm[2]);
        }
        double fac = matrix[perm[2] * 3 + 1] / matrix[perm[1] * 3 + 1];
        matrix[perm[2] * 3 + 1] -= matrix[perm[1] * 3 + 1] * fac;
        matrix[perm[2] * 3 + 2] -= matrix[perm[1] * 3 + 2] * fac;
        b[perm[2]] -= b[perm[1]] * fac;

        std::array<double, 3> result;
        result[2] = b[perm[2]] / matrix[3 * perm[2] + 2];
        result[1] = (b[perm[1]] - matrix[3 * perm[1] + 2] * result[2]) / matrix[3 * perm[1] + 1];
        result[0] = (b[perm[0]] - matrix[3 * perm[0] + 2] * result[2] - matrix[3 * perm[0] + 1] * result[1]) / matrix[3 * perm[0]];

        // std::array<double, 3> err;
        // for (int i = 0; i < 3; ++i) {
        //     err[i] = matrix_orig[i * 3] * result[0] +
        //              matrix_orig[i * 3 + 1] * result[1] +
        //              matrix_orig[i * 3 + 2] * result[2] - b_orig[i];
        // }

        // jassert(! (isnan(result[0]) || isnan(result[1]) || isnan(result[2]) ||
        //            dot(err, err) > 1e-4),
        //         "linear algebra error");
        return result;
    }

    /// Determines barymetric weights of point in relation to four edges of tetrahedron
    static std::array<double, 4> getBarycentricTetrahedronWeights(
        const Point3 &x1, const Point3 &x2, const Point3 &x3, const Point3 &x4, const Point3 &x)
    {
        std::array<double, 9> matrix({{x1[0] - x4[0], x2[0] - x4[0], x3[0] - x4[0],
                                       x1[1] - x4[1], x2[1] - x4[1], x3[1] - x4[1],
                                       x1[2] - x4[2], x2[2] - x4[2], x3[2] - x4[2]}});
        std::array<double, 3> b({{x[0] - x4[0], x[1] - x4[1], x[2] - x4[2]}});
        std::array<double, 3> lambda = solveThreeByThreeSystem(matrix, b);

        return {{lambda[0], lambda[1], lambda[2],
                 1. - lambda[0] - lambda[1] - lambda[2]}};
    }

    IndiceWeights interpolate(const DelaunayTriangulation3 &T, const Point3 &point) const
    {
        IndiceWeights result(4);
        if (T.dimension() == 3)
        { // skip empty triangulations that crash in the
            DelaunayTriangulation3::Locate_type lt;
            int li = -1, lj = -1;
            auto cell_handle = T.locate(point, lt, li, lj, this->start3_);
            if (lt == DelaunayTriangulation3::Locate_type::CELL ||
                lt == DelaunayTriangulation3::Locate_type::EDGE ||
                lt == DelaunayTriangulation3::Locate_type::FACET ||
                lt == DelaunayTriangulation3::Locate_type::VERTEX)
            {
                std::vector<Point3> points;
                std::vector<int> infos;
                for (int i = 0; i < 4; ++i)
                {
                    assert(cell_handle->vertex(i)->info() >= 0);

                    points.push_back(cell_handle->vertex(i)->point());
                    infos.push_back(cell_handle->vertex(i)->info());
                }
                std::array<double, 4> lambda = getBarycentricTetrahedronWeights(
                    points[0], points[1], points[2], points[3], point);

                for (int i = 0; i < 4; ++i)
                {
                    result[i].first = infos[i];
                    result[i].second = lambda[i];
                }

                this->start3_ = cell_handle;
            }
            else if (lt == DelaunayTriangulation3::Locate_type::OUTSIDE_CONVEX_HULL)
            {
                std::cerr << "WARNING: point outside convex hull " << li << "\n";
                for (int i = 0; i < 4; ++i)
                {
                    if (!T.is_infinite(cell_handle->vertex(i)))
                    {
                        std::cerr << i << " " << cell_handle->vertex(i)->info() << "\n";
                    }
                    else
                    {
                        std::cerr << i << " infinite\n";
                    }
                }
                result.clear();
                auto vertex = T.nearest_vertex(point, cell_handle);
                result.emplace_back(vertex->info(), 1.);
            }
            else
            {
                std::cerr << "cgal error? " << lt << " " << li << " " << lt << "\n";
                result.clear();
            }
        }
        return result;
    }
};

/// Facilitates CGAL functionality to implement a nearest neighbour interpolation
class InterpolationNearestNeighbour : public InterpolationStrategy
{

    mutable VertexHandle3 vertex_;

public:
    InterpolationNearestNeighbour() {}

    IndiceWeights interpolate(
        const DelaunayTriangulation3 &T, const Point3 &point) const
    {
        IndiceWeights result;
        if (T.dimension() == 3)
        { // skip empty triangulations that crash
            this->vertex_ = T.nearest_vertex(point,
                                             this->vertex_ != VertexHandle3() ? this->vertex_->cell() : CellHandle3());
            assert(this->vertex_->info() >= 0);
            result.push_back(std::make_pair(this->vertex_->info(), 1.));
        }
        else
        {
            std::cerr << "empty triangulation!\n";
        }
        return result;
    }
};

/// Facilitates CGAL functionality to implement a nearest neighbour interpolation
class InterpolationWeightedDistance : public InterpolationStrategy
{

    mutable CellHandle3 start3_;

public:
    InterpolationWeightedDistance() {}

    IndiceWeights interpolate(
        const DelaunayTriangulation3 &T, const Point3 &point) const
    {
        IndiceWeights result(4);
        if (T.dimension() == 3)
        { // skip empty triangulations that crash in the
            DelaunayTriangulation3::Locate_type lt;
            int li = -1, lj = -1;
            auto cell_handle = T.locate(point, lt, li, lj, this->start3_);
            if (lt == DelaunayTriangulation3::Locate_type::CELL ||
                lt == DelaunayTriangulation3::Locate_type::EDGE ||
                lt == DelaunayTriangulation3::Locate_type::FACET ||
                lt == DelaunayTriangulation3::Locate_type::VERTEX ||
                lt == DelaunayTriangulation3::Locate_type::OUTSIDE_CONVEX_HULL)
            {
                std::vector<int> infos;
                std::vector<double> inv_distances;
                double normalize = 0;
                for (int i = 0; i < 4; ++i)
                {
                    if (!T.is_infinite(cell_handle->vertex(i)))
                    {
                        infos.push_back(cell_handle->vertex(i)->info());
                        double distance = 0;
                        for (int j = 0; j < 3; ++j)
                        {
                            double fnord = cell_handle->vertex(i)->point()[j] - point[j];
                            distance += fnord * fnord;
                        }
                        inv_distances.push_back(1. / std::max(1e-12, sqrt(distance)));
                        normalize += inv_distances.back();
                    }
                }
                for (size_t i = 0; i < infos.size(); ++i)
                {
                    result[i].first = infos[i];
                    result[i].second = inv_distances[i] / normalize;
                }
                this->start3_ = cell_handle;
            }
            else
            {
                std::cerr << "cgal error? " << lt << " " << li << " " << lj << "\n";
                result.clear();
            }
        }
        return result;
    }
};

InterpolationStrategy *InterpolationStrategy::create(
    const std::string strategy)
{
    if (strategy == "NATURALNEIGHBOUR")
    {
        return new InterpolationNaturalNeighbour;
    }
    else if (strategy == "BARYCENTRIC")
    {
        return new InterpolationBarycentric;
    }
    else if (strategy == "NEAREST")
    {
        return new InterpolationNearestNeighbour;
    }
    else if (strategy == "DISTANCE")
    {
        return new InterpolationWeightedDistance;
    }
    else
    {
        jassert(false, "Invalid strategy!" << strategy);
    }
}

// *********************************************************************************
//

/// Identifier and algorithms for interpolation in a specific segment of the
/// Earth surface.
class Patch
{
    double lon_min_;
    double lon_max_;
    double lat_min_;
    double lat_max_;

    /// Pointer to interpolation strategy
    std::unique_ptr<InterpolationStrategy> interp_;
    /// Pointer to a coordinate transformation from geo-coordinates to
    /// interpolation coordinates.
    std::unique_ptr<CoordTransform> trafo_;
    /// Pointer to Delaunay triangulation of patch points.
    std::unique_ptr<DelaunayTriangulation3> triang_;

public:
    inline bool isInitialized() const
    {
        return this->triang_.get() != 0;
    }

    inline bool isInPatch(double lon, double lat) const
    {
        return ((this->lat_min_ <= lat) && (lat <= this->lat_max_) &&
                (this->lon_min_ <= lon) && (lon <= this->lon_max_));
    }

    /// Lazily initializes a patch.
    void initPatch(const int size, const double *alts, const double *lons, const double *lats)
    {
        jassert(this->triang_.get() == 0,
                "Patch not initialized?!" << this->lon_min_ << " " << this->lon_max_ << " " << this->lat_min_ << " " << this->lat_max_);
        this->triang_.reset(new DelaunayTriangulation3);

        double lat_extend = std::abs(this->lat_max_ - this->lat_min_) * 0.25;
        double lon_extend = std::abs(this->lon_max_ - this->lon_min_) * 0.25;

        for (int ip = 0; ip < size; ++ip)
        {
            double lat = lats[ip];
            if ((this->lat_min_ - lat_extend <= lat) &&
                (lat <= this->lat_max_ + lat_extend))
            {
                double lon = fixLongitude(lons[ip]);
                for (int add = -360; add <= 360; add += 360)
                {
                    if ((this->trafo_->wrap() || add == 0) &&
                        (this->lon_min_ - lon_extend <= lon + add) &&
                        (lon + add <= this->lon_max_ + lon_extend))
                    {
                        auto point = this->trafo_->convert(Point3(lon + add, lat, alts[ip]));
                        auto handle = this->triang_->insert(point);
                        handle->info() = ip;
                    }
                }
            }
        }
    }

    /// Returns indice of irregular grid and weigths for interpolation
    IndiceWeights
    getPointsAndWeights(const Point3 &point) const
    {
        return this->interp_->interpolate(*(this->triang_),
                                          this->trafo_->convert(point));
    }

    Patch(double lomi, double loma, double lami, double lama,
          InterpolationStrategy *interp,
          CoordTransform *traf) : lon_min_(lomi), lon_max_(loma),
                                  lat_min_(lami), lat_max_(lama),
                                  interp_(interp), trafo_(traf)
    {
        if (this->lat_min_ > this->lat_max_)
            std::swap(this->lat_min_, this->lat_max_);
    }
};

// *********************************************************************************
//

/// 3D-Interpolation class that segments earth into patches that are each
/// interpolated separately using a stereographic interpolation routine.
/// Triangulation is performed lazily.
class Interpolate3D
{

    const int size_;
    const double *alts_;
    const double *lons_;
    const double *lats_;
    const std::vector<double> scales_;
    /// Disjoint list of patches the union of which shall cover the Earth.
    mutable std::vector<Patch> patches;

    void addPatch(
        const std::string &strategy,
        double lon_min, double lon_max, double lon_center,
        double lat_min, double lat_max, double lat_center, bool wrap)
    {
        this->patches.emplace_back(
            lon_min, lon_max, lat_min, lat_max,
            InterpolationStrategy::create(strategy),
            (this->scales_.size() == 1) ? static_cast<CoordTransform *>(
                                              new CoordTransformStereographic(lon_center, lat_center, wrap, this->scales_[0]))
                                        : static_cast<CoordTransform *>(
                                              new CoordTransformStereographicLevels(lon_center, lat_center, wrap, this->scales_)));
    }

public:
    /// constructor
    /// n ist length of alts, lons and lats arrays
    Interpolate3D(
        const int n,
        const double *alts,
        const double *lons,
        const double *lats,
        const int m,
        const double *scales,
        const std::string &strategy) : size_(n),
                                       alts_(alts),
                                       lons_(lons),
                                       lats_(lats),
                                       scales_(scales, scales + m),
                                       patches()
    {
        // initializes patches on southern and northern hemisphere.
        // Three stripes on each hemisphere with different longitundinal extent.
        for (int sgn = -1; sgn <= 1; sgn += 2)
        {
            const static double lat_delta = 36;
            const static double lon_delta_1 = 36;
            const static double lon_delta_2 = 45;
            for (double lon_min = -180; lon_min < 179; lon_min += lon_delta_1)
            {
                this->addPatch(
                    strategy, lon_min, lon_min + lon_delta_1, lon_min + (lon_delta_1 / 2.),
                    0, sgn * lat_delta, sgn * lat_delta / 2., true);
            }
            for (double lon_min = -180; lon_min < 179; lon_min += lon_delta_2)
            {
                this->addPatch(
                    strategy, lon_min, lon_min + lon_delta_2, lon_min + (lon_delta_2 / 2.),
                    sgn * lat_delta, sgn * 2 * lat_delta, sgn * 1.5 * lat_delta, true);
            }
            this->addPatch(
                strategy, -180, 180, 0, sgn * 2 * lat_delta, sgn * 90, sgn * 90, false);
        }
    }

    IndiceWeights interpolate(
        double alt, double lon, double lat) const
    {
        lon = fixLongitude(lon);
        lat = fixLatitude(lat);
        for (auto &patch : this->patches)
        {
            if (patch.isInPatch(lon, lat))
            {
                if (!patch.isInitialized())
                {
                    patch.initPatch(this->size_, this->alts_, this->lons_, this->lats_);
                }
                return patch.getPointsAndWeights(Point3(lon, lat, alt));
            }
        }
        jassert(false, "patch not found!?" << lon << " " << lat);
    }
};

/// Finds interval containing the supplied point. The lower index of the
/// interval is returned. In case the element is not contained in the array the
/// closest interval is supplied.
template <typename U, typename T>
inline U locateInterval(const T *const xxs, const U n, size_t step, const double x)
{
    register U ilo = 0;
    register U ihi = n - 1;
    register U i;
    while (std::isnan(xxs[ilo * step]) && ilo < ihi)
    {
        ++ilo;
    }
    while (std::isnan(xxs[ihi * step]) && ilo < ihi)
    {
        --ihi;
    }

    if (xxs[ilo] < xxs[ihi * step])
    {
        while (ihi > ilo + 1)
        {
            i = (ihi + ilo) >> 1;
            if (xxs[i * step] > x)
            {
                ihi = i;
            }
            else
            {
                ilo = i;
            }
        }
    }
    else
    {
        while (ihi > ilo + 1)
        {
            i = (ihi + ilo) >> 1;
            if (xxs[i * step] < x)
            {
                ihi = i;
            }
            else
            {
                ilo = i;
            }
        }
    }
    return ilo;
}

void compute_weight_matrix_c(
    double *coords, size_t n0, size_t n1, size_t n2, size_t n3,
    size_t coord_idx,
    double *level, size_t nlev,
    float *data, int *indices, int *indptr, size_t rows, size_t cols,
    bool allow_extrapolate)
{
    size_t istrides[] = {n1 * n2 * n3, n2 * n3, n3, 1};
    size_t istrides_base[] = {n1 * n2 * n3, n2 * n3, n3, 1};
    istrides_base[coord_idx] = 0;
    size_t ins[] = {n0, n1, n2, n3};
    size_t ons[] = {n0, n1, n2, n3};
    ons[coord_idx] = nlev;
    size_t ostrides[] = {ons[1] * ons[2] * ons[3], ons[2] * ons[3], ons[3], 1};

    size_t icstride = istrides[coord_idx];
    size_t icn = ins[coord_idx];
    size_t entry = 0;
    std::cerr << "a" << allow_extrapolate << "a\n";
    for (size_t i0 = 0; i0 < ons[0]; ++i0)
    {
        for (size_t i1 = 0; i1 < ons[1]; ++i1)
        {
            for (size_t i2 = 0; i2 < ons[2]; ++i2)
            {
                for (size_t i3 = 0; i3 < ons[3]; ++i3)
                {
                    size_t idx[] = {i0, i1, i2, i3};
                    size_t ilev = idx[coord_idx];
                    idx[coord_idx] = 0;
                    size_t row = i0 * ostrides[0] + i1 * ostrides[1] + i2 * ostrides[2] + i3 * ostrides[3];
                    double *coords_base = coords + idx[0] * istrides[0] + idx[1] * istrides[1] + idx[2] * istrides[2] + idx[3] * istrides[3];
                    size_t iclo = locateInterval(coords_base, icn, icstride, level[ilev]);
                    double whi = (coords_base[iclo * icstride] - level[ilev]) /
                                 (coords_base[iclo * icstride] - coords_base[(iclo + 1) * icstride]);
                    indptr[row] = entry;
                    indptr[row + 1] = entry + 2;
                    if ((((whi < 0) || (whi > 1)) && (!allow_extrapolate)) || std::isnan(whi))
                    {
                        data[entry] = std::numeric_limits<double>::quiet_NaN();
                        indices[entry] = 0;
                        entry += 1;
                        data[entry] = std::numeric_limits<double>::quiet_NaN();
                        indices[entry] = 0;
                        entry += 1;
                    }
                    else
                    {
                        double wlo = 1. - whi;
                        idx[coord_idx] = iclo;
                        size_t col = idx[0] * istrides[0] + idx[1] * istrides[1] + idx[2] * istrides[2] + idx[3] * istrides[3];
                        data[entry] = wlo;
                        indices[entry] = col;
                        entry += 1;
                        data[entry] = whi;
                        indices[entry] = col + icstride;
                        entry += 1;
                    }
                }
            }
        }
    }
}
