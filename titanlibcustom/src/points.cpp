#include "titanlibcustom.h"

using namespace titanlibcustom;

titanlibcustom::Points::Points() {
    // TODO: Deal with the empty case. Can't really find nearest neighbours then

}
titanlibcustom::Points::Points(vec lats, vec lons, vec elevs, vec lafs, CoordinateType type) {
    int N = lats.size();
    if(lons.size() != N)
        throw std::invalid_argument("Cannot create points with unequal lat and lon sizes");
    if(elevs.size() != 0 && elevs.size() != N)
        throw std::invalid_argument("'elevs' must either be size 0 or the same size at lats/lons");
    if(lafs.size() != 0 && lafs.size() != N)
        throw std::invalid_argument("'lafs' must either be size 0 or the same size at lats/lons");
    mLats = lats;
    mLons = lons;
    mElevs = elevs;
    mLafs = lafs;
    KDTree tree = KDTree(lats, lons, type);
    mTree = tree;
    if(mElevs.size() != N) {
        mElevs.clear();
        mElevs.resize(N, titanlibcustom::MV);
    }
    if(mLafs.size() != N) {
        mLafs.clear();
        mLafs.resize(N, titanlibcustom::MV);
    }
}
titanlibcustom::Points::Points(KDTree tree, vec elevs, vec lafs) {
    mElevs = elevs;
    mLafs = lafs;
    mTree = tree;
    mLats = tree.get_lats();
    mLons = tree.get_lons();
}

int titanlibcustom::Points::get_num_neighbours(float lat, float lon, float radius, bool include_match) const {
    return mTree.get_num_neighbours(lat, lon, radius, include_match);
}

ivec titanlibcustom::Points::get_neighbours_with_distance(float lat, float lon, float radius, vec& distances, bool include_match) const {
    return mTree.get_neighbours_with_distance(lat, lon, radius, distances, include_match);
}

ivec titanlibcustom::Points::get_neighbours(float lat, float lon, float radius, bool include_match) const {
    return mTree.get_neighbours(lat, lon, radius, include_match);
}

ivec titanlibcustom::Points::get_closest_neighbours(float lat, float lon, int num, bool include_match) const {
    return mTree.get_closest_neighbours(lat, lon, num, include_match);
}
int titanlibcustom::Points::get_nearest_neighbour(float lat, float lon, bool include_match) const {
    ivec I = get_closest_neighbours(lat, lon, 1, include_match);
    if(I.size() > 0)
        return I[0];
    else
        return -1;
}
vec titanlibcustom::Points::get_lats() const {
    return mLats;
}
vec titanlibcustom::Points::get_lons() const {
    return mLons;
}
vec titanlibcustom::Points::get_elevs() const {
    return mElevs;
}
vec titanlibcustom::Points::get_lafs() const {
    return mLafs;
}
int titanlibcustom::Points::size() const {
    return mLats.size();
}
titanlibcustom::Points& titanlibcustom::Points::operator=(titanlibcustom::Points other) {
    std::swap(mLats, other.mLats);
    std::swap(mLons, other.mLons);
    std::swap(mElevs, other.mElevs);
    std::swap(mLafs, other.mLafs);
    std::swap(mTree, other.mTree);
    return *this;
}
titanlibcustom::Points::Points(const titanlibcustom::Points& other) {
    mLats = other.mLats;
    mLons = other.mLons;
    mElevs = other.mElevs;
    mLafs = other.mLafs;
    mTree = KDTree(mLats, mLons, mTree.get_coordinate_type());
}
CoordinateType titanlibcustom::Points::get_coordinate_type() const {
    return mTree.get_coordinate_type();
}
