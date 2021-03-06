#include "titanlibcustom.h"

using namespace titanlibcustom;

ivec titanlibcustom::metadata_check(const Points& points, bool check_lat, bool check_lon, bool check_elev, bool check_laf) {
    const int s = points.size();
    const vec& lats = points.get_lats();
    const vec& lons = points.get_lons();
    const vec& elevs = points.get_elevs();
    const vec& lafs = points.get_lafs();

    if( lats.size() != s || lons.size() != s || elevs.size() != s || lafs.size() != s ) {
        throw std::runtime_error("Dimension mismatch");
    }

    ivec flags(s, 0);
    for(int i = 0; i < s; i++) {
        if(check_lat && !titanlibcustom::is_valid(lats[i]))
            flags[i] = 1;
        if(check_lon && !titanlibcustom::is_valid(lons[i]))
            flags[i] = 1;
        if(check_elev && !titanlibcustom::is_valid(elevs[i]))
            flags[i] = 1;
        if(check_laf && !titanlibcustom::is_valid(lafs[i]))
            flags[i] = 1;
    }
    return flags;
}
