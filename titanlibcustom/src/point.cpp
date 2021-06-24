#include "titanlibcustom.h"

using namespace titanlibcustom;

titanlibcustom::Point::Point(float lat, float lon, float elev, float laf, CoordinateType type) {
    this->lat = lat;
    this->lon = lon;
    this->elev = elev;
    this->laf = laf;
    this->type = type;
}
