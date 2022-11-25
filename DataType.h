#ifndef DATATYPE_H_INCLUDED
#define DATATYPE_H_INCLUDED


typedef double  (*MapN2R)(unsigned);    // Ptr2Fun like f: N --> R
typedef double  (*MapR2R)(double, ...); // Ptr2Fun like f: R x ... x R --> R

typedef struct  {double x; double y;} Point2d;
typedef Point2d (*MapP2P)(Point2d p, ...); // Ptr2Fun like f: Point2d --> Poiint2d

#endif // DATATYPE_H_INCLUDED
