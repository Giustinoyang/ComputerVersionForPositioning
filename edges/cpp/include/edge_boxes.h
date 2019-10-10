#include <vector>

// trivial array class encapsulating pointer arrays
template <class T> class Array
{
public:
  Array() { _h=_w=0; _x=0; _free=0; }
  virtual ~Array() { clear(); }
  void clear() { if(_free) delete [] _x; _h=_w=0; _x=0; _free=0; }
  void init(int h, int w) { clear(); _h=h; _w=w; _x=new T[h*w](); _free=1; }
  T& val(size_t c, size_t r) { return _x[c*_h+r]; }
  int _h, _w; T *_x; bool _free;
};

// convenient typedefs
typedef std::vector<float> vectorf;
typedef std::vector<int> vectori;
typedef Array<float> arrayf;
typedef Array<int> arrayi;

// bounding box data structures and routines
typedef struct { int c, r, w, h; float s; } Box;
typedef std::vector<Box> Boxes;
bool boxesCompare( const Box &a, const Box &b );
float boxesOverlap( Box &a, Box &b );
void boxesNms( Boxes &boxes, float thr, float eta, int maxBoxes );

////////////////////////////////////////////////////////////////////////////////

// main class for generating edge boxes
class EdgeBoxGenerator
{
public:
  // method parameters (must be manually set)
  float _alpha, _beta, _eta, _minScore; int _maxBoxes;
  float _edgeMinMag, _edgeMergeThr, _clusterMinMag;
  float _maxAspectRatio, _minBoxArea, _gamma, _kappa;

  // main external routine (set parameters first)
  void generate( Boxes &boxes, arrayf &E, arrayf &O, arrayf &V );

private:
  // edge segment information (see clusterEdges)
  int h, w;                         // image dimensions
  int _segCnt;                      // total segment count
  arrayi _segIds;                   // segment ids (-1/0 means no segment)
  vectorf _segMag;                  // segment edge magnitude sums
  vectori _segR, _segC;             // segment lower-right pixel
  std::vector<vectorf> _segAff;          // segment affinities
  std::vector<vectori> _segAffIdx;       // segment neighbors

  // data structures for efficiency (see prepDataStructs)
  arrayf _segIImg, _magIImg; arrayi _hIdxImg, _vIdxImg;
  std::vector<vectori> _hIdxs, _vIdxs; vectorf _scaleNorm;
  float _scStep, _arStep, _rcStepRatio;

  // data structures for efficiency (see scoreBox)
  arrayf _sWts; arrayi _sDone, _sMap, _sIds; int _sId;

  // helper routines
  void clusterEdges( arrayf &E, arrayf &O, arrayf &V );
  void prepDataStructs( arrayf &E );
  void scoreAllBoxes( Boxes &boxes );
  void scoreBox( Box &box );
  void refineBox( Box &box );
  void drawBox( Box &box, arrayf &E, arrayf &V );
};
