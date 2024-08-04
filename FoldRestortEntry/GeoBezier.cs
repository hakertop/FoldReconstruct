using GeologicalEntity;
using OSGeo.OGR;
using System;
using System.Collections.Generic;
using System.Geometry;
using System.Numerics;
using System.Text;
using GdalLib;
using GeoAPI.Geometries;
using NetTopologySuite.Operation.Buffer;
using System.Linq;
using GeoAPI.Operation.Buffer;

namespace FoldRestortEntry
{
    /// <summary>
    /// Bezier object for the stratum
    /// </summary>
    public class GeoBezier
    {

        /// <summary>
        /// Bezier Object
        /// </summary>
        public int _geoId;
        public string _geoBezierName;
        public Bezier _geoBezier;

        /// <summary>
        /// Geological elements
        /// </summary>
        public OccurrencePoint _geoLeftOccurrence;//the attitude at the left side
        public OccurrencePoint _geoRightOccurrence;//the attitude at the right side
        public double _interlimbangle;//interlimb angle
        public double _azimuthAngleofSectionLine;//the azimuth of the cross-section line (degree)

        /// <summary>
        /// 2D geometric objects
        /// </summary>
        public Vector2 _geoLeftPoint;//control point at the left side
        public Vector2 _geoRightPoint;//control point at the right side
        public Vector2 _geointersectPoint;//The intersection point of the two side rays is used as the control point for the Bezier curve
        public Vector2 _zeroTangentPoint;//Point where the tangent direction is horizontal
        public Geometry _2DBezierGeometry;//2D Bezier curve
        public Geometry _ParaellelGeometry;//isopach curve
        //"Shift distance: When converting the 3D non-core section line to 2D, it needs to be shifted a certain distance to the left relative to the core section line
        public float _geoMoveDistance;


        /// <summary>
        /// 3D geometric objects
        /// </summary>
        public Geometry _3DBezierGeometry;//3D Bezier curve
        public Geometry _topPointsIn3DBezierGeometry;//the toppest point at the 3D Bezier curve

        /// <summary>
        /// file save path
        /// </summary>
        public string _geofilesavepath;

        



       /// <summary>
       /// Constructor 1
       /// </summary>
       /// <param name="_beziername">name</param>
       /// <param name="_pLeft"></param>
       /// <param name="_pRight"></param>
       /// <param name="_movedistance"></param>
       /// <param name="_filepath"></param>
        public GeoBezier(string _beziername, OccurrencePoint _pLeft, OccurrencePoint _pRight, float _movedistance,string _filepath)
        {
            _geoBezierName = _beziername;
            _geoLeftOccurrence = _pLeft;
            _geoRightOccurrence = _pRight;
            _geoMoveDistance = _movedistance;
            _geofilesavepath = _filepath;
        }


        /// <summary>
        /// Constructor 2
        /// </summary>
        /// <param name="_beziername"></param>
        /// <param name="_pLeftOcc"></param>
        /// <param name="_pRightOcc"></param>
        /// <param name="_pLeft"></param>
        /// <param name="_pRight"></param>
        /// <param name="_pMiddlePoint"></param>
        /// <param name="_movedistance"></param>
        /// <param name="_filepath"></param>
        public GeoBezier(string _beziername, OccurrencePoint _pLeftOcc, OccurrencePoint _pRightOcc,Vector2 _pLeft, Vector2 _pRight, Vector2 _pMiddlePoint, float _movedistance, string _filepath)
        {
            
            _geoBezierName = _beziername;
            _geoLeftPoint = _pLeft;
            _geoRightPoint = _pRight;
            _geoLeftOccurrence = _pLeftOcc;
            _geoRightOccurrence = _pRightOcc;
            _geoMoveDistance = _movedistance;
            _geofilesavepath = _filepath; 
            _geointersectPoint = _pMiddlePoint;
        }

        #region Function library

        /// <summary>
        /// Construct a Bezier curve based on three control points
        /// </summary>
        /// <param name="_intervals">Specify the interval size between points when outputting the geometric curve</param>
        public void DrawSectionWithThreeControlPoints(float _intervals)
        {          
            _geoBezier = new Bezier(_geoLeftPoint, _geointersectPoint, this._geoRightPoint);          
            _azimuthAngleofSectionLine = CalAzimuthbytwoPoints(_geoLeftOccurrence.GetX(0), _geoLeftOccurrence.GetY(0), _geoRightOccurrence.GetX(0), _geoRightOccurrence.GetY(0), false);          
            CreateMiddleGeometry(_intervals);
        }

        /// <summary>
        ///  Construct a Bezier curve based on two attitude points
        /// </summary>
        /// <param name="_intervals">Specify the interval distance between points when outputting the geometric curve</param>
        public void DrawSectionByBezierFunction(float _intervals)
        {
            //1、Calculate the azimuth (in degrees) of the line connecting two attitude points
            _azimuthAngleofSectionLine = CalAzimuthbytwoPoints(_geoLeftOccurrence.GetX(0), _geoLeftOccurrence.GetY(0), _geoRightOccurrence.GetX(0), _geoRightOccurrence.GetY(0),false);

            //2、Calculate the apparent dip for each attitude point
            double viewDipLeft = Math.Round(CalculateViewDipAngle(_geoLeftOccurrence.dipAngle, _geoLeftOccurrence.tendency, _azimuthAngleofSectionLine, 0), 4);
            double viewDipRight = Math.Round(CalculateViewDipAngle(_geoRightOccurrence.dipAngle, _geoRightOccurrence.tendency, _azimuthAngleofSectionLine, 0), 4);

            //3 Calculate the distance between two attitude points
            double dis_Left_Right = Math.Sqrt(Math.Pow(_geoLeftOccurrence.GetX(0) - _geoRightOccurrence.GetX(0), 2) + Math.Pow(_geoRightOccurrence.GetY(0) - _geoLeftOccurrence.GetY(0), 2));

            //4、Project the two attitude points onto a 2D plane
            this._geoLeftPoint = new Vector2(0.0f - _geoMoveDistance, (float)_geoLeftOccurrence.GetZ(0));
            this._geoRightPoint = new Vector2((float)dis_Left_Right - _geoMoveDistance, (float)_geoRightOccurrence.GetZ(0));

            //5 Find the intersection point of two rays.
            //5.1 Condition for intersection: The two apparent dips can only intersect if their signs are opposite
            if (viewDipLeft * viewDipRight <= 0)
            {
                //5.2 calculate the intersection point
                _geointersectPoint = CalculateIntersection(this._geoLeftPoint, viewDipLeft, this._geoRightPoint, viewDipRight);

                // 5.2.1 calculate the interlimb angle
                CaculateInterlimbAngle();
                Console.WriteLine("interlimb angle：" + this._interlimbangle);

                //5.3 Construct a quadratic Bézier curve
                _geoBezier = new Bezier(_geoLeftPoint, _geointersectPoint, _geoRightPoint);
                if(_geoBezier==null)
                    Console.WriteLine("Failed to generate the Bézier curve" + _geoBezierName);
                
                //5.4 output
                CreateMiddleGeometry(_intervals);
            }
            else
            {
                Console.WriteLine("An error occurred; please check the attitudes on both sides!");
                Console.WriteLine(this._geoBezierName);
                
                Console.WriteLine("Left-side attitude:"+ _geoLeftOccurrence.dipAngle+" " + _geoLeftOccurrence.strike + " " +_geoLeftOccurrence.tendency);
                Console.WriteLine("Right-side attitude:"+ _geoRightOccurrence.dipAngle+" " + _geoRightOccurrence.strike + " " + _geoRightOccurrence.tendency);
            }
        }

        ///<summary>
        /// Convert true dip to apparent dip
        /// </summary>
        /// <param name="dip">true dip (degree)</param>
        /// <param name="dipDirection"></param>
        /// <param name="azimuth"></param>
        /// <param name="viewAngle"></param>
        /// <returns>apparent dip</returns>
        public double CalculateViewDipAngle(double dip, double dipDirection, double azimuth, double viewAngle)
        {
            
            double dipRad = dip * Math.PI / 180.0;
            
            double dipDirectionRad = dipDirection * Math.PI / 180.0;
           
            double azimuthRad = azimuth * Math.PI / 180.0;
            
            double viewAngleRad = viewAngle * Math.PI / 180.0;

           
            double cosTerm = Math.Tan(dipRad - viewAngleRad) * Math.Cos(dipDirectionRad - azimuthRad);
            double viewDipAngleRad = Math.Atan(cosTerm);
            double viewDipAngle = viewDipAngleRad * 180.0 / Math.PI;

            return viewDipAngle;
        }


        /// <summary>
        /// Calculate the azimuth of a line using the coordinates of two endpoints, with the result in degrees or radians
        /// </summary>
        /// <param name="_x1"></param>
        /// <param name="_y1"></param>
        /// <param name="_x2"></param>
        /// <param name="_y2"></param>
        /// <returns></returns>
        public static double CalAzimuthbytwoPoints(double _x1, double _y1, double _x2, double _y2,bool _isArc)
        {
            double deltaX = _x2 - _x1;
            double deltaY = _y2 - _y1;
            double angle = Math.Atan2(deltaX, deltaY);

            if (angle < 0)
            {
                angle += (2 * Math.PI);
            }

            //Output radians or degrees
            if (_isArc)
                return angle;
            else
                return angle * 180 / Math.PI;
        }

        /// <summary>
        /// Calculate the intersection point of the rays corresponding to the attitude points.
        /// </summary>
        /// <param name="p1">Left-side attitude point</param>
        /// <param name="bearing1">the apparent angle of the left-side attitude point</param>
        /// <param name="p2">Reft-side attitude point</param>
        /// <param name="bearing2">the apparent angle of the reft-side attitude point</param>
        /// <returns></returns>
        public Vector2 CalculateIntersection(Vector2 p1, double bearing1, Vector2 p2, double bearing2)
        {
            // Convert bearings to radians
            double theta1 = -(bearing1 * Math.PI / 180);
            double theta2 = -(bearing2 * Math.PI / 180);

            // Convert points to Cartesian coordinates
            float x1 = p1.X;
            float y1 = p1.Y;
            float x2 = p2.X;
            float y2 = p2.Y;

            // Calculate slopes of lines
            double m1 = Math.Tan(theta1);
            double m2 = Math.Tan(theta2);

            // Calculate y-intercepts of lines
            double b1 = y1 - m1 * x1;
            double b2 = y2 - m2 * x2;

            // Calculate intersection point
            double x = (b2 - b1) / (m1 - m2);
            double y = m1 * x + b1;

            return new Vector2((float)x, (float)y);
        }


        /// <summary>
        /// Output the Bezier curves (2D and 3D)
        /// </summary>
        /// <param name="_bezier"></param>
        public void CreateMiddleGeometry(float _stepSize)
        {
            
            _2DBezierGeometry = new Geometry(wkbGeometryType.wkbLineString);
            _3DBezierGeometry = new Geometry(wkbGeometryType.wkbLineString);

            
            int count = 0;
            for (float t = 0; t <= 1; t += _stepSize)
            {
                Vector2 point = _geoBezier.Position(t);
                _2DBezierGeometry.AddPoint_2D(point.X, point.Y);
                count++;
            }
           
            for (int t = 0; t < count; t++)
            {
                if (t == 0)
                    _3DBezierGeometry.AddPoint(_geoLeftOccurrence.GetX(0), _geoLeftOccurrence.GetY(0), _geoLeftOccurrence.GetZ(0));
                if (t == count)
                    _3DBezierGeometry.AddPoint(_geoRightOccurrence.GetX(0), _geoRightOccurrence.GetY(0), _geoRightOccurrence.GetZ(0));

                _3DBezierGeometry.AddPoint(_geoLeftOccurrence.GetX(0) + (_geoRightOccurrence.GetX(0) - _geoLeftOccurrence.GetX(0)) * t / (count - 1), _geoLeftOccurrence.GetY(0) +
                    (_geoRightOccurrence.GetY(0) - _geoLeftOccurrence.GetY(0)) * t / (count - 1), _2DBezierGeometry.GetY(t));
            }
        }



        /// <summary>
        /// Find the point where the tangent is horizontal (in 2D)
        /// </summary>
        /// <param name="_stepSize">step length</param>
        /// <param name="_thresold">error thresold</param>
        public Geometry CalucateZeroTangent(float _stepSize, float _thresold)
        {
            Geometry pt = new Geometry(wkbGeometryType.wkbPoint);

            float tnt = 0.0f;
            while(pt.GetX(0)==0.0)
            {
                for (float t = 0.0f; t <= 1.0f; t += _stepSize)
                {
                    Vector2 pTangent = _geoBezier.Tangent(t);

                    float pSlope = pTangent.Y / pTangent.X;

                    if (Math.Abs(pSlope) < (_thresold+ tnt))
                    {
                        _zeroTangentPoint = _geoBezier.Position(t);

                        pt.AddPoint_2D(_zeroTangentPoint.X, _zeroTangentPoint.Y);
                    }
                }
                tnt = tnt + 0.01f;

            }
            return pt;
        }


        /// <summary>
        /// Calculate interlimb angle
        /// </summary>
        /// <returns></returns>
        public double CaculateInterlimbAngle()
        {
            Vector2 pAB = Vector2.Subtract (_geoLeftPoint,_geointersectPoint);
            Vector2 pAC = Vector2.Subtract (_geoRightPoint,_geointersectPoint);

            float cosTheta = (Vector2.Dot(pAB, pAC) / (pAB.Length() * pAC.Length()));
            _interlimbangle = Math.Acos(cosTheta) * 180 / Math.PI;

            return _interlimbangle;
        }

        /// <summary>
        /// Get the toppest point on the 3D stratigraphic paleo-boundary
        /// </summary>
        /// <returns></returns>
        public Geometry Get3DTopestPointOnPolyline()
        {
            if(this._3DBezierGeometry==null)
            {
                Console.WriteLine("3D stratigraphic paleo=boundary has not been generated！！！");
            }

            Dictionary<Geometry,double> pPointsOnBezier = new Dictionary<Geometry, double>(this._3DBezierGeometry.GetPointCount());

            for(int i=0;i< this._3DBezierGeometry.GetPointCount();i++)
            {
                Geometry pt = new Geometry(wkbGeometryType.wkbPoint);
                pt.AddPoint(this._3DBezierGeometry.GetX(i), this._3DBezierGeometry.GetY(i), this._3DBezierGeometry.GetZ(i));

                pPointsOnBezier.Add(pt,pt.GetZ(0));
            }
            List<double> pDistances = pPointsOnBezier.Values.ToList();
            pDistances.Sort((a, b) => b.CompareTo(a));


            this._topPointsIn3DBezierGeometry = new Geometry(wkbGeometryType.wkbPoint);
            
            foreach (var vp in pPointsOnBezier.Keys)
            {
                if(pPointsOnBezier[vp]== pDistances[0])
                    this._topPointsIn3DBezierGeometry = vp;
            }

            return this._topPointsIn3DBezierGeometry;
        }


        /// <summary>
        /// Generate a curve to construct an isopach map
        /// </summary>
        /// <param name="_pLine"></param>
        /// <param name="_distance"></param>
        /// <returns></returns>
        public void ConstructParallelGeometry(double _distance)
        {
            //1 Create a buffer for curve recording
            _ParaellelGeometry = new Geometry(wkbGeometryType.wkbLineString);

            //2 Convert lines from GDAL to Geometries in NetTopologySuite
            Coordinate[] coords = new Coordinate[_2DBezierGeometry.GetPointCount()];
            for (int i = 0; i < _2DBezierGeometry.GetPointCount(); i++)
            {
                var coord = new Coordinate(_2DBezierGeometry.GetX(i), _2DBezierGeometry.GetY(i));
                coords[i] = coord;
            }
            var line = new NetTopologySuite.Geometries.LineString(coords);

            //3 buffer type
            var bufferParams = new BufferParameters
            {
                EndCapStyle = EndCapStyle.Flat,
                IsSingleSided = true
            };

            //4 Specify the buffer distance and generate a single-sided buffer
            var polygon = line.Buffer(_distance, bufferParams) as NetTopologySuite.Geometries.Polygon;

            //5 Convert a polygon back to a line
            var buffer = polygon.ExteriorRing;

            //6 Output the coordinates of the single-sided buffer; if only the outer line is needed, filtering is required
            foreach (var coordinate in buffer.Coordinates)
            {
                if (coords.Contains(coordinate))
                    continue;

                _ParaellelGeometry.AddPoint_2D(coordinate.X, coordinate.Y);
            }
        }

        #endregion

    }
}
