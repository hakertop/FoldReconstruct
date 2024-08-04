using GdalLib;
using GeoCommon;
using OSGeo.OGR;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;

namespace GdalLib
{
   public static class GeometrySolve
    {

      
        public static List<Geometry> GetAllPointsFromPolygon(Geometry _polygon, int _dimension)
        {
            int geocount = _polygon.GetGeometryCount();

            List<Geometry> pPointColl = new List<Geometry>(geocount);

            for (int i = 0; i < geocount; i++)
            {
                Geometry geo_i = _polygon.GetGeometryRef(i);

                for (int j = 0; j < geo_i.GetPointCount(); j++)
                {
                    Geometry pt = new Geometry(wkbGeometryType.wkbPoint);

                    if (_dimension == 2)
                    {
                        pt.AddPoint_2D(geo_i.GetX(j), geo_i.GetY(j));
                    }


                    if (_dimension == 3)
                    {
                        pt.AddPoint(geo_i.GetX(j), geo_i.GetY(j), geo_i.GetZ(j));
                    }

                    
                    bool isExist = false;
                    foreach(var vp in pPointColl)
                    {
                        if (vp.GetX(0) == pt.GetX(0) && vp.GetY(0) == pt.GetY(0))
                            isExist = true;
                    }
                    if (!isExist)
                        pPointColl.Add(pt);
                }
            }

            return pPointColl;
        }

       
        public static List<double> GetInterpolatePointSclale(Geometry pCurve)
        {

            int pointcount = pCurve.GetPointCount();

            List<double> scales = new List<double>(pointcount - 2);

           
            Geometry pReferLine = new Geometry(wkbGeometryType.wkbLineString);

            
            Geometry pt2 = new Geometry(wkbGeometryType.wkbPoint);
            pt2.AddPoint_2D(pCurve.GetX(0), pCurve.GetY(0));

            
            Geometry pt3 = new Geometry(wkbGeometryType.wkbPoint);
            pt3.AddPoint_2D(pCurve.GetX(pointcount - 1), pCurve.GetY(pointcount - 1));

            pReferLine.AddPoint_2D(pCurve.GetX(0), pCurve.GetY(0));
            for (int i = 1; i < pointcount - 1; i++)
            {
                Geometry pt1 = new Geometry(wkbGeometryType.wkbPoint);
                pt1.AddPoint_2D(pCurve.GetX(i), pCurve.GetY(i));
                
                Geometry pft = GeometrySolve.GetFootPoint(pt1, pt2, pt3);
                pReferLine.AddPoint_2D(pft.GetX(0), pft.GetY(0));
            }
            pReferLine.AddPoint_2D(pCurve.GetX(pointcount - 1), pCurve.GetY(pointcount - 1));

            double _cumulation = 0;

            for (int i = 1; i < pReferLine.GetPointCount() - 1; i++)
            {
                Geometry pSegment = new Geometry(wkbGeometryType.wkbLineString);
                pSegment.AddPoint_2D(pReferLine.GetX(i), pReferLine.GetY(i));
                pSegment.AddPoint_2D(pReferLine.GetX(i - 1), pReferLine.GetY(i - 1));

                _cumulation = _cumulation + pSegment.Length();
                
                double sca = _cumulation / pReferLine.Length();
                scales.Add(sca);
            }

            return scales;

            //pReferLine.ExportSimpleGeometryToShapfile(_workspacePath, "changelinefoot");
        }


      
        public static float GetAngle(Geometry cen, Geometry first, Geometry second)
        {
            double dx1, dx2, dy1, dy2;
            float angle;

            dx1 = first.GetX(0) - cen.GetX(0);
            dy1 = first.GetY(0) - cen.GetY(0);

            dx2 = second.GetX(0) - cen.GetX(0);
            dy2 = second.GetY(0) - cen.GetY(0);

            float c = (float)Math.Sqrt(dx1 * dx1 + dy1 * dy1) * (float)Math.Sqrt(dx2 * dx2 + dy2 * dy2);

            if (c == 0) return -1;

           
            double bs = Math.Acos((dx1 * dx2 + dy1 * dy2) / c);

            if (!Double.IsNaN(bs))
                angle = (float)(Math.Acos((dx1 * dx2 + dy1 * dy2) / c) * 180 / Math.PI);
            else
                angle = 180f;
          
            return angle;
           
        }


       
        public static Geometry GetFootPoint(Geometry pG1,Geometry pG2,Geometry pG3)
        {
            Geometry _ptFootPoint = new Geometry(wkbGeometryType.wkbPoint);

            Point p1 = new Point(pG1.GetX(0), pG1.GetY(0));
            Point p2 = new Point(pG2.GetX(0), pG2.GetY(0));
            Point p3 = new Point(pG3.GetX(0), pG3.GetY(0));

            Vector p2p3 = p3 - p2;
            Vector p2p1 = p2 - p1;
            Vector BD = p2p3 * p2p1 / ((p2.X - p3.X) * (p2.X - p3.X) + (p2.Y - p3.Y) * (p2.Y - p3.Y)) * p2p3;//BD:P2与目标的构成的向量
            Point p4 = p2 - BD;

            _ptFootPoint.AddPoint_2D(p4.X,p4.Y);

            return _ptFootPoint;

        }

        
        public static bool IsSameOrder(Geometry pLine1,Geometry pLine2)
        {
            if (!pLine1.Intersect(pLine2))
            {
                Geometry _pLine1 = new Geometry(wkbGeometryType.wkbLineString);
                _pLine1.AddPoint_2D(pLine1.GetX(0), pLine1.GetY(0));
                _pLine1.AddPoint_2D(pLine2.GetX(pLine2.GetPointCount()-1), pLine2.GetY(pLine2.GetPointCount() - 1));

                Geometry _pLine2 = new Geometry(wkbGeometryType.wkbLineString);
                _pLine2.AddPoint_2D(pLine2.GetX(0), pLine2.GetY(0));
                _pLine2.AddPoint_2D(pLine1.GetX(pLine1.GetPointCount() - 1), pLine1.GetY(pLine1.GetPointCount() - 1));

                if (_pLine1.Intersect(_pLine2))
                    return true;
                else
                    return false;
            }

            return false;
        }

       
        public static double GetEvelationH(double x1,double y1,double x2,double y2,DEMRaster _raster)
        {
            double z1 = _raster.GetElevation(x1, y1);
            double z2 = _raster.GetElevation(x2, y2);
            double detaZ = z1 - z2;
            return detaZ;
        }


      
        public static double GetDistance2D(double x1, double y1, double x2, double y2)
        {
            return Math.Sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
        }

      
        public static Geometry SwapLineOrder(Geometry _PLine)
        {
            Geometry pSwapLine = new Geometry(wkbGeometryType.wkbLineString);
            for(int i= _PLine.GetPointCount()-1;i>=0;i--)
            {
                pSwapLine.AddPoint_2D(_PLine.GetX(i),_PLine.GetY(i));
            }

            return pSwapLine;
        }


      
        public static Geometry Get3DIntersectPoint(Geometry _pLine1, Geometry _pLine2, DEMRaster _raster)
        {

            Geometry pt = new Geometry(wkbGeometryType.wkbPoint);

            if (!_pLine1.Intersect(_pLine2))
                return null;

            Geometry _intersectPoint = _pLine1.Intersection(_pLine2);

            if(_intersectPoint.GetGeometryType()==wkbGeometryType.wkbPoint)
            {
                pt.AddPoint(_intersectPoint.GetX(0), _intersectPoint.GetY(0), _raster.GetElevation(_intersectPoint.GetX(0), _intersectPoint.GetY(0)));

            }

            return pt;
        }

       


       
        public static Geometry InsertPointsToLine(Geometry _pLine,double _interval)
        {
            Geometry _pNewGeometry = new Geometry(wkbGeometryType.wkbLineString);

            int pointcount = _pLine.GetPointCount();

            for (int i = 0; i < _pLine.GetPointCount()-1; i++)
            {
                _pNewGeometry.AddPoint_2D(_pLine.GetX(i), _pLine.GetY(i));

                
                double dist = Math.Sqrt(Math.Pow(_pLine.GetX(i)- _pLine.GetX(i+1),2)+ Math.Pow(_pLine.GetY(i) - _pLine.GetY(i + 1), 2));
                
                int pCount =Convert.ToInt16(dist / _interval);
                
                for(int j=0;j< pCount;j++)
                {
                    double detaX = _pLine.GetX(i + 1) - _pLine.GetX(i);
                    double detaY = _pLine.GetY(i + 1) - _pLine.GetY(i);
                    double _tx = _pLine.GetX(i) + (j + 1) * detaX / (pCount + 1);
                    double _ty = _pLine.GetY(i) + (j + 1) * detaY / (pCount + 1) ;

                    _pNewGeometry.AddPoint_2D(_tx, _ty);
                }

            }
            _pNewGeometry.AddPoint_2D(_pLine.GetX(pointcount-1), _pLine.GetY(pointcount-1));

            return _pNewGeometry;
        }

       
        public static Geometry InsertPointsToLine(Geometry _pLine, int _pCount)
        {
            Geometry _pNewGeometry = new Geometry(wkbGeometryType.wkbLineString);

            int pointcount = _pLine.GetPointCount();

            for (int i = 0; i < _pLine.GetPointCount() - 1; i++)
            {
                _pNewGeometry.AddPoint_2D(_pLine.GetX(i), _pLine.GetY(i));

               
                double dist = Math.Sqrt(Math.Pow(_pLine.GetX(i) - _pLine.GetX(i + 1), 2) + Math.Pow(_pLine.GetY(i) - _pLine.GetY(i + 1), 2));
                
                for (int j = 0; j < _pCount; j++)
                {
                    double detaX = _pLine.GetX(i + 1) - _pLine.GetX(i);
                    double detaY = _pLine.GetY(i + 1) - _pLine.GetY(i);
                    double _tx = _pLine.GetX(i) + (j + 1) * detaX / (_pCount + 1);
                    double _ty = _pLine.GetY(i) + (j + 1) * detaY / (_pCount + 1);
                    _pNewGeometry.AddPoint_2D(_tx, _ty);
                }

            }
            _pNewGeometry.AddPoint_2D(_pLine.GetX(pointcount - 1), _pLine.GetY(pointcount - 1));

            return _pNewGeometry;
        }

       
        public static  Geometry Create3DLine( Geometry _pLine, DEMRaster _raster)
        {
             
                List<Geometry> pLinePoints = new List<Geometry>();

                Geometry _p3DLine = new Geometry(wkbGeometryType.wkbLineString);

                for (int i = 0; i < _pLine.GetPointCount(); i++)
                {
                    Geometry pt = new Geometry(wkbGeometryType.wkbPoint);
                    pt.AddPoint(_pLine.GetX(i), _pLine.GetY(i), _raster.GetElevation(_pLine.GetX(i), _pLine.GetY(i)));

                    
                    bool isExist = IsExistRepeatPoint(pLinePoints, pt);
                    if (!isExist)
                    {
                        pLinePoints.Add(pt);
                    }
                    pLinePoints.TrimExcess();
                }
                for (int i = 0; i < pLinePoints.Count; i++)
                {
                    _p3DLine.AddPoint(pLinePoints[i].GetX(0), pLinePoints[i].GetY(0), pLinePoints[i].GetZ(0));
                }

            return _p3DLine;
        }

        
        public static bool IsExistRepeatPoint(List<Geometry> _pBoundaryPoints, Geometry _pCheckPoint)
        {
            bool isExist = false;
            foreach (var vt in _pBoundaryPoints)
            {
                if (Math.Abs(_pCheckPoint.GetX(0) - vt.GetX(0)) < 0.01 && Math.Abs(_pCheckPoint.GetY(0) - vt.GetY(0)) < 0.01)
                {
                    isExist = true;
                    break;
                }
            }

            return isExist;
        }

       
        public static Geometry ConverteMultiLineToSingleLine(Geometry pMultiLines)
        {
            Geometry pSingleLine = new Geometry(wkbGeometryType.wkbLineString);
            Eage3D pEage = new Eage3D();
            for(int i=0;i<pMultiLines.GetGeometryCount();i++)
            {
                for(int j=0;j< pMultiLines.GetGeometryRef(i).GetPointCount();j++)
                {
                    pEage.AddVertex(pMultiLines.GetGeometryRef(i).GetX(j), pMultiLines.GetGeometryRef(i).GetY(j));
                }
            }

            foreach(var vt in pEage.vertexList)
            {
                pSingleLine.AddPoint_2D(vt.x,vt.y);
            }

            pSingleLine.AddPoint_2D(pEage.vertexList[0].x, pEage.vertexList[0].y);

            return pSingleLine;
        }


    }
}
