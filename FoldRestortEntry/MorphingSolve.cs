using GeoCommon;
using OSGeo.OGR;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using GdalLib;

using GeoAPI.Geometries;

namespace FoldRestortEntry
{
    public  class MorphingSolve
    {



        /// <summary>
        /// Generate transition Paleo-boundaries by Morphing technique
        /// </summary>
        /// <param name="SRC">The starting boundary</param>
        /// <param name="DEST">The target boundary</param>
        /// <param name="leftControl">The first constraint boundary</param>
        /// <param name="middleControl">intermediate constraint boundaries</param>
        /// <param name="rightControl">The last constraint boundary</param>
        /// <param name="_raster">DEM</param>
        /// <returns>transition Paleo-boundaries</returns>
        public static List<Geometry> Create3DTransitionPaleoBoundariesByMorphing(Geometry SRC, Geometry DEST, Geometry leftControl, Geometry middleControl, Geometry rightControl,DEMRaster _raster)
        {
            //1 Convert all points within the OGR Geometry to GeoAPI Coordinates for easier subsequent usage
            List<Coordinate> pCoorSRC = new List<Coordinate>(SRC.GetPointCount());
            List<Coordinate> pCoorDEST = new List<Coordinate>(DEST.GetPointCount());

            //2 In theory, the number of points in SRC should be consistent with the number of points in DEST
            if (SRC.GetPointCount() == DEST.GetPointCount())
            {
                for (int i = 0; i < SRC.GetPointCount(); i++)
                {
                    pCoorSRC.Add(new Coordinate(SRC.GetX(i), SRC.GetY(i), SRC.GetZ(i)));
                    pCoorDEST.Add(new Coordinate(DEST.GetX(i), DEST.GetY(i), DEST.GetZ(i)));

                }
            }
            else
                Console.WriteLine("The number of points on the starting boundary and the end boundary is inconsistent, please check!!！");

            //3 Convert all points within the OGR Geometry to GeoAPI Coordinates for easier subsequent usage
            List<Coordinate> pCoorleftControl = new List<Coordinate>(leftControl.GetPointCount());
            List<Coordinate> pCoorMiddlecontrol = new List<Coordinate>(middleControl.GetPointCount());
            List<Coordinate> pCoorRightControl = new List<Coordinate>(rightControl.GetPointCount());

            //4 Theoretically, the number of points on the control boundaries should be consistent.
            if ((leftControl.GetPointCount() == middleControl.GetPointCount()) && (middleControl.GetPointCount() == rightControl.GetPointCount()) && (leftControl.GetPointCount() == rightControl.GetPointCount()))
            {
                for (int i = 0; i < leftControl.GetPointCount(); i++)
                {
                    pCoorleftControl.Add(new Coordinate(leftControl.GetX(i), leftControl.GetY(i), leftControl.GetZ(i)));
                    pCoorMiddlecontrol.Add(new Coordinate(middleControl.GetX(i), middleControl.GetY(i), middleControl.GetZ(i)));
                    pCoorRightControl.Add(new Coordinate(rightControl.GetX(i), rightControl.GetY(i), rightControl.GetZ(i)));

                }
            }
            else
                Console.WriteLine("The number of geometric points on the control boundaries is inconsistent, please check!!");

            //5 Calculate the interpolated points based on the corresponding points on SRC and DEST.
            List<List<Coordinate>> pMultiInterpolate = new List<List<Coordinate>>(pCoorSRC.Count-2);
            for (int i = 1; i < pCoorSRC.Count-1; i++)
            {
                double detaX = pCoorDEST[i].X - pCoorSRC[i].X;
                double detaY = pCoorDEST[i].Y - pCoorSRC[i].Y;
                double detaZ = pCoorDEST[i].Z - pCoorSRC[i].Z;

                List<Coordinate> pInterpolatePoints = new List<Coordinate>(pCoorleftControl.Count);

                List<Geometry> plins = new List<Geometry>();

                double detaZleft = pCoorleftControl[pCoorleftControl.Count - 1].Z - pCoorleftControl[0].Z;
                double detaZMiddle = pCoorMiddlecontrol[pCoorMiddlecontrol.Count - 1].Z - pCoorMiddlecontrol[0].Z;
                double detaZRight = pCoorRightControl[pCoorRightControl.Count - 1].Z - pCoorRightControl[0].Z;

                for (int j = 1; j < pCoorleftControl.Count - 1; j++)
                {
                    double xl = pCoorSRC[i].X + ((j * detaX )/ (pCoorleftControl.Count - 1));
                    double yl = pCoorSRC[i].Y + (j* detaY) / (pCoorleftControl.Count - 1);
                    double zl = pCoorSRC[i].Z + (detaZ*j) / (pCoorleftControl.Count - 1);

                    double leftzl= pCoorleftControl[0].Z + (detaZleft * j) / (pCoorleftControl.Count - 1);
                    double middelzl= pCoorMiddlecontrol[0].Z + (detaZMiddle * j) / (pCoorleftControl.Count - 1);
                    double rightzl= pCoorRightControl[0].Z + (detaZRight * j) / (pCoorleftControl.Count - 1);

                    Coordinate p3DInterpolatePoint = new Coordinate(xl,yl,zl);

                    double detaZToLeftControl = pCoorleftControl[j].Z- leftzl;
                    double detaZToMiddleControl = pCoorMiddlecontrol[j].Z- middelzl;
                    double detaZToRightControl = pCoorRightControl[j].Z- rightzl;

                    double distanceToLeftControl = p3DInterpolatePoint.Distance(pCoorleftControl[j]);
                    double distanceToMiddleControl = p3DInterpolatePoint.Distance(pCoorMiddlecontrol[j]);
                    double distanceToRightControl = p3DInterpolatePoint.Distance(pCoorRightControl[j]);


                    if (distanceToLeftControl == 0.0)
                    {
                        p3DInterpolatePoint.Z = pCoorleftControl[j].Z;
                    }
                    else if (distanceToMiddleControl == 0.0)
                    {
                        p3DInterpolatePoint.Z = pCoorMiddlecontrol[j].Z;
                    }
                    else if (distanceToRightControl == 0.0)
                    {
                        p3DInterpolatePoint.Z = pCoorRightControl[j].Z;
                    }
                    //Morphing technique
                    else
                    {
                        double z1w = 1 / distanceToLeftControl;

                        double z2w = 1 / distanceToMiddleControl;
                        double z3w = 1 / distanceToRightControl;
                   
                        double ztolal = z1w + z2w + z3w;

                        double zleft = z1w / ztolal;
                        detaZToLeftControl = zleft * detaZToLeftControl;

                        double zMiddle = z2w / ztolal;
                        detaZToMiddleControl = zMiddle * detaZToMiddleControl;

                        double zRight = z3w / ztolal;
                        detaZToRightControl = zRight * detaZToRightControl;

                        //Compare with the elevation on the original DEM
                        double elevationOfDEM =  _raster.GetElevation(p3DInterpolatePoint.X, p3DInterpolatePoint.Y);

                        if ((zl + (detaZToLeftControl + detaZToMiddleControl + detaZToRightControl)) < elevationOfDEM)
                            p3DInterpolatePoint.Z = elevationOfDEM;
                        else
                            p3DInterpolatePoint.Z = zl + (detaZToLeftControl + detaZToMiddleControl + detaZToRightControl);
                    }
                    pInterpolatePoints.Add(p3DInterpolatePoint);

                }

                pMultiInterpolate.Add(pInterpolatePoints);
            }


            //6 Generate transsition paleo-boundaries
            List<Geometry> pTranLines = new List<Geometry>(pCoorleftControl.Count - 2);
            for (int i = 0; i < pCoorleftControl.Count - 2; i++)
            {
                List<Coordinate> pTranLinePoints = new List<Coordinate>();
                foreach (var vp in pMultiInterpolate)
                {
                    pTranLinePoints.Add(vp[i]);
                }

                pTranLines.Add(ConvertToGeometryLine(pTranLinePoints));

            }


            return pTranLines;
        }

        /// <summary>
        /// Convert the continuous Coordinates into a 3D Geometry line
        /// </summary>
        /// <param name="_pLineCoordintes"></param>
        /// <param name="_raster"></param>
        /// <returns></returns>
        public static Geometry ConvertToGeometryLine(List<Coordinate> _pLineCoordintes)
        {
            Geometry pLine = new Geometry(wkbGeometryType.wkbLineString);
            for (int i = 0; i < _pLineCoordintes.Count; i++)
            {
                pLine.AddPoint(_pLineCoordintes[i].X, _pLineCoordintes[i].Y, _pLineCoordintes[i].Z);
            }

            return pLine;
        }
    }
}
