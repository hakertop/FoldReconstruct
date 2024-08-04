﻿using System;
using System.Collections.Generic;
using System.Text;
using GeoAPI.Geometries;
using GeoCommon;
using GeologicalEntity;
using OSGeo.OGR;

namespace FoldRestortEntry
{
    public class FoldStratum:Stratum
    {
        /// <summary>
        /// Bezier curves
        /// </summary>
        public List<GeoBezier> beziers { get; set; }

        /// <summary>
        /// Cross-Section Lines
        /// </summary>
        public List<Geometry> pSectionLines { get; set; }

        /// <summary>
        /// Hinge point
        /// </summary>
        public List<Geometry> pHingePoints { get; set; }

        /// <summary>
        /// Nodal point on one side of the structural trend (2D)
        /// </summary>
        public Geometry pOneNode { get; set; }
        /// <summary>
        /// Nodal point on the opposite side of the structural trend (2D)
        /// </summary>
        public Geometry pAnotherNode { get; set; }

        /// <summary>
        /// Nodal point on one side of the structural trend (3D)
        /// </summary>
        public Geometry p3DOneNode { get; set; }
        /// <summary>
        /// Nodal point on the opposite side of the structural trend (3D)
        /// </summary>
        public Geometry p3DAnotherNode { get; set; }


        /// <summary>
        /// Outer boundary points (2D)
        /// </summary>
        public List<Geometry> pOutBoundaryPoints { get; set; }

        /// <summary>
        /// Outer boundary points (3D)
        /// </summary>
        public List<IGeometry> p3DOutBoundaryPoints { get; set; }

        /// <summary>
        /// Triangulated mesh of the upper surface
        /// </summary>
        public TriangleNet.Mesh pMeshOfUpTopSurface { get; set; }

        /// <summary>
        /// All 3D points on the upper surface
        /// </summary>
        public List<IGeometry> p3DPointsOfUpTopSurface { get; set; }

        /// <summary>
        /// Triangulated mesh of the lower surface
        /// </summary>
        public TriangleNet.Mesh pMeshOfBottomSurface { get; set; }


        /// <summary>
        /// Constraint points
        /// </summary>
        public List<GeoAPI.Geometries.IGeometry> constrainPoints { get; set; }


        /// <summary>
        /// Transition curve generated by morphing technique
        /// </summary>
        public List<Geometry> transitionLines { get; set; }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="_sid">stratigraphic id</param>
        /// <param name="_scode">stratigraphic code</param>
        /// <param name="_spolygon">stratigraphic polygon</param>
        public FoldStratum(int _sid, string _scode, Polygon3D _spolygon) : base(_sid, _scode, _spolygon)
        {
            beziers = new List<GeoBezier>();
            pSectionLines = new List<Geometry>();
            pHingePoints = new List<Geometry>();
            
            this.pOutBoundaryPoints = ExtractPointsOnOutBoundary(_spolygon.GeoPolygon);
        }


        /// <summary>
        /// Extract points from the outer boundary of the stratigraphic geometric surface
        /// </summary>
        /// <param name="_polygon">stratigraphic geometric surface</param>
        /// <returns></returns>
        public List<Geometry> ExtractPointsOnOutBoundary(Geometry _polygon)
        {
            Geometry pBoundaryPolyline = _polygon.GetBoundary();
            if(pBoundaryPolyline.GetGeometryType()!=wkbGeometryType.wkbLineString)
            {
                pBoundaryPolyline = pBoundaryPolyline.GetGeometryRef(0);
            }
            List<Geometry> pBoundaryPoints = new List<Geometry>(pBoundaryPolyline.GetPointCount());
            for (int i = 0;i< pBoundaryPolyline.GetPointCount();i++)
            {
                Geometry pt = new Geometry(wkbGeometryType.wkbPoint);
                pt.AddPoint_2D(pBoundaryPolyline.GetX(i), pBoundaryPolyline.GetY(i));
                pBoundaryPoints.Add(pt);
            }

            return pBoundaryPoints;
        }
    }
}
