using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows;

using GeoCommon;
using GeologicalEntity;
using GdalLib;

using System.Geometry;
using System.Numerics;

using MathNet.Numerics.Interpolation;

using System.Drawing;
using OSGeo.OGR;
using OSGeo.OSR;
using OSGeo.GDAL;

using NetTopologySuite;
using NetTopologySuite.Operation.Buffer;
using GeoAPI.Geometries;
using GeoAPI.Operation.Buffer;

using g4;



namespace FoldRestortEntry
{

    /// <summary>
    /// Main class for fold paleo-surface reconstruction
    /// </summary>
    public class FPRMain
    {
        /// <summary>
        /// Files save path
        /// </summary>
        private string filesavepath;

        public FPRMain(string _filesavepath)
        {
            this.filesavepath = _filesavepath;
        }




        /// <summary>
        /// Fold paleo-surface reconstruction for non-Dome fold
        /// </summary>
        /// <param name="_rasterD">Raster data</param>
        /// <param name="_foldtype">Fold type: 0 represents a fully closed fold, 1 represents a partially closed fold</param>
        /// <param name="_pStrataFile">Fold stratigraphic surface data (ESRI SHP)</param>
        /// <param name="_pOccurenceFile">Stratigraphic attitudes (ESRI SHP)</param>
        /// <param name="_pSectionLines">Cross-section Lines (ESRI SHP)</param>
        /// <param name="_pExtremeNodes">Extremity points at both ends of the fold (ESRI SHP)</param>
        public void FoldRestort(DEMRaster _rasterD, int _foldtype, string _pStrataFile, string _pOccurenceFile, List<Geometry> _pSectionLines, List<Geometry> _pExtremeNodes)
        {
            #region 1 Read stratigraphic surface data
            Layer shpLayer = GdalLib.LayerHelp.GetLayerByLayerName(_pStrataFile);

            int strataCount = (int)shpLayer.GetFeatureCount(0);
            List<FoldStratum> foldStrata = new List<FoldStratum>(strataCount);

            for (int i = 0; i < strataCount; i++)
            {
                string stratumName = shpLayer.GetFeature(i).GetFieldAsString("GeoCode");
                Geometry pg = shpLayer.GetFeature(i).GetGeometryRef();

                Polygon3D p3d = new Polygon3D(pg, stratumName);
                FoldStratum pStratum = new FoldStratum(i, stratumName, p3d);
                foldStrata.Add(pStratum);
            }

            #endregion

            #region 2 Read attitudes(or orientation) data and save them to an R-tree spatial index for efficient querying

            //2.1 find the position of the orientation point
            Layer occurrenceLayer = GdalLib.LayerHelp.GetLayerByLayerName(_pOccurenceFile);

            //2.2 Create a geometry factory object
            var geometryFactory = new NetTopologySuite.Geometries.GeometryFactory();


            //2.3 Create a R-tree spatial index
            var rtree = new NetTopologySuite.Index.Strtree.STRtree<OccurrencePoint>();

            strataCount = (int)occurrenceLayer.GetFeatureCount(0);
            List<OccurrencePoint> OccurrencePoints = new List<OccurrencePoint>(strataCount);

            for (int i = 0; i < strataCount; i++)
            {
                double dipAngle = occurrenceLayer.GetFeature(i).GetFieldAsDouble("DipAngle");
                double tendency = occurrenceLayer.GetFeature(i).GetFieldAsDouble("Tendency");
                double strike = occurrenceLayer.GetFeature(i).GetFieldAsDouble("Strike");
                Geometry pg = occurrenceLayer.GetFeature(i).GetGeometryRef();
                double pZ = _rasterD.GetElevation(pg.GetX(0), pg.GetY(0));
                OccurrencePoint oneOccurrence = new OccurrencePoint(wkbGeometryType.wkbPoint);
                oneOccurrence.AddPoint(pg.GetX(0), pg.GetY(0), pZ);
                oneOccurrence.dipAngle = dipAngle;
                oneOccurrence.tendency = tendency;
                oneOccurrence.strike = strike;
                OccurrencePoints.Add(oneOccurrence);

                // Save all attitude points to R-tree 
                var oc = geometryFactory.CreatePoint(new Coordinate(oneOccurrence.GetX(0), oneOccurrence.GetY(0)));
                rtree.Insert(oc.EnvelopeInternal, oneOccurrence);
            }
            #endregion

            #region 3 Determine the encompassing or outermost stratigraphic layer by comparing the areas of the bounding rectangles of all stratigraphic layers.
            
            //3.1 Obtain the bounding rectangle for each stratigraphic layer
            Dictionary<OSGeo.OGR.Envelope, FoldStratum> pEnvelopetoStratum = new Dictionary<OSGeo.OGR.Envelope, FoldStratum>();
            for (int i = 0; i < foldStrata.Count; i++)
            {
                OSGeo.OGR.Envelope pe = new OSGeo.OGR.Envelope();
                foldStrata[i].SPolygon.GeoPolygon.GetEnvelope(pe);
                pEnvelopetoStratum.Add(pe, foldStrata[i]);
            }

            //3.2 Calculate the area of each bounding rectangle
            Dictionary<double, OSGeo.OGR.Envelope> pEnvelopeArea = new Dictionary<double, OSGeo.OGR.Envelope>();
            List<OSGeo.OGR.Envelope> pEnvelope = pEnvelopetoStratum.Keys.ToList();
            for (int i = 0; i < pEnvelopetoStratum.Keys.ToList().Count; i++)
            {
                Geometry envelopePolygon = Geometry.CreateFromWkt(string.Format("POLYGON(({0} {1},{2} {1},{2} {3},{0} {3},{0} {1}))", pEnvelope[i].MinX, pEnvelope[i].MinY, pEnvelope[i].MaxX, pEnvelope[i].MaxY));
                pEnvelopeArea.Add(envelopePolygon.Area(), pEnvelope[i]);
            }
            //3.3 Sort by area to obtain the largest area
            List<double> pAreas = pEnvelopeArea.Keys.ToList();
            pAreas.Sort((a, b) => b.CompareTo(a));
            //3.4 Identify the stratigraphic layers corresponding to the largest bounding rectangles
            FoldStratum pOutestStratum = pEnvelopetoStratum[pEnvelopeArea[pAreas[0]]];
            //3.5 Rearrange the foldStrata in ascending order based on area
            List<FoldStratum> pNewOrderStrata = new List<FoldStratum>(foldStrata.Count);
            for (int i = 0; i < pAreas.Count; i++)
            {
                pNewOrderStrata.Insert(0, pEnvelopetoStratum[pEnvelopeArea[pAreas[i]]]);
            }

            #endregion

            #region 4-5 Solve the Cross-sectiion Lines
            /**
             * Calculate the bounding rectangle (not the minimum bounding rectangle) of the outermost stratigraphic layer and export it as a shapefile (SHP). 
             * The long axis of this rectangle can represent the fold's structural direction. 
             * Sample points along the two long sides of the bounding rectangle to construct a cross-section line perpendicular to the fold's structural direction
            */
            List<Geometry> pLines = null;// Cross-section Lines
            GeoAPI.Geometries.IGeometry env = null;//Minimum bounding rectangle
            int[] pOrderofWidth = new int[4];//Record the bounding rectangle's long and short side information, with [0]-[1] representing the long side and [2]-[3] representing the short side.


            // If you aleardy have some cross-section lines in advance
            if (_pSectionLines != null)
            {
                pLines = _pSectionLines;

                //To ensure the program runs correctly at step 8.2.1, it is still necessary to generate a rectangle here
                ConstructSectionLinesfromGeometryLine(pOutestStratum.SPolygon.GeoPolygon.GetGeometryRef(0), -1, 40, ref pOrderofWidth, ref env);

                //Export the rectangle and cross-section line for comparison and manual editing
                env.ExportSimpleGeometryToShapfileByNet(this.filesavepath, "env.shp");
            }
            else
            {
                //Specify cross-section line spacing or the number of cross-section lines
                pLines = ConstructSectionLinesfromGeometryLine(pOutestStratum.SPolygon.GeoPolygon.GetGeometryRef(0), -1, 39, ref pOrderofWidth, ref env);

                
                env.ExportSimpleGeometryToShapfileByNet(this.filesavepath, "env.shp");
                pLines.ExportGeometryToShapfile(this.filesavepath, "sections");
            }
            #endregion

            #region 6 Construct 3D stratigraphic paleo-boundaries based on cross-section lines.

            //6.1 Record the stratigraphic layers intersected by each cross-section line
            Dictionary<Geometry, List<FoldStratum>> pStratas = new Dictionary<Geometry, List<FoldStratum>>();
            for (int i = 0; i < pLines.Count; i++)
            {
                List<FoldStratum> subStrata = new List<FoldStratum>();
                for (int j = 0; j < pNewOrderStrata.Count; j++)
                {
                    if (pLines[i].Intersect(pNewOrderStrata[j].SPolygon.GeoPolygon))
                        subStrata.Add(pNewOrderStrata[j]);
                }
                pStratas.Add(pLines[i], subStrata);
            }

            //6.2 Calculate the intersection points between each cross-section line and the stratigraphic layers, and construct Bezier curves
            Dictionary<Geometry, List<GeoBezier>> pBeziersOnStrata = GetBeziersOnStrata(pStratas, rtree, geometryFactory, _rasterD);

            //6.3 Store in memory (this is the output before optimization)
            for (int i = 0; i < pNewOrderStrata.Count; i++)
            {
                List<Geometry> p3DBezierGeometrys = new List<Geometry>();
                List<Geometry> p2DBezierGeometrys = new List<Geometry>();

                for (int j = 0; j < pNewOrderStrata[i].beziers.Count; j++)
                {
                    GeoBezier _pBeGeo = pNewOrderStrata[i].beziers[j];
                    _pBeGeo.DrawSectionByBezierFunction(0.01f);
                    _pBeGeo.CalucateZeroTangent(0.01f, 0.01f);

                    p2DBezierGeometrys.Add(_pBeGeo._2DBezierGeometry);
                    p3DBezierGeometrys.Add(_pBeGeo._3DBezierGeometry);

                }
                p2DBezierGeometrys.ExportGeometryToShapfile(this.filesavepath, pNewOrderStrata[i].SCode + "2D");
                p3DBezierGeometrys.ExportGeometryToShapfile(this.filesavepath, pNewOrderStrata[i].SCode + "3D");
            }
            #endregion

            #region 7 Optimize 3D stratigraphic paleo-boundaries

            //7.1 Optimize the Bézier curve in cases where intersections occur
            OptimizeBeziers(ref pStratas, 6);

            //7.2 Store in memory (this is the output after optimization)
            for (int i = 0; i < pNewOrderStrata.Count; i++)
            {
                List<Geometry> p3DBezierGeometrys = new List<Geometry>();
                List<Geometry> p2DBezierGeometrys = new List<Geometry>();

                for (int j = 0; j < pNewOrderStrata[i].beziers.Count; j++)
                {
                    GeoBezier _pBeGeo = pNewOrderStrata[i].beziers[j];
                    //_pBeGeo.DrawSectionByBezierFunction(0.01f);
                    _pBeGeo.CalucateZeroTangent(0.01f, 0.01f);

                    p2DBezierGeometrys.Add(_pBeGeo._2DBezierGeometry);
                    p3DBezierGeometrys.Add(_pBeGeo._3DBezierGeometry);

                }
                p2DBezierGeometrys.ExportGeometryToShapfile(this.filesavepath, pNewOrderStrata[i].SCode + "2DMBezier");
                p3DBezierGeometrys.ExportGeometryToShapfile(this.filesavepath, pNewOrderStrata[i].SCode + "3DMBezier");
            }
            #endregion

            #region 8 Morphing interpolation

            //8.1 Reverse the stratigraphic layer collection to sort from outermost to innermost
            pNewOrderStrata.Reverse();

            //8.2 Obtain the endpoints on both sides of the structural direction for each stratigraphic layer

            //8.2.1 First, obtain and record the endpoints on both sides of the outermost stratigraphic layer.
            FoldStratum subStratum = pNewOrderStrata[0];
            if (_pExtremeNodes != null && _pExtremeNodes.Count == 2)
            {
                GetNodesBothStratum(ref subStratum, pLines[0], _pExtremeNodes);
            }
            else
            {
                GetNodesBothStratum(ref subStratum, env, pOrderofWidth);
            }

            //8.2.2 For all stratigraphic layers except the outermost,
            //select the nearest two nodes to pLeftNode and pRightNode as the endpoints along the structural direction
            for (int i = 1; i < pNewOrderStrata.Count; i++)
            {
                FoldStratum vs = pNewOrderStrata[i];
                if (vs != subStratum)
                    GetNodesBothStratum(ref vs, ref subStratum);
            }

            //8.2.3 Convert 2D endpoints to 3D
            List<Geometry> pnodes = new List<Geometry>(pNewOrderStrata.Count * 2);
            foreach (var vStratum in pNewOrderStrata)
            {
                Geometry p3dNode1 = new Geometry(wkbGeometryType.wkbPoint);
                p3dNode1.AddPoint(vStratum.pOneNode.GetX(0), vStratum.pOneNode.GetY(0), _rasterD.GetElevation(vStratum.pOneNode.GetX(0), vStratum.pOneNode.GetY(0)));
                vStratum.p3DOneNode = p3dNode1;

                Geometry p3dNode2 = new Geometry(wkbGeometryType.wkbPoint);
                p3dNode2.AddPoint(vStratum.pAnotherNode.GetX(0), vStratum.pAnotherNode.GetY(0), _rasterD.GetElevation(vStratum.pAnotherNode.GetX(0), vStratum.pAnotherNode.GetY(0)));
                vStratum.p3DAnotherNode = p3dNode2;

                pnodes.Add(vStratum.pOneNode);
                pnodes.Add(vStratum.pAnotherNode);
            }
            pnodes.ExportGeometryToShapfile(this.filesavepath, "pnodes");

            //8.3 Obtain the hinge points for each 3D stratigarphic paleo-boundaries
            //(including the nodes on both sides of the stratigraphic layer and the highest point of each Bezier curve)
            for (int i = 0; i < pNewOrderStrata.Count; i++)
            {
                pNewOrderStrata[i].pHingePoints.Add(pNewOrderStrata[i].p3DOneNode);
                foreach (var pStratumBezier in pNewOrderStrata[i].beziers)
                {
                    //pStratumBezier.DrawSectionByBezierFunction(0.01f);
                    pStratumBezier.Get3DTopestPointOnPolyline();
                    pNewOrderStrata[i].pHingePoints.Add(pStratumBezier._topPointsIn3DBezierGeometry);
                }
                pNewOrderStrata[i].pHingePoints.Add(pNewOrderStrata[i].p3DAnotherNode);

                pNewOrderStrata[i].pHingePoints.ExportGeometryToShapfile(this.filesavepath, pNewOrderStrata[i].SCode + "hige");
            }

            //8.4 (1) Segment the outer boundary of the stratigraphic layer into an equal number of opposite sides,
            //(2) then smoothly interpolate the hinge line,
            //and (3) use morphing technology to interpolate the transition curves.
            List<Geometry> pLinesToModifyDEM = new List<Geometry>();
            MorphingIntepolateTransmitionLines(_rasterD, geometryFactory, ref pNewOrderStrata, ref pLinesToModifyDEM, _foldtype);
            #endregion

            #region 9 Obtain all generated 3D points for interpolation to create a DEM, focusing on the 3D points of the outermost layer

            List<Geometry> pAllPoints = new List<Geometry>();

            //9.1 First, obtain the nodes on the 3D profile line
            for (int i = 0; i < pNewOrderStrata[0].beziers.Count; i++)
            {
                pLinesToModifyDEM.Add(pNewOrderStrata[0].beziers[i]._3DBezierGeometry);
            }

            //9.2 Then, obtain the points on the transition line generated by Morphing interpolation
            for (int i = 0; i < pLinesToModifyDEM.Count; i++)
            {
                for (int k = 0; k < pLinesToModifyDEM[i].GetPointCount(); k++)
                {
                    Geometry pt = new Geometry(wkbGeometryType.wkbPoint);

                    double dProjX = pLinesToModifyDEM[i].GetX(k);
                    double dProjY = pLinesToModifyDEM[i].GetY(k);
                    double dProjz = pLinesToModifyDEM[i].GetZ(k);

                    pt.AddPoint(dProjX, dProjY, dProjz);

                    pAllPoints.Add(pt);
                }
            }

            //9.3 Finally, obtain the points on the outermost boundary of the stratigraphic layer
            foreach (var vp in pNewOrderStrata[0].pOutBoundaryPoints)
            {
                Geometry pt = new Geometry(wkbGeometryType.wkbPoint);
                pt.AddPoint(vp.GetX(0), vp.GetY(0), _rasterD.GetElevation(vp.GetX(0), vp.GetY(0)));
                pAllPoints.Add(pt);
            }

            pAllPoints.ExportGeometryToShapfile(this.filesavepath, "elevationPoints");
            #endregion

            #region 10 Output all of the points on the paleo-boundaries, and reconstruct the paleo-surfaces of different strata
            
            for (int k = 0; k < pNewOrderStrata.Count; k++)
            {
                FoldStratum foldStratum = pNewOrderStrata[k];

                //(1) Add the points on the surface
                List<IGeometry> _topSurfacePoints = new List<IGeometry>();

                for (int j = 0; j < foldStratum.beziers.Count; j++)
                {
                    GeoBezier gb = foldStratum.beziers[j];

                    
                    for (int i = 0; i < gb._3DBezierGeometry.GetPointCount() - 1; i++)
                    {
                        Coordinate subpc = new Coordinate(gb._3DBezierGeometry.GetX(i), gb._3DBezierGeometry.GetY(i), gb._3DBezierGeometry.GetZ(i));

                        IGeometry pc = new NetTopologySuite.Geometries.Point(subpc);

                        if (!IsExistInPointSet(_topSurfacePoints, pc))
                            _topSurfacePoints.Add(pc);
                    }

                }

                for (int j = 0; j < foldStratum.transitionLines.Count; j++)
                {
                  
                    for (int i = 0; i < foldStratum.transitionLines[j].GetPointCount() - 1; i++)
                    {
                        Coordinate subpc = new Coordinate(foldStratum.transitionLines[j].GetX(i), foldStratum.transitionLines[j].GetY(i), foldStratum.transitionLines[j].GetZ(i));

                        IGeometry pc = new NetTopologySuite.Geometries.Point(subpc);

                        if (!IsExistInPointSet(_topSurfacePoints, pc))
                            _topSurfacePoints.Add(pc);
                    }

                }

                //(2) Add boundary points
                List<IGeometry> _boundaryPoints = new List<IGeometry>();
                for (int i = 0; i < foldStratum.pOutBoundaryPoints.Count; i++)
                {
                    double _x = foldStratum.pOutBoundaryPoints[i].GetX(0);
                    double _y = foldStratum.pOutBoundaryPoints[i].GetY(0);
                    double _z = _rasterD.GetElevation(_x, _y);
                    Coordinate subpt = new Coordinate(_x, _y, _z);
                    IGeometry pt = new NetTopologySuite.Geometries.Point(subpt);

                    //if (!IsExistInPointSet(_boundaryPoints, pt))
                    _boundaryPoints.Add(pt);

                    if (!IsExistInPointSet(_topSurfacePoints, pt))
                        _topSurfacePoints.Add(pt);
                }
                foldStratum.p3DOutBoundaryPoints = _boundaryPoints;
                foldStratum.p3DPointsOfUpTopSurface = _topSurfacePoints;

                //(3) Export to the filefolder
                _topSurfacePoints.ExportGeometryToShapfileByNet(this.filesavepath, foldStratum.SCode + "topsurfacepoints");
                _boundaryPoints.ExportGeometryToShapfileByNet(this.filesavepath, foldStratum.SCode + "boundaryPoints");

                //(4) Create tri_mesh
                TriangleNet.Mesh polygonMeshOfTopSurface = CreateDelaunaryTri(_topSurfacePoints, new List<List<IGeometry>>() { _boundaryPoints });

                //(5) export 3D model using the Class DMesh3
                DMesh3 dmesh = new DMesh3();

                List<TriangleNet.Geometry.Vertex> pvers = polygonMeshOfTopSurface.Vertices.ToList();
                foreach (var vt in pvers)
                {
                    double zvalue = GetElevationFromPoint(vt, _topSurfacePoints);

                    Vector3d v3d = new Vector3d(vt.X, vt.Y, zvalue);
                    dmesh.AppendVertex(v3d);
                }
                foreach (var vtri in polygonMeshOfTopSurface.Triangles)
                {
                    dmesh.AppendTriangle(vtri.GetVertexID(0), vtri.GetVertexID(1), vtri.GetVertexID(2));
                }

                //(6) Output
                StandardMeshWriter.WriteMesh(this.filesavepath + "\\" + foldStratum.SCode + ".obj", dmesh, WriteOptions.Defaults);
            }

            #endregion
        }


        /// <summary>
        /// For Domes
        /// </summary>
        /// <param name="_rasterD"></param>
        public void DomeTest(string _shpFile,string _occurrenceFile, DEMRaster _rasterD,double _intervalangle, int _sectionNumber)
        {
            #region 1 Read stratigraphic surface data

            Layer shpLayer = GdalLib.LayerHelp.GetLayerByLayerName(_shpFile);       

            int strataCount = (int)shpLayer.GetFeatureCount(0);
            List<FoldStratum> domeStrata = new List<FoldStratum>(strataCount);
         
            for (int i = 0; i < strataCount; i++)
            {
                string stratumName = shpLayer.GetFeature(i).GetFieldAsString("GeoCode");
                Geometry pg = shpLayer.GetFeature(i).GetGeometryRef();

                Polygon3D p3d = new Polygon3D(pg, stratumName);
                FoldStratum pStratum = new FoldStratum(i, stratumName, p3d);
                domeStrata.Add(pStratum);
            }
            #endregion

            #region 2  Read attitudes(or orientation) data and save them to an R-tree spatial index for efficient querying

            //2.1 find the position of the orientation point
            Layer occurrenceLayer = GdalLib.LayerHelp.GetLayerByLayerName(_occurrenceFile);

            //2.2 Create a geometry factory object
            var geometryFactory = new NetTopologySuite.Geometries.GeometryFactory();


            //2.3 Create a R tree spatial index
            var rtree = new NetTopologySuite.Index.Strtree.STRtree<OccurrencePoint>();

            strataCount = (int)occurrenceLayer.GetFeatureCount(0);
            List<OccurrencePoint> OccurrencePoints = new List<OccurrencePoint>(strataCount);

            for (int i = 0; i < strataCount; i++)
            {
                double dipAngle = occurrenceLayer.GetFeature(i).GetFieldAsDouble("DipAngle");
                double tendency = occurrenceLayer.GetFeature(i).GetFieldAsDouble("Tendency");
                double strike = occurrenceLayer.GetFeature(i).GetFieldAsDouble("Strike");
                Geometry pg = occurrenceLayer.GetFeature(i).GetGeometryRef();
                double pZ = _rasterD.GetElevation(pg.GetX(0), pg.GetY(0));
                OccurrencePoint oneOccurrence = new OccurrencePoint(wkbGeometryType.wkbPoint);
                oneOccurrence.AddPoint(pg.GetX(0), pg.GetY(0), pZ);
                oneOccurrence.dipAngle = dipAngle;
                oneOccurrence.tendency = tendency;
                oneOccurrence.strike = strike;
                OccurrencePoints.Add(oneOccurrence);

                // Save all attitude points to R-tree 
                var oc = geometryFactory.CreatePoint(new Coordinate(oneOccurrence.GetX(0), oneOccurrence.GetY(0)));
                rtree.Insert(oc.EnvelopeInternal, oneOccurrence);
            }
            #endregion

            #region 3 Determine the encompassing or outermost stratigraphic layer by comparing the areas of the bounding rectangles of all stratigraphic layers

            //3.1 Obtain the bounding rectangle for each stratigraphic layer
            Dictionary<OSGeo.OGR.Envelope, FoldStratum> pEnvelopetoStratum = new Dictionary<OSGeo.OGR.Envelope, FoldStratum>();
            for (int i = 0; i < domeStrata.Count; i++)
            {
                OSGeo.OGR.Envelope pe = new OSGeo.OGR.Envelope();
                domeStrata[i].SPolygon.GeoPolygon.GetEnvelope(pe);
                pEnvelopetoStratum.Add(pe, domeStrata[i]);
            }

            //3.2 Calculate the area of each bounding rectangle
            Dictionary<double, OSGeo.OGR.Envelope> pEnvelopeArea = new Dictionary<double, OSGeo.OGR.Envelope>();
            List<OSGeo.OGR.Envelope> pEnvelope = pEnvelopetoStratum.Keys.ToList();
            for (int i = 0; i < pEnvelopetoStratum.Keys.ToList().Count; i++)
            {
                Geometry envelopePolygon = Geometry.CreateFromWkt(string.Format("POLYGON(({0} {1},{2} {1},{2} {3},{0} {3},{0} {1}))", pEnvelope[i].MinX, pEnvelope[i].MinY, pEnvelope[i].MaxX, pEnvelope[i].MaxY));
                pEnvelopeArea.Add(envelopePolygon.Area(), pEnvelope[i]);
            }

            //3.3 Sort by area to obtain the largest area
            List<double> pAreas = pEnvelopeArea.Keys.ToList();
            pAreas.Sort((a, b) => b.CompareTo(a));

            //3.4 Identify the stratigraphic layers corresponding to the largest bounding rectangles
            FoldStratum pOutestStratum = pEnvelopetoStratum[pEnvelopeArea[pAreas[0]]];

            //3.5 Rearrange the foldStrata in ascending order based on area
            List<FoldStratum> pNewOrderStrata = new List<FoldStratum>(domeStrata.Count);
            for (int i = 0; i < pAreas.Count; i++)
            {
                pNewOrderStrata.Insert(0, pEnvelopetoStratum[pEnvelopeArea[pAreas[i]]]);
            }

            #endregion

            #region 4 Create the envelope
            /**
             * Calculate the bounding rectangle (not the minimum bounding rectangle) of the outermost stratigraphic layer and export it as a shapefile (SHP). 
             * The long axis of this rectangle can represent the fold's structural direction.
             */

            int cou = pOutestStratum.SPolygon.GeoPolygon.GetBoundary().GetGeometryRef(0).GetPointCount();
            Coordinate[] corrds = new Coordinate[cou];
            for (int i = 0; i < cou; i++)
            {
                corrds[i] = new Coordinate(pOutestStratum.SPolygon.GeoPolygon.GetGeometryRef(0).GetX(i), pOutestStratum.SPolygon.GeoPolygon.GetGeometryRef(0).GetY(i));
            }
            NetTopologySuite.Geometries.Polygon pNetPolygon = new NetTopologySuite.Geometries.Polygon(new NetTopologySuite.Geometries.LinearRing(corrds));
            GeoAPI.Geometries.IGeometry env = NetTopologySuite.Algorithm.MinimumDiameter.GetMinimumRectangle(pNetPolygon);
            env.ExportSimpleGeometryToShapfileByNet(this.filesavepath, "StrataEnvelope.shp");
            #endregion

            #region 5 Construct the cross-section lines

            /**
             * Approach: Use the centroid of the innermost stratigraphic layer as the starting point for the rays. 
             * Begin from the due north direction and set multiple azimuth angles in a clockwise direction to construct profile lines. 
             * The endpoint of each ray is the intersection with the minimum bounding rectangle.
             */
            
            Geometry pCentreid = pNewOrderStrata[0].SPolygon.GeoPolygon.Centroid();
            Coordinate pCenterPoint = new Coordinate(pCentreid.GetX(0), pCentreid.GetY(0));
            NetTopologySuite.Geometries.LineString pBoundaryLine = new NetTopologySuite.Geometries.LineString(env.Coordinates);
            List<NetTopologySuite.Geometries.LineString> pSectionLines = CreateSectionLines(pCenterPoint, _intervalangle, _sectionNumber, pBoundaryLine);
            List<Geometry> pOutSectionLines = new List<Geometry>();
            foreach(var vp in pSectionLines)
            {
                pOutSectionLines.Add(ConvertNetTLineToOGRLine(vp,false));
            }
            pOutSectionLines.ExportGeometryToShapfile(this.filesavepath, "StrataSectionLines");
            #endregion

            #region 6 Construct 3D stratigraphic lines based on the cross-section lines

            //6.1 Record the stratigraphic layers intersected by each cross-section line
            Dictionary<Geometry, List<FoldStratum>> pStratas = new Dictionary<Geometry, List<FoldStratum>>();
            for (int i = 0; i < pOutSectionLines.Count; i++)
            {
                List<FoldStratum> subStrata = new List<FoldStratum>();
                for (int j = 0; j < pNewOrderStrata.Count; j++)
                {
                    if (pOutSectionLines[i].Intersect(pNewOrderStrata[j].SPolygon.GeoPolygon))
                        subStrata.Add(pNewOrderStrata[j]);
                }
                pStratas.Add(pOutSectionLines[i], subStrata);
            }

            //6.2 Calculate the intersection points between each cross-section line and the stratigraphic layers, and construct Bezier curves
            GetBeziersOnDomeStrata(pStratas, rtree, geometryFactory, _rasterD);

            //6.3  Store in memory (this is the output before optimization)
            for (int i = 0; i < pNewOrderStrata.Count; i++)
            {
                List<Geometry> p3DBezierGeometrys = new List<Geometry>();
                List<Geometry> p2DBezierGeometrys = new List<Geometry>();

                for (int j = 0; j < pNewOrderStrata[i].beziers.Count; j++)
                {
                    GeoBezier _pBeGeo = pNewOrderStrata[i].beziers[j];
                    _pBeGeo.DrawSectionByBezierFunction(0.01f);                              
                    _pBeGeo.CalucateZeroTangent(0.01f, 0.01f);

                    p2DBezierGeometrys.Add(_pBeGeo._2DBezierGeometry);
                    p3DBezierGeometrys.Add(_pBeGeo._3DBezierGeometry);

                }
                p2DBezierGeometrys.ExportGeometryToShapfile(this.filesavepath, pNewOrderStrata[i].SCode + "2D");
                p3DBezierGeometrys.ExportGeometryToShapfile(this.filesavepath, pNewOrderStrata[i].SCode + "3D");
            }
            #endregion

            #region 7 Optimize 3D stratigraphic paleo-boundaries

            //7.1 Optimize the Bézier curve in cases where intersections occur
            OptimizeBeziers(ref pStratas,2);

            //7.2 Store in memory (this is the output after optimization)
            for (int i = 0; i < pNewOrderStrata.Count; i++)
            {
                List<Geometry> p3DBezierGeometrys = new List<Geometry>();
                List<Geometry> p2DBezierGeometrys = new List<Geometry>();

                for (int j = 0; j < pNewOrderStrata[i].beziers.Count; j++)
                {
                    GeoBezier _pBeGeo = pNewOrderStrata[i].beziers[j];
                    _pBeGeo.CalucateZeroTangent(0.01f, 0.01f);
                    _pBeGeo.Get3DTopestPointOnPolyline();

                    p2DBezierGeometrys.Add(_pBeGeo._2DBezierGeometry);
                    p3DBezierGeometrys.Add(_pBeGeo._3DBezierGeometry);

                }
                p2DBezierGeometrys.ExportGeometryToShapfile(this.filesavepath, pNewOrderStrata[i].SCode + "2DMBezier");
                p3DBezierGeometrys.ExportGeometryToShapfile(this.filesavepath, pNewOrderStrata[i].SCode + "3DMBezier");
            }


            //7.3 The elevation and position of the highest points on each Bézier curve within the same stratigraphic layer vary;
            //determine these using the least squares method

            //(1) Determine the z-value of the highest point using the least squares method
            for (int k = 0; k < pNewOrderStrata.Count; k++)
            {
                FoldStratum foldStratum = pNewOrderStrata[k];
                List<double> _3DHighestZOnBeziers = new List<double>();
                for (int i = 0; i < foldStratum.beziers.Count; i++)
                {
                    if (foldStratum.beziers[i]._topPointsIn3DBezierGeometry != null)
                    {
                        _3DHighestZOnBeziers.Add(foldStratum.beziers[i]._topPointsIn3DBezierGeometry.GetZ(0));
                    }
                }
                double averageZValue = _3DHighestZOnBeziers.Average();

                //(2)At the point where the tangent is horizontal (the highest point on the 3D Bézier curve),
                //split the profile line and construct a new Bézier curve

                List<GeoBezier> pNewGeoBeziers = new List<GeoBezier>();

                for (int i = 0; i < foldStratum.beziers.Count; i++)
                {
                    if (foldStratum.beziers[i]._topPointsIn3DBezierGeometry != null)
                    {

                        //At the highest point of the 3D Bézier curve, create two new virtual attitude points with opposite dip directions and a dip angle of 0

                        OccurrencePoint newoccptToLeft = new OccurrencePoint(wkbGeometryType.wkbPoint);
                        newoccptToLeft.AddPoint(pCentreid.GetX(0), pCentreid.GetY(0), averageZValue);
                        newoccptToLeft.dipAngle = 0.0;
                        newoccptToLeft.tendency = foldStratum.beziers[i]._azimuthAngleofSectionLine - 180.0;
                        newoccptToLeft.strike = newoccptToLeft.tendency + 90.0;
                        if (newoccptToLeft.strike > 360.0)
                            newoccptToLeft.strike = newoccptToLeft.strike - 360.0;

                        OccurrencePoint newoccptToRight = new OccurrencePoint(wkbGeometryType.wkbPoint);
                        newoccptToRight.AddPoint(pCentreid.GetX(0), pCentreid.GetY(0), averageZValue);
                        newoccptToRight.dipAngle = 0.0;
                        newoccptToRight.tendency = foldStratum.beziers[i]._azimuthAngleofSectionLine;
                        newoccptToRight.strike = newoccptToLeft.tendency + 90.0;
                        if (newoccptToRight.strike > 360.0)
                            newoccptToRight.strike = newoccptToRight.strike - 360.0;

                        
                        GeoBezier pGeoBezierLeft = new GeoBezier(foldStratum.beziers[i]._geoBezierName + "Left", foldStratum.beziers[i]._geoLeftOccurrence, newoccptToLeft, foldStratum.beziers[i]._geoMoveDistance, this.filesavepath);
                        GeoBezier pGeoBezierRight = new GeoBezier(foldStratum.beziers[i]._geoBezierName + "Left", newoccptToRight, foldStratum.beziers[i]._geoRightOccurrence, foldStratum.beziers[i]._geoMoveDistance, this.filesavepath);

                        
                        pNewGeoBeziers.Add(pGeoBezierLeft);
                        pNewGeoBeziers.Add(pGeoBezierRight);
                    }
                }
                //(3) Replace the original Bézier curve
                foldStratum.beziers = pNewGeoBeziers;

                //(4）Output
                List<Geometry> ptest3DBezierGeometrys = new List<Geometry>();
                List<Geometry> ptest2DBezierGeometrys = new List<Geometry>();

                for (int j = 0; j < foldStratum.beziers.Count; j++)
                {
                    GeoBezier _pBeGeo = foldStratum.beziers[j];
                    _pBeGeo.DrawSectionByBezierFunction(0.01f);
                    _pBeGeo.CalucateZeroTangent(0.01f, 0.01f);
                    _pBeGeo.Get3DTopestPointOnPolyline();

                    ptest2DBezierGeometrys.Add(_pBeGeo._2DBezierGeometry);
                    ptest3DBezierGeometrys.Add(_pBeGeo._3DBezierGeometry);

                }
                ptest2DBezierGeometrys.ExportGeometryToShapfile(this.filesavepath, foldStratum.SCode + "new2DMBezier");
                ptest3DBezierGeometrys.ExportGeometryToShapfile(this.filesavepath, foldStratum.SCode + "new3DMBezier");

            }

            #endregion

            #region 8 Construct the 3D surface of the dome

            for (int k = 0; k < pNewOrderStrata.Count; k++)
            {
                FoldStratum foldStratum = pNewOrderStrata[k];

                #region 8.1 Construct the top surface triangulated mesh (TIN) for each stratigraphic layer.

                //(1) Add points on the top surface
                List<IGeometry> _topSurfacePoints = new List<IGeometry>();
                
                for (int j = 0; j < foldStratum.beziers.Count; j++)
                {
                    GeoBezier gb = foldStratum.beziers[j];

                    for (int i=0;i< gb._3DBezierGeometry.GetPointCount()-1;i++)
                    {
                        Coordinate subpc = new Coordinate(gb._3DBezierGeometry.GetX(i), gb._3DBezierGeometry.GetY(i), gb._3DBezierGeometry.GetZ(i));

                        IGeometry pc = new NetTopologySuite.Geometries.Point(subpc);

                        if (!IsExistInPointSet(_topSurfacePoints,pc))
                            _topSurfacePoints.Add(pc);
                    }
                    
                }

                //(2) Add boundary constrain
                List<IGeometry> _boundaryPoints = new List<IGeometry>();
                for (int i=0;i< foldStratum.pOutBoundaryPoints.Count;i++)
                {
                    double _x = foldStratum.pOutBoundaryPoints[i].GetX(0);
                    double _y = foldStratum.pOutBoundaryPoints[i].GetY(0);
                    double _z = _rasterD.GetElevation(_x, _y);
                    Coordinate subpt = new Coordinate(_x, _y, _z);
                    IGeometry pt = new NetTopologySuite.Geometries.Point(subpt);

                    //if (!IsExistInPointSet(_boundaryPoints, pt))
                        _boundaryPoints.Add(pt);

                    if (!IsExistInPointSet(_topSurfacePoints, pt))
                        _topSurfacePoints.Add(pt);
                }
                foldStratum.p3DOutBoundaryPoints = _boundaryPoints;
                foldStratum.p3DPointsOfUpTopSurface = _topSurfacePoints;

                //(3)  
                _topSurfacePoints.ExportGeometryToShapfileByNet(this.filesavepath, foldStratum.SCode + "topsurfacepoints");
                _boundaryPoints.ExportGeometryToShapfileByNet(this.filesavepath, foldStratum.SCode + "boundaryPoints");

                //(4) generate triangulated irregular network
                TriangleNet.Mesh polygonMeshOfTopSurface = CreateDelaunaryTri(_topSurfacePoints, new List<List<IGeometry>>() { _boundaryPoints });
                foldStratum.pMeshOfUpTopSurface = polygonMeshOfTopSurface;

                polygonMeshOfTopSurface.ExportTriMeshToShapfile(this.filesavepath, foldStratum.SCode + "Triangles");

                //(5) export 3D model type using the Class DMesh3
                DMesh3 dmesh = new DMesh3();

                List<TriangleNet.Geometry.Vertex> pvers = polygonMeshOfTopSurface.Vertices.ToList();
                foreach (var vt in pvers)
                {
                    double zvalue = GetElevationFromPoint(vt, _topSurfacePoints);

                    Vector3d v3d = new Vector3d(vt.X, vt.Y, zvalue);
                    dmesh.AppendVertex(v3d);
                }
                foreach (var vtri in polygonMeshOfTopSurface.Triangles)
                {
                    dmesh.AppendTriangle(vtri.GetVertexID(0), vtri.GetVertexID(1), vtri.GetVertexID(2));
                }

                //(6) Output
                StandardMeshWriter.WriteMesh(this.filesavepath + "\\" + foldStratum.SCode + ".obj", dmesh, WriteOptions.Defaults);

                #endregion

            }

            #endregion
     
           
        }


        #region Preprocessing functions
        /// <summary>
        /// Primarily for cases where the initially generated profile lines are suboptimal, regenerate the profile lines using a rectangle; 
        /// needs to be rewritten (2023-05-27)）
        /// </summary>
        /// <param name="_pShpFileOfPolygon">geometric polygons filepath</param>
        /// <param name="_sampledistance">sample distance</param>
        public void CreateSectionLines(string _pShpFileOfRectangle,double _sampledistance)
        {

            //todo

            Layer shpLayer = GdalLib.LayerHelp.GetLayerByLayerName(_pShpFileOfRectangle);

            int strataCount = (int)shpLayer.GetFeatureCount(0);

            for (int i = 0; i < strataCount; i++)
            {
                Geometry pg = shpLayer.GetFeature(i).GetGeometryRef();

                GeoAPI.Geometries.IGeometry env = null;
                int[] pOrderofWidth = new int[4];
               
                List<Geometry> pLines = ConstructSectionLinesfromGeometryLine(pg.GetGeometryRef(0), _sampledistance, -1, ref pOrderofWidth, ref env);

                env.ExportSimpleGeometryToShapfileByNet(this.filesavepath, "env"+i+ ".shp");
                pLines.ExportGeometryToShapfile(this.filesavepath, "sections"+i);
            }
        }

        /// <summary>
        /// Primarily for cases where the initially generated profile lines are suboptimal, regenerate the profile lines using a rectangle;
        /// </summary>
        /// <param name="_pShpFileOfPolygon">geometric polygons filepath</param>
        /// <param name="_sampledistance">sample distance</param>
        /// <param name="_filepath">file save path</param>
        /// <param name="_shpname">filename</param>
        public void CreateSectionLines(string _pShpFileOfPolygon, double _sampledistance, string _filepath, string _shpname)
        {
            Layer shpLayer = GdalLib.LayerHelp.GetLayerByLayerName(_pShpFileOfPolygon);

            int strataCount = (int)shpLayer.GetFeatureCount(0);

            for (int i = 0; i < strataCount; i++)
            {
                Geometry pg = shpLayer.GetFeature(i).GetGeometryRef();

                GeoAPI.Geometries.IGeometry env = null;//Minimum bounding rectangle (MBR)

                //Record the information of the bounding rectangle’s long and short sides: [0]-[1] represent the long side, and [2]-[3] represent the short side.
                int[] pOrderofWidth = new int[4];

                //Specify profile line spacing or specify the number of profile lines
                List<Geometry> pLines = ConstructSectionLinesfromGeometryLine(pg.GetGeometryRef(0), _sampledistance, -1, ref pOrderofWidth, ref env);

                //Output the rectangle and cross-section lines for observation
                env.ExportSimpleGeometryToShapfileByNet(_filepath, _shpname+"env" + i + ".shp");
                pLines.ExportGeometryToShapfile(_filepath, _shpname+"sections" + i);
            }
        }

        /// <summary>
        /// Primarily for cases where the initially generated cross-section lines are suboptimal, 
        /// regenerate the profile lines using the structural orientation lines of the stratigraphic surface
        /// </summary>
        /// <param name="_pShpFileOfLine">Line representing the strike of the stratigraphy</param>
        /// <param name="_sampledistance">Sample distance</param>
        /// <param name="_verticalDistance">The extended distance of the cross-section line</param>
        public void CreateSectionLines(string _pShpFileOfLine, double _sampledistance, double _verticalDistance, string _filesavepath, string _filename)
        {
            Layer shpLayer = GdalLib.LayerHelp.GetLayerByLayerName(_pShpFileOfLine);

            int strataCount = (int)shpLayer.GetFeatureCount(0);

            for (int i = 0; i < strataCount; i++)
            {
                Geometry pg = shpLayer.GetFeature(i).GetGeometryRef();

                
                List<Geometry> pLines = ConstructSectionLines(pg, _sampledistance, _verticalDistance);

                
                pLines.ExportGeometryToShapfile(_filesavepath, _filename + "SectionLines" + i);
            }
        }

        /// <summary>
        /// Obtain the pre-prepared cross-section line data
        /// </summary>
        /// <param name="_sectionLineShpFile"></param>
        /// <returns></returns>
        public List<Geometry> GetGeometrysfromFile(string _sectionLineShpFile)
        {
            List<Geometry> pAllSections = new List<Geometry>();

            Layer shpLayer = GdalLib.LayerHelp.GetLayerByLayerName(_sectionLineShpFile);

            int strataCount = (int)shpLayer.GetFeatureCount(0);

            for (int i = 0; i < strataCount; i++)
            {
                Geometry pg = shpLayer.GetFeature(i).GetGeometryRef();
                pAllSections.Add(pg);
            }

            return pAllSections;
        }

        #endregion

        #region Functions required in this class

        #region some common functions

        /// <summary>
        /// Given a direction vector, distance, and starting coordinates, calculate the endpoint coordinates
        /// </summary>
        /// <param name="_startPointX"></param>
        /// <param name="_startPointY"></param>
        /// <param name="_directionX"></param>
        /// <param name="_directionY"></param>
        /// <param name="endPoint"></param>
        public void CalculateEndPoint(double _startPointX, double _startPointY,double _directionX,double _directionY,double _distance ,ref double [] _endPoint)
        {
            //length
            double magnitude = Math.Sqrt(_directionX * _directionX + _directionY * _directionY);

            // Normalize the direction vector
            double normalizedDirectionX = _directionX / magnitude;
            double normalizedDirectionY = _directionY / magnitude;

            // Calculate endpoint coordinate
            double endX = _startPointX + normalizedDirectionX * _distance;
            double endY = _startPointY + normalizedDirectionY * _distance;

            _endPoint[0] = endX;
            _endPoint[1] = endY;
        }

        /// <summary>
        /// Find the point in a set of points that is closest to a given line
        /// </summary>
        /// <param name="_points">point set</param>
        /// <param name="_pLine">Line</param>
        /// <returns></returns>
        public Geometry FindNearestPointToLine(List<Geometry> _points, Geometry _pLine)
        {
            Dictionary<double, Geometry> pPointsDistance = new Dictionary<double, Geometry>();
            for (int i = 0; i < _points.Count; i++)
            {
                if (pPointsDistance.Keys.Contains(_points[i].Distance(_pLine)))
                    continue;

                pPointsDistance.Add(_points[i].Distance(_pLine), _points[i]);
            }

            List<double> pDistances = pPointsDistance.Keys.ToList();
            pDistances.Sort((a, b) => a.CompareTo(b));

            return pPointsDistance[pDistances[0]];
        }

        /// <summary>
        /// Find the point in a set of points that is closest to a target point
        /// </summary>
        /// <param name="targetPoint">target point</param>
        /// <param name="points">point set</param>
        /// <returns></returns>
        public Coordinate FindNearestPoint(Coordinate targetPoint, List<Coordinate> points)
        {
            Coordinate nearestPoint = null;
            double minDistance = double.MaxValue;

            foreach (var point in points)
            {
                double distance = targetPoint.Distance(point);

                if (distance < minDistance)
                {
                    minDistance = distance;
                    nearestPoint = point;
                }
            }

            return nearestPoint;
        }

        /// <summary>
        /// Find the point in a set of points that is closest to a target point
        /// </summary>
        /// <param name="targetPoint">target point</param>
        /// <param name="points">point set</param>
        /// <returns></returns>
        public Geometry FindNearestPoint(Geometry targetPoint, List<Geometry> points)
        {
            Geometry nearestPoint = null;
            double minDistance = double.MaxValue;

            foreach (var point in points)
            {
                double distance = targetPoint.Distance(point);

                if (distance < minDistance)
                {
                    minDistance = distance;
                    nearestPoint = point;
                }
            }

            return nearestPoint;
        }

        /// <summary>
        /// Convert points from GDAL Geometry to GeoAPI points
        /// </summary>
        /// <param name="_point">Geometry point</param>
        /// <param name="_3D">Is 3D?</param>
        /// <returns></returns>
        public Coordinate ConvertGeometryToCoordinate(Geometry _point, bool _3D)
        {
            if (_3D)
                return new Coordinate(_point.GetX(0), _point.GetY(0), _point.GetZ(0));
            else
                return new Coordinate(_point.GetX(0), _point.GetY(0));
        }

        /// <summary>
        /// Convert points from GeoAPI points to GDAL Geometry 
        /// </summary>
        /// <param name="_point">Coordinate</param>
        /// <param name="_3D">Is 3D?</param>
        /// <returns></returns>
        public Geometry ConvertCoordinateToGeometry(Coordinate _point, bool _3D)
        {
            Geometry pt = new Geometry(wkbGeometryType.wkbPoint);
            if (_3D)
                pt.AddPoint(_point.X, _point.Y, _point.Z); 
            else
                pt.AddPoint_2D(_point.X, _point.Y);

            return pt;
        }

        /// <summary>
        /// Convert a sequence of 2D Coordinates into a 3D Geometry line
        /// </summary>
        /// <param name="_pLineCoordintes">point set</param>
        /// <param name="_raster">DEM</param>
        /// <returns></returns>
        public Geometry ConvertToGeometryLine(List<Coordinate> _pLineCoordintes, DEMRaster _raster)
        {
            Geometry pLine = new Geometry(wkbGeometryType.wkbLineString);
            for (int i = 0; i < _pLineCoordintes.Count; i++)
            {

                pLine.AddPoint(_pLineCoordintes[i].X, _pLineCoordintes[i].Y, _raster.GetElevation(_pLineCoordintes[i].X, _pLineCoordintes[i].Y));
            }

            return pLine;
        }

        /// <summary>
        /// Convert a sequence of Coordinates into a Geometry line
        /// </summary>
        /// <param name="_pLineCoordintes">Coordinate 点</param>
        /// <param name="_haveZvalue">Does a Z-value already exist?</param>
        /// <returns></returns>
        public Geometry ConvertToGeometryLine(List<Coordinate> _pLineCoordintes,bool _haveZvalue)
        {
            Geometry pLine = new Geometry(wkbGeometryType.wkbLineString);

            if (_haveZvalue)
            {
                for (int i = 0; i < _pLineCoordintes.Count; i++)
                {
                    pLine.AddPoint(_pLineCoordintes[i].X, _pLineCoordintes[i].Y, _pLineCoordintes[i].Z);
                }
            }
            else
            {
                for (int i = 0; i < _pLineCoordintes.Count; i++)
                {
                    pLine.AddPoint_2D(_pLineCoordintes[i].X, _pLineCoordintes[i].Y);
                }
            }

            return pLine;
        }



        /// <summary>
        /// Convert a line from NetTopologySuite to OGR
        /// </summary>
        /// <param name="_pNetTLine">line in NetTopologySuite</param>
        /// <param name="_zvalue">Does a Z-value already exist?</param>
        /// <returns></returns>
        public Geometry ConvertNetTLineToOGRLine(NetTopologySuite.Geometries.LineString _pNetTLine, bool _zvalue)
        {
            Geometry pORGLine = new Geometry(wkbGeometryType.wkbLineString);

            Coordinate[] linePoints = _pNetTLine.Coordinates;
            
            for(int i=0;i< linePoints.Length;i++)
            {
                if(_zvalue)
                {
                    pORGLine.AddPoint(linePoints[i].X, linePoints[i].Y, linePoints[i].Z);
                }
                else
                {
                    pORGLine.AddPoint_2D(linePoints[i].X, linePoints[i].Y);
                }
            }

            return pORGLine;
        }


        /// <summary>
        /// Convert a line from OGR to NetTopologySuite
        /// </summary>
        /// <param name="_pOGRLine">line in ogr</param>
        /// <param name="_zvalue">Does a Z-value already exist?</param>
        /// <returns></returns>
        public NetTopologySuite.Geometries.LineString ConvertOGRLineNetTLine(Geometry _pOGRLine, bool _zvalue)
        {
            Coordinate[] linePoints = new Coordinate[_pOGRLine.GetPointCount()];

            for (int i = 0; i < _pOGRLine.GetPointCount(); i++)
            {
                if (_zvalue)
                {
                    linePoints[i] = new Coordinate(_pOGRLine.GetX(i), _pOGRLine.GetY(i), _pOGRLine.GetZ(i));                  
                }
                else
                {
                    linePoints[i] = new Coordinate(_pOGRLine.GetX(i), _pOGRLine.GetY(i));
                }
            }

            NetTopologySuite.Geometries.LineString pNetTLine = new NetTopologySuite.Geometries.LineString(linePoints);

            return pNetTLine;
        }


        /// <summary>
        /// Check if a point is within the point set by using the absolute difference in coordinates
        /// </summary>
        /// <param name="_points">point set</param>
        /// <param name="_pc">the point need be checked</param>
        /// <returns></returns>
        public bool IsExistInPointSet(List<IGeometry> _points,IGeometry _pc)
        {
            foreach(var vp in _points)
            {
                if (Math.Abs(_pc.Coordinate.X - vp.Coordinate.X) <= 0.001 && Math.Abs(_pc.Coordinate.Y- vp.Coordinate.Y) <= 0.001 && Math.Abs(_pc.Coordinate.Z - vp.Coordinate.Z) <= 0.001)
                    return true;               
            }
            return false;
        }

        /// <summary>
        /// Rotate a geometry surface
        /// </summary>
        /// <param name="_point">Boundary point set of a geometry surface</param>
        /// <param name="_centrePoint">Centroid of a geometry surface</param>
        /// <param name="_angle">Rotation angle (degrees)</param>
        /// <returns></returns>
        private Coordinate RotatePoint(Coordinate _point, IPoint _centrePoint, double _angle)
        {
            double ix = _point.X;
            double iy = _point.Y;
            double cx = _centrePoint.X;
            double cy = _centrePoint.Y;
            
            double radians = _angle * Math.PI / 180;

            double ox, oy;
            ox = (ix - cx) * Math.Cos(radians) - (iy - cy) * Math.Sin(radians) + cx;   //旋转公式
            oy = (ix - cx) * Math.Sin(radians) + (iy - cy) * Math.Cos(radians) + cy;


            return new Coordinate(ox, oy);
        }

        #endregion

        #region Functions related to constructing a triangulated mesh

        /// <summary>
        /// Get the elevation of a point, used only during triangulated mesh construction
        /// </summary>
        /// <param name="_vt"></param>
        /// <param name="_points"></param>
        /// <returns></returns>
        public double GetElevationFromPoint(TriangleNet.Geometry.Vertex _vt, List<IGeometry> _points)
        {
            foreach(var vp in _points)
            {
                if (Math.Abs(_vt.X - vp.Coordinate.X) < 0.001 && Math.Abs(_vt.Y - vp.Coordinate.Y) < 0.001)
                    return vp.Coordinate.Z;
            }

            return -1;
        }

        /// <summary>
        /// Rebuild the triangulated mesh (2D)
        /// </summary>
        /// <param name="_pAllPoints"></param>
        /// <param name="_pBoundaryConstrainPoints"></param>
        /// <returns></returns>
        public TriangleNet.Mesh CreateDelaunaryTri(List<IGeometry> _pAllPoints, List<List<IGeometry>> _pBoundaryConstrainPoints)
        {
            //1. constrain options
            var options = new TriangleNet.Meshing.ConstraintOptions();
            options.SegmentSplitting = 2;
            options.ConformingDelaunay = true;

            //2. constrcut IPolygon
            TriangleNet.Geometry.IPolygon data = new TriangleNet.Geometry.Polygon();
            foreach (var vt in _pAllPoints)
            {
                TriangleNet.Geometry.Vertex pt = new TriangleNet.Geometry.Vertex();
                pt.X = vt.Coordinate.X;
                pt.Y = vt.Coordinate.Y;

                data.Add(pt);
            }

            foreach (var vps in _pBoundaryConstrainPoints)
            {
                List<TriangleNet.Geometry.Vertex> pContour = new List<TriangleNet.Geometry.Vertex>();

                foreach (var vt in vps)
                {
                    TriangleNet.Geometry.Vertex pt = new TriangleNet.Geometry.Vertex();
                    pt.X = vt.Coordinate.X;
                    pt.Y = vt.Coordinate.Y;
                    pContour.Add(pt);
                }

                TriangleNet.Geometry.Contour pNewCon = new TriangleNet.Geometry.Contour(pContour);

                if (_pBoundaryConstrainPoints[0] == vps)
                    data.Add(pNewCon, false);
                else
                    data.Add(pNewCon, true);
            }

            //quality options
            var quality = new TriangleNet.Meshing.QualityOptions();

            TriangleNet.Mesh mesh = null;
            if (data != null)
            {
                mesh = (TriangleNet.Mesh)TriangleNet.Geometry.ExtensionMethods.Triangulate(data, options);
            }

            return mesh;
        }

        /// <summary>
        /// Insert a specified number of geometric points within a geometry surface using the minimum bounding rectangle method
        /// </summary>
        /// <param name="_pPolygonOfStratum">geometric surface</param>
        /// <param name="_interstep">step length</param>
        /// <returns></returns>
        public List<IGeometry> DisperseRoi(Geometry _pPolygonOfStratum, double _interstep)
        {           
            List<IGeometry> pInterPoints = new List<IGeometry>();

            //Get the bounding rectangle
            OSGeo.OGR.Envelope pEnv = new OSGeo.OGR.Envelope();
            _pPolygonOfStratum.GetEnvelope(pEnv);

            //Four extremum values
            double xMin = pEnv.MinX;
            double xMax = pEnv.MaxX;
            double yMin = pEnv.MinY;
            double yMax = pEnv.MaxY;

            //side length
            int xSize = Convert.ToInt32(Math.Round((xMax - xMin) / _interstep));
            int ySize = Convert.ToInt32(Math.Round((yMax - yMin) / _interstep));

            //Select points located within the study area
            for (int i = 0; i <= xSize; i++)
            {
                for (int j = 0; j <= ySize; j++)
                {
                    Geometry point = new Geometry(wkbGeometryType.wkbPoint);
                    point.AddPoint_2D(xMin + i * _interstep, yMin + j * _interstep);

                    IGeometry pt = new NetTopologySuite.Geometries.Point(new Coordinate(xMin + i * _interstep, yMin + j * _interstep));

                    if (_pPolygonOfStratum.Contains(point))
                        pInterPoints.Add(pt);// The point is within the study area
                }
            }

            return pInterPoints;
        }

        #endregion


        #region 4-5 Find the minimum bounding rectangle of the outer stratigraphic layer and construct parallel profile lines within it

        /// <summary>
        /// Construct cross-section lines using the minimum bounding rectangle of a geometry surface
        /// </summary>
        /// <param name="pBoundaryLine">Outer boundary of a geometry surface</param>
        /// <param name="_distance">Specify the distance, if false set a negative value</param>
        /// <param name="_sectionLineNumber">Specify the number of the cross-section lines, if false set a negative value</param>
        /// <param name="_pLengthBoundaryofEnv">Record the information of the rectangle's long and short sides</param>
        /// <param name="_env">Minimum Bounding Rectangle</param>
        /// <returns></returns>
        public List<Geometry> ConstructSectionLinesfromGeometryLine(Geometry pBoundaryLine,double _distance,int _sectionLineNumber, ref int[] _pLengthBoundaryofEnv, ref GeoAPI.Geometries.IGeometry _env)
        {
            int cou = pBoundaryLine.GetPointCount();
            Coordinate[] corrds = new Coordinate[cou];
            for (int i = 0; i < cou; i++)
            {
                corrds[i] = new Coordinate(pBoundaryLine.GetX(i), pBoundaryLine.GetY(i));
            }
            NetTopologySuite.Geometries.Polygon pNetPolygon = new NetTopologySuite.Geometries.Polygon(new NetTopologySuite.Geometries.LinearRing(corrds));
            _env = NetTopologySuite.Algorithm.MinimumDiameter.GetMinimumRectangle(pNetPolygon);
                      
            int[] pOrderofWidth = new int[4];
            List<Geometry> pLines = null;
            
            
            if (_distance>0)
            pLines = ConstructSectionLines(_env, _distance, ref _pLengthBoundaryofEnv);
            
            
            if(_sectionLineNumber>0)
                pLines = ConstructSectionLines(_env, _sectionLineNumber, ref _pLengthBoundaryofEnv);

            return pLines;
        }

        /// <summary>
        /// 5 Construct cross-section lines within the bounding rectangle along the structural direction
        /// </summary>
        /// <param name="_env">the bounding rectangle</param>
        /// <param name="sectlionLineNumber">the number of cross-section lines</param>
        /// <param name="pLengthBoundaryofEnv">record the information of the long side</param>
        /// <returns></returns>
        public List<Geometry> ConstructSectionLines(GeoAPI.Geometries.IGeometry _env,int sectlionLineNumber,ref int[] _pLengthBoundaryofEnv)
        {
            List<Geometry> pLines = new List<Geometry>(sectlionLineNumber);
            if (_env.Coordinates[0].Distance(_env.Coordinates[1]) > _env.Coordinates[1].Distance(_env.Coordinates[2]))
            {
                
                Geometry pLine1 = new Geometry(wkbGeometryType.wkbLineString);
                
                pLine1.AddPoint_2D(_env.Coordinates[0].X, _env.Coordinates[0].Y);
                pLine1.AddPoint_2D(_env.Coordinates[1].X, _env.Coordinates[1].Y);
                Geometry newLine1 = GeometrySolve.InsertPointsToLine(pLine1, sectlionLineNumber);

                Geometry pLine2 = new Geometry(wkbGeometryType.wkbLineString);
                pLine2.AddPoint_2D(_env.Coordinates[3].X, _env.Coordinates[3].Y);
                pLine2.AddPoint_2D(_env.Coordinates[2].X, _env.Coordinates[2].Y);
                Geometry newLine2 = GeometrySolve.InsertPointsToLine(pLine2, sectlionLineNumber);

                for (int i = 1; i < newLine1.GetPointCount() - 1; i++)
                {
                    Geometry subLine = new Geometry(wkbGeometryType.wkbLineString);
                    subLine.AddPoint_2D(newLine1.GetX(i), newLine1.GetY(i));
                    subLine.AddPoint_2D(newLine2.GetX(i), newLine2.GetY(i));

                    pLines.Add(subLine);
                }
                
                _pLengthBoundaryofEnv[0] = 0;
                _pLengthBoundaryofEnv[1] = 1;
                _pLengthBoundaryofEnv[2] = 2;
                _pLengthBoundaryofEnv[3] = 3;

            }
            else
            {
                Geometry pLine1 = new Geometry(wkbGeometryType.wkbLineString);
                pLine1.AddPoint_2D(_env.Coordinates[0].X, _env.Coordinates[0].Y);
                pLine1.AddPoint_2D(_env.Coordinates[3].X, _env.Coordinates[3].Y);
                Geometry newLine1 = GeometrySolve.InsertPointsToLine(pLine1, sectlionLineNumber);

                Geometry pLine2 = new Geometry(wkbGeometryType.wkbLineString);
                pLine2.AddPoint_2D(_env.Coordinates[1].X, _env.Coordinates[1].Y);
                pLine2.AddPoint_2D(_env.Coordinates[2].X, _env.Coordinates[2].Y);
                Geometry newLine2 = GeometrySolve.InsertPointsToLine(pLine2, sectlionLineNumber);

                for (int i = 1; i < newLine1.GetPointCount() - 1; i++)
                {
                    Geometry subLine = new Geometry(wkbGeometryType.wkbLineString);
                    subLine.AddPoint_2D(newLine1.GetX(i), newLine1.GetY(i));
                    subLine.AddPoint_2D(newLine2.GetX(i), newLine2.GetY(i));
                    pLines.Add(subLine);
                }

                
                _pLengthBoundaryofEnv[0] = 0;
                _pLengthBoundaryofEnv[1] = 3;
                _pLengthBoundaryofEnv[2] = 1;
                _pLengthBoundaryofEnv[3] = 2;
            }
            return pLines;
        }


        /// <summary>
        /// 5 Construct cross-section lines within the bounding rectangle along the structural direction
        /// </summary>
        /// <param name="_env">the bounding rectangle</param>
        /// <param name="_distance">distance</param>
        /// <param name="pLengthBoundaryofEnv">record the information of the long side</param>
        /// <returns></returns>
        public List<Geometry> ConstructSectionLines(GeoAPI.Geometries.IGeometry _env, double _distance, ref int[] _pLengthBoundaryofEnv)
        {
            List<Geometry> pLines = new List<Geometry>();
            if (_env.Coordinates[0].Distance(_env.Coordinates[1]) > _env.Coordinates[1].Distance(_env.Coordinates[2]))
            {

                Geometry pLine1 = new Geometry(wkbGeometryType.wkbLineString);
                pLine1.AddPoint_2D(_env.Coordinates[0].X, _env.Coordinates[0].Y);
                pLine1.AddPoint_2D(_env.Coordinates[1].X, _env.Coordinates[1].Y);
                //剖面线数量
                int sectlionLineNumber = (int)(pLine1.Length() / _distance);
                Geometry newLine1 = GeometrySolve.InsertPointsToLine(pLine1, sectlionLineNumber);

                Geometry pLine2 = new Geometry(wkbGeometryType.wkbLineString);
                pLine2.AddPoint_2D(_env.Coordinates[3].X, _env.Coordinates[3].Y);
                pLine2.AddPoint_2D(_env.Coordinates[2].X, _env.Coordinates[2].Y);
                Geometry newLine2 = GeometrySolve.InsertPointsToLine(pLine2, sectlionLineNumber);

                for (int i = 1; i < newLine1.GetPointCount() - 1; i++)
                {
                    Geometry subLine = new Geometry(wkbGeometryType.wkbLineString);
                    subLine.AddPoint_2D(newLine1.GetX(i), newLine1.GetY(i));
                    subLine.AddPoint_2D(newLine2.GetX(i), newLine2.GetY(i));

                    pLines.Add(subLine);
                }
                
                _pLengthBoundaryofEnv[0] = 0;
                _pLengthBoundaryofEnv[1] = 1;
                _pLengthBoundaryofEnv[2] = 2;
                _pLengthBoundaryofEnv[3] = 3;

            }
            else
            {
                Geometry pLine1 = new Geometry(wkbGeometryType.wkbLineString);
                pLine1.AddPoint_2D(_env.Coordinates[0].X, _env.Coordinates[0].Y);
                pLine1.AddPoint_2D(_env.Coordinates[3].X, _env.Coordinates[3].Y);
                
                int sectlionLineNumber = (int)(pLine1.Length() / _distance);
                Geometry newLine1 = GeometrySolve.InsertPointsToLine(pLine1, sectlionLineNumber);

                Geometry pLine2 = new Geometry(wkbGeometryType.wkbLineString);
                pLine2.AddPoint_2D(_env.Coordinates[1].X, _env.Coordinates[1].Y);
                pLine2.AddPoint_2D(_env.Coordinates[2].X, _env.Coordinates[2].Y);
                Geometry newLine2 = GeometrySolve.InsertPointsToLine(pLine2, sectlionLineNumber);

                for (int i = 1; i < newLine1.GetPointCount() - 1; i++)
                {
                    Geometry subLine = new Geometry(wkbGeometryType.wkbLineString);
                    subLine.AddPoint_2D(newLine1.GetX(i), newLine1.GetY(i));
                    subLine.AddPoint_2D(newLine2.GetX(i), newLine2.GetY(i));
                    pLines.Add(subLine);
                }

                
                _pLengthBoundaryofEnv[0] = 0;
                _pLengthBoundaryofEnv[1] = 3;
                _pLengthBoundaryofEnv[2] = 1;
                _pLengthBoundaryofEnv[3] = 2;
            }
            return pLines;
        }

        /// <summary>
        /// 5 Construct cross-section lines along the structural direction on the line segment
        /// </summary>
        /// <param name="_pline">multilines</param>
        /// <param name="_sampledistance">sample distance</param>
        /// <param name="_verticalDistance">vertical distance</param>
        /// <returns></returns>
        public List<Geometry> ConstructSectionLines(Geometry _pline, double _sampledistance, double _verticalDistance)
        {
            
            List<Geometry> pSectionLines = new List<Geometry>();

            for(int i=0;i< _pline.GetPointCount()-1;i++)
            {
                Geometry pLine1 = new Geometry(wkbGeometryType.wkbLineString);
                pLine1.AddPoint_2D(_pline.GetX(i), _pline.GetY(i));
                pLine1.AddPoint_2D(_pline.GetX(i+1), _pline.GetY(i+1));

                int sectlionLineNumber = (int)(pLine1.Length() / _sampledistance);
                Geometry newLine1 = GeometrySolve.InsertPointsToLine(pLine1, sectlionLineNumber);

                
                for(int j=1;j< newLine1.GetPointCount()-1;j++)
                {

                    double Vx = _pline.GetX(i) - _pline.GetX(i + 1);
                    double Vy = _pline.GetY(i) - _pline.GetY(i + 1);

                    double magnitude = Math.Sqrt(Vx * Vx + Vy * Vy);

                    
                    double B1x = newLine1.GetX(j) + (_verticalDistance * Vy) / magnitude;
                    double B1y = newLine1.GetY(j) - (_verticalDistance * Vx) / magnitude;

                    Vx = _pline.GetX(i + 1) - _pline.GetX(i);
                    Vy = _pline.GetY(i + 1) - _pline.GetY(i);

                    
                    double B2x = newLine1.GetX(j) + (_verticalDistance * Vy) / magnitude;
                    double B2y = newLine1.GetY(j) - (_verticalDistance * Vx) / magnitude;

                    
                    Geometry pSectionLine = new Geometry(wkbGeometryType.wkbLineString);
                    pSectionLine.AddPoint_2D(B1x, B1y);
                    pSectionLine.AddPoint_2D(B2x, B2y);

                    pSectionLines.Add(pSectionLine);
                }



                
                if((i+2)<= (_pline.GetPointCount() - 1))
                {
                    
                    
                    double x1 = _pline.GetX(i);
                    double y1 = _pline.GetY(i);

                    
                    double x2 = _pline.GetX(i + 1);
                    double y2 = _pline.GetY(i + 1);

                    
                    double x3 = _pline.GetX(i + 2);
                    double y3 = _pline.GetY(i + 2);

                    
                    Vector2 v21 = new Vector2((float)(x1-x2), (float)(y1-y2));
                    Vector2 v32 = new Vector2((float)(x3 - x2), (float)(y3 - y2));
                    
                    Vector2 NorV21 = Vector2.Normalize(v21);
                    Vector2 NorV32 = Vector2.Normalize(v32);

                    Vector2 bisVector = new Vector2((NorV21.X+ NorV32.X) /2, (NorV21.Y + NorV32.Y) / 2);

                    
                    double[] endpoint1 = new double[2];
                    CalculateEndPoint(x2, y2, bisVector.X, bisVector.Y, _verticalDistance, ref endpoint1);


                    
                    v21 = new Vector2((float)(x2 -x1), (float)(y2 -y1));
                     v32 = new Vector2((float)(x2-x3 ), (float)(y2-y3 ));
                    
                     NorV21 = Vector2.Normalize(v21);
                     NorV32 = Vector2.Normalize(v32);
                    
                    bisVector = new Vector2((NorV21.X + NorV32.X) / 2,(NorV21.Y + NorV32.Y) / 2);

                    
                    double[] endpoint2 = new double[2];
                    CalculateEndPoint(x2, y2, bisVector.X, bisVector.Y, _verticalDistance, ref endpoint2);

                    
                    double Vx = _pline.GetX(i + 1) - _pline.GetX(i);
                    double Vy = _pline.GetY(i + 1) - _pline.GetY(i);

                    double Px = pSectionLines[pSectionLines.Count - 1].GetX(0) - _pline.GetX(i);
                    double Py = pSectionLines[pSectionLines.Count - 1].GetY(0) - _pline.GetY(i);

                    double Mx = endpoint1[0] - _pline.GetX(i);
                    double My = endpoint1[1] - _pline.GetY(i);

                    if((Vx* Py- Vy* Px)*(Vx * My- Vy * Mx)>0)
                    {
                        
                        Geometry pSectionLine = new Geometry(wkbGeometryType.wkbLineString);
                        pSectionLine.AddPoint_2D(endpoint1[0], endpoint1[1]);
                        pSectionLine.AddPoint_2D(endpoint2[0], endpoint2[1]);

                        pSectionLines.Add(pSectionLine);
                    }
                    else
                    {
                        
                        Geometry pSectionLine = new Geometry(wkbGeometryType.wkbLineString);
                        pSectionLine.AddPoint_2D(endpoint2[0], endpoint2[1]);
                        pSectionLine.AddPoint_2D(endpoint1[0], endpoint1[1]);

                        pSectionLines.Add(pSectionLine);
                    }

                    
                }
             
            }

            return pSectionLines;
        }

        #endregion

        #region 6.Construct 3D stratigraphic paleo-boundaries based on cross-section lines
        /// <summary>
        /// 6.2 Calculate the stratigraphic Bezier curve corresponding to each cross-section line
        /// </summary>
        /// <param name="_pStratas">Each cross-section line and the stratigraphy it intersects</param>
        /// <returns></returns>
        public Dictionary<Geometry, List<GeoBezier>> GetBeziersOnStrata(Dictionary<Geometry, List<FoldStratum>> _pStratas, NetTopologySuite.Index.Strtree.STRtree<OccurrencePoint> _rtree,
            NetTopologySuite.Geometries.GeometryFactory _geometryFactory, DEMRaster _rasterd)
        {
            Dictionary<Geometry, List<GeoBezier>> pBeziersOnStrata = new Dictionary<Geometry, List<GeoBezier>>();
            
            List<Geometry> pOccurencePoints = new List<Geometry>();

            List<Geometry> pGeoLines = _pStratas.Keys.ToList<Geometry>();
            for (int i = 0; i < pGeoLines.Count; i++)
            {
               
                Geometry distancePoint = new Geometry(wkbGeometryType.wkbPoint);

                List<GeoBezier> strataGeoBeziers = new List<GeoBezier>(_pStratas[pGeoLines[i]].Count);

                for (int j = 0; j < _pStratas[pGeoLines[i]].Count; j++)
                {
                    Geometry pBoundaryPolyline = _pStratas[pGeoLines[i]][j].SPolygon.GeoPolygon.GetBoundary();
                    //pBoundaryPolyline.ExportSimpleGeometryToShapfile(this.filesavepath, "boundary");

                   
                    if (pBoundaryPolyline.GetGeometryType() == wkbGeometryType.wkbMultiLineString)
                    {
                        pBoundaryPolyline = pBoundaryPolyline.GetGeometryRef(0);
                        //pBoundaryPolyline.ExportSimpleGeometryToShapfile(this.filesavepath, "boundary1");
                    }

                    
                    if (pGeoLines[i].Intersect(pBoundaryPolyline))
                    {
                        
                        Geometry pIntersectPoints = pGeoLines[i].Intersection(pBoundaryPolyline);

                        if (pIntersectPoints.GetGeometryCount() == 2)
                        {
                            
                            Geometry pt1 = pIntersectPoints.GetGeometryRef(0);
                            Geometry pt2 = pIntersectPoints.GetGeometryRef(1);

                            
                            OccurrencePoint pt1_Occurence = FindNearestOccurrence(pt1, _rtree, _geometryFactory, _rasterd);
                            OccurrencePoint pt2_Occurence = FindNearestOccurrence(pt2, _rtree, _geometryFactory, _rasterd);
                            if (pt1_Occurence.dipAngle == pt2_Occurence.dipAngle && pt1_Occurence.tendency == pt2_Occurence.tendency
                                && pt1_Occurence.strike == pt2_Occurence.strike)
                            {
                                Console.WriteLine("There is an error as the two points on the same cross-section line have consistent structural attitudes. Please check!");
                                Console.WriteLine("The issue is with the cross-section line:" + i + "stratum:" + _pStratas[pGeoLines[i]][j].SCode);
                            }

                            
                            
                            Geometry pLineStartPoint = new Geometry(wkbGeometryType.wkbPoint);
                            pLineStartPoint.AddPoint_2D(pGeoLines[i].GetX(0), pGeoLines[i].GetY(0));
                            GeoBezier _pBeGeo = null;
                            if (pt1.Distance(pLineStartPoint) < pt2.Distance(pLineStartPoint))
                            {
                                
                                if (j == 0)
                                    distancePoint.AddPoint_2D(pt1.GetX(0), pt1.GetY(0));

                                float _movedistance = (float)pt1.Distance(distancePoint);

                                _pBeGeo = new GeoBezier(_pStratas[pGeoLines[i]][j].SCode + i, pt1_Occurence, pt2_Occurence, _movedistance, this.filesavepath);

                                pOccurencePoints.Add(pt1);
                                pOccurencePoints.Add(pt2);
                            }
                            else
                            {
                                
                                if (j == 0)
                                    distancePoint.AddPoint_2D(pt2.GetX(0), pt2.GetY(0));

                                float _movedistance = (float)pt2.Distance(distancePoint);

                                _pBeGeo = new GeoBezier(_pStratas[pGeoLines[i]][j].SCode + i, pt2_Occurence, pt1_Occurence, _movedistance, this.filesavepath);
                                pOccurencePoints.Add(pt2);
                                pOccurencePoints.Add(pt1);
                            }
                            _pStratas[pGeoLines[i]][j].beziers.Add(_pBeGeo);
                            _pStratas[pGeoLines[i]][j].pSectionLines.Add(pGeoLines[i]);

                            
                            strataGeoBeziers.Add(_pBeGeo);
                        }
                    }
                }
                pBeziersOnStrata.Add(pGeoLines[i], strataGeoBeziers);
            }
            pOccurencePoints.ExportGeometryToShapfile(this.filesavepath, "OCCURECE");
            return pBeziersOnStrata;
        }


        /// <summary>
        /// Find the nearest structural attitude point to a given geometric point, 
        /// and assign the dip, strike, and dip direction of that structural attitude point to the geometric point
        /// </summary>
        /// <param name="pPoint">geometric point</param>
        /// <param name="rtree"></param>
        /// <param name="geometryFactory"></param>
        /// <returns></returns>
        public OccurrencePoint FindNearestOccurrence(Geometry pPoint,
            NetTopologySuite.Index.Strtree.STRtree<OccurrencePoint> rtree,
            NetTopologySuite.Geometries.GeometryFactory geometryFactory, DEMRaster _raster)
        {
            var gt = geometryFactory.CreatePoint(new Coordinate(pPoint.GetX(0), pPoint.GetY(0)));
            OccurrencePoint occpt = new OccurrencePoint(wkbGeometryType.wkbPoint);
            occpt.AddPoint_2D(pPoint.GetX(0), pPoint.GetY(0));
            var distanceMeasure = new PointDistance();
            OccurrencePoint nearestPoint1 = rtree.NearestNeighbour(gt.EnvelopeInternal, occpt, distanceMeasure);

            OccurrencePoint newoccpt = new OccurrencePoint(wkbGeometryType.wkbPoint);
            newoccpt.AddPoint(pPoint.GetX(0), pPoint.GetY(0), _raster.GetElevation(pPoint.GetX(0), pPoint.GetY(0)));
            newoccpt.dipAngle = nearestPoint1.dipAngle;
            newoccpt.tendency = nearestPoint1.tendency;
            newoccpt.strike = nearestPoint1.strike;

            return newoccpt;
        }


        #endregion

        #region 7 Functions for optimizing stratigraphic paleo-boundaries.

        /// <summary>
        /// 7.1 Optimize the unreasonable situations where the Bezier curves of the outer and inner stratigraphy intersect or 
        /// have incorrect order due to elevation or attitude differences
        /// </summary>
        /// <param name="pStratas">All strata</param>
        ///<param name="n">scale</param>
        public void OptimizeBeziers(ref Dictionary<Geometry, List<FoldStratum>> pStratas,int n)
        {
            int countK = 0;
            foreach (var vg in pStratas.Keys)
            {
                countK++;

                if (pStratas[vg].Count > 1)
                {
                    List<FoldStratum> pFoldStrata = pStratas[vg];
                    for (int i = 0; i < pFoldStrata.Count - 1; i++)
                    {
                        FoldStratum pInStratum = pFoldStrata[i];
                        FoldStratum pOutStratum = pFoldStrata[i + 1];

                        int pInIndex = pInStratum.pSectionLines.FindIndex(x => x.Equals(vg));
                        int pOutIndex = pOutStratum.pSectionLines.FindIndex(x => x.Equals(vg));

                        //Calculate stratigraphic thickness on the left side
                        double _distanceLeft = pOutStratum.beziers[pOutIndex]._geoLeftOccurrence.Distance(pInStratum.beziers[pInIndex]._geoLeftOccurrence);
                        Vector2 normalLeft = pInStratum.beziers[pInIndex]._geoBezier.Normal(0.0f);
                        double angleLeft = GeoBezier.CalAzimuthbytwoPoints(0.0, 0.0, normalLeft.X, normalLeft.Y,true);
                        double thicknessLeft = _distanceLeft * Math.Cos(Math.Abs(angleLeft - (3 * Math.PI / 2)));

                        //Calculate stratigraphic thickness on the right side
                        double _distanceRight = pOutStratum.beziers[pOutIndex]._geoRightOccurrence.Distance(pInStratum.beziers[pInIndex]._geoRightOccurrence);
                        Vector2 normalRight = pOutStratum.beziers[pOutIndex]._geoBezier.Normal(1.0f);
                        double angleRight = GeoBezier.CalAzimuthbytwoPoints(0.0, 0.0, normalRight.X, normalRight.Y,true);
                        double thicknessRight = _distanceRight * Math.Cos(Math.Abs(angleRight - (Math.PI / 2)));

                        //calculate the average thickness of the stratum
                        float _distance = (float)((thicknessLeft + thicknessRight) / n);

                        GeoBezier pInBezier = pInStratum.beziers[pInIndex];
                        GeoBezier pOutBezier = pOutStratum.beziers[pOutIndex];

                        //Three scenarios requiring modification of the Bezier curve for the external stratigraphy
                        //（1）Two Bezier curves intersect within the 2D cross-section
                        //（2）Although they do not intersect, the highest point of the outer stratigraphy’s Bezier curve is lower than the highest point of the inner stratigraphy’s Bezier curve
                        //（3）Even if neither of the above conditions occurs, a vertical line extended downward from the highest point of the inner stratigraphy’s Bezier curve intersects with the outer stratigraphy’s Bezier curve

                        //Obtain the vertical line at the point where the tangent of the inner stratigraphy is zero, to determine if the inner stratigraphy’s Bezier curve is positioned above the outer stratigraphy’s Bezier curve
                        Geometry pVerticalLine = new Geometry(wkbGeometryType.wkbLineString);
                        pVerticalLine.AddPoint_2D(pInStratum.beziers[pInIndex]._zeroTangentPoint.X, pInStratum.beziers[pInIndex]._zeroTangentPoint.Y);
                        pVerticalLine.AddPoint_2D(pInStratum.beziers[pInIndex]._zeroTangentPoint.X, 0.0);

                        //move interval distance (m)
                        float internla = 0.0f;

                        //Set a flag to indicate that only Bezier curves exhibiting any of the above three conditions require processing
                        int isNeedModify = -1;                      

                        while (pInStratum.beziers[pInIndex]._2DBezierGeometry.Intersects(pOutStratum.beziers[pOutIndex]._2DBezierGeometry) ||
                            pOutStratum.beziers[pOutIndex]._zeroTangentPoint.Y < pInStratum.beziers[pInIndex]._zeroTangentPoint.Y ||
                            pVerticalLine.Intersects(pOutStratum.beziers[pOutIndex]._2DBezierGeometry)
                            )
                        {
                            if (pInStratum.beziers[pInIndex]._2DBezierGeometry.Intersects(pOutStratum.beziers[pOutIndex]._2DBezierGeometry))
                                Console.WriteLine("Scenario 1: Intersect" + pInStratum.SCode + " " + pOutStratum.SCode);
                            else if (pOutStratum.beziers[pOutIndex]._zeroTangentPoint.Y < pInStratum.beziers[pInIndex]._zeroTangentPoint.Y)
                                Console.WriteLine("Scenario 2：The tangent point of the outer stratigraphy is located below the inner stratigraphy" + pInStratum.SCode + " " + pOutStratum.SCode);
                            else
                                Console.WriteLine("Scenario 3：Located below" + pInStratum.SCode + " " + pOutStratum.SCode);

                            Console.WriteLine(countK+" "+ pInStratum.SCode+" " + pInStratum.beziers[pInIndex]._zeroTangentPoint.X+" "+pInStratum.beziers[pInIndex]._zeroTangentPoint.Y);
                            Console.WriteLine(countK+" "+ pOutStratum.SCode+" " + pOutStratum.beziers[pOutIndex]._zeroTangentPoint.X+" "+ pOutStratum.beziers[pOutIndex]._zeroTangentPoint.Y);

                            isNeedModify = 0;
                            pOutStratum.beziers[pOutIndex] = ChangeBezierByControlPoint(pInStratum.beziers[pInIndex], pOutStratum.beziers[pOutIndex], internla++, pOutStratum.beziers[pOutIndex]._geoMoveDistance, 0.01f);

                        }

                        if (isNeedModify == 0)
                        {

                            while (Math.Abs(pOutStratum.beziers[pOutIndex]._zeroTangentPoint.Y-pInStratum.beziers[pInIndex]._zeroTangentPoint.Y - _distance) > 10)
                            {
                                double ggg =GetVerticalDistanceOfTwoTangentPoints(pInStratum.beziers[pInIndex]._zeroTangentPoint, pOutStratum.beziers[pOutIndex]._zeroTangentPoint);
                                double ssss = Math.Abs(GetVerticalDistanceOfTwoTangentPoints(pInStratum.beziers[pInIndex]._zeroTangentPoint, pOutStratum.beziers[pOutIndex]._zeroTangentPoint) - _distance);
                                Console.WriteLine(ssss+" "+ggg);
                                pOutStratum.beziers[pOutIndex] = ChangeBezierByControlPoint(pInStratum.beziers[pInIndex], pOutStratum.beziers[pOutIndex], internla++, pOutStratum.beziers[pOutIndex]._geoMoveDistance, 0.01f);
                            }
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Obtain a new Bezier curve by adjusting the control points
        /// </summary>
        /// <param name="_pBezierGeo1">The GeoBezier of the inner strata</param>
        /// <param name="_pBezierGeo2">The GeoBezier of the outer strata</param>
        /// <param name="_distance">The distance that needs to be moved along the extension line</param>
        /// <param name="_movedistance">The distance the Bezier curve needs to be shifted to the left</param>
        /// <param name="_interval">Specify the interval between points when outputting 2D or 3D data</param>
        /// <returns></returns>
        public GeoBezier ChangeBezierByControlPoint(GeoBezier _pBezierGeo1, GeoBezier _pBezierGeo2, float _distance,
            float _movedistance, float _interval)
        {
            Vector2 controlPoint = GetExtendPointByTwopoints(_pBezierGeo1._zeroTangentPoint, _pBezierGeo1._geointersectPoint, _distance);

            
            GeoBezier _pBezierGeo = new GeoBezier(_pBezierGeo2._geoBezierName, _pBezierGeo2._geoLeftOccurrence, _pBezierGeo2._geoRightOccurrence, _pBezierGeo2._geoLeftPoint, _pBezierGeo2._geoRightPoint, controlPoint, _movedistance, this.filesavepath);
            _pBezierGeo.DrawSectionWithThreeControlPoints(_interval);

            _pBezierGeo.CalucateZeroTangent(0.01f, 0.01f);
            return _pBezierGeo;
        }

        /// <summary>
        /// Calculate the vertical distance between two tangent points
        /// </summary>
        /// <param name="_pInTangentPoint">the tangent point of the outer statum</param>
        /// <param name="_pOutTangentPoint">the tangent point of the inner statum</param>
        /// <returns></returns>
        public double GetVerticalDistanceOfTwoTangentPoints(Vector2 _pInTangentPoint,Vector2 _pOutTangentPoint)
        {
            
            double distanceofTangentPoints = Vector2.Distance(_pInTangentPoint, _pOutTangentPoint);
          
            double angleofTangePoints = GeoBezier.CalAzimuthbytwoPoints((double)_pOutTangentPoint.X, (double)_pOutTangentPoint.Y, (double)_pInTangentPoint.X, (double)_pInTangentPoint.Y,true);
            
            if(angleofTangePoints>=Math.PI&&angleofTangePoints<(3*Math.PI/2))
            {
                return distanceofTangentPoints * Math.Cos(angleofTangePoints - Math.PI);
            }
            else if(angleofTangePoints <= Math.PI && angleofTangePoints > ( Math.PI / 2))
            {
                return distanceofTangentPoints * Math.Cos(Math.PI-angleofTangePoints);
            }
            else
            {
                Console.WriteLine("The relative position of the tangent points on the inner and outer Bezier curves is incorrect. Please check！！！");
                return -1.0;
            }
        }

        /// <summary>
        /// Obtain the point on the extended line
        /// </summary>
        /// <param name="_pointA"></param>
        /// <param name="_pointB"></param>
        /// <param name="_distance"></param>
        /// <returns></returns>
        public Vector2 GetExtendPointByTwopoints(Vector2 _pointA, Vector2 _pointB, float _distance)
        {
            float distAB = Vector2.Distance(_pointA, _pointB);

            float unitABx = (_pointB.X - _pointA.X) / distAB;

            float unitABy = (_pointB.Y - _pointA.Y) / distAB;

            float Cx = _pointB.X + unitABx * _distance;
            float Cy = _pointB.Y + unitABy * _distance;

            return new Vector2(Cx, Cy);
        }

        /// <summary>
        /// Create lines of equal inclination, given the tangent vectors at each point on the lower Bezier curve
        /// </summary>
        /// <param name="_pBottomBezier"></param>
        /// <param name="_pUpBezier"></param>
        /// <param name="_stepDistance"></param>
        /// <returns></returns>
        public List<Geometry> CreateDipIsogons(Bezier _pBottomBezier, Bezier _pUpBezier, float _stepSize)
        {
            
            List<Geometry> pDipIsogons = new List<Geometry>();

            for (float i = 0; i <= 1; i += 0.01f)
            {
                
                Vector2 pBottomtangent = _pBottomBezier.Tangent(i);
                Vector2 pBottomPoint = _pBottomBezier.Position(i);

                
                float t = 0.0f; 

                // Loop until a _point with the desired pBottomtangent is found
                while (t <= 1.0f)
                {
                    
                    Vector2 direction = _pUpBezier.Tangent(t);

                   
                    float dotProduct = pBottomtangent.X * direction.X + pBottomtangent.Y * direction.Y;


                    
                    if (Math.Abs((dotProduct / (pBottomtangent.Length() * direction.Length()) - 1.0f)) < 0.001f)
                    {
                        
                        Vector2 pUpPoint = _pUpBezier.Position(t);

                        Geometry pDipIsogon = new Geometry(wkbGeometryType.wkbLineString);
                        pDipIsogon.AddPoint_2D(pBottomPoint.X, pBottomPoint.Y);
                        pDipIsogon.AddPoint_2D(pUpPoint.X, pUpPoint.Y);

                        pDipIsogons.Add(pDipIsogon);
                      
                        break;
                    }

                    
                    t += _stepSize;
                }
            }
            return pDipIsogons;
        }

        #endregion

        #region 8 Morphing插值需要的函数

        /// <summary>
        /// 8.4 Functions required for Morphing interpolation
        /// </summary>
        /// <param name="_rasterD"></param>
        /// <param name="geometryFactory"></param>
        /// <param name="pNewOrderStrata"></param>
        /// <param name="pLinesToModifyDEM"></param>
        /// <param name="foldtype">Fold type，0 represents a completely closed fold，1 represents a uncomletely closed fold</param>
        public void MorphingIntepolateTransmitionLines(DEMRaster _rasterD, NetTopologySuite.Geometries.GeometryFactory geometryFactory, ref List<FoldStratum> pNewOrderStrata,ref List<Geometry> pLinesToModifyDEM,int foldtype)
        {
            
            Dictionary<string, List<Geometry>> pTransmitionLineDict = new Dictionary<string, List<Geometry>>(pNewOrderStrata.Count);//所有过渡曲线

            for (int i = 0; i < pNewOrderStrata.Count; i++)
            {
                //8.4.1 Obtain the outer boundary of the stratigraphy—there may be inner and outer rings or holes within the surface, so directly select the outer boundary
                Geometry pBoundaryPolyline = pNewOrderStrata[i].SPolygon.GeoPolygon.GetBoundary();
                wkbGeometryType ps = pBoundaryPolyline.GetGeometryType();
                if (pBoundaryPolyline.GetGeometryType() == wkbGeometryType.wkbMultiLineString)
                {
                    pBoundaryPolyline = pBoundaryPolyline.GetGeometryRef(0);
                }


                //8.4.2 Extract the points of the outer boundary
                Coordinate[] pOutwardLinePoints = new Coordinate[pBoundaryPolyline.GetPointCount()];
                for (int j = 0; j < pBoundaryPolyline.GetPointCount(); j++)
                {
                    pOutwardLinePoints[j] = new Coordinate(pBoundaryPolyline.GetX(j), pBoundaryPolyline.GetY(j));
                }

                var outBoundaryLine = geometryFactory.CreateLineString(pOutwardLinePoints);

                //8.4.3 Obtain the truncation point on one side
                Coordinate[] pSpiltPointsOneEage = new Coordinate[pNewOrderStrata[i].beziers.Count + 2];
                pSpiltPointsOneEage[0] = new Coordinate(pNewOrderStrata[i].pOneNode.GetX(0), pNewOrderStrata[i].pOneNode.GetY(0));
                for (int j = 0; j < pNewOrderStrata[i].beziers.Count; j++)
                {
                    pSpiltPointsOneEage[j + 1] = new Coordinate(pNewOrderStrata[i].beziers[j]._geoLeftOccurrence.GetX(0), pNewOrderStrata[i].beziers[j]._geoLeftOccurrence.GetY(0));
                }
                pSpiltPointsOneEage[pSpiltPointsOneEage.Length - 1] = new Coordinate(pNewOrderStrata[i].pAnotherNode.GetX(0), pNewOrderStrata[i].pAnotherNode.GetY(0));

              

                ConvertToGeometryLine(pSpiltPointsOneEage.ToList(), false).ExportSimpleGeometryToShapfile(this.filesavepath, "oneeagesplitpoint"+i);

                //8.4.4 Obtain the truncation point on the other side
                Coordinate[] pSpiltPointsAnotherEage = new Coordinate[pNewOrderStrata[i].beziers.Count + 2];
                pSpiltPointsAnotherEage[0] = new Coordinate(pNewOrderStrata[i].pOneNode.GetX(0), pNewOrderStrata[i].pOneNode.GetY(0));
                for (int j = 0; j < pNewOrderStrata[i].beziers.Count; j++)
                {
                    pSpiltPointsAnotherEage[j + 1] = new Coordinate(pNewOrderStrata[i].beziers[j]._geoRightOccurrence.GetX(0), pNewOrderStrata[i].beziers[j]._geoRightOccurrence.GetY(0));
                }
                pSpiltPointsAnotherEage[pSpiltPointsAnotherEage.Length - 1] = new Coordinate(pNewOrderStrata[i].pAnotherNode.GetX(0), pNewOrderStrata[i].pAnotherNode.GetY(0));
                ConvertToGeometryLine(pSpiltPointsAnotherEage.ToList(), false).ExportSimpleGeometryToShapfile(this.filesavepath, "anothereagesplitpoint" + i);

                //8.4.5 Identify a new starting point to modify the order of points in the sequential geometric boundary
                List<Coordinate> pNeworderOutwardLinePoints = new List<Coordinate>(pBoundaryPolyline.GetPointCount());
                int flag = -1;
                for (int j = 0; j < pOutwardLinePoints.Length; j++)
                {
                    if ((pOutwardLinePoints[j].X == pSpiltPointsOneEage[0].X) && (pOutwardLinePoints[j].Y == pSpiltPointsOneEage[0].Y))
                        flag = j;
                }
                for (int j = flag; j < pOutwardLinePoints.Length; j++)
                {
                    pNeworderOutwardLinePoints.Add(pOutwardLinePoints[j]);
                }
                for (int j = 1; j <= flag; j++)
                {
                    pNeworderOutwardLinePoints.Add(pOutwardLinePoints[j]);
                }

                ConvertToGeometryLine(pNeworderOutwardLinePoints, false).ExportSimpleGeometryToShapfile(this.filesavepath, "TestSpline1");

                //8.4.6 Insert the intersection points of the cross-section line and the outer boundary line into the set
                InterPointsToOrderPoints(pNewOrderStrata[i], pNewOrderStrata[i].pSectionLines, ref pNeworderOutwardLinePoints);

                //8.4.7 Obtain the new segments and convert them to Geometry (OGR) for output
                Dictionary<int, List<Coordinate>> pOneEageSegements = GetOrderSegements(pSpiltPointsOneEage, pNeworderOutwardLinePoints);
                Dictionary<int, List<Coordinate>> pAnotherEageSegements = GetOrderSegements(pSpiltPointsAnotherEage, pNeworderOutwardLinePoints);

                //8.4.8 Obtain the maximum value of each pair of opposite sides
                List<int> pMaxNodesNumber = new List<int>(pOneEageSegements.Count);
                foreach (var num in pOneEageSegements.Keys)
                {
                    int maxcount = Math.Max(pOneEageSegements[num].Count, pAnotherEageSegements[num].Count);
                    pMaxNodesNumber.Add(maxcount - 2);

                    if (pOneEageSegements[num].Count > pAnotherEageSegements[num].Count)
                    {
                        List<Coordinate> pNeedAddNodes = pAnotherEageSegements[num];
                        AddPointsToLine(ref pNeedAddNodes, pOneEageSegements[num].Count - pAnotherEageSegements[num].Count);
                    }
                    else if (pOneEageSegements[num].Count < pAnotherEageSegements[num].Count)
                    {
                        List<Coordinate> pNeedAddNodes = pOneEageSegements[num];
                        AddPointsToLine(ref pNeedAddNodes, pAnotherEageSegements[num].Count - pOneEageSegements[num].Count);
                    }
                }

                //8.4.9 Perform Lagrange interpolation on the three-dimensional control points
                List<Geometry> pAllInterpolates = new List<Geometry>();
                Dictionary<int, List<Coordinate>> pMiddleSegments = new Dictionary<int, List<Coordinate>>();
                pMiddleSegments = LagrangeInterpolate1(pNewOrderStrata[i].pHingePoints, pMaxNodesNumber, ref pAllInterpolates);
                Geometry pHingeLine = new Geometry(wkbGeometryType.wkbLineString);
                foreach (var vt in pAllInterpolates)
                    pHingeLine.AddPoint(vt.GetX(0), vt.GetY(0), vt.GetZ(0));
                pHingeLine.ExportSimpleGeometryToShapfile(this.filesavepath, "higelines" + i);

                //8.4.10 Perform morphing interpolation

                //(1) Starting point   
                int firstNodesCount = pNewOrderStrata[i].beziers[0]._3DBezierGeometry.GetPointCount();
                Geometry pStartLine = new Geometry(wkbGeometryType.wkbLineString);
                for (int k = 0; k < firstNodesCount; k++)
                {
                    pStartLine.AddPoint(pNewOrderStrata[i].p3DOneNode.GetX(0), pNewOrderStrata[i].p3DOneNode.GetY(0), pNewOrderStrata[i].p3DOneNode.GetZ(0));
                }

                //(2) ending point
                Geometry pEndLine = new Geometry(wkbGeometryType.wkbLineString);
                for (int k = 0; k < firstNodesCount; k++)
                {
                    pEndLine.AddPoint(pNewOrderStrata[i].p3DAnotherNode.GetX(0), pNewOrderStrata[i].p3DAnotherNode.GetY(0), pNewOrderStrata[i].p3DAnotherNode.GetZ(0));
                }

                //(3) Interpolate transition paleo-boundaries
                List<Geometry> pTransitionLines = new List<Geometry>();

                List<Geometry> pTransitlionLinesNoBoundarypPoints = new List<Geometry>();

                for (int k = 0; k < pOneEageSegements.Count; k++)
                {
                    Geometry p3DOneEageSegement = ConvertToGeometryLine(pOneEageSegements[k], _rasterD);
                    Geometry p3DMiddleEageSegement = ConvertToGeometryLine(pMiddleSegments[k], true);
                    Geometry p3DAnotherEageSegement = ConvertToGeometryLine(pAnotherEageSegements[k], _rasterD);
                  
                    if (k == 0)
                    {
                        List<Geometry> pTransLines = MorphingSolve.Create3DTransitionPaleoBoundariesByMorphing(pStartLine, pNewOrderStrata[i].beziers[0]._3DBezierGeometry,
                            p3DOneEageSegement, p3DMiddleEageSegement, p3DAnotherEageSegement,_rasterD);
                      
                        if (i == 0)
                        {                           
                            pLinesToModifyDEM.AddRange(pTransLines);
                        }

                        pTransitionLines.AddRange(pTransLines);
                       
                    }
                    else if (k == pOneEageSegements.Count - 1)
                    {
                        List<Geometry> pTransLines = MorphingSolve.Create3DTransitionPaleoBoundariesByMorphing(pNewOrderStrata[i].beziers[pNewOrderStrata[i].beziers.Count - 1]._3DBezierGeometry, pEndLine,
                           p3DOneEageSegement, p3DMiddleEageSegement, p3DAnotherEageSegement,_rasterD);

                                           
                        if (i == 0)
                        {
                          
                                pLinesToModifyDEM.AddRange(pTransLines);
                        }

                        pTransitionLines.AddRange(pTransLines);
                    }
                    else
                    {
                        List<Geometry> pTransLines = MorphingSolve.Create3DTransitionPaleoBoundariesByMorphing(pNewOrderStrata[i].beziers[k - 1]._3DBezierGeometry, pNewOrderStrata[i].beziers[k]._3DBezierGeometry,
                           p3DOneEageSegement, p3DMiddleEageSegement, p3DAnotherEageSegement, _rasterD);

                        
                        if (i == 0)
                            pLinesToModifyDEM.AddRange(pTransLines);


                        pTransitionLines.AddRange(pTransLines);
                        pTransitlionLinesNoBoundarypPoints.AddRange(pTransLines);
                    }
                }

                pNewOrderStrata[i].transitionLines = pTransitlionLinesNoBoundarypPoints;
                pTransitionLines.ExportGeometryToShapfile(this.filesavepath, "TransmitLines" + i);
            }
        }
        /// <summary>
        /// Obtain the intersection points of the cross-section line and the stratigraphic outer boundary segments
        /// In step 6, there is an operation to find an intersection point, but it is unclear between which two points it is located
        /// To facilitate accurate insertion into the existing boundary points, it is necessary to retrieve the points again
        /// </summary>
        /// <param name="_sectionLines"></param>
        /// <param name="_pBoundaryPoints"></param>
        public void InterPointsToOrderPoints(FoldStratum _stratum, List<Geometry> _sectionLines, ref List<Coordinate> _pBoundaryPoints)
        {
            Dictionary<int, Coordinate> flags = new Dictionary<int, Coordinate>(_sectionLines.Count * 2);
            for (int i = 0; i < _sectionLines.Count; i++)
            {
                for (int j = 0; j < _pBoundaryPoints.Count - 1; j++)
                {
                    Geometry pline = new Geometry(wkbGeometryType.wkbLineString);
                    pline.AddPoint_2D(_pBoundaryPoints[j].X, _pBoundaryPoints[j].Y);
                    pline.AddPoint_2D(_pBoundaryPoints[j + 1].X, _pBoundaryPoints[j + 1].Y);

                    if (pline.Intersects(_sectionLines[i]))
                    {
                        Geometry pIntersectPoint = pline.Intersection(_sectionLines[i]);
                        Coordinate PS = new Coordinate(pIntersectPoint.GetX(0), pIntersectPoint.GetY(0));

                        flags.Add(j, PS);
                    }

                }
            }

            List<Coordinate> pIntes = flags.Values.ToList();

            int count = 0;
            for (int i = 0; i < pIntes.Count; i++)
            {
                for (int j = 0; j < _stratum.beziers.Count; j++)
                {
                    Coordinate pLeft = new Coordinate(_stratum.beziers[j]._geoLeftOccurrence.GetX(0), _stratum.beziers[j]._geoLeftOccurrence.GetY(0));
                    Coordinate pRight = new Coordinate(_stratum.beziers[j]._geoRightOccurrence.GetX(0), _stratum.beziers[j]._geoRightOccurrence.GetY(0));


                    if (pIntes[i].X == pLeft.X && pIntes[i].Y == pLeft.Y)
                    {
                        count++;
                    }
                    if (pIntes[i].X == pRight.X && pIntes[i].Y == pRight.Y)
                    {
                        count++;
                    }
                }
            }

            List<int> flagint = flags.Keys.ToList();
            flagint.Sort((a, b) => a.CompareTo(b));

            for (int i = 0; i < flagint.Count; i++)
            {
                _pBoundaryPoints.Insert(flagint[i] + i + 1, flags[flagint[i]]);
            }
        }



        /// <summary>
        /// Perform smooth interpolation on the hinge line
        /// </summary>
        /// <param name="pPoints">points on the hinge line</param>
        /// <param name="numPoints">The number of points to be inserted</param>
        /// <param name="interpolatedPoints">Output the interpolated points</param>
        /// <returns></returns>
        public Dictionary<int, List<Coordinate>> LagrangeInterpolate1(List<Geometry> pPoints, List<int> numPoints, ref List<Geometry> interpolatedPoints)
        {
            
            Dictionary<int, List<Coordinate>> pInterpolatePoints = new Dictionary<int, List<Coordinate>>();

            double[] xValues = new double[pPoints.Count];
            double[] yValues = new double[pPoints.Count];          
            double[] zValues = new double[pPoints.Count];

            
            double[] sumDistanceValues = new double[pPoints.Count];

            for (int i = 0; i < pPoints.Count; i++)
            {
                xValues[i] = pPoints[i].GetX(0);
                yValues[i] = pPoints[i].GetY(0);

                if(i==0)               
                    sumDistanceValues[0] = 0.0;
                else
                {
                    double _distance = pPoints[i].Distance(pPoints[i-1]);
                    sumDistanceValues[i] = _distance + sumDistanceValues[i - 1];
                }
                
                zValues[i] = pPoints[i].GetZ(0);
            }

            
            var interpolation = CubicSpline.InterpolateNaturalSorted(sumDistanceValues, zValues);

           
            double sumdistance = 0.0;

            for (int j = 0; j < pPoints.Count - 1; j++)
            {
                
                interpolatedPoints.Add(pPoints[j]);

                
                List<Coordinate> subinterpolatedPoints = new List<Coordinate>(numPoints[j]);

                double start_X = pPoints[j].GetX(0);
                double end_X = pPoints[j + 1].GetX(0);

                double start_Y = pPoints[j].GetY(0);
                double end_Y = pPoints[j + 1].GetY(0);

                Coordinate pStartNode = new Coordinate(start_X, start_Y);
                Coordinate pEndNode = new Coordinate(end_X, end_Y);

                subinterpolatedPoints.Add(new Coordinate(start_X, start_Y, pPoints[j].GetZ(0)));

                
                Coordinate[] pInterPoints = new Coordinate[numPoints[j]];
                for (int i = 1; i <= numPoints[j]; i++)
                {

                    double dataX = start_X + ((end_X - start_X) * i) / (numPoints[j] + 1);
                    double datay = start_Y + ((end_Y - start_Y) * i) / (numPoints[j] + 1);

                    pInterPoints[i - 1] = new Coordinate(dataX, datay);
                }

                
                for (int i = 0; i < pInterPoints.Length; i++)
                {
                    if (i == 0&&j==0)
                    {
                        sumdistance = pInterPoints[i].Distance(pStartNode);
                    }
                    else if(i == 0 && j > 0)
                    {
                        sumdistance =sumdistance+ pInterPoints[i].Distance(pStartNode);
                    }
                    else if(i== pInterPoints.Length-1)
                    {
                        sumdistance = sumdistance + pInterPoints[i].Distance(pInterPoints[i - 1]);

                        sumdistance = sumdistance + pEndNode.Distance(pInterPoints[i]);
                    }
                    else
                    {
                        sumdistance = sumdistance + pInterPoints[i].Distance(pInterPoints[i - 1]);
                    }

                    
                    double value = interpolation.Interpolate(sumdistance);

                    
                    Geometry pt = new Geometry(wkbGeometryType.wkbPoint);
                    pt.AddPoint(pInterPoints[i].X, pInterPoints[i].Y, value);

                    interpolatedPoints.Add(pt);
                    subinterpolatedPoints.Add(new Coordinate(pInterPoints[i].X, pInterPoints[i].Y, value)); ;
                }

                subinterpolatedPoints.Add(new Coordinate(end_X, end_Y, pPoints[j + 1].GetZ(0)));
                interpolatedPoints.Add(pPoints[j + 1]);
                pInterpolatePoints.Add(j, subinterpolatedPoints);
            }

            return pInterpolatePoints;
        }

        /// <summary>
        /// Obtain the endpoints on both sides in the structural direction for all stratigraphic layers except the outermost one
        /// </summary>
        /// <param name="_pNewOrderStrata"></param>
        public void GetNodesBothStratum(ref FoldStratum vStratum, ref FoldStratum subStratum)
        {
            
            Geometry pBounPolyline = vStratum.SPolygon.GeoPolygon.GetBoundary();

            
            if (pBounPolyline.GetGeometryType() == wkbGeometryType.wkbMultiLineString)
            {
                pBounPolyline = pBounPolyline.GetGeometryRef(0);
            }


            List<Coordinate> pCoordInPolygon = new List<Coordinate>(pBounPolyline.GetPointCount());
            int count = pBounPolyline.GetPointCount();
            for (int i = 0; i < pBounPolyline.GetPointCount(); i++)
            {
                
                pCoordInPolygon.Add(new Coordinate(pBounPolyline.GetX(i), pBounPolyline.GetY(i)));
            }
            Coordinate pLeftNearestPoint = FindNearestPoint(new Coordinate(subStratum.pOneNode.GetX(0), subStratum.pOneNode.GetY(0)), pCoordInPolygon);
            Geometry pNewLeftNode = new Geometry(wkbGeometryType.wkbPoint);
            pNewLeftNode.AddPoint_2D(pLeftNearestPoint.X, pLeftNearestPoint.Y);
            vStratum.pOneNode = pNewLeftNode;

            double s1 = subStratum.pAnotherNode.GetX(0);
            double s2 = subStratum.pAnotherNode.GetY(0);

            Coordinate pRightNearestPoint = FindNearestPoint(new Coordinate(subStratum.pAnotherNode.GetX(0), subStratum.pAnotherNode.GetY(0)), pCoordInPolygon);
            Geometry pNewRightNode = new Geometry(wkbGeometryType.wkbPoint);
            pNewRightNode.AddPoint_2D(pRightNearestPoint.X, pRightNearestPoint.Y);
            vStratum.pAnotherNode = pNewRightNode;

        }

        /// <summary>
        /// Obtain the endpoints on both sides in the structural direction of the outermost stratigraphic layer
        /// </summary>
        /// <param name="_pNewOrderStrata"></param>
        /// <param name="_env"></param>
        /// <param name="_pOrderofWidth"></param>
        public void GetNodesBothStratum(ref FoldStratum subStratum, IGeometry _env, int[] _pOrderofWidth)
        {
            Geometry pBoundaryPolyline = subStratum.SPolygon.GeoPolygon.GetBoundary();
            wkbGeometryType ps = pBoundaryPolyline.GetGeometryType();

            
            if (pBoundaryPolyline.GetGeometryType() == wkbGeometryType.wkbMultiLineString)
            {
                pBoundaryPolyline = pBoundaryPolyline.GetGeometryRef(0);
            }

            
            Geometry pLeftLine = new Geometry(wkbGeometryType.wkbLineString);
            pLeftLine.AddPoint_2D(_env.Coordinates[_pOrderofWidth[0]].X, _env.Coordinates[_pOrderofWidth[0]].Y);
            pLeftLine.AddPoint_2D(_env.Coordinates[_pOrderofWidth[2]].X, _env.Coordinates[_pOrderofWidth[2]].Y);

            
            Geometry pRightLine = new Geometry(wkbGeometryType.wkbLineString);
            pRightLine.AddPoint_2D(_env.Coordinates[_pOrderofWidth[1]].X, _env.Coordinates[_pOrderofWidth[1]].Y);
            pRightLine.AddPoint_2D(_env.Coordinates[_pOrderofWidth[3]].X, _env.Coordinates[_pOrderofWidth[3]].Y);

            if (pLeftLine.Intersect(pBoundaryPolyline))
            {
                Geometry pIntersectPoints = pLeftLine.Intersection(pBoundaryPolyline);
                if (pIntersectPoints.GetGeometryType() == wkbGeometryType.wkbPoint)
                {
                    subStratum.pOneNode = pIntersectPoints;
                }
                else
                {
                    subStratum.pOneNode = FindNearestPointToLine(subStratum.pOutBoundaryPoints, pLeftLine);
                }
            }
            else
            {
                subStratum.pOneNode = FindNearestPointToLine(subStratum.pOutBoundaryPoints, pLeftLine);
            }

            if (pRightLine.Intersects(pBoundaryPolyline))
            {
                Geometry pIntersectPoints = pLeftLine.Intersection(pBoundaryPolyline);
                wkbGeometryType MPS = pIntersectPoints.GetGeometryType();
                if (pIntersectPoints.GetGeometryType() == wkbGeometryType.wkbPoint)
                {
                    subStratum.pAnotherNode = pIntersectPoints;
                }
                // If there is more than one intersection point, use the closest point method
                else
                {
                    subStratum.pAnotherNode = FindNearestPointToLine(subStratum.pOutBoundaryPoints, pRightLine);
                }
            }
            else
            {
                subStratum.pAnotherNode = FindNearestPointToLine(subStratum.pOutBoundaryPoints, pRightLine);
            }
        }

        /// <summary>
        /// Obtain the endpoints on both sides in the structural direction of the outermost stratigraphic layer
        /// </summary>
        /// <param name="subStratum"></param>
        /// <param name="_pFirstSectionLine"></param>
        /// <param name="_extremeNodes"></param>
        public void GetNodesBothStratum(ref FoldStratum subStratum, Geometry _pFirstSectionLine, List<Geometry> _extremeNodes)
        {
            Geometry pBoundaryPolyline = subStratum.SPolygon.GeoPolygon.GetBoundary();
            wkbGeometryType ps = pBoundaryPolyline.GetGeometryType();

            
            if (pBoundaryPolyline.GetGeometryType() == wkbGeometryType.wkbMultiLineString)
            {
                pBoundaryPolyline = pBoundaryPolyline.GetGeometryRef(0);
            }

            List<Geometry> pBoundaryPoints = new List<Geometry>(pBoundaryPolyline.GetPointCount());
            for(int i=0;i< pBoundaryPolyline.GetPointCount(); i++)
            {
                Geometry pt = new Geometry(wkbGeometryType.wkbPoint);
                pt.AddPoint_2D(pBoundaryPolyline.GetX(i), pBoundaryPolyline.GetY(i));

                pBoundaryPoints.Add(pt);
            }

            Geometry pOne = FindNearestPoint(_extremeNodes[0], pBoundaryPoints);
            Geometry pTwo = FindNearestPoint(_extremeNodes[1], pBoundaryPoints);    

            double distance1 = pOne.Distance(_pFirstSectionLine);
            double distance2 = pTwo.Distance(_pFirstSectionLine);


            if (distance1< distance2)
            {
                subStratum.pOneNode = pOne;
                subStratum.pAnotherNode = pTwo;
            }
            else
            {
                subStratum.pOneNode = pTwo;
                subStratum.pAnotherNode = pOne;
            }
        }

        /// <summary>
        /// Add a specific number of points between two points
        /// </summary>
        /// <param name="_pNodes"></param>
        /// <param name="_pointcount"></param>
        public void AddPointsToLine(ref List<Coordinate> _pNodes, int _pointcount)
        {
            int fl = (int)_pNodes.Count / (_pointcount + 1);

            for (int i = 1; i <= _pointcount; i++)
            {
                double x = (_pNodes[fl * i].X + _pNodes[fl * i + 1].X) / 2;
                double y = (_pNodes[fl * i].Y + _pNodes[fl * i + 1].Y) / 2;
                _pNodes.Insert(fl, new Coordinate(x, y));
            }

        }


        /// <summary>
        /// Return the segmented line segments
        /// </summary>
        /// <param name="_oneEage"></param>
        /// <param name="_allNodes"></param>
        /// <returns></returns>
        public Dictionary<int, List<Coordinate>> GetOrderSegements(Coordinate[] _oneEage, List<Coordinate> _allNodes)
        {
            int[] posistions = new int[_oneEage.Length];

            for (int j = 0; j < _oneEage.Length; j++)
            {
                for (int i = 0; i < _allNodes.Count - 1; i++)
                {
                    if (_oneEage[j].X == _allNodes[i].X && _oneEage[j].Y == _allNodes[i].Y)
                    {
                        posistions[j] = i;

                        i = _allNodes.Count;
                    }
                }
            }

            Dictionary<int, List<Coordinate>> pDifferentSegment = new Dictionary<int, List<Coordinate>>(posistions.Length - 1);

            if (posistions[1] < posistions[2])
            {
                for (int i = 0; i < posistions.Length - 1; i++)
                {
                    List<Coordinate> pSegmentPoints = new List<Coordinate>();
                    for (int j = posistions[i]; j <= posistions[i + 1]; j++)
                    {
                        pSegmentPoints.Add(_allNodes[j]);
                    }

                    pDifferentSegment.Add(i, pSegmentPoints);
                }
            }
            else
            {
                posistions[0] = _allNodes.Count - 1;
                for (int i = 0; i < posistions.Length - 1; i++)
                {

                    List<Coordinate> pSegmentPoints = new List<Coordinate>();
                    for (int j = posistions[i]; j >= posistions[i + 1]; j--)
                    {
                        pSegmentPoints.Add(_allNodes[j]);
                    }

                    pDifferentSegment.Add(i, pSegmentPoints);
                }
            }
            return pDifferentSegment;
        }

        #endregion

        #region Dome: 5. Construct cross-section lines

        /// <summary>
        /// Create cross-section lines based on the given interval angle
        /// </summary>
        /// <param name="_centroid">Centroid</param>
        /// <param name="_angleInterval">Interval angle</param>
        /// <param name="_envlopeLine"></param>
        /// <returns></returns>
        public List<NetTopologySuite.Geometries.LineString> CreateSectionLines(Coordinate _centroid,double _angleInterval,int _count,NetTopologySuite.Geometries.LineString _envlopeLine)
        {

            List<NetTopologySuite.Geometries.LineString> pSectionLines = new List<NetTopologySuite.Geometries.LineString>();

            
            //Get the number of cross-section lines
            if(_angleInterval < 0&& _count<0)
            {
                Console.WriteLine("When reconstructing the dome, unable to obtain information such as the angle or number of cross-section lines needed for construction");
            }
            int count = 0;
            if(_angleInterval<0)
            {
                count = _count;
            }
            else
            {
                count = (int)(180.0 / _angleInterval);
            }

            double _distance = _envlopeLine.Coordinates[0].Distance(_envlopeLine.Coordinates[2]);

            for (int i=0;i<count;i++)
            {
                double arcAngle = i * _angleInterval * Math.PI / 180.0;
                double px = _centroid.X + _distance * Math.Sin(arcAngle);
                double py = _centroid.Y + _distance * Math.Cos(arcAngle);
                Coordinate point = new Coordinate(px, py);
                NetTopologySuite.Geometries.LineString pOneLine = new NetTopologySuite.Geometries.LineString(new Coordinate[] { _centroid, point });
                IGeometry pt1 = pOneLine.Intersection(_envlopeLine);


                arcAngle = arcAngle+Math.PI;
                px = _centroid.X + _distance * Math.Sin(arcAngle);
                py = _centroid.Y + _distance * Math.Cos(arcAngle);
                point = new Coordinate(px, py);
                pOneLine = new NetTopologySuite.Geometries.LineString(new Coordinate[] { _centroid, point });
                IGeometry pt2 = pOneLine.Intersection(_envlopeLine);

                NetTopologySuite.Geometries.LineString pSectionLine = new NetTopologySuite.Geometries.LineString(new Coordinate[] { pt1.Coordinate, _centroid,pt2.Coordinate });

                pSectionLines.Add(pSectionLine);
            }

            return pSectionLines;
        }

        #endregion

        #region Dome: 6 Construct three-dimensional stratigraphic lines based on cross-section lines
        /// <summary>
        /// 6.2 Calculate the Bezier curve for each stratigraphic layer at the cross-section line
        /// </summary>
        /// <param name="_pStratas">Each cross-section line and the stratigraphic layers it intersects</param>
        /// <returns></returns>
        public void GetBeziersOnDomeStrata(Dictionary<Geometry, List<FoldStratum>> _pStratas, NetTopologySuite.Index.Strtree.STRtree<OccurrencePoint> _rtree,
            NetTopologySuite.Geometries.GeometryFactory _geometryFactory, DEMRaster _rasterd)
        {

            
            List<Geometry> pOccurencePoints = new List<Geometry>();

            
            List<Geometry> pGeoLines = _pStratas.Keys.ToList<Geometry>();

            for (int i = 0; i < pGeoLines.Count; i++)
            {
                // (1) Obtain a reference point
                // This reference point is crucial because, when constructing the 2D Bezier curve, the left control point is positioned at the origin.
                // To ensure that the left side of the Bezier curve for layers other than the core layer is not at the origin, adjust by moving it a certain distance.
                
                Geometry distancePoint = new Geometry(wkbGeometryType.wkbPoint);

                for (int j = 0; j < _pStratas[pGeoLines[i]].Count; j++)
                {
                    Geometry pBoundaryPolyline = _pStratas[pGeoLines[i]][j].SPolygon.GeoPolygon.GetBoundary();
                    //pBoundaryPolyline.ExportSimpleGeometryToShapfile(this.filesavepath, "boundary");

                    //(2) Obtain the outer boundary line

                    // Due to the presence of holes in some stratigraphic surfaces, both inner and outer boundary rings exist.
                    // However, for reconstructing the paleotopography of folded structures, only the outer boundary line is needed
                    if (pBoundaryPolyline.GetGeometryType() == wkbGeometryType.wkbMultiLineString)
                    {
                        pBoundaryPolyline = pBoundaryPolyline.GetGeometryRef(0);
                        //pBoundaryPolyline.ExportSimpleGeometryToShapfile(this.filesavepath, "boundary1");
                    }

                   
                    // (3) First, check if the cross-section line intersects with the outer boundary of the stratigraphy, then find the two intersection points and assign the nearest structural attitude.
                    // (3-1) The previously generated cross-section line contains three points, with the middle point being the centroid of the innermost stratigraphic layer.
                    // To facilitate consistent ordering of structural attitude points later and to avoid multiple intersections caused by irregular outer boundaries, the cross-section line can be split at the midpoint, dividing it into two segments.
                   
                    Geometry pUpLine = new Geometry(wkbGeometryType.wkbLineString); //segment 1
                    pUpLine.AddPoint_2D(pGeoLines[i].GetX(1), pGeoLines[i].GetY(1));
                    pUpLine.AddPoint_2D(pGeoLines[i].GetX(0), pGeoLines[i].GetY(0));

                    Coordinate pCentreid = new Coordinate(pGeoLines[i].GetX(1), pGeoLines[i].GetY(1));

                    Geometry pBottomLine = new Geometry(wkbGeometryType.wkbLineString); //segment 2
                    pBottomLine.AddPoint_2D(pGeoLines[i].GetX(1), pGeoLines[i].GetY(1));
                    pBottomLine.AddPoint_2D(pGeoLines[i].GetX(2), pGeoLines[i].GetY(2));


                    //(3-2) Calculate the two intersection points
                    Geometry pt1 = null;
                    Geometry pt2 = null;

                    //segment
                    if (pUpLine.Intersect(pBoundaryPolyline))
                    {
                        //Get the intersect point
                        Geometry pIntersectPoints = pUpLine.Intersection(pBoundaryPolyline);
                        if (pIntersectPoints.GetGeometryType() == wkbGeometryType.wkbPoint)
                            pt1 = pIntersectPoints;
                        else if (pIntersectPoints.GetGeometryType() == wkbGeometryType.wkbMultiPoint)
                        {
                            //If multiple points exist, select the one closest to the centroid
                            List<Coordinate> pCoordinates = new List<Coordinate>(pIntersectPoints.GetGeometryCount());
                            for (int k = 0; k < pIntersectPoints.GetGeometryCount(); k++)
                            {
                                Geometry pt = pIntersectPoints.GetGeometryRef(k);
                                Coordinate pc = ConvertGeometryToCoordinate(pt, false);
                                pCoordinates.Add(pc);
                            }
                            Coordinate pNearestPoint = FindNearestPoint(pCentreid, pCoordinates);
                            pt1 = ConvertCoordinateToGeometry(pNearestPoint, false);
                        }
                        else
                            Console.WriteLine("Cross-section" + i + "Intersected with the stratigraphic outer boundary, but the result is not a point. Please check！");
                    }
                    else
                        Console.WriteLine("Cross-section" + i + "does not intersect with the stratigraphic outer boundary. Please check！");
                    
                    //segment 2
                    if (pBottomLine.Intersect(pBoundaryPolyline))
                    {
                        
                        Geometry pIntersectPoints = pBottomLine.Intersection(pBoundaryPolyline);
                        if (pIntersectPoints.GetGeometryType() == wkbGeometryType.wkbPoint)
                            pt2 = pIntersectPoints;
                        else if (pIntersectPoints.GetGeometryType() == wkbGeometryType.wkbMultiPoint)
                        {
                            
                            List<Coordinate> pCoordinates = new List<Coordinate>(pIntersectPoints.GetGeometryCount());
                            for (int k = 0; k < pIntersectPoints.GetGeometryCount(); k++)
                            {
                                Geometry pt = pIntersectPoints.GetGeometryRef(k);
                                Coordinate pc = ConvertGeometryToCoordinate(pt, false);
                                pCoordinates.Add(pc);
                            }
                            Coordinate pNearestPoint = FindNearestPoint(pCentreid, pCoordinates);
                            pt2 = ConvertCoordinateToGeometry(pNearestPoint, false);
                        }
                        else
                            Console.WriteLine("Cross-section" + i + "Intersected with the stratigraphic outer boundary, but the result is not a point. Please check！");
                    }
                    else
                        Console.WriteLine("Cross-section" + i + "does not intersect with the stratigraphic outer boundary. Please check！");


                    //(3-3) Obtain the structural attitude information for the two points, including dip angle, dip direction, and strike

                    OccurrencePoint pt1_Occurence = FindNearestOccurrence(pt1, _rtree, _geometryFactory, _rasterd);
                    OccurrencePoint pt2_Occurence = FindNearestOccurrence(pt2, _rtree, _geometryFactory, _rasterd);
                    if (pt1_Occurence.dipAngle == pt2_Occurence.dipAngle && pt1_Occurence.tendency == pt2_Occurence.tendency
                        && pt1_Occurence.strike == pt2_Occurence.strike)
                    {
                        Console.WriteLine("An error occurred as the structural attitudes at two points on the same cross-section line are identical. Please check!");
                        Console.WriteLine("The issue lies with the cross-section line:" + i + "stratum:" + _pStratas[pGeoLines[i]][j].SCode);
                    }

                    //(3-4)Get the reference point on the innermost statum
                    if (j == 0)
                        distancePoint.AddPoint_2D(pt1.GetX(0), pt1.GetY(0));

                    //(3-5)calculate the distance
                    float _movedistance = (float)pt1.Distance(distancePoint);

                    //(3-6)Generate the Bezier curve
                    GeoBezier _pBeGeo = new GeoBezier(_pStratas[pGeoLines[i]][j].SCode + i, pt1_Occurence, pt2_Occurence, _movedistance, this.filesavepath);

                      
                    _pStratas[pGeoLines[i]][j].beziers.Add(_pBeGeo);
                    _pStratas[pGeoLines[i]][j].pSectionLines.Add(pGeoLines[i]);

                    pOccurencePoints.Add(pt1);
                    pOccurencePoints.Add(pt2);
                }
            }
            pOccurencePoints.ExportGeometryToShapfile(this.filesavepath, "occurencepoints");

        }
        #endregion

        #endregion


    }


}
