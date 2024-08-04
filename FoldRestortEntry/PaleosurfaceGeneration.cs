using GeoAPI.Geometries;
using System;
using System.Collections.Generic;
using System.Text;

using GdalLib;
using g4;
using System.Linq;

namespace FoldRestortEntry
{
    public class PaleosurfaceGeneration
    {

        /// <summary>
        /// Reconstruct the stratigraphic paleo-surfaces of different strata
        /// </summary>
        /// <param name="_fileSurfacePoints">3D points on the surface</param>
        /// <param name="_fileBoundaryPoints">Boundary points</param>
        /// <param name="_saveFile">save filepath</param>
        /// <param name="_fileName">file name</param>
        /// <param name="_filetype">3d type (.obj, .fpx,...)</param>
        public void ReconstructStratigraphicPaleosurface(string _fileSurfacePoints, string _fileBoundaryPoints,
            string _saveFile, string _fileName,string _filetype)
        {
            // (1) obtain points from the filefolder
            List<IGeometry> _topSurfacePoints = NetClassExtensionMethods.GetAllGeometry(_fileSurfacePoints);

            List<IGeometry> _boundaryPoints = NetClassExtensionMethods.GetAllGeometry(_fileSurfacePoints);

            // (2) create tri_mesh
            TriangleNet.Mesh polygonMeshOfTopSurface = CreateDelaunaryTri(_topSurfacePoints, new List<List<IGeometry>>() { _boundaryPoints });


            // (3) export 3D file type using the Class DMesh3
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

            // (4) Output
            StandardMeshWriter.WriteMesh(_saveFile + "\\" + _fileName + _filetype, dmesh, WriteOptions.Defaults);
        }




        /// <summary>
        /// Reconstruct the triangulated mesh by TriangleNet
        /// </summary>
        /// <param name="_pAllPoints">All points on the stratigarphic surface</param>
        /// <param name="_pBoundaryConstrainPoints">Boundary constrain points ensure correct geometric topology</param>
        /// <returns></returns>
        public TriangleNet.Mesh CreateDelaunaryTri(List<IGeometry> _pAllPoints, List<List<IGeometry>> _pBoundaryConstrainPoints)
        {
            //1. Constrain options
            var options = new TriangleNet.Meshing.ConstraintOptions();
            options.SegmentSplitting = 2;
            options.ConformingDelaunay = true;

            //2. Construct IPolgon
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

            //Quality options
            var quality = new TriangleNet.Meshing.QualityOptions();

            TriangleNet.Mesh mesh = null;
            if (data != null)
            {
                mesh = (TriangleNet.Mesh)TriangleNet.Geometry.ExtensionMethods.Triangulate(data, options);
            }

            return mesh;
        }


        /// <summary>
        /// Get the elevation of a point. This function only be used when constructing the triangulated mesh
        /// </summary>
        /// <param name="_vt"></param>
        /// <param name="_points"></param>
        /// <returns>return -1 when don't get the elevation</returns>
        public double GetElevationFromPoint(TriangleNet.Geometry.Vertex _vt, List<IGeometry> _points)
        {
            foreach (var vp in _points)
            {
                if (Math.Abs(_vt.X - vp.Coordinate.X) < 0.001 && Math.Abs(_vt.Y - vp.Coordinate.Y) < 0.001)
                    return vp.Coordinate.Z;
            }

            return -1;
        }
    }
}
