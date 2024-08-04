using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using GdalLib;
using GeoCommon;
using GeologicalEntity;
using OSGeo.OGR;

namespace FoldRestortEntry
{
    /// <summary>
    /// voronoi subdivision
    /// </summary>
    public class StrataSubdivision
    {
        /// <summary>
        /// file save path
        /// </summary>
        public string filesavepath;

        public StrataSubdivision(string _filesavepath)
        {
            this.filesavepath = _filesavepath;
        }

        /// <summary>
        /// Main function
        /// </summary>
        /// <param name="_file"></param>
        public void DataEnty(string _file)
        {
            
            Layer pLayer = LayerHelp.GetLayerByLayerName(_file);

            int fieldIndex = pLayer.GetLayerDefn().GetFieldIndex("NAME");

            Feature feature = null;

            
            List<Vertex3D> pPointsInPolygon = new List<Vertex3D>();

            
            while ((feature = pLayer.GetNextFeature()) != null)
            {
                string p1Value = feature.GetFieldAsString(fieldIndex);

                
                if (p1Value == "P2")
                {
                    Geometry coreGeo = feature.GetGeometryRef();
                    wkbGeometryType ps = coreGeo.GetGeometryType();
                    List<Geometry> pAllPoints = GeometrySolve.GetAllPointsFromPolygon(coreGeo, 2);

                    pAllPoints.ExportGeometryToShapfile(this.filesavepath, "Point1");

                    foreach (var pt in pAllPoints)
                    {
                        Vertex3D vd = new Vertex3D(pt.GetX(0), pt.GetY(0));
                        vd.belongToPolygon = p1Value;
                        pPointsInPolygon.Add(vd);
                    }
                }
               
                if (p1Value == "C1")
                {
                    Geometry coreGeo = feature.GetGeometryRef();
                    wkbGeometryType ps = coreGeo.GetGeometryType();
                    List<Geometry> pAllPoints = GeometrySolve.GetAllPointsFromPolygon(coreGeo, 2);
                    foreach (var pt in pAllPoints)
                    {
                        Vertex3D vd = new Vertex3D(pt.GetX(0), pt.GetY(0));
                        vd.belongToPolygon = p1Value;
                        pPointsInPolygon.Add(vd);
                    }
                }

                if (p1Value == "G2")
                {
                    Geometry coreGeo = feature.GetGeometryRef();
                    wkbGeometryType ps = coreGeo.GetGeometryType();
                    List<Geometry> pAllPoints = GeometrySolve.GetAllPointsFromPolygon(coreGeo, 2);
                    foreach (var pt in pAllPoints)
                    {
                        Vertex3D vd = new Vertex3D(pt.GetX(0), pt.GetY(0));
                        vd.belongToPolygon = p1Value;
                        pPointsInPolygon.Add(vd);
                    }
                }
            }

            
            string[] objectNames = new string[3]{ "P2", "C1","G2"};
            CreateVoronoi(pPointsInPolygon, objectNames);
        }

     
        public void CreateVoronoi(List<Vertex3D> _nodes, string[] _objectNames)
        {
            
            List<TriangleNet.Geometry.Vertex> pVertexList = GetListVetexs(_nodes);
            //pVertexList.ExportTrianglePointToShapefile(this.filesavepath, "pointBO");

            
            TriangleNet.Mesh mesh = GetTriMesh(pVertexList);

            
            TriMesh3D trimesh = new TriMesh3D();
            foreach (var item in mesh.Triangles)
            {
                TriangleNet.Geometry.Vertex p1, p2, p3;
                TriangleNet.Topology.Triangle pTinTri = item;
                p1 = item.GetVertex(0);
                p2 = item.GetVertex(1);
                p3 = item.GetVertex(2);
                
                trimesh.AddTriangle(p1.X, p1.Y, 0, p2.X, p2.Y, 0, p3.X, p3.Y, 0);
            }
            //trimesh.ExportTriMeshToShapfile(this.filesavepath, "delaunaryBO");

           
            TriangleNet.Voronoi.Legacy.SimpleVoronoi sv = new TriangleNet.Voronoi.Legacy.SimpleVoronoi(mesh);


            
            ICollection<TriangleNet.Voronoi.Legacy.VoronoiRegion> pRegion = sv.Regions;
            List<TriangleNet.Geometry.Vertex> pVertexs = new List<TriangleNet.Geometry.Vertex>();
            foreach (var mV in mesh.Vertices)
            {
                pVertexs.Add(mV);
            }

            
            List<List<TriangleNet.Voronoi.Legacy.VoronoiRegion>> pClassficationVoronoi = ClassificationVoronoiRegion(_objectNames, pRegion, pVertexs);

            
            Dictionary<string, Geometry> pDictVoronoi = GetMultiPolygonMergedVoronois(pClassficationVoronoi, _objectNames);


           
            Dictionary<string, Geometry> pSection = new Dictionary<string, Geometry>();
            for(int i=0;i< pDictVoronoi.Keys.Count-1;i++)
            {
                for(int j=1;j< pDictVoronoi.Keys.Count;j++)
                {
                    if(pDictVoronoi[pDictVoronoi.Keys.ToList()[i]].Intersect(pDictVoronoi[pDictVoronoi.Keys.ToList()[j]]))
                    {
                        Geometry gLine = pDictVoronoi[pDictVoronoi.Keys.ToList()[i]].Intersection(pDictVoronoi[pDictVoronoi.Keys.ToList()[j]]);
                        pSection.Add(pDictVoronoi.Keys.ToList()[i]+ pDictVoronoi.Keys.ToList()[j],gLine);
                    }
                }
            }
            pSection.Values.ToList().ExportGeometryToShapfile(this.filesavepath, "intersectlines");

        }

        
        public Geometry GetNewSection(Geometry pOriginalGeo, List<List<Vertex3D>> pOriginalVertex, double height)
        {
            List<Vertex3D> pNewOriginalVertex = new List<Vertex3D>();

            foreach (var vts in pOriginalVertex)
            {
                foreach (var vt in vts)
                {
                    pNewOriginalVertex.Add(vt);
                }
            }

            int count = pOriginalGeo.GetGeometryCount();

            Geometry pLine = pOriginalGeo.GetGeometryRef(0);
            Geometry pRing = new Geometry(wkbGeometryType.wkbLinearRing);
            for (int i = 0; i < pLine.GetPointCount(); i++)
            {
                bool flag = true;
                for (int j = 0; j < pNewOriginalVertex.Count; j++)
                {
                    if (Math.Round(pLine.GetX(i), 3) == Math.Round(pNewOriginalVertex[j].x, 3) && Math.Round(pLine.GetY(i), 3) == Math.Round(pNewOriginalVertex[j].y, 3))
                    {
                        pRing.AddPoint(pLine.GetX(i), pLine.GetY(i), pNewOriginalVertex[j].z);
                        flag = false;
                    }
                }

                if (flag)
                {
                    pRing.AddPoint(pLine.GetX(i), pLine.GetY(i), height);
                }
            }

            Geometry pPolygon = new Geometry(wkbGeometryType.wkbPolygon);
            pPolygon.AddGeometry(pRing);
            return pPolygon;

        }


     
        public List<List<TriangleNet.Voronoi.Legacy.VoronoiRegion>> ClassificationVoronoiRegion(string[] _objectNames,
            ICollection<TriangleNet.Voronoi.Legacy.VoronoiRegion> pRegion
            , List<TriangleNet.Geometry.Vertex> pVertexList)
        {
            List<List<TriangleNet.Voronoi.Legacy.VoronoiRegion>> pListVoronoi = new List<List<TriangleNet.Voronoi.Legacy.VoronoiRegion>>();
            foreach (var objectName in _objectNames)
            {
                List<TriangleNet.Voronoi.Legacy.VoronoiRegion> voronoiReg = new List<TriangleNet.Voronoi.Legacy.VoronoiRegion>();
                foreach (var vregion in pRegion)
                {
                    TriangleNet.Geometry.Vertex vt = pVertexList[vregion.ID];
                    if (vt.NAME == objectName)
                    {
                        voronoiReg.Add(vregion);
                    }
                }
                pListVoronoi.Add(voronoiReg);
            }
            return pListVoronoi;
        }

        /// <summary>
        /// combine the classfied Voronoi polygons into a polygon
        /// </summary>
        /// <param name="pClassficationVoronoi"></param>
        /// <param name="_objectNames"></param>
        /// <returns></returns>
        public Dictionary<string, Geometry> GetMultiPolygonMergedVoronois(List<List<TriangleNet.Voronoi.Legacy.VoronoiRegion>> pClassficationVoronoi,string[] _objectNames)
        {
            Dictionary<string, Geometry> pDictVoronoi = new Dictionary<string, Geometry>();
            for (int i = 0; i < pClassficationVoronoi.Count; i++)
            {
                List<TriangleNet.Voronoi.Legacy.VoronoiRegion> pvor = pClassficationVoronoi[i];

                //obtain all Voronoi faces
                List<Geometry> polygonlist = new List<Geometry>();
                for (int j = 0; j < pvor.Count; j++)
                {
                    Geometry linev = new Geometry(wkbGeometryType.wkbLinearRing);
                    List<TriangleNet.Geometry.Point> pvs = pvor[j].Vertices as List<TriangleNet.Geometry.Point>;
                    for (int k = 0; k < pvs.Count; k++)
                    {
                        linev.AddPoint_2D(pvs[k].X, pvs[k].Y);
                    }
                    linev.AddPoint_2D(pvs[0].X, pvs[0].Y);
                    Geometry polygonv = new Geometry(wkbGeometryType.wkbPolygon);
                    polygonv.AddGeometry(linev);
                    polygonlist.Add(polygonv);
                }

                polygonlist.ExportGeometryToShapfile(this.filesavepath, "vorvint" + i);

                // combine all triangle faces
                Geometry gp = new Geometry(wkbGeometryType.wkbPolygon);
                for (int j = 0; j < polygonlist.Count; j++)
                {
                    double pArea = polygonlist[j].GetArea();
                    gp = gp.Union(polygonlist[j]);
                }

                // output the data
                List<Geometry> gs = new List<Geometry>();
                gs.Add(gp);
                gs.ExportGeometryToShapfile(this.filesavepath, "combinge" + i);

                pDictVoronoi.Add(_objectNames[i], gp);
            }
            return pDictVoronoi;
        }

        /// <summary>
        /// return vertex set
        /// </summary>
        /// <param name="pEages"></param>
        /// <returns></returns>
        private List<TriangleNet.Geometry.Vertex> GetListVetexs(List<Vertex3D> _nodes)
        {
            List<TriangleNet.Geometry.Vertex> pVertexList = new List<TriangleNet.Geometry.Vertex>();
            int i = 0;
            foreach (var vt in _nodes)
            {
                TriangleNet.Geometry.Vertex vr = new TriangleNet.Geometry.Vertex();
                vr.ID = i;
                vr.X = vt.x;
                vr.Y = vt.y;
                vr.NAME = vt.belongToPolygon;
                pVertexList.Add(vr);

                i++;
            }
            
            return pVertexList;
        }


        /// <summary>
        /// triangulation
        /// </summary>
        /// <param name="pA"></param>
        /// <param name="pB"></param>
        /// <returns></returns>
        private TriangleNet.Mesh GetTriMesh(List<TriangleNet.Geometry.Vertex> pA)
        {
            #region trianglation
          
           
            var options = new TriangleNet.Meshing.ConstraintOptions();
            options.SegmentSplitting = 1;
            options.ConformingDelaunay = false;
            options.Convex = false;

            
            var quality = new TriangleNet.Meshing.QualityOptions();
            TriangleNet.Geometry.IPolygon input = GetPolygon(pA);
            //TriangleNet.Geometry.Contour con = GetContourByTriangle(pA);
            //input.Add(con, false);


            TriangleNet.Mesh mesh = null;
            if (input != null)
            {
                mesh = (TriangleNet.Mesh)TriangleNet.Geometry.ExtensionMethods.Triangulate(input, options);

            }

            return mesh;
            #endregion

        }

        /// <summary>
        /// Greate the Contour object 
        /// </summary>
        /// <param name="pA"></param>
        /// <param name="pB"></param>
        /// <returns></returns>
        public TriangleNet.Geometry.Contour GetContourByTriangle(List<TriangleNet.Geometry.Vertex> pA)
        {
            List<TriangleNet.Geometry.Vertex> pv = new List<TriangleNet.Geometry.Vertex>();

            foreach (var vt in pA)
            {
                pv.Add(vt);
            }
            TriangleNet.Geometry.Contour pNewCon = new TriangleNet.Geometry.Contour(pv);
            return pNewCon;
        }

        /// <summary>
        /// create POlygon
        /// </summary>
        /// <param name="drillList"></param>
        /// <returns></returns>
        private TriangleNet.Geometry.IPolygon GetPolygon(List<TriangleNet.Geometry.Vertex> pA)
        {
            TriangleNet.Geometry.IPolygon data = new TriangleNet.Geometry.Polygon();

            foreach (var vt in pA)
            {
                TriangleNet.Geometry.Vertex triVertex = new TriangleNet.Geometry.Vertex(vt.X, vt.Y);
                triVertex.NAME = vt.NAME;
                data.Add(triVertex);

            }
            return data;
        }
    }
}
