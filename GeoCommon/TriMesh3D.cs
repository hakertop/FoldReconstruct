using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace GeoCommon
{
    public class TriMesh3D
    {
        /// <summary>
        /// ID
        /// </summary>
        public int id;

        /// <summary>
        /// name
        /// </summary>
        public string name;
        /// <summary>
        /// 
        /// </summary>
        public IList<Vertex3D> VertexList;

        /// <summary>
        /// 
        /// </summary>
        public IList<Triangle3D> triangleList;

        /// <summary>
        /// 
        /// </summary>
        public TriMesh3D()
        {
            name = "untitled";

            VertexList = new List<Vertex3D>();

            triangleList = new List<Triangle3D>();
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="z"></param>
        public int AddVertex3D(double x, double y, double z)
        {

            
            Vertex3D geoVertex3D = new Vertex3D();
            geoVertex3D.id = VertexList.Count;
            geoVertex3D.x = Math.Round(x, 2);
            geoVertex3D.y = Math.Round(y, 2);
            geoVertex3D.z = Math.Round(z, 2);

            foreach (var Vertex3D in VertexList)
            {
                if (Math.Abs(Vertex3D.x - geoVertex3D.x) < 0.01 &&
                    Math.Abs(Vertex3D.y - geoVertex3D.y) < 0.01 &&
                    Math.Abs(Vertex3D.z - geoVertex3D.z) < 0.01)
                {
                    return Vertex3D.id;
                }
            }

            VertexList.Add(geoVertex3D);

            return geoVertex3D.id;
        }

        /// <summary>
        /// add a triangle into the triangle mesh  -way 1
        /// </summary>
        /// <param name="id"></param>
        /// <param name="v0"></param>
        /// <param name="v1"></param>
        /// <param name="v2"></param>
        public void AddTriangle(int id, int v0, int v1, int v2)
        {
            Triangle3D geoTri = new Triangle3D();

            geoTri.id = id;
            geoTri.v0 = v0;
            geoTri.v1 = v1;
            geoTri.v2 = v2;

            triangleList.Add(geoTri);
        }

        /// <summary>
        /// add a triangle into the triangle mesh  -way 2
        /// </summary>
        /// <param name="x0"></param>
        /// <param name="y0"></param>
        /// <param name="z0"></param>
        /// <param name="x1"></param>
        /// <param name="y1"></param>
        /// <param name="z1"></param>
        /// <param name="x2"></param>
        /// <param name="y2"></param>
        /// <param name="z2"></param>
        public void AddTriangle(double x0, double y0, double z0,
            double x1, double y1, double z1, double x2, double y2, double z2)
        {
            int gv0 = AddVertex3D(x0, y0, z0);
            int gv1 = AddVertex3D(x1, y1, z1);
            int gv2 = AddVertex3D(x2, y2, z2);

            Triangle3D geoTri = new Triangle3D();
            geoTri.id = triangleList.Count;
            geoTri.v0 = gv0;
            geoTri.v1 = gv1;
            geoTri.v2 = gv2;
            triangleList.Add(geoTri);
        }

        /// <summary>
        /// get a triangle
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public List<Triangle3D> getTriangleByCoord(double x, double y)
        {
            List<Triangle3D> triBedrock = new List<Triangle3D>();

            for (int i = 0; i < triangleList.Count; i++)
            {
                Triangle3D tri = triangleList[i];

                if ((Math.Round(VertexList[tri.v0].x, 2) == Math.Round(x, 2) && Math.Round(VertexList[tri.v0].y, 2) == Math.Round(y, 2))
                        || (Math.Round(VertexList[tri.v1].x, 2) == Math.Round(x, 2) && Math.Round(VertexList[tri.v1].y, 2) == Math.Round(y, 2))
                        || (Math.Round(VertexList[tri.v2].x, 2) == Math.Round(x, 2) && Math.Round(VertexList[tri.v2].y, 2) == Math.Round(y, 2)))
                    triBedrock.Add(tri);
            }

            return triBedrock;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="_workspacePath"></param>
        /// <param name="fileName"></param>
        public void ExportTriMesh(string _workspacePath, string fileName)
        {
            StreamWriter streamWriter = new StreamWriter(_workspacePath + "\\" + fileName + ".txt");
            streamWriter.WriteLine(VertexList.Count.ToString() + " " + triangleList.Count.ToString());
            for (int i = 0; i < VertexList.Count; i++)
            {
                Vertex3D v = VertexList[i];
                streamWriter.WriteLine(v.x.ToString() + " " + v.y.ToString() + " " + v.z.ToString());
            }
            for (int i = 0; i < triangleList.Count; i++)
            {
                Triangle3D tri = triangleList[i];
                streamWriter.WriteLine("3 " + tri.v0.ToString() + " " + tri.v1.ToString() + " " + tri.v2.ToString());
            }

            streamWriter.Close();
        }
    }
}
