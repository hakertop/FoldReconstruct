using OSGeo.OGR;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GeoCommon
{
    /// <summary>
    /// 
    /// </summary>
    public class Eage3D
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
        /// vertex set
        /// </summary>
        public List<Vertex3D> vertexList;

        /// <summary>
        /// shp 
        /// </summary>
        public Geometry GeoPolyline { get; set; }


        /// <summary>
        /// 
        /// </summary>
        public Eage3D()
        {
            
            vertexList = new List<Vertex3D>();
            Geometry GeoPolyline = new Geometry(wkbGeometryType.wkbLineString);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="pGeoPolyline"></param>
        /// <param name="pPolylineName"></param>
        public Eage3D(Geometry pGeoPolyline)
        {
            
            vertexList = new List<Vertex3D>();
            GeoPolyline = pGeoPolyline;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="modelVertex"></param>
        /// <returns></returns>
        public int AddVertex(Vertex3D polygonVertex)
        {
            
            Vertex3D geoVertex = new Vertex3D();
            geoVertex.name = polygonVertex.name;
            geoVertex.id = vertexList.Count;

            geoVertex.x = polygonVertex.x;
            geoVertex.y = polygonVertex.y;
          

            foreach (var vertex in vertexList)
            {
                if (Math.Abs(vertex.x - geoVertex.x) < 0.01 &&
                    Math.Abs(vertex.y - geoVertex.y) < 0.01)
                {
                    return vertex.id;
                }
            }

            vertexList.Add(geoVertex);

            return geoVertex.id;
        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="z"></param>
        public int AddVertex(double x, double y)
        {
            Vertex3D geoVertex = new Vertex3D();
            geoVertex.id = vertexList.Count;
           

            geoVertex.x = x;
            geoVertex.y = y;

            foreach (var vertex in vertexList)
            {
                if (Math.Abs(vertex.x - geoVertex.x) < 0.01 &&
                    Math.Abs(vertex.y - geoVertex.y) < 0.01)
                {
                    return vertex.id;
                }
            }

            vertexList.Add(geoVertex);

            return geoVertex.id;
        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="pEageA"></param>
        /// <param name="pEageB"></param>
        public Vertex3D Get3DCentralPoint()
        {
            double distX = 0;
            double distY = 0;

            foreach (var item in this.vertexList)
            {
                distX = distX + item.x;
                distY = distY + item.y;
            }

            int count = this.vertexList.Count;
            Vertex3D vt = new Vertex3D(distX / count, distY / count);

            return vt;
        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="pEageA"></param>
        /// <param name="pEageB"></param>
        public Vertex3D Get2DCentralPoint()
        {
            double distX = 0;
            double distY = 0;

            foreach (var item in this.vertexList)
            {
                distX = distX + item.x;
                distY = distY + item.y;
            }

            int count = this.vertexList.Count;
            Vertex3D vt = new Vertex3D(distX / count, distY / count);

            return vt;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public double GetLength()
        {
            return GeoPolyline.Length();
        }

    }
}
