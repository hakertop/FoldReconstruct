using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using OSGeo.OGR;


namespace GeoCommon
{
    /// <summary>
    /// 
    /// </summary>
    public class Polygon3D
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
        /// 
        /// </summary>
        public List<Eage3D> eageList;

        /// <summary>
        /// max thickness
        /// </summary>
        public double SMaxT { get; set; }

        /// <summary>
        /// min thickness
        /// </summary>
        public double SMinT { get; set; }


        /// <summary>
        /// 
        /// </summary>
        public Geometry GeoPolygon { get; set; }

        /// <summary>
        /// 
        /// </summary>
        public Polygon3D()
        {
            name = "untitled";

            vertexList = new List<Vertex3D>();

            eageList = new List<Eage3D>();

            GeoPolygon = new Geometry(wkbGeometryType.wkbPolygon);
        }

        /// <summary>
        /// 
        /// </summary>
        public Polygon3D(Geometry pGeoPolygon,String pPolygonName)
        {
            name = pPolygonName;

            vertexList = new List<Vertex3D>();

            eageList = new List<Eage3D>();

            GeoPolygon = pGeoPolygon;
        }



        /// <summary>
        /// add vertex
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="z"></param>
        public int AddVertex(Vertex3D modelVertex)
        {
            
            Vertex3D geoVertex = new Vertex3D();
            geoVertex.id = vertexList.Count;
            
            geoVertex.x = modelVertex.x;
            geoVertex.y = modelVertex.y;

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
        /// add vertex
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="z"></param>
        public int AddVertex(double x, double y)
        {

            
            Vertex3D geoVertex = new Vertex3D();
            geoVertex.id = vertexList.Count;
            geoVertex.x = x ;
            geoVertex.y = y ;

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
        /// add eage
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="z"></param>
        public void AddEage(Eage3D modelEage)
        {
            eageList.Add(modelEage);
        }

        /// <summary>
        /// add eage
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="z"></param>
        public void AddEage(List<Vertex3D> listVertex)
        {
            Eage3D eage = new Eage3D();
            eage.id = eageList.Count;
     
            eageList.Add(eage);
        }

        /// <summary>
        /// get the polygon area
        /// </summary>
        /// <returns></returns>
       public float GetArea()
        {
           return (float)GeoPolygon.GetArea();
        }



    }
}
