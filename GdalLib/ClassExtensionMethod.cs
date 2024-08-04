using OSGeo.OGR;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using GeoCommon;

using GeoAPI.IO;
using GeoAPI.Geometries;
using NetTopologySuite;

namespace GdalLib
{
    public static class ClassExtensionMethod
    {


      
        public static void ExportTriMeshToShapfile(this TriMesh3D _triMesh, string _workSpacePath, string _fileName)
        {
            
            string pszDriverName = "ESRI Shapefile";
            
            OSGeo.OGR.Driver poDriver = OSGeo.OGR.Ogr.GetDriverByName(pszDriverName);
            if (poDriver == null)
                throw new Exception("Driver Error");

            
            OSGeo.GDAL.Gdal.SetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
            
            OSGeo.GDAL.Gdal.SetConfigOption("SHAPE_ENCODING", "");


            
            OSGeo.OGR.DataSource poDS;
            poDS = poDriver.CreateDataSource(_workSpacePath + "\\" + _fileName + ".shp", null);
            if (poDS == null)
                throw new Exception("DataSource Creation Error");

            
            Layer poLayer = poDS.CreateLayer(_fileName, null, wkbGeometryType.wkbPolygon, null);
            if (poLayer == null)
                throw new Exception("Layer Creation Failed");

            FieldDefn oFieldID = new FieldDefn("FieldID", FieldType.OFTInteger);
            poLayer.CreateField(oFieldID, 1);

           
            Feature poFeature = new Feature(poLayer.GetLayerDefn());

            for (int i = 0; i < _triMesh.triangleList.Count; i++)
            {
                Triangle3D tri = _triMesh.triangleList[i];
                Vertex3D vTri0 = _triMesh.VertexList[tri.v0];
                Vertex3D vTri1 = _triMesh.VertexList[tri.v1];
                Vertex3D vTri2 = _triMesh.VertexList[tri.v2];


                string polygonTri = "POLYGON((" + vTri0.x + " " + vTri0.y + "," + vTri1.x + " " + vTri1.y + "," + vTri2.x + " " + vTri2.y + "," + vTri0.x + " " + vTri0.y + "))";
                poFeature.SetField(0, i);
                Geometry polyTri = Geometry.CreateFromWkt(polygonTri);
                poFeature.SetGeometry(polyTri);

                poLayer.CreateFeature(poFeature);

            }
            poDS.Dispose();
            poLayer.Dispose();
            poFeature.Dispose();
        }


        public static void ExportTriMeshToShapfile(this TriangleNet.Mesh polygonMesh, string _workSpacePath, string _fileName)
        {
            List<IGeometry> pgeoers = new List<IGeometry>();
            foreach (var vtri in polygonMesh.Triangles)
            {
                TriangleNet.Geometry.Vertex vt1 = vtri.GetVertex(0);
                TriangleNet.Geometry.Vertex vt2 = vtri.GetVertex(1);
                TriangleNet.Geometry.Vertex vt3 = vtri.GetVertex(2);

                NetTopologySuite.Geometries.LinearRing plr = new NetTopologySuite.Geometries.LinearRing(new Coordinate[]
                {
            new Coordinate(vt1.X, vt1.Y),
            new Coordinate(vt2.X, vt2.Y),
            new Coordinate(vt3.X, vt3.Y),
            new Coordinate(vt1.X, vt1.Y)});
                NetTopologySuite.Geometries.Polygon polg = new NetTopologySuite.Geometries.Polygon(plr);
                pgeoers.Add(polg);
            }
            pgeoers.ExportGeometryToShapfileByNet(_workSpacePath, _fileName);
        }


       
        public static void ExportTrianglePointToShapefile(this List<TriangleNet.Geometry.Vertex> PNewEage, string pSavePath, string pName)
        {
            List<Geometry> ListPT = new List<Geometry>();
            for (int i = 0; i < PNewEage.Count; i++)
            {
                Geometry pt = new Geometry(wkbGeometryType.wkbPoint);
                pt.AddPoint_2D(PNewEage[i].X, PNewEage[i].Y);
                ListPT.Add(pt);
            }

            ListPT.ExportGeometryToShapfile(pSavePath, pName);
        }

       
        public static void ExportSimpleGeometryToShapfile(this Geometry geometryCollection, string _workSpacePath, string _fileName)
        {
            
            string pszDriverName = "ESRI Shapefile";
            
            OSGeo.OGR.Driver poDriver = OSGeo.OGR.Ogr.GetDriverByName(pszDriverName);
            if (poDriver == null)
                throw new Exception("Driver Error");

            
            OSGeo.GDAL.Gdal.SetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
            
            OSGeo.GDAL.Gdal.SetConfigOption("SHAPE_ENCODING", "");


           
            OSGeo.OGR.DataSource poDS;
            poDS = poDriver.CreateDataSource(_workSpacePath + "\\" + _fileName + ".shp", null);//如果原始文件夹内有该要素数据，会覆盖。
            if (poDS == null)
                throw new Exception("DataSource Creation Error");


            
            Layer poLayer = poDS.CreateLayer(_fileName, null, geometryCollection.GetGeometryType(), null);
            if (poLayer == null)
                throw new Exception("Layer Creation Failed");



            FieldDefn oFieldID = new FieldDefn("FieldID", FieldType.OFTInteger);
            poLayer.CreateField(oFieldID, 1);

            
            Feature poFeature = new Feature(poLayer.GetLayerDefn());


            poFeature.SetField(0, 0);

            poFeature.SetGeometry(geometryCollection);

            poLayer.CreateFeature(poFeature);


            poDS.Dispose();
            poLayer.Dispose();
            poFeature.Dispose();
        }




        
        public static void ExportGeometryToShapfile(this List<Geometry> geometryCollection, string _workSpacePath, string _fileName)
        {
            if (geometryCollection.Count <= 0)
                return;

            
            string pszDriverName = "ESRI Shapefile";
            
            OSGeo.OGR.Driver poDriver = OSGeo.OGR.Ogr.GetDriverByName(pszDriverName);
            if (poDriver == null)
                throw new Exception("Driver Error");


           
            OSGeo.OGR.DataSource poDS;
            poDS = poDriver.CreateDataSource(_workSpacePath + "\\" + _fileName + ".shp", null);
            if (poDS == null)
                throw new Exception("DataSource Creation Error");

       

            Layer poLayer = poDS.CreateLayer(_fileName, null, geometryCollection[0].GetGeometryType(), null);
            if (poLayer == null)
                throw new Exception("Layer Creation Failed");

            FieldDefn oFieldID = new FieldDefn("FieldID", FieldType.OFTInteger);
            poLayer.CreateField(oFieldID, 1);

         
            Feature poFeature = new Feature(poLayer.GetLayerDefn());

            for (int i = 0; i < geometryCollection.Count; i++)
            {

                poFeature.SetField(0, i);

                poFeature.SetGeometry(geometryCollection[i]);

                poLayer.CreateFeature(poFeature);

            }
            poDS.Dispose();
            poLayer.Dispose();
            poFeature.Dispose();
        }

      
        public static void Export3DPoints(this List<Geometry> geometryCollection, string _workSpacePath, string _fileName)
        {
            if (geometryCollection.Count <= 0)
                return;

          
            string pszDriverName = "ESRI Shapefile";
            
            OSGeo.OGR.Driver poDriver = OSGeo.OGR.Ogr.GetDriverByName(pszDriverName);
            if (poDriver == null)
                throw new Exception("Driver Error");


            OSGeo.OGR.DataSource poDS;
            poDS = poDriver.CreateDataSource(_workSpacePath + "\\" + _fileName + ".shp", null);
            if (poDS == null)
                throw new Exception("DataSource Creation Error");


            Layer poLayer = poDS.CreateLayer(_fileName, null, wkbGeometryType.wkbPoint, null);
            if (poLayer == null)
                throw new Exception("Layer Creation Failed");

            FieldDefn oFieldID = new FieldDefn("FieldID", FieldType.OFTInteger);
            poLayer.CreateField(oFieldID, 1);

            FieldDefn oFieldID_1 = new FieldDefn("z", FieldType.OFTReal);
            poLayer.CreateField(oFieldID_1, 1);

           
            Feature poFeature = new Feature(poLayer.GetLayerDefn());

            for (int i = 0; i < geometryCollection.Count; i++)
            {
                Geometry pt = new Geometry(wkbGeometryType.wkbPoint);
                pt.AddPoint_2D(geometryCollection[i].GetX(0), geometryCollection[i].GetY(0));


                poFeature.SetField(0, i);
                poFeature.SetField(1, geometryCollection[i].GetZ(0));

                poFeature.SetGeometry(pt);

                poLayer.CreateFeature(poFeature);

            }
            poDS.Dispose();
            poLayer.Dispose();
            poFeature.Dispose();
        }

      
        public static void ExportTriMeshToShapfile(string _workSpacePath, string _fileName)
        {           
            string pszDriverName = "ESRI Shapefile";
           
            OSGeo.OGR.Driver poDriver = OSGeo.OGR.Ogr.GetDriverByName(pszDriverName);
            if (poDriver == null)
                throw new Exception("Driver Error");

            
            OSGeo.GDAL.Gdal.SetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
            
            OSGeo.GDAL.Gdal.SetConfigOption("SHAPE_ENCODING", "");


           
            OSGeo.OGR.DataSource poDS;
            poDS = poDriver.CreateDataSource(_workSpacePath + "\\" + _fileName + ".shp", null);//如果原始文件夹内有该要素数据，会覆盖。
            if (poDS == null)
                throw new Exception("DataSource Creation Error");
            
            Layer poLayer = poDS.CreateLayer(_fileName, null, wkbGeometryType.wkbPoint, null);
            if (poLayer == null)
                throw new Exception("Layer Creation Failed");

           
            FieldDefn oFieldID = new FieldDefn("FieldID", FieldType.OFTInteger);
            poLayer.CreateField(oFieldID, 1);
          
            Feature poFeature = new Feature(poLayer.GetLayerDefn());
            for (int i = 0; i <= 5; i++)
            {
                Geometry pt = new Geometry(wkbGeometryType.wkbPoint);
                pt.AddPoint_2D(i * 3.14, i * 3.14);

               
                poFeature.SetField(0, i);
                poFeature.SetGeometry(pt);

                poLayer.CreateFeature(poFeature);
            }
            
            poDS.Dispose();
            poLayer.Dispose();
            poFeature.Dispose();
        }

    }
}
