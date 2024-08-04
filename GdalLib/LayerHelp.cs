using OSGeo.OGR;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GdalLib
{
   public static class LayerHelp
    {
        /// <summary>
        /// Getting feature based on feature paths and feature names
        /// </summary>
        /// <param name="filePath"></param>
        /// <param name="layerName"></param>
        /// <returns></returns>
        public static Layer GetLayerByLayerName(string filePath, string layerName)
        {
            DataSource ds = Ogr.Open(filePath, 0);
            if (ds == null)
            {
                throw new Exception("can not open" + filePath);
            }

            OSGeo.OGR.Driver drv = ds.GetDriver();
            if (drv == null)
            {
                throw new Exception("can not get the drive，please check！");
            }

            Layer drillLayer = ds.GetLayerByName(layerName);
            if (drillLayer == null)
                throw new Exception("fail to get the feature");

            return drillLayer;

            
        }


        /// <summary>
        /// Getting feature based on feature path
        /// </summary>
        /// <param name="filePath"></param>
        /// <param name="layerName"></param>
        /// <returns></returns>
        public static Layer GetLayerByLayerName(string pFilePath)
        {
            string filePath = System.IO.Path.GetDirectoryName(pFilePath);
            string layerName = System.IO.Path.GetFileNameWithoutExtension(pFilePath);
            DataSource ds = Ogr.Open(filePath, 0);
            if (ds == null)
            {
                throw new Exception("can not open" + filePath);
            }

            OSGeo.OGR.Driver drv = ds.GetDriver();
            if (drv == null)
            {
                throw new Exception("can not get drive，please check！");
            }

            Layer drillLayer = ds.GetLayerByName(layerName);
            if (drillLayer == null)
                throw new Exception("fail to get feature");

            return drillLayer;
        }


        /// <summary>
        /// Get all the geometry in the feature class
        /// </summary>
        /// <param name="pLayFile"></param>
        /// <returns></returns>
        public static Dictionary<int, Geometry> GetGeometryList(string pLayFile)
        {
            Dictionary<int, Geometry> pGeometry = new Dictionary<int, Geometry>();

            Layer pContourLayer = LayerHelp.GetLayerByLayerName(System.IO.Path.GetDirectoryName(pLayFile), System.IO.Path.GetFileNameWithoutExtension(pLayFile));
            int PFeatureCount = (int)pContourLayer.GetFeatureCount(0);
            for (int i = 0; i < PFeatureCount; i++)
            {
                //int id = pContourLayer.GetFeature(i).GetFieldAsInteger("Id");
                Geometry pg = pContourLayer.GetFeature(i).GetGeometryRef();
               
                pGeometry.Add(i, pg);
            }
            return pGeometry;
        }


        /// <summary>
        /// Form a separate file for each feature in the feature class
        /// </summary>
        /// <param name="_shpPath"></param>
        /// <param name="_savePath"></param>
        public static void MultiPartToSinglePart(string _shpPath,string _savePath,string _fieldName)
        {
            Layer pContourLayer = LayerHelp.GetLayerByLayerName(System.IO.Path.GetDirectoryName(_shpPath), System.IO.Path.GetFileNameWithoutExtension(_shpPath));

            int PFeatureCount = (int)pContourLayer.GetFeatureCount(0);
            for (int i = 0; i < PFeatureCount; i++)
            {
                int id = pContourLayer.GetFeature(i).GetFieldAsInteger(_fieldName);
                Geometry pg = pContourLayer.GetFeature(i).GetGeometryRef();

                ExportSimpleGeometryToShapfile(pg, _savePath, "Right"+id.ToString());
            }
        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="_triMesh"></param>
        /// <param name="_workSpacePath"></param>
        /// <param name="_fileName"></param>
        public static void ExportSimpleGeometryToShapfile( Geometry geometryCollection, string _workSpacePath, string _fileName)
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


            Layer poLayer = poDS.CreateLayer(_fileName, null, wkbGeometryType.wkbTINZ, null);
            if (poLayer == null)
                throw new Exception("Layer Creation Failed");

            FieldDefn oFieldID = new FieldDefn("FieldID", FieldType.OFTInteger);
            poLayer.CreateField(oFieldID, 1);

            Feature poFeature = new Feature(poLayer.GetLayerDefn());

            poFeature.SetField(0, 0);

            Geometry kkk =geometryCollection.GetGeometryRef(0);
            wkbGeometryType ptt = kkk.GetGeometryType();

            poFeature.SetGeometry(geometryCollection);

            poLayer.CreateFeature(poFeature);


            poDS.Dispose();
            poLayer.Dispose();
            poFeature.Dispose();
        }
    }
}
