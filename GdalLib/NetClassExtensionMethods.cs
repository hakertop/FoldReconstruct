using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using GeoAPI.Geometries;
using NetTopologySuite.Features;
using NetTopologySuite.IO;

namespace GdalLib
{
    public static class NetClassExtensionMethods
    {
        /// <summary>
        /// Exporting individual geometries to element classes
        /// </summary>
        /// <param name="_triMesh"></param>
        /// <param name="_workSpacePath"></param>
        /// <param name="_fileName"></param>
        public static void ExportSimpleGeometryToShapfileByNet(this IGeometry geometry, string _workSpacePath, string _fileName)
        {
            
            var feature = new Feature(geometry, new AttributesTable(new Dictionary<string, object>
        {
            {"id", 1},
            {"name", "p1"}
        }));

            var vt = new System.Collections.ObjectModel.Collection<IFeature>();
            vt.Add(feature as IFeature);
            var featureCollection = new FeatureCollection(vt);
            
            var writer = new ShapefileDataWriter(_workSpacePath + "\\" + _fileName);
            writer.Header = ShapefileDataWriter.GetHeader(featureCollection.Features[0], 1);
            writer.Write(featureCollection.Features);

            
        }

        /// <summary>
        /// Exporting multiple Geometries to an element class
        /// </summary>
        /// <param name="_triMesh"></param>
        /// <param name="_workSpacePath"></param>
        /// <param name="_fileName"></param>
        public static void ExportGeometryToShapfileByNet(this List<IGeometry> geometryCollection, string _workSpacePath, string _fileName)
        {
            if (geometryCollection.Count <= 0)
                return;
                        
            var features = new List<IFeature>();
           
            DbaseFieldDescriptor field1 = new DbaseFieldDescriptor();
            field1.Name = "id";
            field1.DbaseType = 'N';
            field1.Length = 10;
            field1.DecimalCount = 0;

            DbaseFieldDescriptor field2 = new DbaseFieldDescriptor();
            field2.Name = "name";
            field2.DbaseType = 'C';
            field2.Length = 20;
            field2.DecimalCount = 0;

            DbaseFieldDescriptor[] fields = new DbaseFieldDescriptor[]{field1,field2};

            
            for(int i=0;i< geometryCollection.Count;i++)
            {
                
                var attributes = new AttributesTable();
                attributes.Add("id", i);
                attributes.Add("name", "m"+i);
                var feature1 = new Feature(geometryCollection[i], attributes);
                features.Add(feature1);
            }


            
            var writer = new ShapefileDataWriter(_workSpacePath + "\\" + _fileName);
            writer.Header = ShapefileDataWriter.GetHeader(fields, features.Count);
            writer.Write(features);
        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="_shpFile"></param>
        /// <param name="fields"></param>
        /// <param name="fieldInformation"></param>
        /// <returns></returns>
        public static List<IGeometry> GetAllGeometry(string _shpFile, List<string> fields, ref List<string[]> fieldInformation)
        {

            var reader = new NetTopologySuite.IO.ShapefileDataReader(_shpFile, NetTopologySuite.Geometries.GeometryFactory.Default);

            int[] flags = new int[fields.Count];

            int count = 0;

            foreach (var vt in fields)
            {
                for (int i = 0; i < reader.DbaseHeader.NumFields; i++)
                {
                    var fieldName = reader.DbaseHeader.Fields[i].Name;

                    if (vt == fieldName)
                    {
                        flags[count] = i;
                        count++;
                        break;
                    }
                }
            }

            List<IGeometry> pnewGeos = new List<IGeometry>();

            List<string> types = new List<string>();

            while (reader.Read())
            {
                var geometry = reader.Geometry;

                string gtype = geometry.GeometryType;

                string[] fieldvalues = new string[flags.Length];

                int ct = 0;
                foreach (var vt in flags)
                {
                    var fieldValue = reader.GetValue(vt + 1).ToString();
                    if (fieldValue.Contains("\0"))
                    {
                        fieldvalues[ct] = "0";
                    }
                    else
                    {
                        fieldvalues[ct] = fieldValue;
                    }

                    ct++;
                }
                fieldInformation.Add(fieldvalues);

                types.Add(gtype);
                pnewGeos.Add(geometry);
            }
           
            reader.Close();
           
            return pnewGeos;

        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="_shpFile"></param>
        /// <returns></returns>
        public static List<IGeometry> GetAllGeometry(string _shpFile)
        {

            var reader = new NetTopologySuite.IO.ShapefileDataReader(_shpFile, NetTopologySuite.Geometries.GeometryFactory.Default);
          
            List<IGeometry> pnewGeos = new List<IGeometry>();      

            while (reader.Read())
            {
                var geometry = reader.Geometry;
                pnewGeos.Add(geometry);
            }

            reader.Close();

            return pnewGeos;

        }

    } 
}
