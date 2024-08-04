using OSGeo.GDAL;
using OSGeo.OGR;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GdalLib
{
    public class DEMRaster
    {
        
        private string _rasterFilePath;

       
        public Dataset _ds;

        public DEMRaster()
        {
            
        }

        public DEMRaster(string filePath)
        {
            _rasterFilePath = filePath;
            
            _ds = Gdal.Open(_rasterFilePath, Access.GA_Update);

                       
        }


        public DEMRaster(string directoryPath, string fileName)
        {
            _rasterFilePath = directoryPath + "\\" + fileName;
            
            _ds = Gdal.Open(_rasterFilePath, Access.GA_Update);

        }


       
        public double GetPixel()
        {
            double[] adfGeoTransform = new double[6];
            _ds.GetGeoTransform(adfGeoTransform);

            return adfGeoTransform[1];
        }

      
        public double GetElevation(double dProjX, double dProjY)
        {
            try
            {
                Band Band = _ds.GetRasterBand(1);
                             
                int width = Band.XSize;
                int height = Band.YSize;

                
                double[] adfGeoTransform = new double[6];
                _ds.GetGeoTransform(adfGeoTransform);
                
                double dTemp = adfGeoTransform[1] * adfGeoTransform[5] - adfGeoTransform[2] * adfGeoTransform[4];
                double dCol = 0.0, dRow = 0.0;
                dCol = (adfGeoTransform[5] * (dProjX - adfGeoTransform[0]) -
                    adfGeoTransform[2] * (dProjY - adfGeoTransform[3])) / dTemp + 0.5;
                dRow = (adfGeoTransform[1] * (dProjY - adfGeoTransform[3]) -
                    adfGeoTransform[4] * (dProjX - adfGeoTransform[0])) / dTemp + 0.5;
                int dc = Convert.ToInt32(dCol);
                int dr = Convert.ToInt32(dRow);


               
                double[] data = new double[1 * 1];
                CPLErr err = Band.ReadRaster(dc, dr, 1, 1, data, 1, 1, 0, 0);
                Band.Dispose();
                double elvate = data[0];
                return elvate;
            }
            catch
            {
                return 0.0;
            }
        }


       
        public void ModifyElevation(List<Geometry> pLines, string _filepath, string _demname)
        {
            OSGeo.GDAL.Driver driver = _ds.GetDriver();

            string projection = _ds.GetProjection(); 

            int width = _ds.RasterXSize;
            int height = _ds.RasterYSize;

            Band Band = _ds.GetRasterBand(1);

            double[] adfGeoTransform = new double[6];
            _ds.GetGeoTransform(adfGeoTransform);

            double dTemp = adfGeoTransform[1] * adfGeoTransform[5] - adfGeoTransform[2] * adfGeoTransform[4];
            double dCol = 0.0, dRow = 0.0;


            for (int i = 0; i < pLines.Count; i++)
            {
                for (int k = 0; k < pLines[i].GetPointCount(); k++)

                {
                    double dProjX = pLines[i].GetX(k);
                    double dProjY = pLines[i].GetY(k);

                    dCol = (adfGeoTransform[5] * (dProjX - adfGeoTransform[0]) -
                        adfGeoTransform[2] * (dProjY - adfGeoTransform[3])) / dTemp + 0.5;
                    dRow = (adfGeoTransform[1] * (dProjY - adfGeoTransform[3]) -
                        adfGeoTransform[4] * (dProjX - adfGeoTransform[0])) / dTemp + 0.5;
                    int dc = Convert.ToInt32(dCol);
                    int dr = Convert.ToInt32(dRow);

                    
                    double[] data = new double[1 * 1];

                    Band.ReadRaster(dc, dr, 1, 1, data, 1, 1, 0, 0);
                    data[0] = pLines[i].GetZ(k);
                    Band.WriteRaster(dc, dr, 1, 1, data, 1, 1, 0, 0);
                }
            }

            float[] newdata = new float[width * height];
            Band.ReadRaster(0, 0, width, height, newdata, width, height, 0, 0);

           for(int i=0;i< newdata.Length;i++)
            {
                if (newdata[i] > 2000)
                    newdata[i] = 0;
            }

            List<Geometry> pAllGeometry = new List<Geometry>();

            for (int row = 0; row < height; row++)
            {
                for (int col = 0; col < width; col++)
                {
                    Geometry pt = new Geometry(wkbGeometryType.wkbPoint);

                    double x = adfGeoTransform[0] + (col) * adfGeoTransform[1] + (row ) * adfGeoTransform[2];
                    double y = adfGeoTransform[3] + (col) * adfGeoTransform[4] + (row ) * adfGeoTransform[5];


                 
                    double[] data = new double[1 * 1];
                    Band.ReadRaster(col, row, 1, 1, data, 1, 1, 0, 0);

                    double z = -1;
                    if (data[0] > 2000)
                    {
                        z = 0.0;
                    }
                    else
                    {
                        z = data[0];
                    }
                    pt.AddPoint(x,y,z);

                    pAllGeometry.Add(pt);
                }
            }

            pAllGeometry.ExportGeometryToShapfile(_filepath,"elevationPoints");


            string outputFilename = _filepath + "\\" + _demname;
            Dataset outputDataset = driver.Create(outputFilename, width, height, 1, DataType.GDT_Float32, new string[] { }); 

            outputDataset.SetGeoTransform(adfGeoTransform); 
            outputDataset.SetProjection(projection); 

            Band outputBand = outputDataset.GetRasterBand(1); 
            outputBand.WriteRaster(0, 0, width, height, newdata, width, height, 0, 0); 

            outputBand.FlushCache(); 
            outputDataset.FlushCache(); 

            outputBand.Dispose(); 
            outputDataset.Dispose(); 

            Band.Dispose();

        }



        /// <summary>
        /// read raster
        /// </summary>
        public static void ReadRaster(string strFile)
        {
            
            Gdal.AllRegister();

            Dataset ds = Gdal.Open(strFile, Access.GA_ReadOnly);

            if (ds == null)
            {
                Console.WriteLine("can not open：" + strFile);
                System.Environment.Exit(-1);
            }

            OSGeo.GDAL.Driver drv = ds.GetDriver();
            if (drv == null)
            {
                Console.WriteLine("can not open：" + strFile);
                System.Environment.Exit(-1);
            }

            Console.WriteLine("RasterCount:" + ds.RasterCount);
            Console.WriteLine("RasterSize:" + ds.RasterXSize + " " + ds.RasterYSize);
        }

       
    }
}
