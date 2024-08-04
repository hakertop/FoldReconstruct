using OSGeo.GDAL;
using OSGeo.OGR;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using System.Windows.Forms;

using FoldRestortEntry;
using NetTopologySuite;
using GdalLib;
using GeologicalEntity;


namespace TectonicReconstructEntry
{
    class Program
    {
        [STAThread]
        static void Main(string[] args)
        {

            Gdal.AllRegister();
            Ogr.RegisterAll();            
            OSGeo.GDAL.Gdal.SetConfigOption("GDAL_FILENAME_IS_UTF8", "YES");          
            OSGeo.GDAL.Gdal.SetConfigOption("SHAPE_ENCODING", "UTF-8");


            // An example of Wulongshan Dome

            // please change the below directory according to your path @"..\TectonicReconstruction\Output"
            FPRMain pSege = new FPRMain(@"E:\PhD\MyEssay\SecondValley\Manuscripts\TectonicReconstruction\Output");

            // please change the below directory according to your path @"..\TectonicReconstruction\data_Dome\WuLongShanDEM2.tif"           
            DEMRaster pRaster = new DEMRaster(@"E:\PhD\MyEssay\SecondValley\Manuscripts\TectonicReconstruction\data_Dome\WuLongShanDEM2.tif");

            // please change the below directory according to your path @"..\TectonicReconstruction\data_Dome\ProjectSubDomeStrata1.shp"
            string shpFile = @"E:\PhD\MyEssay\SecondValley\Manuscripts\TectonicReconstruction\data_Dome\ProjectSubDomeStrata1.shp";

            // please change the below directory according to your path @"..\TectonicReconstruction\data_Dome\ProjectSubDomeAttitudes1.shp"
            string occurrenceFile = @"E:\PhD\MyEssay\SecondValley\Manuscripts\TectonicReconstruction\data_Dome\ProjectSubDomeAttitudes1.shp";
            
            // run this example
            pSege.DomeTest(shpFile, occurrenceFile, pRaster, 2.0, -1);

            MessageBox.Show("Run Successfully");

        }
    }
}
