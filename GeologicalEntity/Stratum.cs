using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using GeoCommon;

using GeoAPI;
using GeoAPI.Geometries;

namespace GeologicalEntity
{
    public enum StratumStyle
    {

    }
   public class Stratum
    {

        /// <summary>
        /// stratum id
        /// </summary>
        public int SID { get; set; }
        /// <summary>
        /// stratum name
        /// </summary>
        public string SName { get; set; }

        /// <summary>
        /// stratum code
        /// </summary>
        public string SCode { get; set; }

        /// <summary>
        /// stratum age
        /// </summary>
        public string SAge { get; set; }

        /// <summary>
        /// stratum style
        /// </summary>
        public StratumStyle SStyle { get; set; }

        
        /// <summary>
        /// geometric attitude (use the Geometry class in library GDAL)
        /// </summary>
        public Polygon3D SPolygon { get; set; }

        /// <summary>
        /// geometric attitude(use the Geometry class in library NetTopologySuite)
        /// </summary>
        public IGeometry NPolygon { get; set; }

        /// <summary>
        /// constructor 1
        /// </summary>
        public Stratum()
        {

        }

        /// <summary>
        /// constructor 2
        /// </summary>
        /// <param name="_sid"></param>
        /// <param name="_scode"></param>
        /// <param name="_spolygon"></param>
        public Stratum(int _sid, string _scode,Polygon3D _spolygon)
        {
            this.SID = _sid;
            this.SCode = _scode;
            this.SPolygon = _spolygon;
        }

        /// <summary>
        /// constructor 3 Using
        /// </summary>
        /// <param name="_sid"></param>
        /// <param name="_scode"></param>
        /// <param name="_npolygon"></param>
        public Stratum(int _sid, string _scode, IGeometry _npolygon)
        {
            this.SID = _sid;
            this.SCode = _scode;
            this.NPolygon = _npolygon;
        }
    }
}
