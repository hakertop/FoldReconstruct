using OSGeo.OGR;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeologicalEntity
{
   public class OccurrencePoint: OSGeo.OGR.Geometry
    {
        /// <summary>
        /// dip angle
        /// </summary>
        public double dipAngle;

        /// <summary>
        /// dip
        /// </summary>
        public double tendency;

        /// <summary>
        /// strike
        /// </summary>
        public double strike;

        /// <summary>
        /// occurrence or attitude
        /// </summary>
        /// <param name="_geometryType">geometric type</param>
        public OccurrencePoint(wkbGeometryType _geometryType):base(_geometryType)
        {
             
        }
    }
}
