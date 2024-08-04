using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GeoCommon
{
    /// <summary>
    /// 2D Vertex
    /// </summary>
    public class Vertex3D
    {
        /// <summary>
        /// name
        /// </summary>
        public string name;

        /// <summary>
        /// 
        /// </summary>
        public string belongToEage;

        /// <summary>
        /// 
        /// </summary>
        public string belongToPolygon;

        /// <summary>
        /// 
        /// </summary>
        public int id;

        /// <summary>
        /// 
        /// </summary>
        public double x;

        /// <summary>
        /// 
        /// </summary>
        public double y;

        /// <summary>
        /// 
        /// </summary>
        public double z;

        public Vertex3D()
        {

        }

        public Vertex3D(double _x, double _y,double _z)
        {
            this.x = _x;
            this.y = _y;
            this.z = _z;
        }

        public Vertex3D(double _x, double _y)
        {
            this.x = _x;
            this.y = _y;
            this.z = 0.0;
        }
    }
}
